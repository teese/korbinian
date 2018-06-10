import ast
import csv
import os
import sys
import time
import unicodedata
from multiprocessing import Pool

import numpy as np
import pandas as pd

import korbinian
import korbinian.utils as utils
# import debugging tools
from korbinian.utils import pn, aaa


def prepare_protein_list(s, pathdict, logging):
    """ Sets up the file locations in the DataFrame containing the list of proteins for analysis.

    Parameters
    ----------
    s : dict
        Dictionary of settings derived from settings excel file.
        Columns "Parameter" and "Value" are converted to key/value in the dictionary, respectively.
    pathdict : dict
        Dictionary of the key paths and files associated with that List number.

    DataFrames
    ----------
    df = from parsed_csv.
        index = accessions
        columns = full_seq, topo, seqlen, M_indices etc

    df_PLS = Protein List Summary, which summarises filtering steps
                                                       v  date
        list_number                                   37  2018
        n_prot_AFTER_dropping_SCAMPI_nonTM_seqences  285  2018
        n_prot_AFTER_dropping_non_trusted_SiPe       265  2018
        n_prot_AFTER_dropping_without_list_TMDs      261  2018
        n_prot_AFTER_dropping_without_n_TMDs         261  2018
        n_prot_AFTER_dropping_with_X_in_seq          248  2018
        n_prot_AFTER_dropping_missing_TM_residues     248  2018

    Saved Files and Figures
    -----------------------
    pathdict["list_summary_csv"] : csv file
        Input CSV file is overwritten at end of function, including the extra file locations.

    """
    logging.info('~~~~~~~~~~~~                     starting prepare_protein_list                      ~~~~~~~~~~~~')

    # create dataframe for protein list summary (df_PLS) to record filtering steps
    df_PLS = pd.DataFrame(columns=["v", "date"])
    # get the timestamp for current time
    t = time.ctime(time.time())
    df_PLS.loc["list_number", :] = (s["list_number"], t)

    # DEPRECATED. AM NOT SURE WHY I WOULD NEED TO OPEN AN OLD LIST SUMMARY.
    # # create or open dataframe for protein list summary
    # if os.path.isfile(pathdict["prot_list_summary_csv"]):
    #     df_PLS = pd.read_csv(pathdict["prot_list_summary_csv"], index_col=0)
    # else:
    #     df_PLS = pd.DataFrame(columns=["v", "date"])

    df = pd.read_csv(pathdict["list_csv"], sep=",", quoting=csv.QUOTE_NONNUMERIC, index_col=0, low_memory=False)

    # add locations of the saved files for each protein to the dataframe
    df = add_file_paths_to_list_csv(df, s)

    # number of initial proteins before filtering
    n_initial_prot = df.shape[0]
    logging.info('\nn_initial_prot: {}'.format(n_initial_prot))

    # # convert to python tuple from stringlist
    df.TM_indices = df.TM_indices.dropna().apply(lambda x : tuple(ast.literal_eval(x)))

    ########################################################################################
    #                                                                                      #
    #     exclusion of proteins not judged to be membranous according to SCAMPI            #
    #                                                                                      #
    ########################################################################################

    # drop proteins that are non-transmembrane according to SCAMPI
    # get accession numbers of transmembrane proteins from SCAMPI output list query.TM_list.txt
    SCAMPI_nonTM_path = pathdict['SCAMPI_nonTM']
    n_prot_AFTER_dropping_SCAMPI_nonTM_seqences = 'SCAMPI nonTMD file not found!'
    if os.path.isfile(SCAMPI_nonTM_path):
        modification_date = time.ctime(os.path.getmtime(SCAMPI_nonTM_path))
        SCAMPI_nonTM_list = []
        with open(SCAMPI_nonTM_path) as source:
            for line in source:
                line = line.strip()
                SCAMPI_nonTM_list.append(line)
        keep_after_SCAMPI = []
        for acc in df.index:
            if not acc in SCAMPI_nonTM_list:
                keep_after_SCAMPI.append(acc)
        df = df.loc[keep_after_SCAMPI,:]
        #df = df.drop(SCAMPI_nonTM_list, axis=0)
        n_prot_AFTER_dropping_SCAMPI_nonTM_seqences = df.shape[0]
    else:
        modification_date = None

    # save and log filtering results
    df_PLS.loc["n_prot_AFTER_dropping_SCAMPI_nonTM_seqences",:] = (n_prot_AFTER_dropping_SCAMPI_nonTM_seqences, t)
    logging.info('n_prot_AFTER_dropping_SCAMPI_nonTM_seqences: {}'.format(n_prot_AFTER_dropping_SCAMPI_nonTM_seqences))
    if modification_date is not None:
        logging.info("modification_date of scampi file: {} ---".format(modification_date))

    ########################################################################################
    #                                                                                      #
    #         exclusion of proteins whose TM01 may in fact be a signal peptide            #
    #                                                                                      #
    ########################################################################################
    # load SignalP file containing all proteins that have signal peptides
    SignalP_SiPe_path = pathdict['SignalP_SiPe_acc']
    # create variable to hold the date of the SignalP analysis
    modification_date_SignalP = None

    if "betabarrel" not in df.columns:
        df["betabarrel"] = False

    if os.path.isfile(SignalP_SiPe_path) and True not in df.betabarrel.tolist():
        ########################################################################################
        #                                                                                      #
        #         TMSEG alpha-helical exclusion of non-trusted signal peptides                 #
        #                                                                                      #
        ########################################################################################
        modification_date_SignalP = time.ctime(os.path.getmtime(SignalP_SiPe_path))
        SignalP_acc_containing_SiPe = korbinian.prot_list.SCAMPI.get_SignalP_SiPe_acc(SignalP_SiPe_path)
        # create column containing bool if signal peptide is predicted by SignalP
        for acc in df.index:
            if acc in SignalP_acc_containing_SiPe:
                df.at[acc, 'SignalP_SiPe'] = True
            else:
                df.at[acc, 'SignalP_SiPe'] = False
        # drop all proteins where Uniprot did not get signal peptide
        TM01_potential_SiPe_acc_list = []

        # repeat filter
        acc_with_topo = df.query("TM_indices != () & number_of_TMDs > 1").index.tolist()

        for acc in acc_with_topo:
            TM01_start = df.at[acc, "TM_indices"][0][0]
            # if UniProt DOESN'T think there is a signal peptide, but SignalP does
            # there is a chance that the "TM01" is really a signal peptide
            if df.at[acc, 'uniprot_SiPe'] == False and df.at[acc, 'SignalP_SiPe'] == True:
                # but it's only a danger if the TM01 is in the first 40 residues
                if TM01_start < 15:
                    TM01_potential_SiPe_acc_list.append(acc)

        df = df.drop(TM01_potential_SiPe_acc_list, axis=0)

    elif True in df.betabarrel.tolist():
        # for old parsed data, add BB_SiPe (can be removed after re-running all parsing)
        if "BB_SiPe" not in df.columns:
            df["BB_SiPe"] = df.SP01_start.notnull()
        ########################################################################################
        #                                                                                      #
        #                betabarrel exclusion of non-trusted signal peptides                   #
        #                                                                                      #
        ########################################################################################
        # first select those where no signal peptide was found
        df_removed = df.loc[df["BB_SiPe"] == False]
        # from these, select those where the TM01_start is in the first 15 residues
        df_removed = df_removed.loc[df_removed["TM01_start"] < 15]
        # make a list of those removed
        TM01_potential_SiPe_acc_list = list(df_removed.index)
        # drop and note change in length of the dataframe[to be replaced later with a bool column]
        df = df.drop(TM01_potential_SiPe_acc_list, axis=0)
    else:
        #n_prot_AFTER_dropping_non_trusted_SiPe = 'SignalP output file not found!'
        TM01_potential_SiPe_acc_list = []
        logging.info('SignalP filter not conducted : SignalP_SiPe_path file not found!'.format(TM01_potential_SiPe_acc_list))

    # save and log filtering results
    n_prot_AFTER_dropping_non_trusted_SiPe = df.shape[0]
    df_PLS.loc["n_prot_AFTER_dropping_non_trusted_SiPe",:] = (n_prot_AFTER_dropping_non_trusted_SiPe, t)
    logging.info('n_prot_AFTER_non_trusted_Signal_Peptides: {}'.format(n_prot_AFTER_dropping_non_trusted_SiPe))
    if len(TM01_potential_SiPe_acc_list) > 0:
        logging.info('TM01_potential_SiPe_acc_list: {}'.format(TM01_potential_SiPe_acc_list))
    if modification_date_SignalP is not None:
        logging.info("modification_date of SignalP file: {} ---".format(modification_date_SignalP))

    # KEPT FOR LEGACY REASONS ONLY. FILES SHOULD ALWAYS BE NAMED BASED ON ACC ONLY.
    df['protein_name'] = df.index.astype(str)

    if s["add_user_subseqs"] == True:
        ########################################################################################
        #                                                                                      #
        #      Add user-selected sequences from csv or excel("SE01", "SE02" etc)               #
        #                                                                                      #
        ########################################################################################
        df = add_user_subseqs(df, pathdict)

    # convert the list_of_TMDs to a python object, if it is a string
    if isinstance(df['list_of_TMDs'].dropna().iat[0], str):
        for colname in ['list_of_TMDs', 'list_of_TMDs_excl_SP']:
            df[colname] = df[colname].replace("[]", np.nan)
            df[colname] = df[colname].replace("nan", np.nan)
            df[colname] = df[colname].dropna().apply(lambda x: ast.literal_eval(x))
        #df['list_of_TMDs'] = df['list_of_TMDs'].dropna().apply(lambda x : ast.literal_eval(x))
        #df['list_of_TMDs_excl_SP'] = df['list_of_TMDs_excl_SP'].dropna().apply(lambda x: ast.literal_eval(x))
    #drop any proteins that do not have a valid list_of_TMDs
    df = df.loc[df['list_of_TMDs'].notnull()]
    df = df.loc[df['list_of_TMDs_excl_SP'].notnull()]

    # save and log filtering results
    n_prot_AFTER_dropping_without_list_TMDs = df.shape[0]
    df_PLS.loc["n_prot_AFTER_dropping_without_list_TMDs", :] = (n_prot_AFTER_dropping_without_list_TMDs, t)
    logging.info('n_prot_AFTER_dropping_without_list_TMDs: {}'.format(n_prot_AFTER_dropping_without_list_TMDs)) # line 107

    # drop any proteins where the topology source does not match the required source of TM definitions
    df = df.loc[df["topol_source"] == s["TM_def"]]
    n_prot_AFTER_filt_topology_source = df.shape[0]
    df_PLS.loc["n_prot_AFTER_filt_topology_source", :] = (n_prot_AFTER_filt_topology_source, t)

    ########################################################################################
    #                                                                                      #
    #                drop based on the accepted number of TMDs for that dataset            #
    #              (based on min_TMDs and max_TMDs in "list" tab of settings excel file)   #
    #                                                                                      #
    ########################################################################################
    min_TMDs = s["min_TMDs"]
    max_TMDs = s["max_TMDs"]

    # drop any proteins without a number of TMDs
    df.dropna(subset=["number_of_TMDs"], inplace=True)
    df_PLS.loc["n_prot_AFTER_dropping_without_n_TMDs",:] = (df.shape[0], t)

    # filter by number of TMDs, based on min and max in lists tab of excel file
    df = df.loc[df["number_of_TMDs_excl_SP"].apply(lambda x: min_TMDs <= x <= max_TMDs)]
    n_prot_AFTER_n_TMDs_cutoff = df.shape[0]
    if n_prot_AFTER_n_TMDs_cutoff == 0:
        raise ValueError("The list {} has no valid proteins after number_of_TMDs_excl_SP cutoff.".format(s["list_number"]))
    logging.info('n_prot_AFTER_n_TMDs_cutoff: {}'.format(n_prot_AFTER_n_TMDs_cutoff)) # line 215

    ########################################################################################
    #                                                                                      #
    #     slice out all the TMD_seq_plus_surr, based on settings (e.g. 10aa each side)     #
    #                                                                                      #
    ########################################################################################

    max_num_TMDs = int(df["number_of_TMDs_excl_SP"].max())
    max_list_TMDs = ["TM{:02d}".format(TMD_Nr) for TMD_Nr in range(1, max_num_TMDs + 1)]

    n_aa_before_tmd = s["n_aa_before_tmd"]
    n_aa_after_tmd = s["n_aa_after_tmd"]

    # calculate TM plus surr for ALL sequences, overwriting if necessary, in case this is changed later
    # currently the loop is run for each TMD, based on the sequence with the most TMDs
    for i in range(1, int(max_num_TMDs) + 1):
        TMD = 'TM%02d' % i
        # get the indices for TMD plus surrounding sequence
        df = korbinian.prot_list.prot_list.get_indices_TMD_plus_surr_for_summary_file(df, TMD, n_aa_before_tmd, n_aa_after_tmd)
        # slice out the TMD_seq_plus_surr for each TMD
        df['%s_seq_plus_surr' % TMD] = df[df['%s_start' % TMD].notnull()].apply(utils.slice_uniprot_TMD_plus_surr_seq, args=(TMD,), axis=1)

    # create column with all TMD residues joined togethor
    seq_cols = ["{}_seq".format(TMD) for TMD in max_list_TMDs]
    df['TMD_seq_joined'] = df.loc[:,seq_cols].fillna("").sum(axis=1)
    df['n_TMD_res'] = df['TMD_seq_joined'].str.len()

    # calculate percentage length of TMD region in full protein                   #
    df['TMD_seq_joined_len'] = df['TMD_seq_joined'].str.len()
    df['perc_TMD'] = df['TMD_seq_joined_len'] / df['seqlen'] * 100

    # list of TMDs excluding TM01
    df["list_of_TMDs_excl_TM01"] = df.loc[df["number_of_TMDs_excl_SP"] >= 2]["list_of_TMDs_excl_SP"].apply(lambda x: x[1:])

    len_cols = ["{}_seqlen".format(TMD) for TMD in max_list_TMDs]
    df['len_TMD_mean'] = df.loc[:, len_cols].mean(axis=1)

    df['last_TMD'] = df['list_of_TMDs_excl_SP'].dropna().apply(lambda x : x[-1])

    df_min_3_TMDs = df.loc[df["number_of_TMDs_excl_SP"] >= 3]
    df["list_of_central_TMDs"] = df_min_3_TMDs["list_of_TMDs_excl_SP"].apply(lambda x : x[1:-1])

    ########################################################################################
    #                                                                                      #
    #      Calculate lipophilicity for all TMDs in list               #
    #                                                                                      #
    ########################################################################################
    for TMD_Nr in range(1, max_num_TMDs + 1):
        TMD = "TM{:02d}".format(TMD_Nr)
        # calculate TMD length
        df["{}_seqlen".format(TMD)] = df["{}_seq".format(TMD)].str.len()
        # calculate lipophilicity
        df["{}_lipo".format(TMD)] = df["{}_seq".format(TMD)].dropna().apply(utils.calc_lipophilicity)
    if "SP01_seq" in df.columns:
        df["SP01_lipo".format(TMD)] = df["SP01_seq".format(TMD)].dropna().apply(utils.calc_lipophilicity)

    lipo_cols = ["{}_lipo".format(TMD) for TMD in max_list_TMDs]
    df['lipo_mean_all_TMDs_mean'] = df.loc[:, lipo_cols].mean(axis=1)

    df["lipo_mean_all_TM_res"] = df['TMD_seq_joined'].apply(utils.calc_lipophilicity)

    ########################################################################################
    #                                                                                      #
    #      Calculate lipophilicity for first, last, and central TMDs in list               #
    #                                                                                      #
    ########################################################################################
    # create an empty column, for single-pass etc that do not have a calculated value
    df['lipo_mean_central_TMDs'] = np.nan
    df['lipo_mean_excl_TM01'] = np.nan

    # add a column that holds all joined TMD sequences
    for n, acc in enumerate(df.index):
        #if n % 20 == 0 and n != 0:
        #    sys.stdout.write('.'), sys.stdout.flush()
        #    if n % 600 == 0:
        #        sys.stdout.write('\n'), sys.stdout.flush()

        list_of_TMDs = df.at[acc, 'list_of_TMDs']
        list_of_TMDs_excl_SP = df.at[acc, "list_of_TMDs_excl_SP"]
        last_TMD = df.at[acc, 'last_TMD']
        # get lipo of last TMD

        df.at[acc, 'lipo_last_TMD'] = df.at[acc, '%s_lipo' % last_TMD]
        df.at[acc, 'last_TMD_seq'] = df.at[acc, '%s_seq' % last_TMD]
        df.at[acc, 'last_TMD_seqlen'] = df.at[acc, '%s_seqlen' % last_TMD]

        if df.at[acc, "number_of_TMDs_excl_SP"] >= 3:
        # calculate mean lipophilicity excluding first and last TM
            list_of_central_TMDs_excl_SP = list_of_TMDs_excl_SP[1:-1]
            lipo_cols_central = ["{}_lipo".format(TMD) for TMD in list_of_central_TMDs_excl_SP]
            lipo_sel_ser_central = df.loc[acc,lipo_cols_central]
            df.at[acc, 'lipo_mean_central_TMDs'] = lipo_sel_ser_central.mean()

        if df.at[acc, "number_of_TMDs_excl_SP"] >= 2:
        # calculate mean lipophilicity excluding first and last TM
            list_of_TMDs_excl_SP_excl_TM01 = list_of_TMDs_excl_SP[1:]
            lipo_cols_excl_TM01 = ["{}_lipo".format(TMD) for TMD in list_of_TMDs_excl_SP_excl_TM01]
            lipo_sel_ser_excl_TM01 = df.loc[acc,lipo_cols_excl_TM01]
            df.at[acc, 'lipo_mean_excl_TM01'] = lipo_sel_ser_excl_TM01.mean()
    #sys.stdout.write("\n")

    # drop any proteins with X in sequence
    df['X_in_seq'] = df["full_seq"].str.contains("X")
    list_acc_X_in_seq = df.loc[df['X_in_seq']].index.tolist()
    df = df.loc[df['X_in_seq'] != True]
    n_prot_AFTER_dropping_with_X_in_seq = df.shape[0]
    df_PLS.loc["n_prot_AFTER_dropping_with_X_in_seq",:] = (n_prot_AFTER_dropping_with_X_in_seq, t)
    logging.info('n_prot_AFTER_dropping_with_X_in_seq: {}'.format(n_prot_AFTER_dropping_with_X_in_seq)) # line 297
    if len(list_acc_X_in_seq) > 0:
        logging.info('list_acc_X_in_seq: {}'.format(list_acc_X_in_seq))

    # drop proteins without any TMD residues
    # this is an old filter for unusual cases in UniProt, such as when TM01_start is listed as "?" or ">30" in UniProt record.
    df = df.loc[df['n_TMD_res'] >= 1]
    n_prot_AFTER_dropping_missing_TM_residues = df.shape[0]
    df_PLS.loc["n_prot_AFTER_dropping_missing_TM_residues",:] = (n_prot_AFTER_dropping_missing_TM_residues, t)
    logging.info('n_prot_AFTER_dropping_missing_TM_residues: {}'.format(n_prot_AFTER_dropping_missing_TM_residues)) # line 297

    ########################################################################################
    #                                                                                      #
    #           check for subcellular location (PM, ER, Golgi) in keywords                 #
    #                                                                                      #
    ########################################################################################
    KW_search_list = ['Cell_membrane', 'Endoplasmic_reticulum', 'Golgi_apparatus']

    if 'uniprot_KW' in df.columns:
        # apply ast.literal_eval to every item in df['uniprot_KW']
        if isinstance(df['uniprot_KW'].iloc[0], str):
            df['uniprot_KW'] = df['uniprot_KW'].apply(lambda x: ast.literal_eval(x))
        for KW in KW_search_list:
            df[KW] = df['uniprot_KW'].apply(korbinian.cons_ratio.keywords.KW_list_contains_any_desired_KW, args=([KW],))
        # check for specific keywords with altered names in final file
        df['GPCR'] = df['uniprot_KW'].apply(korbinian.cons_ratio.keywords.KW_list_contains_any_desired_KW, args=(['G-protein coupled receptor'],))
        df['olfactory_receptor'] = df['prot_descr'].apply(korbinian.cons_ratio.keywords.KW_list_contains_any_desired_KW, args=(['Olfactory receptor'],))

    else:
        df['GPCR'] = False
        df['olfactory_receptor'] = False
        for KW in KW_search_list:
            df[KW] = False

    ########################################################################################
    #                                                                                      #
    #             count signal peptides and save a new column "SiPe_in_dataset"            #
    #                                                                                      #
    ########################################################################################
    # number of proteins with signal_peptides
    n_SP_ser = df["number_of_TMDs"] - df["number_of_TMDs_excl_SP"]
    n_prot_with_SP = n_SP_ser.sum()
    logging.info("\nnumber of proteins with signal peptides = {}".format(n_prot_with_SP))
    if "SiPe" in s["regions"] and n_prot_with_SP >= 1:
        df["SiPe_in_dataset"] = True
    else:
        df["SiPe_in_dataset"] = False

    # save to a csv
    df.to_csv(pathdict["list_csv"], sep=",", quoting=csv.QUOTE_NONNUMERIC)
    df_PLS.to_csv(pathdict["prot_list_summary_csv"])

    ########################################################################################
    #                                                                                      #
    #                      calculate random TM and nonTM identity                          #
    #                                                                                      #
    ########################################################################################
    # only calculate the rand_TM and rand_nonTM using the mathematical method
    calc_randTM_and_randnonTM(s, pathdict, logging)
    logging.info('~~~~~~~~~~~~                     finished prepare_protein_list                      ~~~~~~~~~~~~')



def get_topology_for_prot_list(s, pathdict, logging):
    """ Gets membrane protein topology from alternative sources, and slices out TMD sequences from full sequence.

    Parameters
    ----------
    s : dict
        Dictionary of settings derived from settings excel file.
        Columns "Parameter" and "Value" are converted to key/value in the dictionary, respectively.
    pathdict : dict
        Dictionary of the key paths and files associated with that List number.
    logging : logging.Logger
        Logger for printing to console and logfile.

    DataFrames
    ----------
    df = from parsed_csv.
        index = accessions
        columns = full_seq, topo, seqlen, M_indices etc

    Saved Files and Figures
    -----------------------
    IMPORTANT: saves as ListXX.csv, which will be overwritten by prepare_protein_list

    """
    logging.info('~~~~~~~~~~~~                 starting get_topology_for_prot_list                   ~~~~~~~~~~~~')

    df = pd.read_csv(pathdict["list_parsed_csv"], sep = ",", quoting = csv.QUOTE_NONNUMERIC, index_col = 0, low_memory=False)
    n_initial_prot = df.shape[0]
    # convert to python tuple from stringlist
    df.TM_indices = df.TM_indices.dropna().apply(lambda x : tuple(ast.literal_eval(x)))

    predictions_dir = os.path.join(s["data_dir"], "predictions")

    ########################################################################################
    #                                                                                      #
    #                     get TM topology from TMSEG output files                          #
    #  (only for cases where the TM_def in the lists tab of the settings file is "TMSEG"   #
    #                                                                                      #
    ########################################################################################
    if "TMSEG" in s['TM_def']:
        list_acc_without_TMSEG = []
        for acc in df.index:
            sys.stdout.write("."), sys.stdout.flush()
            TMSEG_file_is_good = True
            full_seq = df.at[acc, "full_seq"]
            TMSEG_fastastyle_file = os.path.join(predictions_dir, acc[0:2], "{}_TMSEG_fastastyle.txt".format(acc))
            TMSEG_file = os.path.join(predictions_dir, acc[0:2], "{}_TMSEG.txt".format(acc))

            if any([os.path.isfile(TMSEG_file), os.path.isfile(TMSEG_fastastyle_file)]):
                if os.path.isfile(TMSEG_file):
                    with open(TMSEG_file) as f:
                        split_record = f.read().strip().split("\n")
                        seq = split_record[-2]
                        if full_seq != seq:
                            logging.warning("\n{} seq from TMSEG_fastastyle.txt does not match original seq with this acc\nseq:{}\nfse:{}".format(acc, seq, full_seq))
                            list_acc_without_TMSEG.append(acc)
                            continue
                        topo = split_record[-1]


                elif os.path.isfile(TMSEG_fastastyle_file):
                    with open(TMSEG_fastastyle_file) as f:
                        split_record = f.read().strip().split("\n")
                        seq = split_record[1]
                        if full_seq != seq:
                            logging.warning("\n{} seq from TMSEG_fastastyle.txt does not match original seq with this acc\nseq:{}\nfse:{}".format(acc, seq, full_seq))
                            list_acc_without_TMSEG.append(acc)
                            continue
                        topo = split_record[2]

                else:
                    raise ValueError("how did we end up here? something must be wrong.")

                # overwrite any uniprot TM indices
                TM_indices = korbinian.prot_list.parse_TMSEG.get_TM_indices_from_TMSEG_topo_str(topo)
                df.at[acc, "TM_indices"] = TM_indices
                df.at[acc, "topol_source"] = "TMSEG"
                number_of_TMDs = len(TM_indices)
                df.at[acc, "number_of_TMDs"] = number_of_TMDs
                df.at[acc, "number_of_TMDs_excl_SP"] = number_of_TMDs
                list_of_TMDs = ["TM{:02d}".format(n) for n in range(1, number_of_TMDs + 1)]
                df.at[acc, "list_of_TMDs"] = list_of_TMDs
                df.at[acc, "list_of_TMDs_excl_SP"] = list_of_TMDs
            else:
                list_acc_without_TMSEG.append(acc)
                logging.warning("{} no TMSEG files found".format(acc))
                df.at[acc, "TM_indices"] = ()
                df.at[acc, "topol_source"] = "noTMSEG"
                df.at[acc, "number_of_TMDs"] = np.nan
                df.at[acc, "list_of_TMDs"] = np.nan
        sys.stdout.write("\n")

        TMSEG_missing_dir = os.path.join(s["data_dir"], "summaries", "{:02d}".format(s["list_number"]), "List{:02d}_predictions".format(s["list_number"]), "TMSEG", "missing")
        if not os.path.isdir(TMSEG_missing_dir):
            os.makedirs(TMSEG_missing_dir)
        for acc in list_acc_without_TMSEG:
            fasta_file = os.path.join(TMSEG_missing_dir, "{}.fasta".format(acc))
            with open(fasta_file, "w") as f:
                f.write(">{}\n{}\n".format(acc, df.at[acc, "full_seq"]))
                sys.stdout.write("-")

    # if SCAMPI is the source of TM topology (rare!)
    # UNCERTAIN IF THIS CODE STILL WORKS
    if s['TM_def'] == "SCAMPI":
        df = korbinian.prot_list.SCAMPI.read_scampi_data(pathdict, s, logging, df)

    if s['TM_def'] == "STR_PDBTM":
        # files are derived from elsewhere. Simply update the topology source.
        df["topol_source"] = "STR_PDBTM"
        df["list_of_TMDs_excl_SP"] = df["list_of_TMDs"]
        df["number_of_TMDs_excl_SP"] = df["number_of_TMDs"]
        df["uniprot_SiPe"] = False

    # accessions with topology info
    acc_with_topo = df.query("TM_indices != () & number_of_TMDs >= 1 & topol_source == '{}'".format(s['TM_def'])).index.tolist()

    ########################################################################################
    #                                                                                      #
    #                          slicing of TMD and nonTMD regions                           #
    #                                                                                      #
    ########################################################################################
    # create a list of [TM01, TM02, TM03, etc.
    long_list_of_TMDs = []
    for i in range(1, 50):
        long_list_of_TMDs.append("TM{:02d}".format(i))

    df = df.loc[acc_with_topo, :]
    n_prot_with_topo = df.shape[0]

    sys.stdout.write("\ntopology available for {}/{} proteins. "
                     "TM definitions derived from {}\n".format(n_prot_with_topo, n_initial_prot, str(s['TM_def'])))
    # save to a csv
    df.to_csv(pathdict["list_csv"], sep=",", quoting=csv.QUOTE_NONNUMERIC)
    logging.info('~~~~~~~~~~~~                 finished get_topology_for_prot_list                   ~~~~~~~~~~~~')


def slice_TMDs_in_prot_list(s, pathdict, logging):
    """ Gets membrane protein topology from alternative sources, and slices out TMD sequences from full sequence.

    Parameters
    ----------
    s : dict
        Dictionary of settings derived from settings excel file.
        Columns "Parameter" and "Value" are converted to key/value in the dictionary, respectively.
    pathdict : dict
        Dictionary of the key paths and files associated with that List number.
    logging : logging.Logger
        Logger for printing to console and logfile.

    DataFrames
    ----------
    df = from parsed_csv.
        index = accessions
        columns = full_seq, topo, seqlen, M_indices etc

    Saved Files and Figures
    -----------------------
    IMPORTANT: saves as ListXX.csv, which will be overwritten by prepare_protein_list

    """
    logging.info('~~~~~~~~~~~~                 starting slice_TMDs_in_prot_list                   ~~~~~~~~~~~~')

    df = pd.read_csv(pathdict["list_csv"], sep = ",", quoting = csv.QUOTE_NONNUMERIC, index_col = 0, low_memory=False)
    n_initial_prot = df.shape[0]
    # convert to python tuple from stringlist
    df.TM_indices = df.TM_indices.dropna().apply(lambda x : tuple(ast.literal_eval(x)))

    # create a list of [TM01, TM02, TM03, etc.  to TM50]
    long_list_of_TMDs = []
    for i in range(1, 50):
        long_list_of_TMDs.append("TM{:02d}".format(i))

    for n, acc in enumerate(df.index):
        # get nested tuple of TMDs
        nested_tup_TMs = df.at[acc, "TM_indices"]
        # slice long list of TMD names to get an appropriate list for that protein [TM01, TM02, TM03, etc.
        len_nested_tup_TMs = len(nested_tup_TMs)
        list_of_TMDs = long_list_of_TMDs[:len_nested_tup_TMs]
        # add that list to the dataframe (could also be added as a stringlist, but that's irritating somehow)
        # df.at[acc, 'list_of_TMDs'] = list_of_TMDs
        df.set_value(acc, "list_of_TMDs", list_of_TMDs)
        # set seq for slicing
        full_seq = df.at[acc, "full_seq"]
        # topo = dft.loc[acc, "Topology"]
        # iterate through all the TMDs of that protein, slicing out the sequences
        for i, TMD in enumerate(list_of_TMDs):
            TMD = list_of_TMDs[i]
            start, end = nested_tup_TMs[i]
            # with UniProt indexing, need to slice with -1, not like python index style
            df.at[acc, "%s_start" % TMD] = start
            df.at[acc, "%s_end" % TMD] = end
            # for python indexing of the TMD rather than uniprot, the start should be minus 1
            python_indexing_tuple = (start - 1, end)
            df.at[acc, "%s_seq" % TMD] = utils.slice_with_listlike(full_seq, python_indexing_tuple)
            df.at[acc, "%s_seqlen" % TMD] = len(df.at[acc, "%s_seq" % TMD])

        if len(nested_tup_TMs) > 0:
            df.at[acc, "nonTMD_seq"] = get_nonTM_seq_using_nested_TM_indices(full_seq, nested_tup_TMs)
        else:
            df.at[acc, "nonTMD_seq"] = np.nan
    # save to a csv
    df.to_csv(pathdict["list_csv"], sep=",", quoting=csv.QUOTE_NONNUMERIC)
    sys.stdout.write("Done!\n")
    logging.info('~~~~~~~~~~~~                 finished slice_TMDs_in_prot_list                   ~~~~~~~~~~~~')

def get_indices_TMD_plus_surr_for_summary_file(dfsumm, TMD, n_aa_before_tmd, n_aa_after_tmd):
    """Takes a summary dataframe (1 row for each protein) and slices out the TMD seqs.

    Returns the dataframe with extra columns, TM01_start, TM01_start_plus_surr_seq, etc

    Parameters
    ----------
    dfsumm : pd.DataFrame
        Input dataframe with sequences for slicing
    TMD : str
        String denoting transmembrane domain number (e.g. "TM01")
    n_aa_before_tmd
    n_aa_after_tmd

    Returns
    -------
    dfsumm : pd.DataFrame
        Returns the dataframe with extra columns, TM01_start, TM01_start_plus_surr_seq, etc
    """
    # instead of integers showing the start or end of the TMD, some people write strings into the
    # UniProt database, such as "<5" or "?"
    # to avoid the bugs that this introduces, it is necessary to convert all strings to np.nan (as floats),
    # using the convert objects function. The numbers can then be converted back from floats to integers.
    dfsumm['%s_start' % TMD] = pd.to_numeric(dfsumm['%s_start' % TMD]).dropna().astype('int64')
    dfsumm['%s_end' % TMD] = pd.to_numeric(dfsumm['%s_end' % TMD]).dropna().astype('int64')
    # determine the position of the start of the surrounding sequence
    dfsumm['%s_start_plus_surr' % TMD] = dfsumm['%s_start' % TMD] - n_aa_before_tmd
    # replace negative values with zero. (slicing method was replaced with lambda function to avoid CopyWithSetting warning)
    dfsumm['%s_start_plus_surr' % TMD] = dfsumm['%s_start_plus_surr' % TMD].apply(lambda x: x if x > 0 else 0)
    dfsumm['%s_end_plus_surr' % TMD] = dfsumm['%s_end' % TMD] + n_aa_after_tmd
    # create a boolean series, describing whether the end_surrounding_seq_in_query is longer than the protein seq
    series_indices_longer_than_prot_seq = dfsumm.apply(find_indices_longer_than_prot_seq, args=(TMD,), axis=1)
    # obtain the indices of proteins in the series
    indices_longer_than_prot_seq = series_indices_longer_than_prot_seq[series_indices_longer_than_prot_seq].index
    # use indices to select the main dataframe, and convert these end_surrounding_seq_in_query values to the seqlen value
    dfsumm.loc[indices_longer_than_prot_seq, '%s_end_plus_surr' % TMD] = dfsumm.loc[indices_longer_than_prot_seq, 'seqlen']

    return dfsumm


def find_indices_longer_than_prot_seq(df, TMD):
    """ Finds indices that are longer than the protein sequence.

    Small lambda-like function used to determine where the TMD plus surrounding exceeds the length of the protein sequence

    Parameters
    ----------
    df : pd.DataFrame
        Dataframe containing the columns for analysis.
    TMD : str
        String denoting transmembrane domain number (e.g. "TM01")
    -------

    """
    return df['%s_end_plus_surr'%TMD] > df['seqlen']

def calc_randTM_and_randnonTM(s, pathdict, logging, seq_len=1000, number_seq=1000, multiprocessing_mode=False):
    """ Calculates the random identity for the TM and nonTM regions of a particular dataset.

    Can either use the mathematical method, or extensive randomisation.

    Parameters
    ----------
    s : dict
        Settings dictionary extracted from excel settings file.
    pathdict : dict
        Dictionary of the key paths and files associated with that List number.
    logging : logging.Logger
        Logger for printing to console and logfile.
    seq_len: int
        length of randomly created sequences. To achieve a more plausible result using randomisation method,
        greater values (> 5000) are recommended. Defalut value: 1000
    number_seq: int
        number of aligned sequences. Larger values are recommended. Default value: 1000
    multiprocessing_mode : bool
        If TRUE, the random identity will be calculated using the randomisation method, taking advantage
        of python multiprocessing capability.

    Saved Files and Figures
    -----------------------
    rand_ident_TM_csv : csv
        Tab-separated CSV file with random identity and input data for TM, saved from output_ser
        index : "random_sequence_identity_output", "input_A", "input_C", etc
        values : 0.124, 0.05 etc
    rand_ident_nonTM_csv : csv
        Tab-separated CSV file with random identity and input data for nonTM, saved from output_ser
        index : "random_sequence_identity_output", "input_A", "input_C", etc
        values : 0.056, 0.05 etc
    """
    logging.info('~~~~~~~~~~~~                    starting calc_randTM_and_randnonTM                  ~~~~~~~~~~~~')
    ########################################################################################
    #                                                                                      #
    #                      calculate random TM and nonTM identity                          #
    #                                                                                      #
    ########################################################################################
    # ensure folder is exists (/summaries/ListXX_rand/ListXX...)
    utils.make_sure_path_exists(pathdict["rand_ident_TM_csv"], isfile=True)
    # calculate AA propensity and random ident for TM
    seq_list_csv_in = pathdict["list_csv"]
    aa_prop_csv_out_TM = pathdict["rand_ident_TM_csv"][:-11] + "aa_prop_TM.csv"
    col_name = "TMD_seq_joined"
    rand_ident_TM_csv = pathdict["rand_ident_TM_csv"]
    ident = 0.0
    # calculate aa propensity for all joined TM sequences
    korbinian.MSA_normalisation.calc_aa_propensity_from_csv_col(seq_list_csv_in, aa_prop_csv_out_TM, col_name)
    if multiprocessing_mode == True:
        d = aa_prop_csv_out_TM, rand_ident_TM_csv, seq_len, number_seq, ident, multiprocessing_mode
        d_list = [d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d]
        with Pool(processes = s["multiprocessing_cores"]) as pool:
            return_statement_list = pool.map(calc_random_aa_ident_multiprocessing, d_list)
        output_ser = return_statement_list[0][1]
        random_aa_identity_list = [return_statement[0] for return_statement in return_statement_list]
        random_aa_identity = np.array(random_aa_identity_list).mean()
        sys.stdout.write("\nrandom_aa_identity : {}\n, random_aa_identity_list : {}".format(random_aa_identity, random_aa_identity_list))
        output_ser["random_sequence_identity_output"] = random_aa_identity
        # save the series as csv file
        output_ser.to_csv(rand_ident_TM_csv, sep="\t")
    else:
        # calculate random aa identity based on aa propensity
        korbinian.MSA_normalisation.calc_random_aa_ident(aa_prop_csv_out_TM, rand_ident_TM_csv)
    # remove temp file with aa propensity (it is also saved in the rand_ident csv)
    #os.remove(aa_prop_csv_out_TM)

    # calculate AA propensity and random ident for nonTM
    aa_prop_csv_out_nonTM = pathdict["rand_ident_TM_csv"][:-11] + "aa_prop_nonTM.csv"
    rand_ident_nonTM_csv = pathdict["rand_ident_nonTM_csv"]
    col_name = "nonTMD_seq"
    # calculate aa propensity for nonTM sequence
    korbinian.MSA_normalisation.calc_aa_propensity_from_csv_col(seq_list_csv_in, aa_prop_csv_out_nonTM, col_name)
    if multiprocessing_mode == True:
        d = aa_prop_csv_out_nonTM, rand_ident_nonTM_csv, seq_len, number_seq, ident, multiprocessing_mode
        d_list = [d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d,d]
        with Pool(processes = s["multiprocessing_cores"]) as pool:
            return_statement_list = pool.map(calc_random_aa_ident_multiprocessing, d_list)
        output_ser = return_statement_list[0][1]
        random_aa_identity_list = [return_statement[0] for return_statement in return_statement_list]
        random_aa_identity = np.array(random_aa_identity_list).mean()
        sys.stdout.write("random_aa_identity, {} \nrandom_aa_identity_list {}".format(random_aa_identity, random_aa_identity_list))
        output_ser["random_sequence_identity_output"] = random_aa_identity
        # save the series as csv file
        output_ser.to_csv(rand_ident_nonTM_csv, sep="\t")
    else:
        # calculate random aa identity based on aa propensity
        korbinian.MSA_normalisation.calc_random_aa_ident(aa_prop_csv_out_nonTM, rand_ident_nonTM_csv)
    # remove temp file with aa propensity (it is also saved in the rand_ident csv)
    #os.remove(aa_prop_csv_out_nonTM)
    logging.info('~~~~~~~~~~~~                    finished calc_randTM_and_randnonTM                  ~~~~~~~~~~~~')

def calc_random_aa_ident_multiprocessing(d):
    """ Runs calc_random_aa_ident using for multiprocessing.

    Takes only a single input value, d, which is required for multiprocessing with Pool.

    Parameters
    ----------
    d : list
        List of input values for calc_random_aa_ident

    Returns
    -------
    Only returns values when multiprocessing_mode = True (no file is saved)
    random_aa_identity : float
        Random amino acid identity, based on the limited AA propensity in the dataset.
    output_ser = pd.Series
        Output pandas series consisting of the random_aa_identity, and the input aa_prop_ser derived from aa_prop_csv_in.
    """
    aa_prop_csv_out, rand_ident_TM_csv, seq_len, number_seq, ident, multiprocessing_mode = d
    random_aa_identity, output_ser = korbinian.MSA_normalisation.calc_random_aa_ident_via_randomisation(aa_prop_csv_out, rand_ident_TM_csv, seq_len=seq_len, number_seq=number_seq, ident=0.0, multiprocessing_mode=multiprocessing_mode)
    return random_aa_identity, output_ser


def get_nonTM_seq_using_nested_TM_indices(s, t):
    """Get the nonTM (EM) region of protein sequence based on nested tuple of TM indices.

    Parameters
    ----------
    s : str
        Full sequence to be sliced
        e.g. 'MKLLRRAWRRRAALGLGTLALCGAALLYLARCAAEPGDPRAMSGRSPPPPAPARAAAFLAVLVA'
    t : tuple
        Nested tuple of TM regions
        e.g. ((11, 30),)
        IMPORTANT: Indexing is using UniProt style (1:3), not python style (0:3)

    Returns
    -------
    nonTMD_seq : str
        Non-TMD sequence (all residues surrounding TMD)
    """
    TM01_start = t[0][0]
    last_TMD_end = t[-1][1]
    nonTMD_first = s[0: int(TM01_start) - 1]
    nonTMD_seq = nonTMD_first
    # only for multipass proteins, generate sequences between TMDs
    if len(t) > 1:
        for TM_Nr in range(len(t) - 1):
            # the TMD is the equivalent item in the list
            TMD = "TM{:02d}".format(TM_Nr + 1)
            # the next TMD, which contains the end index, is the next item in the list
            next_TMD = "TM{:02d}".format(TM_Nr + 2)

            TM_start = t[TM_Nr][0]
            TM_end = t[TM_Nr][1]
            next_TM_start = t[TM_Nr + 1][0]
            next_TM_end = t[TM_Nr + 1][1]

            between_TM_and_TMplus1 = s[TM_end: next_TM_start - 1]

            nonTMD_seq += between_TM_and_TMplus1

    # sequence from last TMD to C-term.
    nonTMD_last = s[last_TMD_end:]
    nonTMD_seq += nonTMD_last
    return nonTMD_seq


def add_file_paths_to_list_csv(df, s):
    """Adds file paths related to individual proteins to the list csv.

    E.g., The location of the BLAST-like homologue xml file.

    Parameters
    ----------
    df : pd.DataFrame
        Dataframe derived from list_csv, consisting of a list of proteins
        index = accessions
        columns = 'TM_indices', 'comments_subcellular_location_uniprot', 'full_seq', etc
    s : dict
        Settings dictionary extracted from excel settings file.

    Returns
    -------
    df : pd.DataFrame
        Modified version of original dataframe.

    Usage
    -----
    df = add_file_paths_to_list_csv(df, s)
    """

    ########################################################################################
    #                                                                                      #
    #                                     setup file paths                                 #
    #                                                                                      #
    ########################################################################################
    df['first_two_letters_of_uniprot_acc'] = df.index.str[0:2]

    simap_dir = s["simap_dir"]
    utils.make_sure_path_exists(simap_dir)
    homol_dir = os.path.join(s["data_dir"], "homol") # D:\Databases\homol
    parsed_dir = os.path.join(homol_dir, "parsed") # D:\Databases\homol\parsed
    sliced_dir = os.path.join(homol_dir, "sliced") # D:\Databases\homol\sliced
    utils.make_sure_path_exists(parsed_dir)
    utils.make_sure_path_exists(sliced_dir)

    ########################################################################################
    #                                     SIMAP directory                                 #
    ########################################################################################
    df['simap_filename_base'] = simap_dir + '/' + df.first_two_letters_of_uniprot_acc + '/' + df.index.astype(str)
    # normalise path to suit operating system
    df['simap_filename_base'] = df['simap_filename_base'].apply(lambda x: os.path.normpath(x))

    # create filenames for simap output
    df['SIMAP_tar'] = df.simap_filename_base + '_SIMAP.tar.gz'
    df['SIMAP_feature_table_XML_path'] = df.simap_filename_base + '_feature_table.xml'
    df['SIMAP_homol_XML_path'] = df.simap_filename_base + '_homologues.xml'
    your_name = unicodedata.normalize('NFKD', s["your_name"][:20]).encode('ascii', 'ignore').decode("utf-8")
    df['SIMAP_download_date_file_path'] = df.simap_filename_base + '--{}--{}.txt'.format(time.strftime("%Y%m%d"), your_name)


    ########################################################################################
    #                                     parsed directory                                 #
    ########################################################################################
    df['homol_parsed'] = parsed_dir + '/' + df.first_two_letters_of_uniprot_acc + '/' + df.index.astype(str) # D:\Databases\homol\parsed\A0\A0A0A1Y489
    df['homol_parsed'] = df['homol_parsed'].apply(lambda x: os.path.normpath(x))
    # ORIG: create filename for csv parsed from homologue XML file, stored temp as file, then zipped and pickled
    df['SIMAP_orig_csv'] = df['homol_parsed'] + '_orig.csv'
    # ORIG: create filename for csv with alignment_pretty, for visual analysis of homologues
    df['SIMAP_align_pretty_csv'] = df['homol_parsed'] + '_align_pretty.csv'
    # ORIG: create filename pickled dataframe
    df['homol_df_orig_pickle'] = df['homol_parsed'] + '_df_orig.pickle'
    # ORIG: create filename for zip that holds the XML parsed to a table (i.e. pandas dataframe, pickled)
    df['homol_df_orig_zip'] = df['homol_parsed'] + '_homol_orig_table.zip'
    # filename for header csv, which holds info at top of SIMAP XML
    df['simap_header_info_csv'] = df['homol_parsed'] + '_simap_header_info.csv'

    ########################################################################################
    #                                                                                      #
    #       create a flexible basename (flexibase) that can distinguish between sliced     #
    #                     files with and without signal peptides                           #
    #                                                                                      #
    ########################################################################################
    TM_definitions_dir = s["TM_def"]
    region_dir = s["regions"]
    # create directory, e.g.D:\Databases\homol\sliced\TMSEG\SiPe_TM\P9\
    df['homol_flexidir'] = sliced_dir + '/' + TM_definitions_dir + '/' + region_dir + '/' + df.first_two_letters_of_uniprot_acc
    # create basename, e.g.D:\Databases\homol\sliced\TMSEG\SiPe_TM\P9\P95834
    df['homol_flexibase'] = df['homol_flexidir'] + '/' + df.index.astype(str) + "_" + TM_definitions_dir + "_" + region_dir
    # normalise paths
    df['homol_flexidir'] = df['homol_flexidir'].apply(lambda x: os.path.normpath(x))
    df['homol_flexibase'] = df['homol_flexibase'].apply(lambda x: os.path.normpath(x))

    # SLICED (FASTA AND AAIMON): create filename for zip that holds the tables (pickled dataframes) for each TMD
    df['fa_cr_sliced_TMDs_zip'] = df['homol_flexibase'] + '_fa_cr_sliced_TMDs.zip'

    # FASTA: create filename for zip that holds the tables (pickled dataframes) for each TMD
    df['fa_TMDs_zip'] = df['homol_flexibase'] + '_fa_TMDs.zip'
    # FASTA: create filename for zip that holds the .fas files
    df['fa_fasta_zip'] = df['homol_flexibase'] + '_fa_fasta.zip'

    # CONS_RATIOS: create filename for zip holding the tables (pickled dataframes with seqs) for each TMD
    df['homol_cr_TMDs_zip'] = df['homol_flexibase'] + '_cr_tables.zip'
    # CONS_RATIOS: create filename for zip holding the ratios
    df['homol_cr_ratios_zip'] = df['homol_flexibase'] + '_cr_ratios.zip'
    # CONS_RATIOS: create filename for zip holding the figures
    df['homol_cr_figs_zip'] = df['homol_flexibase'] + '_cr_figs.zip'

    # GAPS: create filename for zip holding the tables (pickled dataframes with seqs) for each TMD
    df['homol_gap_tables_zip'] = df['homol_flexibase'] + '_gap_tables.zip'
    # GAPS: create filename for zip holding the gap numbers (pickled dataframes) for each TMD
    df['homol_gap_nums_zip'] = df['homol_flexibase'] + '_gap_nums.zip'
    # GAPS: create filename for zip folding the figures (pickled dataframes) for each TMD
    df['homol_gap_figs_zip'] = df['homol_flexibase'] + '_gap_figs.zip'

    # FASTAGAP: create filename for zip that holds the .fas files
    df['fastagap_zip'] = df['homol_flexibase'] + '_fastagap.zip'
    df['fastagap_pos_arrays_zip'] = df['homol_flexibase'] + '_fastagap_pos_arrays.zip'
    df['fastagap_base'] = df['homol_flexibase'] + '_homol_seq_plus_surr_'

    # FASTA: create basal name for fasta file with the TMD seqs (eg A0A1F4_EYS_DROME_homol_seq_ + TM01.fas)
    df['fasta_file_BASENAME'] = df.index.astype(str) + '_homol_seq_'
    df['fasta_file_BASENAMEPATH'] = df.homol_flexibase + '_homol_seq_'
    # FASTA: name the fasta file with surrounding seq (eg A0A1F4_EYS_DROME_homol_seq_plus_surr_ + TM01.fas)
    df['fasta_file_plus_surr_BASENAME'] = df.index.astype(str) + '_homol_seq_plus_surr_'
    df['fasta_file_plus_surr_BASENAMEPATH'] = df.homol_flexibase + '_homol_seq_plus_surr_'

    # create a basename for the output histograms
    df['AAIMON_hist_path_prefix'] = df.homol_flexibase + '_AAIMON_hist'
    df['norm_scatter_path_prefix'] = df.homol_flexibase + '_norm_scatter'

    ########################################################################################
    #                                                                                      #
    #           Old stuff, with one single outputfile for different functions              #
    #                                                                                      #
    ########################################################################################
    # df['SIMAP_csv_from_XML'] = df.index.astype(str) + '.csv'
    # df['SIMAP_orig_csv'] = df.simap_filename_base + '.csv'
    # df['SIMAP_csv_from_XML_tarfile'] = df.simap_filename_base + '.csv.tar.gz'
    # df['SIMAP_csv_analysed'] = df.index.astype(str) + '_analysed.csv'
    # df['SIMAP_csv_analysed_path'] = df.simap_filename_base + '_analysed.csv'
    # df['output_tarfile'] = df.index.astype(str) + '_outputfiles.tar.gz'
    # df['output_tarfile_path'] = df.simap_filename_base + '_outputfiles.tar.gz'
    df['csv_SIMAP_homologues_kept_for_statistical_analysis'] = df.simap_filename_base + '_homologues_for_stat.csv'
    # name the second fasta file (eg. A0T0U2_PSBE_THAPS_simap_TMD_seq_homol_&_surrounding.fas)
    df['fast_homol_kept_stat_analysis'] = df.simap_filename_base + '_simap_TMD_seq_kept_stat_analysis.fas'
    df['csv_file_av_cons_ratios_hits'] = df.simap_filename_base + '_cons_ratios.csv'
    df['csv_file_av_cons_ratios_hits_BASENAME'] = df.index.astype(str) + '_cons_ratios_'
    df['csv_file_av_cons_ratios_hits_BASENAMEPATH'] = df.simap_filename_base + '_cons_ratios_'

    return df

def add_user_subseqs(df, pathdict):
    """Rarely used method to add user-derived subsequences for analysis.

    Can be used to compare evolution of a TMD to a particular protein section (e.g. signal peptide, GPI anchor)

    User sequences are labelled "SE01_seq" and so on.

    They are added to the list of TMDs, which then looks like ["TM01", "SE01"]

    Parameters
    ----------
    df : pd.DataFrame
        Dataframe derived from list_csv, consisting of a list of proteins
        index = accessions
        columns = 'TM_indices', 'comments_subcellular_location_uniprot', 'full_seq', etc
    pathdict : dict
        Dictionary of the key paths and files associated with that List number.

    Returns
    -------
    df : pd.DataFrame
        Updated dataframe.

    Usage
    -----
    if s["add_user_subseqs"] == True:
        df = add_user_subseqs(df, pathdict)
    """

    ########################################################################################
    #                                                                                      #
    #      Add user-selected sequences from csv or excel("SE01", "SE02" etc)               #
    #                                                                                      #
    ########################################################################################

    """determine if the user sequences are submitted as excel or CSV
    Either way, they should have the following format:

                SE01_seq 	                            SE02_seq
    uniprot_acc
    A5HEI4 	    LLLSLAFMEALTIYGLVVALVLLFA 	            NaN
    A5U127 	    VDLAVAVVIGTAFTALVTKFTDSIITPLI 	        NaN
    A8EVM5 	    NaN 	                            YAWVFFIPFIFV
    B0R2U4 	    NaN 	                            LGSLFTVIAADIGMCVTGLA
    B0SR19 	    ISRNMYIMFFLGVVLWFVYGI 	                 NaN
    """
    if os.path.isfile(pathdict["list_user_subseqs_csv"]):
        user_subseqs_file = pathdict["list_user_subseqs_csv"]
    elif os.path.isfile(pathdict["list_user_subseqs_xlsx"]):
        user_subseqs_file = pathdict["list_user_subseqs_xlsx"]
    else:
        raise FileNotFoundError("add_user_TMDs is marked as True, but we can't find a file containing the "
                                "user sequences to add to the list of TMDs (e.g.{})".format(pathdict["list_user_subseqs_xlsx"]))

    # load as a pandas dataframe
    if user_subseqs_file[-4:] == ".csv":
        df_SE = pd.read_csv(user_subseqs_file, index_col=0)
    elif user_subseqs_file[-4:] == "xlsx":
        df_SE = pd.read_excel(user_subseqs_file, index_col=0)

    # df_SE with the selected sequences should be cropped to only include the proteins in the original list (df)
    # first find the common indices between the two dataframes
    common = set(df.index).intersection(set(df_SE.index))
    # reindex the dataframe to drop any proteins not in the original list
    df_SE = df_SE.reindex(index=common)

    # create a series of lists of SEs ([SE01, SE02] etc) for each protein
    nested_list_of_SEs = []
    for row in df_SE.index:
        # get list of SE_seqs (SE01seq, SE02seq, etc)
        list_of_SEs_with_seq = list(df_SE.loc[row, :].dropna().index)
        # drop the _seq (SE01, SE02, etc)
        list_of_SEs = [s[:-4] for s in list_of_SEs_with_seq]
        nested_list_of_SEs.append(list_of_SEs)
    # convert nested list to pandas series
    list_of_SEs_ser = pd.Series(nested_list_of_SEs, index=df_SE.index)

    # SHOULD NO LONGER BE A STRING; SHOULD ALL BE PYTHON LIST FORMAT
    # if type(df["list_of_TMDs"][0]) == str:
    #     df["list_of_TMDs"] = df["list_of_TMDs"].str.strip("'[]'").str.split("', '")

    # append the list of SEs (selected sequences) to the list of TMDs e.g. [TM01, TM02, SE01]
    df["list_of_TMDs"] = df["list_of_TMDs"] + list_of_SEs_ser
    # add the sequences (SE01_seq etc) to the main dataframe
    df = pd.concat([df, df_SE], axis=1)

    return df

def you_want_to_see_some_deprecated_stuff():
    # code may be useful for special cases, like signal peptides
    pass

    # DEPRECATED. TMSEG results now found from individual files.
    # # drop all TMSEG nonTM proteins
    # drop_TMSEG_nonTM = False
    # n_prot_AFTER_dropping_TMSEG_nonTM_proteins = 'TMSEG output file not found!'
    # if os.path.isfile(pathdict['TMSEG_nonTM']):
    #     modification_date_TMSEG_nonTM = time.ctime(os.path.getmtime(pathdict['TMSEG_nonTM']))
    #     TMSEG_nonTM_list = []
    #     with open(pathdict['TMSEG_nonTM']) as source:
    #         for line in source:
    #             line = line.strip()
    #             TMSEG_nonTM_list.append(line)
    #
    #     keep_after_TMSEG = []
    #     for acc in df.index:
    #         if not acc in TMSEG_nonTM_list:
    #             keep_after_TMSEG.append(acc)
    #
    #     df = df.loc[keep_after_TMSEG, :]
    #     n_prot_AFTER_dropping_TMSEG_nonTM_proteins = df.shape[0]
    # else:
    #     modification_date_TMSEG_nonTM = None
    #
    # df_PLS.loc["n_prot_AFTER_dropping_TMSEG_nonTM_proteins",:] = (n_prot_AFTER_dropping_TMSEG_nonTM_proteins, t)

    # DEPRECATED signal peptide method
    # # load PrediSi file containing all proteins that have signal peptides
    # PrediSi_inpath = os.path.join('%s_SCAMPI' % pathdict["base_filename_summaries"], 'SiPe_PrediSi.txt')
    # n_prot_AFTER_dropping_non_trusted_SiPe_PrediSi = 'PrediSi output file not found!'
    # if os.path.isfile(pathdict['TMSEG_nonTM']):
    #     modification_date_PrediSi = time.ctime(os.path.getmtime(pathdict['TMSEG_nonTM']))
    #     PrediSi_SiPe_list = korbinian.cons_ratio.SCAMPI.get_PrediSi_SiPe_acc(PrediSi_inpath)
    #     # create column containing bool if signal peptide is predicted by PrediSi
    #     for acc in df.index:
    #         if acc in PrediSi_SiPe_list:
    #             df.loc[acc, 'PrediSi_SiPe'] = True
    #         else:
    #             df.loc[acc, 'PrediSi_SiPe'] = False
    #     # drop all proteins where Uniprot did not get signal peptide
    #     drop = []
    #     for acc in df.index:
    #         if df.loc[acc, 'uniprot_SiPe'] == False and df.loc[acc, 'PrediSi_SiPe'] == True:
    #             drop.append(acc)
    #     df = df.drop(drop, axis=0)
    #     n_prot_AFTER_dropping_non_trusted_SiPe_PrediSi = df.shape[0]
    # else:
    #    modification_date_PrediSi = None



    # if len(list_of_TMDs) == 1:
    #     # there is only one TMD, add it as the TMD_seq_joined and as the lipo mean
    #     df.at[acc, 'TMD_seq_joined'] = df.at[acc, 'TM01_seq']
    #
    #     df.at[acc, 'lipo_mean_all_TMDs_mean'

    # elif len(list_of_TMDs) > 1:
    ########################################################################################
    #                    create joined string with all TM sequences                        #
    #                       measure lipophilicity all TM sequences                         #
    #                             get length of every TMD                                  #
    ########################################################################################
    # TMD_seq_joined = ''
    # lipo_list = []
    # TMD_seqlen_list = []


    # for TMD in list_of_TMDs_excl_SP:
    #     if "{}_seq".format(TMD) not in df.columns:
    #         list_acc_missing_TM_indices.append(acc)
    #         # skip this TMD. Add protein to list to be dropped.
    #         continue
    #     seq = df.loc[acc, '%s_seq' % TMD]
    #     # get length of TMD sequence
    #     seqlen = len(seq)
    #     df.loc[acc, '%s_seqlen' %TMD] = seqlen
    #     TMD_seqlen_list.append(seqlen)
    #     if type(seq) is float:
    #         list_acc_missing_TM_indices.append(acc)
    #         # skip this TMD. Add protein to list to be dropped.
    #         continue
    #     # add each TM sequence to the end of the growing string
    #     TMD_seq_joined += seq
    #
    #     # calculate lipophilicity, add new col for each TMD, add to list in order to calc mean for all TMDs
    #     lipo = utils.calc_lipophilicity(seq)
    #     df.loc[acc, '%s_lipo' % TMD] = lipo
    #     if TMD != 'SP01':
    #         lipo_list.append(lipo)

    # df.loc[acc, 'TMD_seq_joined'] = TMD_seq_joined
    # # calc the mean lipophilicity
    # # note this is the mean of each TMD separately, not the lipo of the joined sequence
    # df.loc[acc, 'lipo_mean_all_TMDs_mean'] = np.array(lipo_list).mean()
    # # calc the mean seqlen of all TMDs
    # df.loc[acc, 'len_TMD_mean'] = np.array(TMD_seqlen_list).mean()
    # get last TMD
    # last_TMD = list_of_TMDs_excl_SP[-1]
    # if 'SP01' in list_of_TMDs:
    #     last_TMD = list_of_TMDs[-2]
    # else:
    #     last_TMD = list_of_TMDs[-1]



    # DEPRECATED
    # # Exclude TMDs that do not have a minimum hydrophobicity
    # lipo_cutoff = s["max_lipo_list"]
    # list_acc_lipo_mean_above_cutoff = list(df['lipo_mean_all_TMDs_mean'].loc[df['lipo_mean_all_TMDs_mean'] > lipo_cutoff].index)
    # # convert current dataframe index to a set
    # index_set = set(df.index)
    # # convert list of acc to drop to a set
    # set_acc_lipo_mean_above_cutoff = set(list_acc_lipo_mean_above_cutoff)
    # # find acc to drop that are actually still in the index by looking for the overlap of both sets
    # acc_lipo_mean_above_cutoff_to_remove = index_set.intersection(set_acc_lipo_mean_above_cutoff)
    #
    # lipo_dropped = list(df.loc[acc_lipo_mean_above_cutoff_to_remove, 'lipo_mean_all_TMDs_mean'])
    #
    # # drop rows
    # df = df.drop(acc_lipo_mean_above_cutoff_to_remove)
    # n_prot_AFTER_lipo_cutoff = df.shape[0]
    # df_PLS.loc["n_prot_AFTER_lipo_cutoff",:] = (n_prot_AFTER_lipo_cutoff, t)
    # logging.info('n_prot_AFTER_lipo_cutoff: {}'.format(n_prot_AFTER_lipo_cutoff))# line 311
    # if list_acc_lipo_mean_above_cutoff != []:
    #     logging.info('list_acc_lipo_mean_above_cutoff: {}\n'.format(list_acc_lipo_mean_above_cutoff))
    #     logging.info('lipo_dropped: {}\n'.format([ '%.3f' % elem for elem in lipo_dropped ]))


    # DEPRECATED. list_of_TMDs_excl_SP should now come directly from UniProt or TMSEG parsing.
    # # filter by number of TMDs (excluding signal peptides)
    # df["list_of_TMDs_excl_SP"] = ""
    # for acc in df.index:
    #     list_of_TMDs = df.loc[acc, "list_of_TMDs"].copy()
    #     # remove SP01 from the list, if it is present
    #     if "SP01" in list_of_TMDs:
    #         list_of_TMDs.remove("SP01")
    #     df.at[acc, "number_of_TMDs_excl_SP"] = len(list_of_TMDs)
    #     # also create a list of TMDs excluding the signal peptide (may already exist, depending on source of proteins)
    #     df.set_value(acc, "list_of_TMDs_excl_SP", list_of_TMDs)
    #
