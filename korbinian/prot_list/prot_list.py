from time import strftime
from multiprocessing import Pool
import ast
import csv
import korbinian
import korbinian.utils as utils
import numpy as np
import os
import pandas as pd
import sys
import time
import unicodedata


def prepare_protein_list(s, pathdict, logging):
    """ Sets up the file locations in the DataFrame containing the list of proteins for analysis.

    Parameters
    ----------
    s : dict
        Dictionary of settings derived from settings excel file.
        Columns "Parameter" and "Value" are converted to key/value in the dictionary, respectively.
    pathdict : dict
        Dictionary of the key paths and files associated with that List number.

    Saved Files and Figures
    -----------------------
    pathdict["list_summary_csv"] : csv file
        Input CSV file is overwritten at end of function, including the extra file locations.

    """
    logging.info('~~~~~~~~~~~~                     starting prepare_protein_list                      ~~~~~~~~~~~~')
    df = pd.read_csv(pathdict["list_parsed_csv"], sep = ",", quoting = csv.QUOTE_NONNUMERIC, index_col = 0, low_memory=False)
    n_initial_prot = df.shape[0]

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
        df = df.drop(SCAMPI_nonTM_list, axis=0)
        n_prot_AFTER_dropping_SCAMPI_nonTM_seqences = df.shape[0]
    else:
        modification_date = None

    # load SignalP file containing all proteins that have signal peptides
    SignalP_SiPe_path = os.path.join('%s_SCAMPI'%pathdict["base_filename_summaries"], 'gff.txt')
    n_prot_AFTER_dropping_non_trusted_SiPe = 'SignalP output file not found!'
    if os.path.isfile(SignalP_SiPe_path):
        modification_date_SignalP = time.ctime(os.path.getmtime(SignalP_SiPe_path))
        SignalP_acc_containing_SiPe = korbinian.cons_ratio.SCAMPI.get_SignalP_SiPe_acc(SignalP_SiPe_path)
        # create column containing bool if signal peptide is predicted by SignalP
        for acc in df.index:
            if acc in SignalP_acc_containing_SiPe:
                df.loc[acc, 'SignalP_SiPe'] = True
            else:
                df.loc[acc, 'SignalP_SiPe'] = False
        # drop all proteins where Uniprot did not get signal peptide
        drop = []
        for acc in df.index:
            # if UniProt DOESN'T think there is a signal peptide, but SignalP does
            # there is a chance that the "TM01" is really a signal peptide
            if df.loc[acc, 'uniprot_SiPe'] == False and df.loc[acc, 'SignalP_SiPe'] == True:
                # but it's only a danger if the TM01 is in the first 40 residues
                if df.loc[acc, 'TM01_start'] < 40:
                    drop.append(acc)
        df = df.drop(drop, axis=0)
        n_prot_AFTER_dropping_non_trusted_SiPe = df.shape[0]
    else:
        modification_date_SignalP = None

    # load PrediSi file containing all proteins that have signal peptides
    PrediSi_inpath = os.path.join('%s_SCAMPI' % pathdict["base_filename_summaries"], 'SiPe_PrediSi.txt')
    n_prot_AFTER_dropping_non_trusted_SiPe_PrediSi = 'PrediSi output file not found!'
    if os.path.isfile(PrediSi_inpath):
        modification_date_PrediSi = time.ctime(os.path.getmtime(PrediSi_inpath))
        PrediSi_SiPe_list = korbinian.cons_ratio.SCAMPI.get_PrediSi_SiPe_acc(PrediSi_inpath)
        # create column containing bool if signal peptide is predicted by PrediSi
        for acc in df.index:
            if acc in PrediSi_SiPe_list:
                df.loc[acc, 'PrediSi_SiPe'] = True
            else:
                df.loc[acc, 'PrediSi_SiPe'] = False
        # drop all proteins where Uniprot did not get signal peptide
        drop = []
        for acc in df.index:
            if df.loc[acc, 'uniprot_SiPe'] == False and df.loc[acc, 'PrediSi_SiPe'] == True:
                drop.append(acc)
        df = df.drop(drop, axis=0)
        n_prot_AFTER_dropping_non_trusted_SiPe_PrediSi = df.shape[0]
    else:
        modification_date_PrediSi = None


    if s['use_scampi_TM_regions']:
        df = korbinian.cons_ratio.SCAMPI.read_scampi_data(pathdict, s, logging, df)

    if "uniprot_entry_name" in df.columns:
        # join the accession and entry name to create a "protein name" for naming files
        df['protein_name'] = df.uniprot_acc + '_' + df.uniprot_entry_name
    else:
        # the list of proteins did not come from UniProt. Simply use the accession to name the files.
        df['protein_name'] = df.index

    # convert the list_of_TMDs to a python object, if it is a string
    if not s['use_scampi_TM_regions']:
        df['list_of_TMDs'] = df['list_of_TMDs'].dropna().apply(lambda x : ast.literal_eval(x))

    if s["add_user_subseqs"] == True:
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

    ########################################################################################
    #                                                                                      #
    #               drop any proteins that do not have a valid list_of_TMDs                #
    #                                                                                      #
    ########################################################################################
    df = df.loc[df['list_of_TMDs'].notnull()]
    n_prot_AFTER_dropping_without_list_TMDs = df.shape[0]

    ########################################################################################
    #                                                                                      #
    #                                     setup file paths                                 #
    #                                                                                      #
    ########################################################################################
    df['first_two_letters_of_uniprot_acc'] = df['protein_name'].str[0:2]

    simap_dir = s["simap_dir"]
    utils.make_sure_path_exists(simap_dir)
    homol_dir = os.path.join(s["data_dir"], "homol")
    utils.make_sure_path_exists(homol_dir)

    df['simap_filename_base'] = simap_dir + '/' + df.first_two_letters_of_uniprot_acc + '/' + df.protein_name
    # normalise path to suit operating system
    df['simap_filename_base'] = df['simap_filename_base'].apply(lambda x: os.path.normpath(x))

    # create the homologue basename, e.g. "D:\Databases\homol\P0\P0A334_KCSA_STRLI"
    # note that because UniProt protein names change, this should at some stage be changed from protein_name to acc alone
    df['homol_base'] = homol_dir + '/' + df.first_two_letters_of_uniprot_acc + '/' + df.protein_name
    # normalise path to suit operating system
    df['homol_base'] = df['homol_base'].apply(lambda x : os.path.normpath(x))

    # create filenames for simap output
    df['SIMAP_tar'] = df.simap_filename_base + '_SIMAP.tar.gz'
    df['SIMAP_feature_table_XML_path'] = df.simap_filename_base + '_feature_table.xml'
    df['SIMAP_homol_XML_path'] = df.simap_filename_base + '_homologues.xml'
    your_name = unicodedata.normalize('NFKD', s["your_name"][:20]).encode('ascii', 'ignore').decode("utf-8")
    df['SIMAP_download_date_file_path'] = df.simap_filename_base + '--{}--{}.txt'.format(strftime("%Y%m%d"), your_name)

    # ORIG: create filename for csv parsed from homologue XML file, stored temp as file, then zipped and pickled
    df['SIMAP_orig_csv'] = df['homol_base'] + '_orig.csv'
    # ORIG: create filename for csv with alignment_pretty, for visual analysis of homologues
    df['SIMAP_align_pretty_csv'] = df['homol_base'] + '_align_pretty.csv'
    # ORIG: create filename pickled dataframe
    df['homol_df_orig_pickle'] = df['homol_base'] + '_df_orig.pickle'
    # ORIG: create filename for zip that holds the XML parsed to a table (i.e. pandas dataframe, pickled)
    df['homol_df_orig_zip'] = df['homol_base'] + '_homol_orig_table.zip'

    # SLICED (FASTA AND AAIMON): create filename for zip that holds the tables (pickled dataframes) for each TMD
    df['fa_cr_sliced_TMDs_zip'] = df['homol_base'] + '_fa_cr_sliced_TMDs.zip'

    # FASTA: create filename for zip that holds the tables (pickled dataframes) for each TMD
    df['fa_TMDs_zip'] = df['homol_base'] + '_fa_TMDs.zip'
    # FASTA: create filename for zip that holds the .fas files
    df['fa_fasta_zip'] = df['homol_base'] + '_fa_fasta.zip'

    # CONS_RATIOS: create filename for zip holding the tables (pickled dataframes with seqs) for each TMD
    df['homol_cr_TMDs_zip'] = df['homol_base'] + '_cr_tables.zip'
    # CONS_RATIOS: create filename for zip holding the ratios
    df['homol_cr_ratios_zip'] = df['homol_base'] + '_cr_ratios.zip'
    # CONS_RATIOS: create filename for zip holding the figures
    df['homol_cr_figs_zip'] = df['homol_base'] + '_cr_figs.zip'

    # GAPS: create filename for zip holding the tables (pickled dataframes with seqs) for each TMD
    df['homol_gap_tables_zip'] = df['homol_base'] + '_gap_tables.zip'
    # GAPS: create filename for zip holding the gap numbers (pickled dataframes) for each TMD
    df['homol_gap_nums_zip'] = df['homol_base'] + '_gap_nums.zip'
    # GAPS: create filename for zip folding the figures (pickled dataframes) for each TMD
    df['homol_gap_figs_zip'] = df['homol_base'] + '_gap_figs.zip'

    # FASTAGAP: create filename for zip that holds the .fas files
    df['fastagap_zip'] = df['homol_base'] + '_fastagap.zip'
    df['fastagap_pos_arrays_zip'] = df['homol_base'] + 'fastagap_pos_arrays.zip'
    df['fastagap_base'] = df['homol_base'] + '_homol_seq_plus_surr_'

    # FASTA: create basal name for fasta file with the TMD seqs (eg A0A1F4_EYS_DROME_homol_seq_ + TM01.fas)
    df['fasta_file_BASENAME'] = df.protein_name + '_homol_seq_'
    df['fasta_file_BASENAMEPATH'] = df.homol_base + '_homol_seq_'
    # FASTA: name the fasta file with surrounding seq (eg A0A1F4_EYS_DROME_homol_seq_plus_surr_ + TM01.fas)
    df['fasta_file_plus_surr_BASENAME'] = df.protein_name + '_homol_seq_plus_surr_'
    df['fasta_file_plus_surr_BASENAMEPATH'] = df.homol_base + '_homol_seq_plus_surr_'

    # create a basename for the output histograms
    df['AAIMON_hist_path_prefix'] = df.homol_base + '_AAIMON_hist'
    df['norm_scatter_path_prefix'] = df.homol_base + '_norm_scatter'

    ########################################################################################
    #                                                                                      #
    #           Old stuff, with one single outputfile for different functions              #
    #                                                                                      #
    ########################################################################################
    # df['SIMAP_csv_from_XML'] = df.protein_name + '.csv'
    # df['SIMAP_orig_csv'] = df.simap_filename_base + '.csv'
    # df['SIMAP_csv_from_XML_tarfile'] = df.simap_filename_base + '.csv.tar.gz'
    # df['SIMAP_csv_analysed'] = df.protein_name + '_analysed.csv'
    # df['SIMAP_csv_analysed_path'] = df.simap_filename_base + '_analysed.csv'
    # df['output_tarfile'] = df.protein_name + '_outputfiles.tar.gz'
    # df['output_tarfile_path'] = df.simap_filename_base + '_outputfiles.tar.gz'
    df['csv_SIMAP_homologues_kept_for_statistical_analysis'] = df.simap_filename_base + '_homologues_for_stat.csv'
    # name the second fasta file (eg. A0T0U2_PSBE_THAPS_simap_TMD_seq_homol_&_surrounding.fas)
    df['fast_homol_kept_stat_analysis'] = df.simap_filename_base + '_simap_TMD_seq_kept_stat_analysis.fas'
    df['csv_file_av_cons_ratios_hits'] = df.simap_filename_base + '_cons_ratios.csv'
    df['csv_file_av_cons_ratios_hits_BASENAME'] = df.protein_name + '_cons_ratios_'
    df['csv_file_av_cons_ratios_hits_BASENAMEPATH'] = df.simap_filename_base + '_cons_ratios_'


    ########################################################################################
    #                                                                                      #
    #                drop based on the accepted number of TMDs for that dataset            #
    #                                                                                      #
    ########################################################################################

    min_TMDs = s["min_TMDs"]
    max_TMDs = s["max_TMDs"]
    analyse_signal_peptides = s["SiPe"]
    if analyse_signal_peptides == True:
        max_TMDs += 1
    df.dropna(subset=["number_of_TMDs"], inplace=True)
    df = df.loc[df["number_of_TMDs"].apply(lambda x: min_TMDs <= x <= max_TMDs)]
    n_prot_AFTER_n_TMDs_cutoff = df.shape[0]
    if n_prot_AFTER_n_TMDs_cutoff == 0:
        raise ValueError("The list {} has no valid proteins after n_TMDs_cutoff.".format(s["list_number"]))

    ########################################################################################
    #                                                                                      #
    #     slice out all the TMD_seq_plus_surr, based on settings (e.g. 10aa each side)     #
    #                                                                                      #
    ########################################################################################

    max_num_TMDs = df["number_of_TMDs"].max()
    n_aa_before_tmd = s["n_aa_before_tmd"]
    n_aa_after_tmd = s["n_aa_after_tmd"]

    #if 'TM01_seq_plus_surr' not in df.columns:
    # calculate TM plus surr for ALL sequences, overwriting if necessary, in case this is changed later
    # currently the loop is run for each TMD, based on the sequence with the most TMDs
    for i in range(1, int(max_num_TMDs) + 1):
        TMD = 'TM%02d' % i
        # get the indices for TMD plus surrounding sequence
        df = korbinian.prot_list.prot_list.get_indices_TMD_plus_surr_for_summary_file(df, TMD, n_aa_before_tmd, n_aa_after_tmd)
        # slice out the TMD_seq_plus_surr for each TMD
        df['%s_seq_plus_surr' % TMD] = df[df['%s_start' % TMD].notnull()].apply(utils.slice_uniprot_TMD_plus_surr_seq, args=(TMD,), axis=1)


    ########################################################################################
    #                                                                                      #
    #       X in sequence, joined TMD and nonTMD seq, lipophilicity calc                   #
    #                                                                                      #
    ########################################################################################
    # indicate that the prepare_protein_list function has been run
    df['prepare_protein_list'] = True

    # add a column that holds all joined TMD sequences, drop proteins with 'X' in full sequence
    logging.info('joining TMD sequences, dropping proteins with "X", calculating lipophilicity:')
    list_acc_X_in_seq = []
    list_acc_missing_TM_indices = []
    for n, acc in enumerate(df.index):
        if n % 20 == 0 and n != 0:
            sys.stdout.write('.'), sys.stdout.flush()
            if n % 600 == 0:
                sys.stdout.write('\n'), sys.stdout.flush()
        ########################################################################################
        #                               X in sequence                                          #
        ########################################################################################
        if 'X' in df.loc[acc, 'full_seq']:
            # remove protein from dataframe if sequence contains "X". Add acc to list.
            list_acc_X_in_seq.append(acc)
            # skip to next protein
            continue

        list_of_TMDs = df.loc[acc, 'list_of_TMDs']
        # if len(list_of_TMDs) == 1:
        #     # there is only one TMD, add it as the TMD_seq_joined and as the lipo mean
        #     df.loc[acc, 'TMD_seq_joined'] = df.loc[acc, 'TM01_seq']
        #
        #     df.loc[acc, 'lipo_mean_all_TMDs'

        # elif len(list_of_TMDs) > 1:
        ########################################################################################
        #                    create joined string with all TM sequences                        #
        #                       measure lipophilicity all TM sequences                         #
        #                             get length of every TMD                                  #
        ########################################################################################
        TMD_seq_joined = ''
        lipo_list = []
        TMD_seqlen_list = []
        for TMD in list_of_TMDs:
            if "{}_seq".format(TMD) not in df.columns:
                list_acc_missing_TM_indices.append(acc)
                # skip this TMD. Add protein to list to be dropped.
                continue
            seq = df.loc[acc, '%s_seq' % TMD]
            # get length of TMD sequence
            seqlen = len(seq)
            df.loc[acc, '%s_seqlen' %TMD] = seqlen
            TMD_seqlen_list.append(seqlen)
            if type(seq) is float:
                list_acc_missing_TM_indices.append(acc)
                # skip this TMD. Add protein to list to be dropped.
                continue
            # add each TM sequence to the end of the growing string
            TMD_seq_joined += seq

            # calculate lipophilicity, add new col for each TMD, add to list in order to calc mean for all TMDs
            lipo = utils.calc_lipophilicity(seq)
            df.loc[acc, '%s_lipo' % TMD] = lipo
            lipo_list.append(lipo)

        df.loc[acc, 'TMD_seq_joined'] = TMD_seq_joined
        # calc the mean lipophilicity
        # note this is the mean of each TMD separately, not the lipo of the joined sequence
        df.loc[acc, 'lipo_mean_all_TMDs'] = np.array(lipo_list).mean()
        # calc the mean seqlen of all TMDs
        df.loc[acc, 'len_TMD_mean'] = np.array(TMD_seqlen_list).mean()
        # get last TMD
        last_TMD = list_of_TMDs[-1]
        df.loc[acc, 'last_TMD'] = last_TMD
        # get lipo of last TMD
        df.loc[acc, 'lipo_last_TMD'] = df.loc[acc, '%s_lipo' % last_TMD]
        seq_last_TMD = df.loc[acc, '%s_seq' % last_TMD]
        df.loc[acc, 'last_TMD_seq'] = seq_last_TMD
        df.loc[acc, 'last_TMD_seqlen'] = len(seq_last_TMD)

        if len(list_of_TMDs) >= 3:
        # calculate mean lipophilicity excluding first and last TM
            list_of_central_TMDs = list_of_TMDs[1:-1]
            lipo_list_central_TMDs = []
            for TMD in list_of_central_TMDs:
                seq = df.loc[acc, '%s_seq' % TMD]
                lipo = utils.calc_lipophilicity(seq)
                lipo_list_central_TMDs.append(lipo)
            df.loc[acc, 'lipo_mean_central_TMDs'] = np.array(lipo_list_central_TMDs).mean()
        elif len(list_of_TMDs) == 1:
            df.loc[acc, 'lipo_mean_central_TMDs'] = df.loc[acc, 'TM01_lipo']
        else:
            df.loc[acc, 'lipo_mean_central_TMDs'] = np.nan

        if len(list_of_TMDs) >= 2:
        # calculate mean lipophilicity excluding first and last TM
            list_of_TMDs_excl_TM01 = list_of_TMDs[1:]
            lipo_list_TMDs_excl_TM01 = []
            for TMD in list_of_TMDs_excl_TM01:
                seq = df.loc[acc, '%s_seq' % TMD]
                lipo = utils.calc_lipophilicity(seq)
                lipo_list_TMDs_excl_TM01.append(lipo)
            df.loc[acc, 'lipo_mean_excl_TM01'] = np.array(lipo_list_TMDs_excl_TM01).mean()
        elif len(list_of_TMDs) == 1:
            df.loc[acc, 'lipo_mean_excl_TM01'] = df.loc[acc, 'TM01_lipo']
        else:
            df.loc[acc, 'lipo_mean_excl_TM01'] = np.nan

    df = df.drop(list_acc_X_in_seq)
    n_prot_AFTER_dropping_with_X_in_seq = df.shape[0]

    # the ones without TMD_seq are actually proteins where the indices were a string, for example
    # TM01_start was "?"
    df = df.drop(list_acc_missing_TM_indices)
    n_prot_AFTER_dropping_missing_TM_indices = df.shape[0]

    lipo_cutoff = s["max_lipo_list"]

    list_acc_lipo_mean_above_cutoff = list(df['lipo_mean_all_TMDs'].loc[df['lipo_mean_all_TMDs'] > lipo_cutoff].index)

    # convert current dataframe index to a set
    index_set = set(df.index)
    # convert list of acc to drop to a set
    set_acc_lipo_mean_above_cutoff = set(list_acc_lipo_mean_above_cutoff)
    # find acc to drop that are actually still in the index by looking for the overlap of both sets
    acc_lipo_mean_above_cutoff_to_remove = index_set.intersection(set_acc_lipo_mean_above_cutoff)

    lipo_dropped = list(df.loc[acc_lipo_mean_above_cutoff_to_remove, 'lipo_mean_all_TMDs'])

    # drop rows
    df = df.drop(acc_lipo_mean_above_cutoff_to_remove)
    n_prot_AFTER_lipo_cutoff = df.shape[0]

    ########################################################################################
    #                                                                                      #
    #           check for subcellular location (PM, ER, Golgi) in keywords                 #
    #                                                                                      #
    ########################################################################################
    if 'uniprot_KW' in df.columns:
        # apply ast.literal_eval to every item in df['uniprot_KW']
        if isinstance(df['uniprot_KW'][0], str):
            df['uniprot_KW'] = df['uniprot_KW'].apply(lambda x: ast.literal_eval(x))
        # check for specific keywords
        df['Cell_membrane'] = df['uniprot_KW'].apply(utils.KW_list_contains_any_desired_KW, args=(['Cell membrane'],))
        df['Endoplasmic_reticulum'] = df['uniprot_KW'].apply(utils.KW_list_contains_any_desired_KW, args=(['Endoplasmic reticulum'],))
        df['Golgi_apparatus'] = df['uniprot_KW'].apply(utils.KW_list_contains_any_desired_KW, args=(['Golgi apparatus'],))

    ########################################################################################
    #                                                                                      #
    #                add percentage length of TMD region in full protein                   #
    #                                                                                      #
    ########################################################################################
    df['TMD_seq_joined_len'] = df['TMD_seq_joined'].str.len()
    df['perc_TMD'] = df['TMD_seq_joined_len'] / df['seqlen'] * 100

    ########################################################################################
    #                                                                                      #
    #                           Print record of dropped proteins                           #
    #                                                                                      #
    ########################################################################################
    logging.info('\nn_initial_prot: {}'.format(n_initial_prot))

    logging.info('n_prot_AFTER_dropping_SCAMPI_nonTM_seqences: {}'.format(n_prot_AFTER_dropping_SCAMPI_nonTM_seqences))
    if modification_date is not None:
        logging.info("modification_date of scampi file: {}".format(modification_date))

    logging.info('n_prot_AFTER_non_trusted_Signal_Peptides: {}'.format(n_prot_AFTER_dropping_non_trusted_SiPe))
    if modification_date_SignalP is not None:
        logging.info("modification_date of SignalP file: {}".format(modification_date_SignalP))

    logging.info('n_prot_AFTER_non_trusted_Signal_Peptides_PrediSi: {}'.format(n_prot_AFTER_dropping_non_trusted_SiPe_PrediSi))
    if modification_date_SignalP is not None:
        logging.info("modification_date of PrediSi file: {}".format(modification_date_PrediSi))

    logging.info('n_prot_AFTER_dropping_without_list_TMDs: {}'.format(n_prot_AFTER_dropping_without_list_TMDs)) # line 107

    logging.info('n_prot_AFTER_n_TMDs_cutoff: {}'.format(n_prot_AFTER_n_TMDs_cutoff)) # line 215

    logging.info('n_prot_AFTER_dropping_with_X_in_seq: {}'.format(n_prot_AFTER_dropping_with_X_in_seq)) # line 297
    if list_acc_X_in_seq != []:
        logging.info('list_acc_X_in_seq: {}'.format(list_acc_X_in_seq))

    logging.info('n_prot_AFTER_dropping_missing_TM_indices: {}'.format(n_prot_AFTER_dropping_missing_TM_indices)) # line 297
    if list_acc_missing_TM_indices != []:
        logging.info('list_acc_missing_TM_indices: {}'.format(list_acc_missing_TM_indices))

    logging.info('n_prot_AFTER_lipo_cutoff: {}'.format(n_prot_AFTER_lipo_cutoff))# line 311
    if list_acc_lipo_mean_above_cutoff != []:
        logging.info('list_acc_lipo_mean_above_cutoff: {}\n'.format(list_acc_lipo_mean_above_cutoff))
        logging.info('lipo_dropped: {}\n'.format([ '%.3f' % elem for elem in lipo_dropped ]))

    ########################################################################################
    #                                                                                      #
    #                                     Save to CSV                                      #
    #                                                                                      #
    ########################################################################################
    # drop any rows that have no data for any proteins (e.g. TM02_start for single-pass proteins)
    df.dropna(axis=1, inplace=True, how="all")
    # save to a csv
    df.to_csv(pathdict["list_csv"], sep=",", quoting=csv.QUOTE_NONNUMERIC)
    ########################################################################################
    #                                                                                      #
    #                      calculate random TM and nonTM identity                          #
    #                                                                                      #
    ########################################################################################
    # only calculate the rand_TM and rand_nonTM using the mathematical method
    calc_randTM_and_randnonTM(s, pathdict, logging)
    logging.info('~~~~~~~~~~~~                     finished prepare_protein_list                      ~~~~~~~~~~~~')

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
