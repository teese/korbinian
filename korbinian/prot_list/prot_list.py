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
    if "uniprot_entry_name" in df.columns:
        # join the accession and entry name to create a "protein name" for naming files
        df['protein_name'] = df.uniprot_acc + '_' + df.uniprot_entry_name
    else:
        # the list of proteins did not come from UniProt. Simply use the accession to name the files.
        df['protein_name'] = df.index

    # convert the list_of_TMDs to a python object, if it is a string
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
    n = 0
    logging.info('joining TMD sequences, dropping proteins with "X", calculating lipophilicity:')
    list_acc_X_in_seq = []
    list_acc_missing_TM_indices = []
    for acc in df.index:
        n += 1
        if n % 20 == 0:
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
        ########################################################################################
        TMD_seq_joined = ''
        lipo_list = []
        for TMD in list_of_TMDs:
            if "{}_seq".format(TMD) not in df.columns:
                list_acc_missing_TM_indices.append(acc)
                # skip this TMD. Add protein to list to be dropped.
                continue
            seq = df.loc[acc, '%s_seq' % TMD]
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
    logging.info('n_initial_prot: {}'.format(n_initial_prot))

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
    logging.info("calculating rand_TM and rand_nonTM...")
    calc_randTM_and_randnonTM(s, pathdict, seq_len=1000, number_seq=1000)

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

def calc_randTM_and_randnonTM(s, pathdict, seq_len, number_seq, multiprocessing_mode=False):
    ########################################################################################
    #                                                                                      #
    #                      calculate random TM and nonTM identity                          #
    #                                                                                      #
    ########################################################################################
    # ensure folder is exists (/summaries/ListXX_rand/ListXX...)
    utils.make_sure_path_exists(pathdict["rand_ident_TM_csv"], isfile=True)
    # calculate AA propensity and random ident for TM
    seq_list_csv_in = pathdict["list_csv"]
    aa_prop_csv_out_TM = pathdict["rand_ident_TM_csv"][:-4] + "aa_prop_temp_TM.csv"
    col_name = "TMD_seq_joined"
    rand_ident_TM_csv = pathdict["rand_ident_TM_csv"]
    ident = 0.7
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
            sys.stdout.write("\nrandom_aa_identity, random_aa_identity_list", random_aa_identity, random_aa_identity_list)
            output_ser["random_sequence_identity_output"] = random_aa_identity
            # save the series as csv file
            output_ser.to_csv(rand_ident_TM_csv, sep="\t")
    else:
        # calculate random aa identity based on aa propensity
        korbinian.MSA_normalisation.calc_random_aa_ident(aa_prop_csv_out_TM, rand_ident_TM_csv, seq_len=seq_len, number_seq=number_seq, ident=ident)
    # remove temp file with aa propensity (it is also saved in the rand_ident csv)
    os.remove(aa_prop_csv_out_TM)

    # calculate AA propensity and random ident for nonTM
    aa_prop_csv_out_nonTM = pathdict["rand_ident_nonTM_csv"][:-4] + "aa_prop_temp_nonTM.csv"
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
            sys.stout.write("random_aa_identity, {} random_aa_identity_list {}".format(random_aa_identity, random_aa_identity_list))
            output_ser["random_sequence_identity_output"] = random_aa_identity
            # save the series as csv file
            output_ser.to_csv(rand_ident_nonTM_csv, sep="\t")
    else:
        # calculate random aa identity based on aa propensity
        korbinian.MSA_normalisation.calc_random_aa_ident(aa_prop_csv_out_nonTM, rand_ident_nonTM_csv, seq_len=seq_len, number_seq=number_seq, ident=ident)
    # remove temp file with aa propensity (it is also saved in the rand_ident csv)
    os.remove(aa_prop_csv_out_nonTM)

def calc_random_aa_ident_multiprocessing(d):

    aa_prop_csv_out, rand_ident_TM_csv, seq_len, number_seq, ident, multiprocessing_mode = d
    random_aa_identity, output_ser = korbinian.MSA_normalisation.calc_random_aa_ident(aa_prop_csv_out, rand_ident_TM_csv, seq_len=seq_len, number_seq=number_seq, ident=0.7, multiprocessing_mode=multiprocessing_mode)
    return random_aa_identity, output_ser
