import csv
import os
import pandas as pd
from time import strftime
import unicodedata
import korbinian.utils as utils

def setup_file_locations_in_df(s, pathdict):
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
    df = pd.read_csv(pathdict["list_summary_csv"], sep = ",", quoting = csv.QUOTE_NONNUMERIC, index_col = 0)
    if "uniprot_entry_name" in df.columns:
        # join the accession and entry name to create a "protein name" for naming files
        df['protein_name'] = df.uniprot_acc + '_' + df.uniprot_entry_name
    else:
        # the list of proteins did not come from UniProt. Simply use the accession to name the files.
        df['protein_name'] = df.index

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

        # append the list of SEs (selected sequences) to the list of TMDs e.g. [TM01, TM02, SE01]
        if type(df["list_of_TMDs"][0]) == str:
            df["list_of_TMDs"] = df["list_of_TMDs"].str.strip("'[]'").str.split("', '")
        df["list_of_TMDs"] = df["list_of_TMDs"] + list_of_SEs_ser
        # add the sequences (SE01_seq etc) to the main dataframe
        df = pd.concat([df, df_SE], axis=1)

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

    # FASTA: create basal name for fasta file with the TMD seqs (eg A0A1F4_EYS_DROME_homol_seq_ + TM01.fas)
    df['fasta_file_BASENAME'] = df.protein_name + '_homol_seq_'
    df['fasta_file_BASENAMEPATH'] = df.homol_base + '_homol_seq_'
    # FASTA: name the fasta file with surrounding seq (eg A0A1F4_EYS_DROME_homol_seq_plus_surr_ + TM01.fas)
    df['fasta_file_plus_surr_BASENAME'] = df.protein_name + '_homol_seq_plus_surr_'
    df['fasta_file_plus_surr_BASENAMEPATH'] = df.homol_base + '_homol_seq_plus_surr_'

    # create a basename for the output histograms
    df['AAIMON_hist_path_prefix'] = df.homol_base + '_AAIMON_hist'

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
    #                                     Save to CSV                                       #
    #                                                                                      #
    ########################################################################################
    # indicate that the setup_file_locations_in_df function has been run
    df['setup_file_locations_in_df'] = True
    # save to a csv
    df.to_csv(pathdict["list_summary_csv"], sep=",", quoting=csv.QUOTE_NONNUMERIC)


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