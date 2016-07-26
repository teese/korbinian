import csv
import os
import pandas as pd
from time import strftime
import korbinian.mtutils as utils

def setup_keyboard_interrupt_and_error_logging(set_, list_number):
    ''' -------Setup keyboard interrupt----------
    '''
    # import arcgisscripting
    import signal
    def ctrlc(sig, frame):
        raise KeyboardInterrupt("CTRL-C!")

    signal.signal(signal.SIGINT, ctrlc)
    '''+++++++++++++++LOGGING++++++++++++++++++'''
    date_string = strftime("%Y%m%d")

    # designate the output logfile
    logfile = os.path.join(set_["logfile_folder"],'List%s_%s_logfile.log' % (list_number, date_string))
    # a file to keep a record of the log settings used for that script
    utils.setup_error_logging(logfile)

def create_settingsdict(excel_file_with_settings):
    sheetnames = ["run_settings", "file_locations", "variables"]
    set_ = {}
    for sheetname in sheetnames:
        # open excel file as pandas dataframe
        dfset = pd.read_excel(excel_file_with_settings, sheetname=sheetname)
        # exclude row with notes, set parameter as index
        dfset = dfset[["parameter", "value"]].dropna()
        dfset.set_index("parameter", inplace=True)
        # convert true-like strings to True, and false-like strings to False
        dfset.value = dfset.value.apply(utils.convert_truelike_to_bool, convert_nontrue=False)
        dfset.value = dfset.value.apply(utils.convert_falselike_to_bool)
        # convert to dictionary
        sheet_as_dict = dfset.to_dict()["value"]
        # join dictionaries together
        set_.update(sheet_as_dict)

    list_paths_to_normalise = ["data_folder", "uniprot_folder", "data_harddrive", "eaSimap_path", "logfile_folder",
                               "summaries_folder", "simap_database_folder", "list_of_uniprot_accessions"]
    # normalise the paths for selected columns, so that they are appropriate for the operating system
    for path in list_paths_to_normalise:
        if path in set_:
            set_[path] = os.path.normpath(set_[path])

    return set_


def create_pathdict(base_filename_summaries):
    pathdict = {}
    pathdict["base_filename_summaries"] = base_filename_summaries
    pathdict["dfout01_uniprotcsv"] = '%s_uniprot.csv' % base_filename_summaries
    pathdict["dfout02_uniprotTcsv"] = '%s_uniprotT.csv' % base_filename_summaries
    pathdict["dfout03_uniprotxlsx"] = '%s_uniprot.xlsx' % base_filename_summaries
    pathdict["dfout04_uniprotcsv_incl_paths"] = '%s_uniprot_incl_paths.csv' % base_filename_summaries
    pathdict["dfout05_simapcsv"] = '%s_simap.csv' % base_filename_summaries
    pathdict["dfout06_simapxlsx"] = '%s_simap.xlsx' % base_filename_summaries
    pathdict["dfout07_simapnonred"] = '%s_simapnonred.csv' % base_filename_summaries
    pathdict["dfout08_simap_AAIMON"] = '%s_simap_AAIMON.csv' % base_filename_summaries
    pathdict["dfout09_simap_AAIMON_02"] = '%s_simap_AAIMON_02.csv' % base_filename_summaries
    pathdict["dfout10_uniprot_gaps"] = '%s_gap_densities.csv' % base_filename_summaries
    pathdict["dfout11_gap_test_out_png"] = '%s_gap_test_out.png' % base_filename_summaries
    pathdict["dfout12"] = 0
    pathdict["dfout13"] = 0
    pathdict["csv_file_with_histogram_data"] = '%s_histogram.csv' % base_filename_summaries
    pathdict["csv_file_with_histogram_data_normalised"] = '%s_histogram_normalised.csv' % base_filename_summaries
    pathdict["csv_file_with_histogram_data_normalised_redundant_removed"] = '%s_histogram_normalised_redundant_removed.csv' % base_filename_summaries
    pathdict["csv_file_with_md5_for_each_query_sequence"] = '%s_query_md5_checksums.csv' % base_filename_summaries
    pathdict["csv_av_cons_ratio_all_proteins"] = '%s_cons_ratios_nonred_av.csv' % base_filename_summaries
    pathdict["csv_std_cons_ratio_all_proteins"] = '%s_cons_ratios_nonred_std.csv' % base_filename_summaries
    pathdict["create_graph_of_gap_density_png"] = '%s_create_graph_of_gap_density.png' % base_filename_summaries
    return pathdict

def setup_file_locations_in_df(df, set_, pathdict):
    # set up a folder to hold the SIMAP BLAST-like output
    # note that at the moment, files cannot be compressed
    simap_data_folder = os.path.join(set_['data_folder'], 'simap')
    if "uniprot_entry_name" in df.columns:
        # join the accession and entry name to create a "protein name" for naming files
        df['A2_protein_name'] = df.A1_uniprot_accession + '_' + df.uniprot_entry_name
    else:
        # the list of proteins did not come from UniProt. Simply use the accession to name the files.
        df['A2_protein_name'] = df.A1_uniprot_accession
    df['first_two_letters_of_uniprot_acc'] = df['A1_uniprot_accession'].apply(lambda x: x[0:2])
    df['simap_filename_base_linuxpath'] = simap_data_folder + '/' + df.first_two_letters_of_uniprot_acc + '/' + df.A2_protein_name
    df['simap_filename_base'] = df['simap_filename_base_linuxpath'].apply(lambda x: os.path.normpath(x))
    df.drop('simap_filename_base_linuxpath', axis=1, inplace=True)

    # create filenames for simap output
    df['SIMAP_tarfile'] = df.simap_filename_base + '_SIMAP.tar.gz'
    df['SIMAP_feature_table_XML_file'] = df.A2_protein_name + '_feature_table.xml'
    df['SIMAP_feature_table_XML_file_path'] = df.simap_filename_base + '_feature_table.xml'
    df['SIMAP_homologues_XML_file'] = df.A2_protein_name + '_homologues.xml'
    df['SIMAP_homologues_XML_file_path'] = df.simap_filename_base + '_homologues.xml'
    df['SIMAP_csv_from_XML'] = df.A2_protein_name + '.csv'
    df['SIMAP_csv_from_XML_path'] = df.simap_filename_base + '.csv'
    df['SIMAP_csv_from_XML_tarfile'] = df.simap_filename_base + '.csv.tar.gz'
    df['SIMAP_csv_analysed'] = df.A2_protein_name + '_analysed.csv'
    df['SIMAP_csv_analysed_path'] = df.simap_filename_base + '_analysed.csv'
    df['output_tarfile'] = df.A2_protein_name + '_outputfiles.tar.gz'
    df['output_tarfile_path'] = df.simap_filename_base + '_outputfiles.tar.gz'
    df['csv_SIMAP_homologues_kept_for_statistical_analysis'] = df.simap_filename_base + '_homologues_for_stat.csv'

    # name the fasta file with the TMD seqs (eg A0A1F4_EYS_DROME_TMD_sequences_of_homologues.fas)
    df['fasta_file_path'] = df.simap_filename_base + '_simap_TMD_seq_homologues.fas'

    # name the second fasta file (eg. A0T0U2_PSBE_THAPS_simap_TMD_seq_homol_&_surrounding.fas)
    df['fasta_file_plus_surr_path'] = df.simap_filename_base + '_simap_TMD_seq_homol_&_surrounding.fas'
    df[
        'fasta_file_with_homologues_kept_for_statistical_analysis'] = df.simap_filename_base + '_simap_TMD_seq_kept_stat_analysis.fas'
    df['csv_file_av_cons_ratios_hits'] = df.simap_filename_base + '_cons_ratios.csv'
    '''
    FOR multiple TMDs, create a BASE from which the TMDs can be numbered
    '''
    df['fasta_file_BASENAME'] = df.A2_protein_name + '_simap_seq_homologues_'
    df['fasta_file_BASENAMEPATH'] = df.simap_filename_base + '_simap_seq_homologues_'

    # name the second fasta file (eg. A0T0U2_PSBE_THAPS_simap_TMD_seq_homol_&_surrounding.fas)
    df['fasta_file_plus_surr_path_BASENAME'] = df.A2_protein_name + '_simap_seq_homol_&_surrounding_'
    df['fasta_file_plus_surr_path_BASENAMEPATH'] = df.simap_filename_base + '_simap_seq_homol_&_surrounding_'

    # create a basename for the output histograms
    df['AAIMON_hist_BASENAME'] = df.A2_protein_name + '_AAIMON_hist'
    df['AAIMON_hist_BASENAMEPATH'] = df.simap_filename_base + '_AAIMON_hist'
    df['AASMON_hist_BASENAME'] = df.A2_protein_name + '_AASMON_hist'
    df['AASMON_hist_BASENAMEPATH'] = df.simap_filename_base + '_AASMON_hist'

    df['csv_file_av_cons_ratios_hits_BASENAME'] = df.A2_protein_name + '_cons_ratios_'
    df['csv_file_av_cons_ratios_hits_BASENAMEPATH'] = df.simap_filename_base + '_cons_ratios_'

    # save to a csv
    df.to_csv(pathdict["dfout04_uniprotcsv_incl_paths"], sep=",", quoting=csv.QUOTE_NONNUMERIC)
    # save to Excel
    writer = pd.ExcelWriter(pathdict["dfout03_uniprotxlsx"])  # engine='xlsxwriter'
    df.to_excel(writer, sheet_name='dfout02')
    df.T.to_excel(writer, sheet_name='dfout01')
    writer.save()
    writer.close()
    return

####################################################################
#               HOPEFULLY NO LONGER NECESSARY
####################################################################
#A## variables are included only to help navigate the document in PyCharm
# A04_setup_df_dtypes = set_["run_A04_setup_df_dtypes"]
# if A04_setup_df_dtypes:
#     logging.info('~~~~~~~~~~~~  starting A04_setup_df_dtypes   ~~~~~~~~~~~~')
#     #test if the dataframe has already been created, otherwise reopen from csv file
#     # if 'df' in globals():
#     #     if isinstance(df, pd.DataFrame):
#     #         logging.info('first protein acc = %s, df already exists, continuing with A04_setup_df_dtypes' % df.iloc[0][0])
#     # else:
#     #     logging.info('df loaded from %s' % pathdict["dfout01_uniprotcsv"])
#     df = pd.read_csv(pathdict["dfout01_uniprotcsv"], sep=",", quoting=csv.QUOTE_NONNUMERIC,index_col=0)
#     '''Create Pandas Dataframe with the protein name and file locations, etc'''
#     original_columns = list(df.columns)
#     columns_added_after_SIMAP_analysis = ['organism_domain', 'protein_name', 'SIMAP_feature_table_XML_file_path',
#                                           'SIMAP_homologues_XML_file_path', 'SIMAP_csv_analysed_path',
#                                           'SIMAP_input_seq_details_dict', 'SIMAP_filter_string',
#                                           'SIMAP_resultSpecification_dict', 'database_details_dict',
#                                           'simap_version', 'SIMAP_total_hits',
#                                           'fasta_file_path',
#                                           'fasta_file_plus_surr_path',
#                                           'query_md5', 'query_length', 'query_selfscore', 'query_sequenceid',
#                                           'total_number_of_simap_hits',
#                                           'query_sequence_from_homologue_XML_file',
#                                           'number_of_hits_in_homologue_XML_file',
#                                           'kept_after_redundancy_check',
#                                           'number_of_valid_hits',
#                                           'df_mem_nonmem_ratios',
#                                           'mean_ratio_ident_mem_to_nonmem',
#                                           'stdev_ratio_ident_mem_to_nonmem',
#                                           'protein_kept_for_statistical_analysis',
#                                           'fasta_file_with_homologues_kept_for_statistical_analysis',
#                                           'csv_file_av_cons_ratios_hits',
#                                           'csv_SIMAP_homologues_kept_for_statistical_analysis'
#                                           ]
#     new_unique_column_list = set(original_columns + columns_added_after_SIMAP_analysis)
#     #add extra columns
#     df = df.reindex(index=df.index, columns=new_unique_column_list)
#     #sort columns
#     df = df.sort_index(axis=1)
#
#     #to avoid trouble with dtypes, change all the new columns to dtype=object
#     for column in columns_added_after_SIMAP_analysis:
#         df[column] = df[column].astype(object)
#
#     '''Useful debugging tool: check for duplicates'''
#     import collections
#     list01 = original_columns + columns_added_after_SIMAP_analysis
#     duplicate_columns = [x for x, y in collections.Counter(list01).items() if y > 1]
#     logging.info('%d duplicate_columns found %s' % (len(duplicate_columns), duplicate_columns))