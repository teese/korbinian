import csv
import os
import pandas as pd
from time import strftime
import korbinian.mtutils as utils
import signal
import unicodedata

def setup_keyboard_interrupt_and_error_logging(set_, list_number):
    ''' -------Setup keyboard interrupt----------
    '''
    # import arcgisscripting

    def ctrlc(sig, frame):
        raise KeyboardInterrupt("CTRL-C!")
    signal.signal(signal.SIGINT, ctrlc)
    '''+++++++++++++++LOGGING++++++++++++++++++'''
    date_string = strftime("%Y%m%d_%H_%M_%S")

    # designate the output logfile
    logfile = os.path.join(set_["logfile_dir"],'List%s_%s_logfile.log' % (list_number, date_string))
    # a file to keep a record of the log settings used for that script
    logging = utils.setup_error_logging(logfile)
    return logging

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

    list_paths_to_normalise = ["data_folder", "uniprot_folder", "data_harddrive", "eaSimap_path", "logfile_dir",
                               "summaries_dir", "simap_database_dir", "list_of_uniprot_accessions"]
    # normalise the paths for selected columns, so that they are appropriate for the operating system
    for path in list_paths_to_normalise:
        if path in set_:
            set_[path] = os.path.normpath(set_[path])

    return set_


def create_pathdict(base_filename_summaries):
    pathdict = {}
    pathdict["base_filename_summaries"] = base_filename_summaries
    # currently the protein list summary (each row is a protein, from uniprot etc) is a csv file
    pathdict["list_summary_csv"] = '%s_summary.csv' % base_filename_summaries
    # add the base path for the sub-sequences (SE01, SE02, etc) added by the user
    pathdict["list_user_subseqs_csv"] = '%s_user_subseqs.csv' % base_filename_summaries
    pathdict["list_user_subseqs_xlsx"] = '%s_user_subseqs.xlsx' % base_filename_summaries

    pathdict["dfout01_uniprotcsv"] = '%s_uniprot.csv' % base_filename_summaries
    pathdict["dfout02_uniprotTcsv"] = '%s_uniprotT.csv' % base_filename_summaries
    pathdict["dfout03_uniprotxlsx"] = '%s_uniprot.xlsx' % base_filename_summaries
    pathdict["dfout04_uniprotcsv_incl_paths"] = '%s_uniprot_incl_paths.csv' % base_filename_summaries
    #pathdict["dfout05_simapcsv"] = '%s_simap.csv' % base_filename_summaries
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

