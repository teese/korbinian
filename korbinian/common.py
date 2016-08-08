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
    date_string = strftime("%Y%m%d")

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

def setup_file_locations_in_df(set_, pathdict):
    df = pd.read_csv(pathdict["list_summary_csv"])
    # set up a folder to hold the SIMAP BLAST-like output
    # note that at the moment, files cannot be compressed
    simap_database_dir = set_['simap_database_dir']
    homol_dir = set_["homol_dir"]

    if "uniprot_entry_name" in df.columns:
        # join the accession and entry name to create a "protein name" for naming files
        df['protein_name'] = df.uniprot_acc + '_' + df.uniprot_entry_name
    else:
        # the list of proteins did not come from UniProt. Simply use the accession to name the files.
        df['protein_name'] = df.uniprot_acc
    df['first_two_letters_of_uniprot_acc'] = df['uniprot_acc'].apply(lambda x: x[0:2])
    df['simap_filename_base'] = simap_database_dir + '/' + df.first_two_letters_of_uniprot_acc + '/' + df.protein_name
    # normalise path to suit operating system
    df['simap_filename_base'] = df['simap_filename_base'].apply(lambda x: os.path.normpath(x))

    df['homol_base'] = homol_dir + '/' + df.first_two_letters_of_uniprot_acc + '/' + df.protein_name
    # normalise path to suit operating system
    df['homol_base'] = df['homol_base'].apply(lambda x : os.path.normpath(x))

    # create filenames for simap output
    df['SIMAP_tar'] = df.simap_filename_base + '_SIMAP.tar.gz'
    df['SIMAP_feature_table_XML_file_path'] = df.simap_filename_base + '_feature_table.xml'
    df['SIMAP_homologues_XML_file_path'] = df.simap_filename_base + '_homologues.xml'
    your_name = unicodedata.normalize('NFKD', set_["your_name"][:20]).encode('ascii', 'ignore').decode("utf-8")
    df['SIMAP_download_date_file_path'] = df.simap_filename_base + '.{}.{}.txt'.format(strftime("%Y%m%d"), your_name)

    df['homol_csv_zip'] = df['homol_base'] + '.csv.zip'
    df['homol_fasta_zip'] = df['homol_base'] + '_fasta.zip'
    df['homol_cons_ratio_zip'] = df['homol_base'] + '_cons_ratio.zip'
    df['homol_gaps_zip'] = df['homol_base'] + '_gaps.zip'
    df['SIMAP_temp_csv_from_XML_path'] = df['homol_base'] + '.csv'

    # name the fasta file with the TMD seqs (eg A0A1F4_EYS_DROME_TMD_sequences_of_homologues.fas)
    '''
    FOR multiple TMDs, create a BASE from which the TMDs can be numbered
    '''
    df['fasta_file_BASENAME'] = df.protein_name + '_homol_seq_'
    df['fasta_file_BASENAMEPATH'] = df.homol_base + '_homol_seq_'

    # name the second fasta file (eg. A0T0U2_PSBE_THAPS_simap_TMD_seq_homol_&_surrounding.fas)
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
    # df['SIMAP_temp_csv_from_XML_path'] = df.simap_filename_base + '.csv'
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
    # indicate that the setup_file_locations_in_df function has been run
    df['setup_file_locations_in_df'] = True
    # save to a csv
    df.to_csv(pathdict["list_summary_csv"], sep=",", quoting=csv.QUOTE_NONNUMERIC)