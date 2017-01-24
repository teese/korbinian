#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Author:         Mark Teese
Created:        Fri Oct 25 17:45:09 2013
Dependencies:   Python 3.x
                numpy
                Bio
                Java
                Javascript for SIMAP accession
                pandas
Purpose:        Analysis of BLAST-like SIMAP data
Credits:        All sections by Mark Teese. The Javascript for SIMAP accession was written by Peter Honigschmidt.
Further Details:For details regarding the java access to SIMAP, see here: E:\Stephis\Projects\Programming\Python\programs\eaSimap_MT_notes.txt

"""
import argparse
import os
import korbinian
import ast

# read the command line arguments
parser = argparse.ArgumentParser()
# add only a single argument, the path to the settings file.
parser.add_argument("-s",  # "-settingsfile",
                    help=r'Full path to your excel settings file.'
                         r'E.g. "C:\Path\to\your\settingsfile.xlsx"')

if __name__ == "__main__":
    print('\nRun korbinian as follows:')
    print(r'python "C:\Path\to\run_korbinian.py" -s "C:\Path\to\your\settingsfile.xlsx"')
    # get the command-line arguments
    args = parser.parse_args()
    # args.s is the excel_settings_file input by the user
    # convert the excel settings file to a python dictionary, s
    s = korbinian.common.create_settingsdict(args.s)

    protein_list_number = ast.literal_eval(str(s['protein_list_number']))
    if type(protein_list_number) is int:
        list_number = protein_list_number
        korbinian.run_korbinian.run_statements(s, list_number)
    else:
        print('\nmultiple lists entered! starting serial analysis of multiple lists!\n\nlists to analyse: {}\n\n'.format(protein_list_number))
        for list_number in protein_list_number:
            korbinian.run_korbinian.run_statements(s, list_number)
            print('\n\n~~~~~~~~~~~~         List {} finished           ~~~~~~~~~~~~\n\n'.format(list_number))

def run_statements(s, list_number):
    # setup error logging
    logging = korbinian.common.setup_keyboard_interrupt_and_error_logging(s, list_number)
    # print the list number describing the protein list
    logging.warning("list_number : {}".format(list_number))

    # set a base folder for the summaries, e.g. "D:\Databases\summaries\05\" for list 05
    base_filename_summaries = os.path.join(s["data_dir"], "summaries", '%02d' % list_number, 'List%02d' % list_number)

    # create dictionary of paths for output files
    # for example the basic pathdict["list_summary_csv"] for list 5 is "D:\Databases\summaries\05\List05_summary.csv"
    pathdict = korbinian.common.create_pathdict(base_filename_summaries)

    ########################################################################################
    #                                                                                      #
    #      prot_list, OMPdb (create a list of proteins from the OMPdb database)            #
    #                                                                                      #
    ########################################################################################

    if s["OMPdb_extract_omp_IDs_from_nr_fasta"]:
        ListXX_OMPdb_nr_fasta = os.path.join(s["data_dir"], "OMPdb", "List{:02d}_OMPdb_nr_fasta.txt".format(list_number))
        ListXX_OMPdb_nr_acc = os.path.join(s["data_dir"], "OMPdb", "List{:02d}_OMPdb_nr_acc.txt".format(list_number))
        korbinian.prot_list.parse_OMPdb.extract_omp_IDs_from_nr_fasta(ListXX_OMPdb_nr_fasta, ListXX_OMPdb_nr_acc, logging)

    if s["OMPdb_parse_OMPdb_all_selected_to_csv"]:
        ListXX_OMPdb_nr_acc = os.path.join(s["data_dir"], "OMPdb", "List{:02d}_OMPdb_nr_acc.txt".format(list_number))
        ListXX_OMPdb_redundant_flatfile = os.path.join(s["data_dir"], "OMPdb", "List{:02d}_OMPdb_redundant_flatfile.flat".format(list_number))
        OMPdb_list_summary_csv = os.path.join(s["data_dir"], "OMPdb", "List{:02d}_OMPdb_summary.csv".format(list_number))
        korbinian.prot_list.parse_OMPdb.parse_OMPdb_all_selected_to_csv(ListXX_OMPdb_nr_acc, ListXX_OMPdb_redundant_flatfile, OMPdb_list_summary_csv, logging)

    if s["OMPdb_get_TM_indices_and_slice"]:
        OMPdb_list_summary_csv = os.path.join(s["data_dir"], "OMPdb", "List{:02d}_OMPdb_summary.csv".format(list_number))
        list_summary_csv = pathdict["list_summary_csv"]
        korbinian.prot_list.parse_OMPdb.get_omp_TM_indices_and_slice_from_summary_table(OMPdb_list_summary_csv, list_summary_csv, logging)

    ########################################################################################
    #                                                                                      #
    #      prot_list, UniProt (create a list of proteins from the UniProt database)        #
    #                                                                                      #
    ########################################################################################

    # define the uniprot directory with selected records
    uniprot_dir_sel = os.path.join(s["data_dir"], 'uniprot', 'selected')
    selected_uniprot_records_flatfile = os.path.join(uniprot_dir_sel, 'List%02d_selected_uniprot_records_flatfile.txt' % list_number)
    excelfile_with_uniprot_accessions = os.path.join(base_filename_summaries, '.xlsx')

    if s["create_nonred_uniprot_flatfile_via_uniref"] == True:
        korbinian.prot_list.uniprot_nonredundant.create_nonred_uniprot_flatfile_via_uniref(s, uniprot_dir_sel, list_number, selected_uniprot_records_flatfile, logging)

    if s["run_parse_large_flatfile_with_list_uniprot_accessions"]:
        input_accession_list_path = os.path.join(s["data_dir"], "uniprot", "selected", "List{:02d}_uniprot_accessions.txt".format(list_number))
        korbinian.prot_list.uniprot_retrieve.parse_large_flatfile_with_list_uniprot_accessions(input_accession_list_path, uniprot_dir_sel, list_number, logging, selected_uniprot_records_flatfile)

    if s["run_retrieve_uniprot_data_for_acc_list_in_xlsx_file"]:
        input_uniprot_flatfile = "needs to be defined if you use this function!"
        korbinian.prot_list.uniprot_retrieve.retrieve_uniprot_data_for_acc_list_in_xlsx_file(excelfile_with_uniprot_accessions, input_uniprot_flatfile, selected_uniprot_records_flatfile, logging)

    if s["run_create_csv_from_uniprot_flatfile"]:
        ''' ~~ DETERMINE START AND STOP INDICES FOR TMD PLUS SURROUNDING SEQ ~~ '''
        n_aa_before_tmd = s["n_aa_before_tmd"]
        n_aa_after_tmd = s["n_aa_after_tmd"]
        list_summary_csv_path = pathdict["list_summary_csv"]
        korbinian.prot_list.uniprot_parse.create_csv_from_uniprot_flatfile(selected_uniprot_records_flatfile, n_aa_before_tmd, n_aa_after_tmd, s['analyse_signal_peptides'], logging, list_summary_csv_path)

    ########################################################################################
    #                                                                                      #
    #                            run_setup_df_file_locations                               #
    #                                                                                      #
    ########################################################################################

    if s["run_setup_df_file_locations"]:
        korbinian.prot_list.prot_list.setup_file_locations_in_df(s, pathdict)

    ########################################################################################
    #                                                                                      #
    #                         run simap download, parse simap                              #
    #                                                                                      #
    ########################################################################################
    if s["run_download_homologues"]:
        korbinian.simap_download.download_homologues_from_simap(pathdict, s, logging)

    if s["run_parse_simap_to_csv"]:
        korbinian.simap_parse.run_parse_simap_to_csv(pathdict, s, logging)

    ########################################################################################
    #                                                                                      #
    #            run_create_fasta, run_calculate_AAIMON_ratios                             #
    #                                                                                      #
    ########################################################################################

    if s["run_slice_TMDs_from_homologues"]:
        korbinian.cons_ratio.slice.run_slice_TMDs_from_homologues(pathdict, s, logging)

    if s["run_create_fasta"]:
        korbinian.fasta.run_create_fasta(pathdict, s, logging)

    if s["run_calculate_AAIMON_ratios"]:
        korbinian.cons_ratio.cons_ratio.run_calculate_AAIMON_ratios(pathdict, s, logging)

    if s['run_filter_truncated_alignments']:
        korbinian.cons_ratio.cons_ratio.throw_out_truncated_sequences(pathdict, s, logging)

    if s["run_gather_AAIMON_ratios"]:
        korbinian.cons_ratio.gather.gather_AAIMON_ratios(pathdict, logging, s)

    ########################################################################################
    #                                                                                      #
    #                             gap density analysis                                     #
    #                                                                                      #
    ########################################################################################

    if s["run_calculate_gap_densities"]:
        korbinian.gap.run_calculate_gap_densities(pathdict, s, logging)

    if s["run_gather_gap_densities"]:
        korbinian.gap.gather_gap_densities(pathdict, s, logging)

    if s["run_create_graph_of_gap_density"]:
        korbinian.gap_figs.create_graph_of_gap_density(pathdict, s, logging)

    if s["run_fastagap_save"]:
        korbinian.fastagap.run_fastagap_save(pathdict, s, logging)

    if s["run_calc_fastagap_densities"]:
        korbinian.fastagap.run_calc_fastagap_densities(pathdict, s, logging)


    ########################################################################################
    #                                                                                      #
    #                 conservation ratio (AAIMON ratio) figures                            #
    #                                                                                      #
    ########################################################################################

    '''+++++++++++++++ Summary figures describing the conservation ratios of proteins in the list ++++++++++++++++++'''
    if s["run_save_figures_describing_proteins_in_list"]:
        korbinian.cons_ratio.figs.save_figures_describing_proteins_in_list(pathdict, s, logging)

    '''+++++++++++++++ Summary figures describing the conservation ratios of proteins in the list ++++++++++++++++++'''
    if s["run_compare_lists"]:
        korbinian.cons_ratio.compare_lists.compare_rel_con_lists(pathdict, s, logging)


    if s["run_keyword_analysis"]:
        korbinian.cons_ratio.keywords.keyword_analysis(pathdict, s, logging)


    if s['send_email_when_finished']:
        korbinian.utils.send_email_when_finished(s, pathdict, list_number)