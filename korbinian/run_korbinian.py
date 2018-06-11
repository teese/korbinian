#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Author:         Mark Teese
Created:        Fri Oct 25 17:45:09 2013
Dependencies:   Python 3.x
                numpy
                pandas
                Bio
                Java (for SIMAP homologue downloads)
Purpose:        Protein Data Science
                Analysis of evolutionary sequences of transmembrane proteins
Credits:        Mark Teese
                Martin Ortner
                Shenger Wang
                Rimma Jenske
License         Released under the permissive MIT license.
"""
import argparse
import ast
import os
import sys
from shutil import copyfile

import numpy as np
import pandas as pd

import korbinian
import korbinian.utils as utils


#Initialize a modified ArgumentParser from the argparse package
class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message + "\n")
        self.print_help()
        sys.exit(2)

# read the command line arguments
parser = MyParser()
# add only a single argument, the path to the settings file.
parser.add_argument("-s",  # "-settingsfile",
                    help=r'Full path to your excel settings file.'
                         r'E.g. "C:\Path\to\your\settingsfile.xlsx"', required=True)

if __name__ == "__main__":
    sys.stdout.write('\nRun korbinian as follows:')
    sys.stdout.write(r'python "C:\Path\to\run_korbinian.py" -s "C:\Path\to\your\settingsfile.xlsx"')
    # get the command-line arguments
    args = parser.parse_args()
    #check if the provided path lead to an existing file -> IF NOT = exit program
    if not (os.path.exists(args.s) and os.path.isfile(args.s) and os.path.getsize(args.s) > 0):
        sys.stderr.write("\n" + "Provided path to the settings file leads to a non existing or empty file" + "\n")
        sys.exit(2)
    # args.s is the excel_settings_file input by the user
    # convert the excel settings file to a python dictionary, s
    s = korbinian.common.create_settingsdict(args.s)

    # to run "compare lists", replace the list_number with "compare"
    if s['protein_list_number'] == 'compare':
        # open the tab containing the list-specific settings as a dataframe
        df_lists_tab = pd.read_excel(s["excel_file_with_settings"], sheetname="lists", index_col=0)
        # add the dataframe "lists" to the existing settings dictionary
        # this adds max_lipo_homol, rand_TM, rand_nonTM, etc to the dictionary
        #s.update(df_lists.to_dict())
        korbinian.cons_ratio.compare_lists.compare_lists(s, df_lists_tab)
    else:
        # if list_number is not "compare", run either a single list, or a list of protein lists
        protein_list_number = ast.literal_eval(str(s['protein_list_number']))
        # run for a single protein list
        if type(protein_list_number) is int:
            # there is only one list to run, so add the provided number to settings and run
            s["list_number"] = protein_list_number
            korbinian.run_korbinian.run_statements(s)
        elif type(protein_list_number) in (list, tuple):
            # run for multiple lists, ONE AFTER ANOTHER
            # note this is different from the compare_lists function, which compares the output for multiple lists, and creates figures
            sys.stdout.write('\nStarting serial analysis of multiple lists. Any python errors will stop all processing. \n\nlists to analyse: {}\n\n'.format(protein_list_number))
            for list_number in protein_list_number:
                # for each number in the list, add to settings and run
                s["list_number"] = list_number
                korbinian.run_korbinian.run_statements(s)
        else:
            raise ValueError("protein_list_number is neither an int nor a list. Panic and check your code and settings file!")

def run_statements(s):
    list_number = s["list_number"]
    # setup error logging
    logging = korbinian.common.setup_keyboard_interrupt_and_error_logging(s, list_number)
    # print the list number describing the protein list
    logging.warning("list_number : {}".format(list_number))

    # open the tab containing the list-specific settings as a dataframe
    df_list_settings = pd.read_excel(s["excel_file_with_settings"], sheetname="lists", index_col=0)
    relevant_row = df_list_settings.loc[list_number, :].to_dict()
    if np.nan in relevant_row.values():
        raise ValueError("The row for List{} in the lists tab of the settings file is missing some values.".format(list_number))
    # add the relevant row (e.g. for List01) to the existing settings dictionary
    # this adds max_lipo_homol, rand_TM, rand_nonTM, etc to the dictionary
    s.update(relevant_row)

    # set a base folder for the summaries, e.g. "D:\Databases\summaries\05\" for list 05
    base_filename_summaries = os.path.join(s["data_dir"], "summaries", '%02d' % list_number, 'List%02d' % list_number)

    # create dictionary of paths for output files
    # for example the basic pathdict["list_csv"] for list 5 is "D:\Databases\summaries\05\List05_summary.csv"
    pathdict = korbinian.common.create_pathdict(base_filename_summaries, s)

    utils.make_sure_path_exists(pathdict["settings_copy_xlsx"], isfile=True)
    # copy the settings file used for the analysis
    copyfile(s["excel_file_with_settings"], pathdict["settings_copy_xlsx"])

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
        OMPdb_list_csv = os.path.join(s["data_dir"], "OMPdb", "List{:02d}_OMPdb_summary.csv".format(list_number))
        korbinian.prot_list.parse_OMPdb.parse_OMPdb_all_selected_to_csv(ListXX_OMPdb_nr_acc, ListXX_OMPdb_redundant_flatfile, OMPdb_list_csv, logging, s)

    if s["OMPdb_get_TM_indices_and_slice"]:
        OMPdb_list_csv = os.path.join(s["data_dir"], "OMPdb", "List{:02d}_OMPdb_summary.csv".format(list_number))
        list_parsed_csv = pathdict["list_parsed_csv"]
        OMPdb_topology_reliability_cutoff = s["OMPdb_topology_reliability_cutoff"]
        korbinian.prot_list.parse_OMPdb.get_omp_TM_indices_and_slice_from_summary_table(OMPdb_list_csv, list_parsed_csv, OMPdb_topology_reliability_cutoff, logging, s)

    ########################################################################################
    #                                                                                      #
    #      prot_list, UniProt (create a list of proteins from the UniProt database)        #
    #                                                                                      #
    ########################################################################################

    # define the uniprot directory with selected records
    uniprot_dir = os.path.join(s["data_dir"], 'uniprot')
    selected_uniprot_records_flatfile = os.path.join(uniprot_dir, 'selected', 'List%02d_selected_uniprot_records_flatfile.txt' % list_number)

    if s["parse_large_flatfile_with_list_uniprot_accessions"]:
        input_accession_list_path = os.path.join(s["data_dir"], "uniprot", "selected", "List{:02d}_uniprot_accessions.txt".format(list_number))
        korbinian.prot_list.uniprot_retrieve.parse_large_flatfile_with_list_uniprot_accessions(s, input_accession_list_path, uniprot_dir, logging, selected_uniprot_records_flatfile)

    if s["retrieve_uniprot_data_for_acc_list_in_xlsx_file"]:
        input_uniprot_flatfile = "function not implemented!"
        excelfile_with_uniprot_accessions = os.path.join(base_filename_summaries, '.xlsx')
        korbinian.prot_list.uniprot_retrieve.retrieve_uniprot_data_for_acc_list_in_xlsx_file(excelfile_with_uniprot_accessions, input_uniprot_flatfile, selected_uniprot_records_flatfile, logging)

    if s["create_nonred_uniprot_flatfile_via_uniref"] == True:
        korbinian.prot_list.uniprot_nonredundant.create_nonred_uniprot_flatfile_via_uniref(s, uniprot_dir, selected_uniprot_records_flatfile, logging)

    if s["parse_flatfile_to_csv"]:
        n_aa_before_tmd = s["n_aa_before_tmd"]
        n_aa_after_tmd = s["n_aa_after_tmd"]
        list_parsed_csv = pathdict["list_parsed_csv"]
        # whether to analyse signal peptides for this dataset
        analyse_sp = True if "SiPe" in s["regions"] else False
        output = korbinian.prot_list.uniprot_parse.parse_flatfile_to_csv(selected_uniprot_records_flatfile, n_aa_before_tmd, n_aa_after_tmd, analyse_sp, logging, list_parsed_csv)
        logging.info(output)

    ########################################################################################
    #                                                                                      #
    #                             Membrane protein filtering                               #
    #                                                                                      #
    ########################################################################################

    if s["TM_filtering_SCAMPI"]:
        korbinian.filtering.scampi.run_filtering(pathdict, s, logging)

    if s["TM_filtering_SignalP"]:
        korbinian.filtering.signalP.run_filtering(pathdict, s, logging)

    if s["TM_filtering_TMSEG_PSI-BLAST_createDatabase"]:
        korbinian.blast.psiblast.create_BLAST_database(pathdict, s, logging)
    if s["TM_filtering_TMSEG_PSI-BLAST"]:
        korbinian.blast.psiblast.run_BLAST(pathdict, s, logging)
    if s["TM_filtering_TMSEG"]:
        korbinian.filtering.tmseg.run_filtering(pathdict, s, logging)

    ########################################################################################
    #                                                                                      #
    #                            prepare_protein_list                                      #
    #                                                                                      #
    ########################################################################################

    if s["prepare_protein_list"]:
        korbinian.prot_list.prot_list.get_topology_for_prot_list(s, pathdict, logging)
        korbinian.prot_list.prot_list.slice_TMDs_in_prot_list(s, pathdict, logging)
        korbinian.prot_list.prot_list.prepare_protein_list(s, pathdict, logging)

    if s['generate_scampi_input_files']:
        korbinian.prot_list.SCAMPI.generate_scampi_input_files(pathdict, s, logging)

    if s['generate_SignalP_input_files']:
        korbinian.prot_list.SCAMPI.generate_SignalP_input_files(pathdict, s, logging)

    ########################################################################################
    #                                                                                      #
    #                         run simap download, parse simap                              #
    #                                                                                      #
    ########################################################################################
    if s["download_homologues"]:
        korbinian.simap_download.download_homologues_from_simap(pathdict, s, logging)

    if s["parse_simap_to_csv"]:
        korbinian.simap_parse.run_parse_simap_to_csv(pathdict, s, logging)

    ########################################################################################
    #                                                                                      #
    #                         run BLAST and parse results                                  #
    #                                                                                      #
    ########################################################################################
    if s["BLASTp_search"]:
        if s["BLAST_modus"] == "online":
            korbinian.blast.blastp.run_BLAST_online(pathdict, s, logging)
        elif s["BLAST_modus"] == "local":
            korbinian.blast.blastp.run_BLAST_local(pathdict, s, logging)

    if s["BLAST_parser"]:
        korbinian.blast.blast_parser.run(pathdict, s, logging)

    ########################################################################################
    #                                                                                         #
    #            run_create_fasta, run_calculate_AAIMON_ratios                             #
    #                                                                                      #
    ########################################################################################

    if s["slice_TMDs_from_homologues"]:
        korbinian.cons_ratio.slice.run_slice_TMDs_from_homologues(pathdict, s, logging)

    if s["create_fasta"]:
        korbinian.fasta.run_create_fasta(pathdict, s, logging)

    if s["calculate_AAIMON_ratios"]:
        korbinian.cons_ratio.cons_ratio.run_calculate_AAIMONs(pathdict, s, logging)

    if s['filter_truncated_alignments']:
        korbinian.cons_ratio.cons_ratio.throw_out_truncated_sequences(pathdict, s, logging)

    if s["gather_AAIMON_ratios"]:
        # reassign pathdict that could have been recreated during gather depending on settings
        pathdict = korbinian.cons_ratio.gather.gather_AAIMONs(pathdict, logging, s)

    ########################################################################################
    #                                                                                      #
    #                             gap density analysis                                     #
    #                                                                                      #
    ########################################################################################

    if s["calculate_gap_densities"]:
        korbinian.gap.run_calculate_gap_densities(pathdict, s, logging)

    if s["gather_gap_densities"]:
        korbinian.gap.gather_gap_densities(pathdict, s, logging)

    if s["create_graph_of_gap_density"]:
        korbinian.gap_figs.create_graph_of_gap_density(pathdict, s, logging)

    if s["save_fastagap"]:
        korbinian.fastagap.save_fastagap(pathdict, s, logging)

    if s["calc_fastagap_densities"]:
        korbinian.fastagap.run_calc_fastagap_densities(pathdict, s, logging)


    ########################################################################################
    #                                                                                      #
    #                 conservation ratio (AAIMON ratio) figures                            #
    #                                                                                      #
    ########################################################################################
    if s["run_keyword_analysis"]:
        output = korbinian.cons_ratio.keywords.keyword_analysis(pathdict, s, logging)
        logging.info(output)

    '''+++++++++++++++ Summary figures describing the conservation ratios of proteins in the list ++++++++++++++++++'''
    if s["save_figures_describing_proteins_in_list"]:
        return_statement = korbinian.cons_ratio.figs.save_figures_describing_proteins_in_list(pathdict, s, logging)
        logging.info(return_statement)

    '''+++++++++++++++ Summary figures describing the conservation ratios of proteins in the list ++++++++++++++++++'''
    # if s["compare_lists"]:
    #     korbinian.cons_ratio.compare_lists_old.compare_rel_con_lists(pathdict, s, logging)

    if "gather_pretty_alignments" in s.keys():
        if s["gather_pretty_alignments"]:
            korbinian.cons_ratio.gather.gather_pretty_alignments(pathdict, logging, s)

    if s['send_email_when_finished']:
        korbinian.utils.send_email_when_finished(s, pathdict)

    sys.stdout.write('\n~~~~~~~~~~~~         List {} finished           ~~~~~~~~~~~~\n'.format(list_number))
