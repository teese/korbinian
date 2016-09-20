#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Author:         Mark Teese
Created:        Fri Oct 25 17:45:09 2013
Dependencies:   utilities file with many small functions "E:/Stephis/Projects/Programming/Python/scripts/mtutils.py"
                Python 3.x
                numpy
                Biopython.
                Java
                Javascript for SIMAP accession
Purpose:        Analysis of BLAST-like SIMAP data
Credits:        All sections by Mark Teese. The Javascript for SIMAP accession was written by Peter Honigschmidt.
Further Details:For details regarding the java access to SIMAP, see here: E:\Stephis\Projects\Programming\Python\programs\eaSimap_MT_notes.txt

PROBLEMS and TASKS
 - check why for Q3KR37_GRM1B_HUMAN homologue 119 the TMD hit and hit_plus_surr_seq do not match! + LLLVISCVICFRYCSVVILSPHHPFHSTPPLFPGTRGLCASACFPTPSALSPHPSPWLGCRVLTPQIRVLLSVPAHRIPPARWHAWVTLSHSRAPSAGRAQGASSALRLLFPPSLVLLVVLNM	DSFHLQSVSKLLLVISCV
 - check why the total_number_of_simap_hits gives the wrong number 
 (df.loc[acc, 'total_number_of_simap_hits'] = query_sequence_node[0].attrib['number_hits'])
 - delete protein_name from output excel file
 - change disallowed words not in desc to true if the description file is empty
 - write average AAIMON for all TMDs into the simap summary csv file
 - filter AAIMON by number of VALID homologues
 - double-check that inplace = True for indexing really doesn't return a copy
 - check if there is a filter for very small nonTMD regions (<20)
 - check if written correctly in new files
 query_md5    query_selfscore    query_sequence_from_homologue_XML_file    query_sequenceid
 - check why csv file with variants gives extra lines for Q14524
~~ DONE - create nonredundant dataset for AAIMON statistics (PREFERABLY WITH BLAST or Uniprot90!) ~~ DONE
 - Nonredundancy test needs to be tested: sp|Q8N349|OR2LD_HUMAN, sp|Q8NH16|OR2L2 together in list 79
   [could be just slight alignment differences]
 - print out proteins with highest and lowest AAIMON for visual analysis
 - check relationship to abundance of each amino acid, and GxxxG motif
 - calculate hydrophobic core / lipid interface conservation as hinted in Mokrab et al.
 - examine beta-barrel proteins and GPI anchors
 - examine TMD/REST correlation to tissues, and cellular location e.g. "cis-Golgi network membrane"
 - include percentage identity of TMD in statistical analysis filters!!! FIX FILTERS!!!
 - if errors in XML, query_md5 = nan, need to drop these from the final dataframe before collecting md5s into a list for redundancy check
 - include analysis of tmd_len/nonTMD_len vs the AAIMON ratio
 - include analysis within multiple cutoffs
 - use numpy to judge how normal the histograms are
 - check why total number of homologues yields extremely high values in final graph
 - remove try/except loops for corrupt files: 
 if isinstance(<path_to_file>, file):
    with open(<path_to_file>) as blabla:
        .......do some crazy stuff....
else:
   sys.stderror("File {0} is corrupt... what an ass...".format(<path_to_file>))

 - check if regex.compile speeds up the regex search (p208, McKinney)
#both stat and fasta filters. NEEDS TO BE UPDATED LATER TO ONLY REFLECT WORKING STAT FILTERS
        df_SIMAP_homologues_kept_for_statistical_analysis = df_SIMAP_homologues_kept_for_statistical_analysis_stat_only[
        df_SIMAP_homologues_kept_for_statistical_analysis_stat_only.match_TMD_added_to_FastA_alignment == True]
 - change matplotlib backend to agg http://stackoverflow.com/questions/21321292/using-matplotlib-when-display-is-undefined
 - specify the dtypes in some of the columns in the csv-file before loading. change np.nan to empty strings, for example. This avoids the DtypeWarning.
 - FIX:no data saved for parameters in simap_homologue_root[0][0][0][0].iter('parameters'):
 - calculate all AAIMAN, not just those within cutoff
 - alter calculations to avaid measring the following (maybe both gappedID and nongappedID cutoffs?)
 
1510917          150 ISLDLYAGALFVHICLGWNFYLSTILTLGITALYTIAGGLAAVIYTDALQ    199
                     |||||||||||||||||||||||||||| : |||||              
64390643         151 ISLDLYAGALFVHICLGWNFYLSTILTLTVAALYTIT-------------    187

1510917          200 TLIMVVGAVILTIKAFDQIGGYGQLEAAYAQAIPSRTIANTTCHLPRTDA    249
                                   |||||||| |||||||||:||| : |||||||| ||
64390643         188 --------------AFDQIGGYEQLEAAYAQAVPSRIVPNTTCHLPRADA    223
 - in rare cases the only alignable sequence is the TMD!! The sequence excluding TMD is therefore '', too small for an analysis, but as long as the program runs successfully this can be filtered out later, at the "calculate conservation" step using pandas dataframes
TM01_perc_sim starts with 0 instead of 1.0. Simply remove column?
 - TM01 plus surr gives a wierd result for O49929 in the uniprot file. Why?
"""
import os
import korbinian
import pandas as pd

def run_korbinian(excel_file_with_settings):

    set_ = korbinian.common.create_settingsdict(excel_file_with_settings)
    list_number = set_["uniprot_list"]

    logging = korbinian.common.setup_keyboard_interrupt_and_error_logging(set_, list_number)
    logging.warning("list_number : {}".format(list_number))

    uniprot_folder_sel = os.path.join(set_["uniprot_dir"], 'selected')
    selected_uniprot_records_flatfile = os.path.join(uniprot_folder_sel,'List%02d_selected_uniprot_records_flatfile.txt' % list_number)

    # set a base folder for the summaries, e.g. "D:\Databases\main\summaries\05\List05 [[_summary.csv]]"
    base_filename_summaries = os.path.join(set_["summaries_dir"], '%02d' % list_number, 'List%02d' % list_number)
    excelfile_with_uniprot_accessions = os.path.join(base_filename_summaries, '.xlsx')

    # create dictionary of paths for output files
    pathdict = korbinian.common.create_pathdict(base_filename_summaries)

    ########################################################################################
    #                                                                                      #
    #      prot_list, OMPdb (create a list of proteins from the OMPdb database)            #
    #                                                                                      #
    ########################################################################################

    if set_["OMPdb_extract_omp_IDs_from_nr_fasta"]:
        korbinian.prot_list.extract_omp_IDs_from_nr_fasta(set_["omp_nr_fasta"], set_["omp_ID_nr_txt"], logging)

    if set_["OMPdb_parse_OMPdb_all_selected_to_csv"]:
        korbinian.prot_list.parse_OMPdb_all_selected_to_csv(set_["omp_ID_nr_txt"], set_["OMPdb_all_flatfile"], set_["OMPdb_summary_nr_csv"], logging)

    if set_["OMPdb_get_omp_TM_indices_and_slice_from_summary_table"]:
        korbinian.prot_list.get_omp_TM_indices_and_slice_from_summary_table(set_["OMPdb_summary_nr_csv"], set_["OMPdb_summary_csv_with_TM_seqs"], logging)

    ########################################################################################
    #                                                                                      #
    #      prot_list, UniProt (create a list of proteins from the UniProt database)        #
    #                                                                                      #
    ########################################################################################

    if set_["create_nonred_uniprot_flatfile_via_uniref"] == True:
        korbinian.prot_list.create_nonred_uniprot_flatfile_via_uniref(set_, uniprot_folder_sel, list_number, selected_uniprot_records_flatfile, logging)

    if set_["run_parse_large_flatfile_with_list_uniprot_accessions"]:
        input_accession_list= set_["list_of_uniprot_accessions"]
        korbinian.prot_list.parse_large_flatfile_with_list_uniprot_accessions(input_accession_list, uniprot_folder_sel, list_number, logging, selected_uniprot_records_flatfile)

    if set_["run_retrieve_uniprot_data_for_acc_list_in_xlsx_file"]:
        korbinian.prot_list.retrieve_uniprot_data_for_acc_list_in_xlsx_file(excelfile_with_uniprot_accessions, logging, selected_uniprot_records_flatfile)

    if set_["run_create_csv_from_uniprot_flatfile"]:
        ''' ~~ DETERMINE START AND STOP INDICES FOR TMD PLUS SURROUNDING SEQ ~~ '''
        n_aa_before_tmd = set_["n_aa_before_tmd"]
        n_aa_after_tmd = set_["n_aa_after_tmd"]
        list_summary_csv_path = pathdict["list_summary_csv"]
        korbinian.prot_list.create_csv_from_uniprot_flatfile(selected_uniprot_records_flatfile, n_aa_before_tmd, n_aa_after_tmd, set_['analyse_signal_peptides'], logging, list_summary_csv_path)
    ########################################################################################
    #                                                                                      #
    #                            run_setup_df_file_locations                               #
    #                                                                                      #
    ########################################################################################

    if set_["run_setup_df_file_locations"]:
        korbinian.prot_list.setup_file_locations_in_df(set_, pathdict)

    ########################################################################################
    #                                                                                      #
    #                         run simap download, parse simap                              #
    #                                                                                      #
    ########################################################################################

    if set_["run_retrieve_simap_feature_table_and_homologues_from_list_in_csv"]:
        korbinian.simap.download_homologues_from_simap(pathdict, set_, logging)

    if set_["run_parse_simap_to_csv"]:
        korbinian.simap.parse_SIMAP_to_csv(pathdict, set_, logging)

    ########################################################################################
    #                                                                                      #
    #            run_create_fasta, run_calculate_AAIMON_ratios                             #
    #                                                                                      #
    ########################################################################################

    if set_["slice_TMDs_from_homologues"]:
        korbinian.cons_ratio.slice_TMDs_from_homologues(pathdict, set_, logging)

    if set_["run_create_fasta"]:
        korbinian.fasta.filter_and_save_fasta(pathdict, set_, logging)

    if set_["run_calculate_AAIMON_ratios"]:
        korbinian.cons_ratio.calculate_AAIMON_ratios(pathdict, set_, logging)

    if set_["run_gather_AAIMON_ratios"]:
        korbinian.cons_ratio.gather_AAIMON_ratios(pathdict, set_, logging)

    ########################################################################################
    #                                                                                      #
    #                             gap density analysis                                     #
    #                                                                                      #
    ########################################################################################

    if set_["run_calculate_gap_densities"]:
        korbinian.gap.calculate_gap_densities(pathdict, set_, logging)

    if set_["run_create_graph_of_gap_density"]:
        korbinian.gap.create_graph_of_gap_density(pathdict, set_, logging)

    ########################################################################################
    #                                                                                      #
    #                 conservation ratio (AAIMON ratio) figures                            #
    #                                                                                      #
    ########################################################################################

    '''+++++++++++++++ Summary figures describing the conservation ratios of proteins in the list ++++++++++++++++++'''
    if set_["run_save_figures_describing_proteins_in_list"]:
        korbinian.cons_ratio.save_figures_describing_proteins_in_list(pathdict, set_, logging)

    '''+++++++++++++++ Summary figures describing the conservation ratios of proteins in the list ++++++++++++++++++'''
    if set_["run_compare_lists"]:
        korbinian.cons_ratio.compare_rel_con_lists(pathdict, set_, logging)

    ########################################################################################
    #                                                                                      #
    #                                     old stuff                                        #
    #                                                                                      #
    ########################################################################################

    '''+++++++++++++++TMD CONSERVATION (OLD)++++++++++++++++++'''
    if set_["old_calculate_TMD_conservation"]:
        korbinian.old.OLD_calculate_TMD_conservation(pathdict, set_, logging)

    '''The gapped identity analysis gives the average TMD mem/nonmem conservation for homologues at different levels of similarity (eg. close homologues vs far homologues)
    This is a funny system. The dictionary keys are the 5000 hits. This would be much more efficient and simple to re-write in pandas style.
    '''
    if set_["old_calculate_TMD_conservation_by_gappedIdentity"]:
        korbinian.old.OLD_calculate_TMD_conservation_by_gappedIdentity(pathdict, set_, logging, list_number)

    if set_["old_run_stat_analysis_sim_ratios_in_dfout05"]:
        korbinian.old.old_run_stat_analysis_sim_ratios_in_dfout05(pathdict, set_, logging)

    if set_["old_fix_dfout05_simapcsv_by_adding_query_md5"]:
        korbinian.old.OLD_fix_dfout05_simapcsv_by_adding_query_md5(pathdict, logging)