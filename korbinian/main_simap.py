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
 - conduct word-cloud analysis of functions for high and low AAIMON
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
"""

# def run_main():
import os

'''
These parameters need to be changed manually by each user, and for each time the program is run:
'''
settings_number = 2
main_folder_harddrive = 'D:\\'
#main_folder_list = [main_folder_harddrive, 'teeselab', 'mark', 'files', 'main']
#main_folder_list = [main_folder_harddrive, 'Schweris', 'Projects', 'Programming', 'Python', 'files', 'MAIN']
main_folder_list = [main_folder_harddrive, 'Databases', 'main']
main_folder = os.path.join(*main_folder_list)
excel_file_with_settings = os.path.join(main_folder, 'settings', 'korbinian_python_run_settings_%02d.xlsx' % settings_number)

'''
The MAIN script starts here
'''
# import publicly available python modules
import csv
import xml.etree.ElementTree as ET
import logging
import pandas as pd

import importlib
import os
import korbinian
from korbinian.uniprot import convert_uniprot_list_to_nonred_ff_via_uniref, parse_large_flatfile_with_list_uniprot_accessions, retrieve_uniprot_data_for_acc_list_in_xlsx_file, create_csv_from_uniprot_flatfile
import korbinian.mtutils as utils
import korbinian.rimma_utility as r_utils
import tlabtools.tools as tools
import tlabtools.figures as figs
# reload the modules, so the interpreter doesn't need to be restarted after they are updated
importlib.reload(utils)
importlib.reload(r_utils)
importlib.reload(tools)
importlib.reload(figs)

settingsdict = korbinian.common.create_settingsdict(excel_file_with_settings)

list_number = settingsdict["run_settings"]["uniprot_list"]

# import private tlab python modules
# the utils file contains a number of necessary scripts, this needs to be visible to the system.
# deprecated in favour of the "setup.py" method promoted by Alex
#if settingsdict["file_locations"]["utils"] not in sys.path:
#    sys.path.append(settingsdict["file_locations"]["utils"])

korbinian.common.setup_keyboard_interrupt_and_error_logging(main_folder, list_number, settings_number)

# Folder, where SIMAP download is stored; requires a lot of harddrive; 
data_folder = settingsdict["file_locations"]["data_folder"]
uniprot_folder = os.path.join(data_folder, 'D_uniprot', 'selected')
uniprot_flatfile_of_selected_records = os.path.join(uniprot_folder,'List%s_selected_uniprot_records_flatfile.txt' % list_number)

base_filename_summaries = os.path.join(main_folder, 'summaries', 'List%02d' % list_number)
#base_filename_input_acc_lists = '%s/input_acc_lists/List%02d'  % (main_folder, list_number)
excelfile_with_uniprot_accessions = os.path.join(base_filename_summaries, '.xlsx')

list_of_uniprot_accessions = 'List%02d_uniprot_accessions.txt' % list_number
# create dictionary of paths for output files
pathdict = korbinian.common.create_pathdict(main_folder, base_filename_summaries, list_of_uniprot_accessions, list_number)

#A## variables are included only to help navigate the document in PyCharm
A00_convert_uniprot_list_to_nonred_ff_via_uniref = settingsdict["run_settings"]["uniprot.inputs.convert_uniprot_list_to_nonred_ff_via_uniref"]
if A00_convert_uniprot_list_to_nonred_ff_via_uniref:
    logging.info('~~~~~~~~~~~~  starting A00_convert_redundant_uniprot_list_to_nonred_ff_via_uniref   ~~~~~~~~~~~~')
    overwrite_selected_ff = settingsdict["variables"]["uniprot.convert_uniprot_list_to_nonred_ff_via_uniref.overwrite_selected_ff"]
    if os.path.isfile(uniprot_flatfile_of_selected_records) == False or overwrite_selected_ff == True:
        convert_uniprot_list_to_nonred_ff_via_uniref(settingsdict, list_number, uniprot_folder, logging,
                                                     uniprot_flatfile_of_selected_records)

A01_parse_large_flatfile_with_list_uniprot_accessions = settingsdict["run_settings"]["uniprot.inputs.parse_large_flatfile_with_list_uniprot_accessions"]
if A01_parse_large_flatfile_with_list_uniprot_accessions:
    parse_large_flatfile_with_list_uniprot_accessions(uniprot_folder, list_number, logging, uniprot_flatfile_of_selected_records)

A02_retrieve_uniprot_data_for_acc_list_in_xlsx_file = settingsdict["run_settings"]["uniprot.inputs.retrieve_uniprot_data_for_acc_list_in_xlsx_file"]
if A02_retrieve_uniprot_data_for_acc_list_in_xlsx_file:
    retrieve_uniprot_data_for_acc_list_in_xlsx_file(excelfile_with_uniprot_accessions, logging, uniprot_flatfile_of_selected_records)

A03_create_csv_from_uniprot_flatfile = settingsdict["run_settings"]["uniprot.parse.create_csv_from_uniprot_flatfile"]
if A03_create_csv_from_uniprot_flatfile:
    create_csv_from_uniprot_flatfile(uniprot_flatfile_of_selected_records, settingsdict, logging, pathdict)

#A## variables are included only to help navigate the document in PyCharm
A05_setup_df_file_locations = settingsdict["run_settings"]["uniprot.parse.A05_setup_df_file_locations"]
if A05_setup_df_file_locations:
    df = pd.read_csv(pathdict["dfout01_uniprotcsv"], sep=",", quoting=csv.QUOTE_NONNUMERIC, index_col=0)
    df = korbinian.common.setup_file_locations_in_df(df, settingsdict, pathdict)

A06_retrieve_simap_feature_table_and_homologues_from_list_in_csv = settingsdict["run_settings"]["simap.download.retrieve_simap_feature_table_and_homologues_from_list_in_csv"]
if A06_retrieve_simap_feature_table_and_homologues_from_list_in_csv:
    korbinian.simap.download_homologues_from_simap(pathdict, settingsdict, logging)

#A## variables are included only to help navigate the document in PyCharm
A07_parse_SIMAP_to_csv = settingsdict["run_settings"]["simap.parse.parse_simap_to_csv"]
if A07_parse_SIMAP_to_csv:
    korbinian.simap.parse_SIMAP_to_csv(pathdict, settingsdict, logging)

#A## variables are included only to help navigate the document in PyCharm
A08_calculate_AAIMON_ratios = settingsdict["run_settings"]["analyse.calculate_AAIMON_ratios"]
if A08_calculate_AAIMON_ratios:
    korbinian.simap.calculate_AAIMON_ratios(pathdict, settingsdict, logging)

#A## variables are included only to help navigate the document in PyCharm
A08a_calculate_gap_densities = settingsdict["run_settings"]["analyse.calculate_gap_densities"]

if A08a_calculate_gap_densities:
    korbinian.gap.calculate_gap_densities(pathdict, settingsdict, logging)

#A## variables are included only to help navigate the document in PyCharm
A08b_create_graph_of_gap_density = settingsdict["run_settings"]["plot.create_graph_of_gap_density"]
if A08b_create_graph_of_gap_density:
    korbinian.gap.create_graph_of_gap_density(pathdict, settingsdict, logging)

'''+++++++++++++++ Summary figures describing the conservation ratios of proteins in the list ++++++++++++++++++'''
A09_save_figures_describing_proteins_in_list = settingsdict["run_settings"]["plot.save_figures_describing_proteins_in_list"]
if A09_save_figures_describing_proteins_in_list:
    korbinian.rel_cons.save_figures_describing_proteins_in_list(pathdict, settingsdict, logging)

'''+++++++++++++++ Summary figures describing the conservation ratios of proteins in the list ++++++++++++++++++'''
A10_compare_lists = settingsdict["run_settings"]["plot.compare_lists"]
if A10_compare_lists:
    korbinian.rel_con.compare_rel_con_lists(pathdict, settingsdict, logging)

'''+++++++++++++++TMD CONSERVATION (OLD)++++++++++++++++++'''

#A## variables are included only to help navigate the document in PyCharm
B01_OLD_calculate_TMD_conservation = settingsdict["run_settings"]["simap.OLD.calculate_TMD_conservation"]
if B01_OLD_calculate_TMD_conservation:
    korbinian.old.OLD_calculate_TMD_conservation(pathdict, settingsdict, logging)

'''The gapped identity analysis gives the average TMD mem/nonmem conservation for homologues at different levels of similarity (eg. close homologues vs far homologues)
This is a funny system. The dictionary keys are the 5000 hits. This would be much more efficient and simple to re-write in pandas style.
'''
#A## variables are included only to help navigate the document in PyCharm
B02_OLD_calculate_TMD_conservation_by_gappedIdentity = settingsdict["run_settings"]["simap.OLD.calculate_TMD_conservation_by_gappedIdentity"]
if B02_OLD_calculate_TMD_conservation_by_gappedIdentity:
    korbinian.old.OLD_calculate_TMD_conservation_by_gappedIdentity(pathdict, settingsdict, logging, list_number)

#A## variables are included only to help navigate the document in PyCharm
A10_conduct_statistical_analysis_of_sim_ratios_saved_in_dfout05_simapcsv = settingsdict["run_settings"]["simap.OLD.conduct_statistical_analysis_of_sim_ratios_saved_in_dfout05_simapcsv"]
if A10_conduct_statistical_analysis_of_sim_ratios_saved_in_dfout05_simapcsv:
    korbinian.old.OLD_conduct_statistical_analysis_of_sim_ratios_saved_in_dfout05_simapcsv(pathdict, settingsdict, logging)

#A## variables are included only to help navigate the document in PyCharm
B03_OLD_fix_dfout05_simapcsv_by_adding_query_md5 = settingsdict["run_settings"]["simap.OLD.fix_dfout05_simapcsv_by_adding_query_md5"]
if B03_OLD_fix_dfout05_simapcsv_by_adding_query_md5:
    korbinian.old.OLD_fix_dfout05_simapcsv_by_adding_query_md5(pathdict, logging)