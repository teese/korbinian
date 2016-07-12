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
settings_number = 8
main_folder_harddrive = '/nas'
main_folder_list = [main_folder_harddrive, 'teeselab', 'students', 'rimma', 'omp']
main_folder = os.path.join(*main_folder_list)
json_file_with_settings = os.path.join(main_folder, 'settings', 'settings_%02d.json' % settings_number)
'''
The MAIN script starts here
'''
import json
import csv
import xml.etree.ElementTree as ET
import numpy as np
from Bio import SwissProt
from Bio import SeqIO
import logging
import tarfile
from time import strftime
import pandas as pd
import re as re
import matplotlib.pyplot as plt
import platform
import itertools
from scipy.stats import ttest_ind
import sys
import utility_rimma as r_utils
import matplotlib.patches as patches
import os
import ast 



#open the settings file
with open(json_file_with_settings, 'r') as f:
    settings = json.load(f)

list_number = settings["run_settings"]["uniprot_list"]

#the utils file contains a number of necessary scripts, this needs to be visible to the system
if settings["file_locations"]["utils"] not in sys.path:
    sys.path.append(settings["file_locations"]["utils"])
import mtutils as utils
# import mtutils, reload(mtutils) might be necessary if the session remembers the old module
# import imp
# imp.reload(utils)

''' -------Setup keyboard interrupt----------
'''
#import arcgisscripting
import signal
def ctrlc(sig, frame):
    raise KeyboardInterrupt("CTRL-C!")
signal.signal(signal.SIGINT, ctrlc)
'''+++++++++++++++LOGGING++++++++++++++++++'''
date_string_with_hour = strftime("%Y-%m-%d-%H-%M")
date_string = strftime("%Y%m%d")

#designate the output logfile
logfile = os.path.join(main_folder, 'logfiles','List%s_Settings%s_%s_logfile.log' % (list_number, settings_number, date_string))
#a file to keep a record of the log settings used for that script
utils.setup_error_logging(logfile)

#if settings['logging_settings']['suppress_error_logging_to_console']:
#    logging.setLevel('WARNING')
#'''+++++++++++++++FILENAMES++++++++++++++++++'''


# Folder, where SIMAP download is stored; requires a lot of harddrive; 
data_folder = settings["file_locations"]["data_folder"]


uniprot_flatfile_of_selected_records = os.path.join(data_folder, 'D_uniprot', 'selected',
                                                    'List%s_selected_uniprot_records_flatfile.txt' % list_number)

base_filename_summaries = os.path.join(main_folder, 'summaries', 'List%02d' % list_number)
#base_filename_input_acc_lists = '%s/input_acc_lists/List%02d'  % (main_folder, list_number)
excelfile_with_uniprot_accessions = os.path.join(base_filename_summaries, '.xlsx')
list_of_uniprot_accessions = 'List%02d_uniprot_accessions.txt' % list_number
list_of_uniprot_accessions_path = os.path.join(main_folder, 'input_acc_lists', list_of_uniprot_accessions)
dfout01_uniprotcsv = '%s_uniprot.csv' % base_filename_summaries
dfout02_uniprotTcsv = '%s_uniprotT.csv' % base_filename_summaries
dfout03_uniprotxlsx = '%s_uniprot.xlsx' % base_filename_summaries
dfout04_uniprotcsv_incl_paths = '%s_uniprot_incl_paths.csv' % base_filename_summaries
dfout05_simapcsv = '%s_simap.csv' % base_filename_summaries
dfout06_simapxlsx = '%s_simap.xlsx' % base_filename_summaries
dfout07_simapnonred = '%s_simapnonred.csv' % base_filename_summaries
dfout08_simap_AAIMON = '%s_simap_AAIMON.csv' % base_filename_summaries
dfout09_simap_AAIMON_02 = '%s_simap_AAIMON_02.csv' % base_filename_summaries
dfout10_uniprot_gaps = '%s_gap_densities.csv' % base_filename_summaries
dfout11_ = 0
dfout12_ = 0
dfout13_ = 0

csv_file_with_histogram_data = os.path.join(main_folder, 'List%02d_histogram.csv' % list_number)
csv_file_with_histogram_data_normalised = os.path.join(main_folder, 'List%02d_histogram_normalised.csv' % list_number)
csv_file_with_histogram_data_normalised_redundant_removed = os.path.join(main_folder,'List%02d_histogram_normalised_redundant_removed.csv' % list_number)
csv_file_with_md5_for_each_query_sequence = os.path.join(main_folder, 'List%02d_query_md5_checksums.csv' % list_number)

csv_av_cons_ratio_all_proteins = '%s_cons_ratios_nonred_av.csv' % base_filename_summaries
csv_std_cons_ratio_all_proteins = '%s_cons_ratios_nonred_std.csv' % base_filename_summaries

A00_convert_uniprot_list_to_nonred_ff_via_uniref = settings["run_settings"]["create_uniprot_flatfile_with_selected_seqs"][
                                               "convert_uniprot_list_to_nonred_ff_via_uniref"]["run"]
overwrite_selected_ff = settings["run_settings"]["create_uniprot_flatfile_with_selected_seqs"]["convert_uniprot_list_to_nonred_ff_via_uniref"]["overwrite_selected_ff"]
if A00_convert_uniprot_list_to_nonred_ff_via_uniref:
    logging.info('~~~~~~~~~~~~  starting A00_convert_uniprot_list_to_nonred_ff_via_uniref   ~~~~~~~~~~~~')
    if os.path.isfile(uniprot_flatfile_of_selected_records) == False or overwrite_selected_ff == True:
        # load uniref cutoff used (typicall 50, for UniRef50)
        uniref_cutoff = settings["run_settings"]["create_uniprot_flatfile_with_selected_seqs"]["convert_uniprot_list_to_nonred_ff_via_uniref"]["uniref_cluster_cutoff"]
        # load filters used when downloading original UniRef cluster from UniProt website    
        uniref_filters = settings["run_settings"]["create_uniprot_flatfile_with_selected_seqs"]["convert_uniprot_list_to_nonred_ff_via_uniref"]["uniref_filters"]
        # use cutoff and filters to select the correct uniref file from the databases folder
        uniref_csv_with_clusters = 'UniRef%s_%s.tab' % (uniref_cutoff, '_'.join(uniref_filters))
        uniref_csv_with_clusters_path = os.path.join(data_folder, 'D_uniref', uniref_csv_with_clusters)
        #create a new dataframe with the uniref csv file, containing the accessions of the reference sequences
        df_uniref = pd.read_table(uniref_csv_with_clusters_path)
        #to simplify, keep only the columns with the uniprot accessions
        #for backwards compatibility, check if the old header is used
        if 'Cluster member(s)' in df_uniref.columns:
            df_uniref.rename(columns={'Cluster member(s)': 'Cluster members'})
        #to simplify, keep only the columns with the uniprot accessions
        df_uniref = df_uniref[['Cluster ID', 'Cluster members']]# 'Organisms'
        df_uniref = df_uniref.set_index('Cluster ID')
        #change the cluster ID to the index
        df_uniref.index.names = ['cluster_ID']
        #convert the list of uniprot accessions in each cluster to a python list format
        df_uniref['all_acc_in_cluster'] = df_uniref['Cluster members'].apply(lambda x : [x.strip() for x in x.split(';')])
        nonred_folder = os.path.join(data_folder, 'D_uniprot', 'redundant')
        selected_uniprot_acc_nonred= 'List%02d_selected_uniprot_acc_nonred.tab' % list_number
        selected_uniprot_acc_nonred_path = os.path.join(nonred_folder, selected_uniprot_acc_nonred)
        if os.path.isfile(selected_uniprot_acc_nonred_path) == False:
            logging.warning('warning, file with nonredundant uniprot acc does not exist : %s' % 
                            selected_uniprot_acc_nonred_path)
        #open up large csv containing the uniprot acc of all single-pass proteins in uniprot
        df_uniprot_acc = pd.read_table(selected_uniprot_acc_nonred_path)
        #set the uniprot acc as the index
        df_uniprot_acc = df_uniprot_acc.set_index('Entry')
        #create an empty column to contain the cluster ID
        df_uniprot_acc['Cluster ID'] = ''
        #create a new dataframe with the uniref csv file, containing the accessions of the reference sequences
        #find the appropriate uniref cluster containing each acc
        acc_counter = 0
        for acc in df_uniprot_acc.index:
            #check if the accession is already a uniref representative (saves searching time!)
            if 'UniRef%i_%s' % (uniref_cutoff, acc) in df_uniref.index:
                df_uniprot_acc.loc[acc, 'Cluster ID'] = 'UniRef50_%s' % acc
            else:
                #if it is not a uniref representative, go through each uniref cluster, checking to see if it is a member
                for acc_cluster in list(df_uniref['all_acc_in_cluster']):
                    if acc in acc_cluster:
                        #if the acc is a member of a uniref cluster, add the cluster name to the original dataframe
                        df_uniprot_acc.loc[acc, 'Cluster ID'] = 'UniRef50_%s' % acc_cluster[0]
            acc_counter += 1
            if acc_counter % 100 == 0:
                logging.info('%i records checked for redundancy' % acc_counter)
        #determine which cluster IDs are in the database only once, as these uniprot entries are already nonredundant
        series_unique_bool = df_uniprot_acc['Cluster ID'].value_counts() == 1
        list_uniref_clusters_matching_only_one_acc =  series_unique_bool.loc[series_unique_bool == True].index
        #now use this list to label the original sequences that are nonredundant
        df_uniprot_acc['nonred'] = df_uniprot_acc['Cluster ID'].apply(
                                   lambda x : x in list_uniref_clusters_matching_only_one_acc)
    
        #obtain the appropriate flatfile, using the filters in the settings
        #note that the flatfile parsing is quite slow, so smaller flatfiles give faster results
        large_uniprot_ff_filters = settings["run_settings"]["create_uniprot_flatfile_with_selected_seqs"][
                                   "convert_uniprot_list_to_nonred_ff_via_uniref"]["large_uniprot_ff_filters"]
        large_uniprot_ff = 'UniProt_FF_%s.txt' % '_'.join(large_uniprot_ff_filters)
        large_uniprot_ff_path = os.path.join(data_folder, 'D_uniprot', large_uniprot_ff)
        
        #create a list of uniprot accessions that are nonredundant
        df_uniprot_acc_nonred = df_uniprot_acc.loc[df_uniprot_acc['nonred'] == True]
        list_nonred_acc = list(df_uniprot_acc_nonred.index)
        
        #create a uniprot flatfile containing only the desired nonredundant accessions
        utils.retrieve_selected_uniprot_records_from_flatfile(list_nonred_acc,
                                                              large_uniprot_ff_path,
                                                              uniprot_flatfile_of_selected_records)
        number_nonredundant_records = len(list_nonred_acc)
        number_total_records = df_uniprot_acc.shape[0]
        logging.info('convert_uniprot_list_to_nonred_ff_via_uniref is finished,\n%i from %i records were nonredundant' % 
                     (number_nonredundant_records,number_total_records))
        logging.info('(%i redundant sequences removed from dataset)' % (number_total_records - number_nonredundant_records))

#A## variables are included only to help navigate the document in PyCharm
A01_parse_large_flatfile_with_list_uniprot_accessions = settings['run_settings']["create_uniprot_flatfile_with_selected_seqs"]["parse_large_flatfile_with_list_uniprot_accessions"]
if A01_parse_large_flatfile_with_list_uniprot_accessions:
    logging.info('~~~~~~~~~~~~  starting A01_parse_large_flatfile_with_list_uniprot_accessions   ~~~~~~~~~~~~')
    #parse_large_flatfile_with_list_uniprot_accessions(list_of_uniprot_accessions, uniprot_flatfile_all_single_pass, uniprot_flatfile_of_selected_records)
    #def parse_large_flatfile_with_list_uniprot_accessions(input_accession_list, input_uniprot_flatfile, output_uniprot_flatfile):
    input_accession_list = list_of_uniprot_accessions
    input_uniprot_flatfile = uniprot_flatfile_all_single_pass
    output_uniprot_flatfile = uniprot_flatfile_of_selected_records
    #from Bio import SeqIO
    #create a list of all the uniprot accessions of the proteins to be selected from the larger uniprot file 
    accession_list = [line.strip() for line in open(input_accession_list, "r")]
    uniprot_index_handle = SeqIO.index(input_uniprot_flatfile, "swiss")
    with open(output_uniprot_flatfile, "wb") as output:
        for acc in accession_list:
            try:
                #add the selected records to the file, but adds a new line after each line! Doesn't affect later conversion to SeqRecord object
                output.write(uniprot_index_handle.get_raw(acc))
            except KeyError:
                logging.info("No SwissProt record found in %s for %s." % (input_uniprot_flatfile, acc))


#A## variables are included only to help navigate the document in PyCharm
A02_retrieve_uniprot_data_for_acc_list_in_xlsx_file = settings["run_settings"]['create_uniprot_flatfile_with_selected_seqs']["retrieve_uniprot_data_for_acc_list_in_xlsx_file"]
if A02_retrieve_uniprot_data_for_acc_list_in_xlsx_file:
    logging.info('~~~~~~~~~~~~  starting A02_retrieve_uniprot_data_for_acc_list_in_xlsx_file   ~~~~~~~~~~~~')
    #take list of acc, search in default uniprot flatfile. If missing, download from uniprot server.
    input_uniprot_flatfile = uniprot_flatfile_all_human_membrane_compressed
    df_uniprot_accessions = pd.read_excel(excelfile_with_uniprot_accessions, sheetname='uniprot_numbers')
    #remove proteins that are marked as 'not included in analysis'
    df_uniprot_accessions = df_uniprot_accessions[df_uniprot_accessions['include_in_analysis'] == True]
    #accession_list = [line.strip() for line in open(input_accession_list, "r")]
    uniprot_index_handle = SeqIO.index(input_uniprot_flatfile, "swiss")
    with open(uniprot_flatfile_of_selected_records, "wb") as output:
        for uniprot_accession in df_uniprot_accessions['A1_uniprot_accession']:
            try:
                # this adds the selected records to the file, but adds a new line after each line!
                #Doesn't affect conversion to SeqRecord object)
                assert isinstance(uniprot_index_handle, object)
                output.write(uniprot_index_handle.get_raw(uniprot_accession))
            except KeyError:
                print("No SwissProt record found in %s for %s." % (input_uniprot_flatfile, uniprot_accession))

#A## variables are included only to help navigate the document in PyCharm
A03_create_csv_from_uniprot_flatfile = settings["run_settings"]["create_csv_from_uniprot_flatfile"]
if A03_create_csv_from_uniprot_flatfile:
    #create_csv_from_uniprot_flatfile(input_file=uniprot_flatfile_of_selected_records, output_file=csv_file_with_uniprot_data)
    ## open uniprot flatfile
    #def create_csv_from_uniprot_flatfile(input_file, output_file):
    #global uniprot_record, uni_dict, record
    #input_file=uniprot_flatfile_of_selected_records
    #output_file=csv_file_with_uniprot_data
    logging.info('~~~~~~~~~~~~  starting A03_create_csv_from_uniprot_flatfile   ~~~~~~~~~~~~')
    uniprot_dict_all_proteins = {}
    with open(uniprot_flatfile_of_selected_records, "r")as f:
        records = SwissProt.parse(f)
        count_of_uniprot_records_processed = 0
        count_of_uniprot_records_added_to_csv = 0
        for record in records:
            #uni_dict = utils.create_dict_of_data_from_uniprot_record(record)
            #create an empty output dictionary to hold the uniprot data for each record
            output_dict = {}
            #extract the subcellular location detail from the (poorly organized and unsorted) uniprot comments section
            comments_dict = {}
            #utils.create_dictionary_of_comments(record, comments_dict)
            try:
                for comment in record.comments:
                    # splits comments based on first ":" symbol, creates a list called split_comment
                    split_comment = comment.strip().split(': ', 1)
                    # several comments have the same name. need to check if it is already in the dictionary
                    if split_comment[0] in comments_dict:
                        # list the different comments, one after another
                        comments_dict[split_comment[0]] += ", %s" % split_comment[1]
                    else:
                        comments_dict[split_comment[0]] = split_comment[1]
                output_dict['comments_subcellular_location_uniprot'] = comments_dict['SUBCELLULAR LOCATION']
            except (AttributeError, KeyError):
                #there are no comments in this uniprot file!
                logging.info('no comments in Uniprot file')
                output_dict['comments_subcellular_location_uniprot'] = ''

            #use regex to search for text describing subcellular locations
            #[ -]? accepts either space, hyphen, or no dividing character
            regex_word_dict = {'multipass': ['multi[ -]?(pass|span)', 'poly[ -]?topic'],
                               'singlepass': ['single[ -]?(pass|span)', 'bi[ -]?topic'],
                               'membrane': ['membran', 'lipid[ -](anchor|bound)'],
                               'typeI': ['type[ -](one|1|I)[ -]membran'],
                               'typeII': ['type[ -](two|2|II)[ -]membran']}
            #comments_subcellular_location_uniprot = 'Membrane; Bitopictype I membrane protein.'
            regex_subcell_loc_dict = {}
            for search_word in regex_word_dict:
                regex_subcell_loc_dict[search_word] = False
                regex_search_list = regex_word_dict[search_word]
                for regex_search_string in regex_search_list:
                    #search for the regex string, ignoring any mismatches in upper or lower case
                    comment_match = re.search(regex_search_string, 
                                              output_dict['comments_subcellular_location_uniprot'],
                                              re.IGNORECASE
                                              )
                    if comment_match:
                        regex_subcell_loc_dict[search_word] = True
                        #else:
                        #print('pattern not found')
                        #regex_subcell_loc_result_dict[search_word] = False
            #the dictionary could also be nested within one column
            #output_dict['regex_subcell_loc_dict'] = regex_subcell_loc_dict
            #add all of the fields to the dictionary
            output_dict.update(regex_subcell_loc_dict)

            #print accession number
            logging.info(record.accessions[0])

            #add data to dictionary
            output_dict['A1_uniprot_accession'] = record.accessions[0]
            output_dict['organism'] = record.organism
            output_dict['uniprot_entry_name'] = record.entry_name
            output_dict['uniprot_gene_name'] = record.gene_name
            output_dict['uniprot_descr'] = record.description
            output_dict['uniprot_seq'] = record.sequence
            output_dict['uniprot_orgclass'] = record.organism_classification
            output_dict['uniprot_all_accessions'] = record.accessions
            output_dict['uniprot_KW'] = record.keywords
            output_dict['uniprot_features'] = record.features
            output_dict['uniprot_seqlen'] = record.sequence_length

            #create a list of all the feature types (signal, transmem, etc)
            list_of_feature_types_in_uniprot_record = []
            for sublist in record.features:
                list_of_feature_types_in_uniprot_record.append(sublist[0])
                #logging.info(sublist)

            #list of the features that we want in the final csv
            desired_features_in_uniprot = ['TRANSMEM', 'VARIANT', 'CONFLICT', 'VAR_SEQ', 'VARSPLIC', 'TOPO_DOM']
            desired_features_in_uniprot_dict = {}

            location_of_tmds_in_feature_list = []
            location_of_non_tmds_in_feature_list = []
				
				
            for feature in desired_features_in_uniprot:
                if feature in list_of_feature_types_in_uniprot_record:
                    #find the features in the feature list.
                    #For polytopic membrane proteins, there will be more than one tmd (labelled "TRANSMEM".
                    location_of_features_in_feature_list = [i for i, x in
                                                            enumerate(list_of_feature_types_in_uniprot_record) if
                                                            x == feature]
                    desired_features_in_uniprot_dict[feature] = location_of_features_in_feature_list
                    if feature == 'TRANSMEM':
                        location_of_tmds_in_feature_list = location_of_features_in_feature_list
                        #sort list to be sure that the "transmem" notation is definitely ordered correctly,
                        #as this order determines the TMD name
                        location_of_tmds_in_feature_list.sort()
                    if feature == 'TOPO_DOM':
                        location_of_non_tmds_in_feature_list = location_of_features_in_feature_list
                        #sort list to be sure that the "transmem" notation is definitely ordered correctly,
                        #as this order determines the TMD name
                        location_of_non_tmds_in_feature_list.sort()	

            #count the number of "TRANSMEM" TMDs listed in the feature-list
            output_dict['number_of_TMDs_in_uniprot_feature_list'] = len(location_of_tmds_in_feature_list)
            
            #information about location of first non-tmd (extracellular or perplasmic/cytoplasmic)
            if len(location_of_non_tmds_in_feature_list)>0:
                output_dict['loc_start']= record.features[location_of_non_tmds_in_feature_list[0]][3]
                output_dict['n_term_ec'] = "Extracellular" in output_dict["loc_start"]
            else: 
                output_dict['loc_start']= np.nan
                output_dict['n_term_ec'] = np.nan
            		
			
            if output_dict['number_of_TMDs_in_uniprot_feature_list'] > 0:
                list_of_TMDs = []
                for TMD_location in location_of_tmds_in_feature_list:
                    #consequtively number the TMDs based on the "TRANSMEM" location in the feature list
                    TMD = 'TM%02d' % (location_of_tmds_in_feature_list.index(TMD_location) + 1)
                    list_of_TMDs.append(TMD)
                    #add the start and stop of each TMD, and the comments
                    output_dict['%s_start' % TMD] = record.features[TMD_location][1]
                    output_dict['%s_end' % TMD] = record.features[TMD_location][2]
                    output_dict['%s_description' % TMD] = record.features[TMD_location][3]
				

                #add the list of TMD names to the dictionary and dataframe
                output_dict['list_of_TMDs'] = list_of_TMDs


                #create a numpy array of any sequence variants are in the TMD region
                list_of_variant_types_in_uniprot = ['VARIANT', 'CONFLICT', 'VARSPLIC', 'VAR_SEQ']
                for TMD in list_of_TMDs:
                    #array_of_all_variants_in_tmd = np.zeros(4)
                    array_of_all_variants_in_tmd = ([])
                    for variant_type in list_of_variant_types_in_uniprot:
                        if variant_type in desired_features_in_uniprot_dict.keys():
                            #if that variant is in the uniprot data for that protein, create a list of the indices showing where that variant is found
                            list_of_variant_locations = list(desired_features_in_uniprot_dict[variant_type])
                            #get the specific start, end and details of that variant
                            for v in range(len(list_of_variant_locations)):
                                #get start
                                start_of_variant_in_seq = record.features[list_of_variant_locations[v]][1]
                                #get end
                                end_of_variant_in_seq = record.features[list_of_variant_locations[v]][2]
                                #get description
                                variant_description = record.features[list_of_variant_locations[v]][3]
                                variant_feature_identifier = record.features[list_of_variant_locations[v]][4]
                                #check if the variant is in the tmd
                                start_of_variant_is_after_start_of_tmd = True if start_of_variant_in_seq > output_dict[
                                    '%s_start' % TMD] else False
                                end_of_variant_is_before_end_of_tmd = True if end_of_variant_in_seq < output_dict[
                                    '%s_end' % TMD] else False
                                variant_is_in_tmd = True if all([start_of_variant_is_after_start_of_tmd,
                                                                 end_of_variant_is_before_end_of_tmd]) else False
                                # if the variants are the tmd region, add to numpy array
                                if variant_is_in_tmd:
                                    #create array of the variant data
                                    variant_array = np.array(
                                        [variant_type, start_of_variant_in_seq, end_of_variant_in_seq,
                                         variant_description, variant_feature_identifier])
                                    if array_of_all_variants_in_tmd != ([]):
                                        #add array with the data for this variant to the array/list for all variants
                                        array_of_all_variants_in_tmd = np.row_stack(
                                            (array_of_all_variants_in_tmd, variant_array))
                                    else:
                                        #if the array is empty, replace the array for all variants with the array for the first variant
                                        array_of_all_variants_in_tmd = variant_array
                    #if there were variants added, convert to string and add them to the output dictionary
                    if array_of_all_variants_in_tmd != ([]):
                        output_dict['%s_seq_variants' % TMD] = str(array_of_all_variants_in_tmd)

            count_of_uniprot_records_processed += 1
            #nest each dictionary containing the data for each protein into a large dictionary that contains all data from all proteins
            uniprot_dict_all_proteins[output_dict['A1_uniprot_accession']] = output_dict

        #convert that nested dict into a pandas dataframe
        df = pd.DataFrame(uniprot_dict_all_proteins).sort_index()
        #count records in dataframe
        count_of_uniprot_records_added_to_csv = len(df.columns)
        #flip rows and columns (transverse)
        df = df.T
        
        ''' ~~ DETERMINE START AND STOP INDICES FOR TMD PLUS SURROUNDING SEQ ~~ '''
        aa_before_tmd = settings["simap_match_filters"]["fasta_files"]["aa_before_tmd"]
        aa_after_tmd = settings["simap_match_filters"]["fasta_files"]["aa_after_tmd"]
        #determine max number of TMD columns that need to be created
        max_num_TMDs = df['number_of_TMDs_in_uniprot_feature_list'].max()
        #currently the loop is run for each TMD, based on the sequence with the most TMDs        
        for i in range(1, max_num_TMDs + 1):
            TMD = 'TM%02d' % i
            TMD_seq_name = '%s_seq' % TMD
            df['start_surrounding_seq_in_query_%s' % TMD] = df['%s_start' % TMD] - aa_before_tmd - 1
            df['start_surrounding_seq_in_query_%s' % TMD][df['start_surrounding_seq_in_query_%s' % TMD] < 0] = 0
            df['end_surrounding_seq_in_query_%s' % TMD] = df['%s_end' % TMD] + aa_after_tmd
            df['end_surrounding_seq_in_query_%s' % TMD][df['end_surrounding_seq_in_query_%s' % TMD] > df['uniprot_seqlen']] = df['uniprot_seqlen']
            #set up empty column for each sequence [[necessary??]]
            df['%s_with_surrounding_seq' % TMD] = ''
            
            
        ''' ~~   SLICE TMDS FROM UNIPROT SEQ    ~~ '''
        #iterate through each TMD, slicing out the relevant sequence. 
        #If there is no TMD, the cells will contain np.nan
        for i in range(1, max_num_TMDs + 1):
            TMD = 'TM%02d' % i
            #slice TMD
            df['%s_seq' % TMD] = df[df['%s_start' % TMD].notnull()].apply(utils.slice_uniprot_TMD_seq, args = (TMD,), axis=1)
            #slice TMD plus surrounding seq
            df['%s_with_surrounding_seq' % TMD] = df[df['%s_start' % TMD].notnull()].apply(utils.slice_uniprot_TMD_plus_surr_seq, args = (TMD,), axis=1)
        #save to a csv
        df.to_csv(dfout01_uniprotcsv, sep=",", quoting=csv.QUOTE_NONNUMERIC)
        #flip rows and columns, and save again as csv for an alternative view
        df.T.to_csv(dfout02_uniprotTcsv, sep=",", quoting=csv.QUOTE_NONNUMERIC)
        #save both the inverted and original data as worksheets in excel
        #first convert python datatypes to strings, as these currently give a TypeError
        df['uniprot_orgclass'] = df['uniprot_orgclass'].astype(str)
        df['organism_domain'] = df.uniprot_orgclass.apply(lambda x: x.strip("'[]").split("', '")[0])
        df['uniprot_all_accessions'] = df['uniprot_all_accessions'].astype(str)
        df['uniprot_KW'] = df['uniprot_KW'].astype(str)
        df['uniprot_features'] = df['uniprot_features'].astype(str)
        df['list_of_TMDs'] = df['list_of_TMDs'].astype(str)
        #save to Excel
        writer = pd.ExcelWriter(dfout03_uniprotxlsx)  #engine='xlsxwriter'
        df.to_excel(writer, sheet_name='dfout02')
        df.T.to_excel(writer, sheet_name='dfout01')
        writer.save()
        writer.close()

    logging.info('A03_create_csv_from_uniprot_flatfile was successful:'
                 '\n\t%i uniprot records processed\n\t%i uniprot records parsed to csv' % (
                 count_of_uniprot_records_processed, count_of_uniprot_records_added_to_csv))

#A## variables are included only to help navigate the document in PyCharm
A04_setup_df_dtypes = settings["run_settings"]["A04_setup_df_dtypes"]
if A04_setup_df_dtypes:
    logging.info('~~~~~~~~~~~~  starting A04_setup_df_dtypes   ~~~~~~~~~~~~')
    #test if the dataframe has already been created, otherwise reopen from csv file
    if 'df' in globals():    
        if isinstance(df, pd.DataFrame):
            logging.info('first protein acc = %s, df already exists, continuing with A04_setup_df_dtypes' % df.iloc[0][0])
    else:
        logging.info('df loaded from %s' % dfout01_uniprotcsv)
        df = pd.read_csv(dfout01_uniprotcsv, sep=",", quoting=csv.QUOTE_NONNUMERIC,index_col=0) 
    '''Create Pandas Dataframe with the protein name and file locations, etc'''
    original_columns = list(df.columns)
    columns_added_after_SIMAP_analysis = ['organism_domain', 'protein_name', 'SIMAP_feature_table_XML_file_path',
                                          'SIMAP_homologues_XML_file_path', 'SIMAP_csv_analysed_path',
                                          'SIMAP_input_seq_details_dict', 'SIMAP_filter_string',
                                          'SIMAP_resultSpecification_dict', 'database_details_dict',
                                          'simap_version', 'SIMAP_total_hits',
                                          'fasta_file_path',
                                          'fasta_file_plus_surr_path',
                                          'query_md5', 'query_length', 'query_selfscore', 'query_sequenceid',
                                          'total_number_of_simap_hits',
                                          'query_sequence_from_homologue_XML_file',
                                          'number_of_hits_in_homologue_XML_file',
                                          'kept_after_redundancy_check',
                                          'number_of_valid_hits',
                                          'df_mem_nonmem_ratios',
                                          'mean_ratio_ident_mem_to_nonmem',
                                          'stdev_ratio_ident_mem_to_nonmem',
                                          'protein_kept_for_statistical_analysis',
                                          'fasta_file_with_homologues_kept_for_statistical_analysis',
                                          'csv_file_av_cons_ratios_hits',
                                          'csv_SIMAP_homologues_kept_for_statistical_analysis'
                                          ]
    new_unique_column_list = set(original_columns + columns_added_after_SIMAP_analysis)
    #add extra columns
    df = df.reindex(index=df.index, columns=new_unique_column_list)
    #sort columns
    df = df.sort_index(axis=1)

    #to avoid trouble with dtypes, change all the new columns to dtype=object
    for column in columns_added_after_SIMAP_analysis:
        df[column] = df[column].astype(object)

    '''Useful debugging tool: check for duplicates'''
    import collections
    list01 = original_columns + columns_added_after_SIMAP_analysis
    duplicate_columns = [x for x, y in collections.Counter(list01).items() if y > 1]
    logging.info('%d duplicate_columns found %s' % (len(duplicate_columns), duplicate_columns))

#A## variables are included only to help navigate the document in PyCharm
A05_setup_df_file_locations = settings["run_settings"]["A05_setup_df_file_locations"]
if A05_setup_df_file_locations:
    logging.info('~~~~~~~~~~~~  starting A05_setup_df_file_locations   ~~~~~~~~~~~~')
    #test if the dataframe has already been created, otherwise reopen from csv file
    if 'df' in globals():    
        if isinstance(df, pd.DataFrame):
            logging.info('first protein acc = %s, df already exists, continuing with A05_setup_df_file_locations' % df.iloc[0][0])
    else:
        logging.info('df loaded from %s' % dfout01_uniprotcsv)
        df = pd.read_csv(dfout01_uniprotcsv, sep=",", quoting=csv.QUOTE_NONNUMERIC,index_col=0) 
    #set up a folder to hold the SIMAP BLAST-like output
    #note that at the moment, files cannot be compressed
    simap_data_folder = os.path.join(settings['file_locations']['data_folder'], 'simap')
    #join the accession and entry name to create a "protein name" for naming files
    df['A2_protein_name'] = df.A1_uniprot_accession
    df['first_two_letters_of_uniprot_acc'] = df['A1_uniprot_accession'].apply(lambda x: x[0:2])
    df['simap_filename_base_linuxpath'] = simap_data_folder + '/' + df.first_two_letters_of_uniprot_acc + '/' + df.A2_protein_name
    df['simap_filename_base'] = df['simap_filename_base_linuxpath'].apply(lambda x: os.path.normpath(x))
    df.drop('simap_filename_base_linuxpath', axis=1, inplace=True)
    
    #create filenames for simap output
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
    
    #name the fasta file with the TMD seqs (eg A0A1F4_EYS_DROME_TMD_sequences_of_homologues.fas)
    df['fasta_file_path'] = df.simap_filename_base + '_simap_TMD_seq_homologues.fas'
    
    #name the second fasta file (eg. A0T0U2_PSBE_THAPS_simap_TMD_seq_homol_&_surrounding.fas)
    df['fasta_file_plus_surr_path'] = df.simap_filename_base + '_simap_TMD_seq_homol_&_surrounding.fas'
    df[ 'fasta_file_with_homologues_kept_for_statistical_analysis'] = df.simap_filename_base + '_simap_TMD_seq_kept_stat_analysis.fas'
    df['csv_file_av_cons_ratios_hits'] = df.simap_filename_base + '_cons_ratios.csv'
    '''
    FOR multiple TMDs, create a BASE from which the TMDs can be numbered
    '''
    df['fasta_file_BASENAME'] = df.A2_protein_name + '_simap_seq_homologues_'
    df['fasta_file_BASENAMEPATH'] = df.simap_filename_base + '_simap_seq_homologues_'

    #name the second fasta file (eg. A0T0U2_PSBE_THAPS_simap_TMD_seq_homol_&_surrounding.fas)
    df['fasta_file_plus_surr_path_BASENAME'] = df.A2_protein_name + '_simap_seq_homol_&_surrounding_'
    df['fasta_file_plus_surr_path_BASENAMEPATH'] = df.simap_filename_base + '_simap_seq_homol_&_surrounding_'

    #create a basename for the output histograms
    df['AAIMON_hist_BASENAME'] = df.A2_protein_name + '_AAIMON_hist'
    df['AAIMON_hist_BASENAMEPATH'] = df.simap_filename_base + '_AAIMON_hist'
    df['AASMON_hist_BASENAME'] = df.A2_protein_name + '_AASMON_hist'
    df['AASMON_hist_BASENAMEPATH'] = df.simap_filename_base + '_AASMON_hist'

    df['csv_file_av_cons_ratios_hits_BASENAME'] = df.A2_protein_name + '_cons_ratios_'
    df['csv_file_av_cons_ratios_hits_BASENAMEPATH'] = df.simap_filename_base + '_cons_ratios_'
    
    #save to a csv
    df.to_csv(dfout04_uniprotcsv_incl_paths, sep=",", quoting=csv.QUOTE_NONNUMERIC)
    #save to Excel
    writer = pd.ExcelWriter(dfout03_uniprotxlsx)  #engine='xlsxwriter'
    df.to_excel(writer, sheet_name='dfout02')
    df.T.to_excel(writer, sheet_name='dfout01')
    writer.save()
    writer.close()


#A## variables are included only to help navigate the document in PyCharm
A06_retrieve_simap_feature_table_and_homologues_from_list_in_csv = settings["run_settings"][
    "retrieve_simap_feature_table_and_homologues_from_list_in_csv"]
#'''+++++++++++++++SIMAP++++++++++++++++++'''
if A06_retrieve_simap_feature_table_and_homologues_from_list_in_csv:
    logging.info('~~~~~~~~~~~~  starting A06_retrieve_simap_feature_table_and_homologues_from_list_in_csv   ~~~~~~~~~~~~')
    #test if the dataframe has already been created, otherwise reopen from csv file
    if 'df' in globals():    
        if isinstance(df, pd.DataFrame):
            logging.info('first protein acc = %s, df already exists' % df.iloc[0][0])
    else:
        logging.info('df loaded from %s' % dfout04_uniprotcsv_incl_paths)
        df = pd.read_csv(dfout04_uniprotcsv_incl_paths, sep=",", quoting=csv.QUOTE_NONNUMERIC,index_col=0)
        
    #def retrieve_simap_feature_table_and_homologues_from_list_in_csv(input_file, list_of_keys, settings):
    '''
    First prepare the csv file from the uniprot record.
    Run this to save the files based on their domain.
    '''
    #global list_of_files_with_feature_tables, list_of_files_with_homologues
    #The SIMAP download settings can be altered as desired, using the json settings file
    max_hits = settings["simap_homologue_download"]["max_hits"]
    max_memory_allocation = settings["simap_homologue_download"]["java_max_RAM_memory_allocated_to_simap_download"]
    database = settings["simap_homologue_download"]["database"]
    taxid = settings["simap_homologue_download"]["taxid"]  # eg.'7227' for Drosophila melanogaster
    extra_search_string = ''

    enough_hard_drive_space = True
    try:
        byteformat = "GB"
        data_harddrive = settings["file_locations"]["data_harddrive"]
        size = utils.get_free_space(data_harddrive, byteformat)
        logging.info('Hard disk remaining space =')
        logging.info(size)
        if size[0] < 5:
            raise utils.HardDriveSpaceException(
                "Hard drive space limit reached, there is only %s %s space left." % (size[0], size[1]))
            enough_hard_drive_space = False
    except utils.HardDriveSpaceException as e:
        logging.warning(e)
    if enough_hard_drive_space:
        #iterate over each uniprot record contained in the dataframe. note that acc = uniprot accession number
        number_of_files_not_found = 0
        for acc in df.index:
            protein_name = df.loc[acc, 'A2_protein_name']
            query_sequence_length = df.loc[acc, 'uniprot_seqlen']
            input_sequence = df.loc[acc, 'uniprot_seq']
            ''' windows has a character limit in the command prompt in theory of 8191 characters,
            but the command line java command seems to cause errors with sequences above 3000 amino acids.
            Assume that this problem only applies to Windows, 
            therefore in Windows systems limit the java string to proteins less than 3000 amino acids.
            The character limit can be adjusted in the settings file
            '''
            if 'Windows' in str(platform.system()):
                if query_sequence_length < settings["simap_homologue_download"]["max_query_sequence_length"]:
                    download_homologues = False
                    if query_sequence_length > settings["simap_homologue_download"]["max_query_sequence_length"]:
                        logging.warning('%s cannot be processed into a java command in windows OS,'
                                        'as the sequence is longer than %i characters (%i). Moving to next sequence' % (
                                        protein_name, settings["simap_homologue_download"]["max_query_sequence_length"],
                                        query_sequence_length))
                else:
                    download_homologues = True
            else:
                download_homologues = True
            if download_homologues == True:
                simap_data_folder = os.path.join(settings['file_locations']['data_folder'], 'simap')
                subfolder = simap_data_folder+ "/" + df.loc[acc, 'first_two_letters_of_uniprot_acc']
                if os.path.isdir(subfolder) == False:
                    os.mkdir(subfolder)
                #check which files exist. This is useful, because it is not possible to open the tarfile as 'a:gz', 
                #therefore you cannot add files to an existing tarfile)
                feature_table_XML_exists, homologues_XML_exists, SIMAP_tarfile_exists = utils.check_tarfile(df, acc)
                if not SIMAP_tarfile_exists:
                    if not feature_table_XML_exists:
                        #download feature table from SIMAP
                        utils.retrieve_simap_feature_table(input_sequence,
                                                           output_file=df.loc[acc, 'SIMAP_feature_table_XML_file_path'])
                    if not homologues_XML_exists:
                        #download homologue file from SIMAP
                        utils.retrieve_simap_homologues(input_sequence,
                                                        output_file=df.loc[acc, 'SIMAP_homologues_XML_file_path'],
                                                        database=database, max_hits=max_hits,
                                                        max_memory_allocation=max_memory_allocation, taxid=taxid,
                                                        extra_search_string=extra_search_string)
                        #now check again if the files exist
                    feature_table_XML_exists, homologues_XML_exists, SIMAP_tarfile_exists = utils.check_tarfile(df, acc)
                    if not homologues_XML_exists:
                        #add one to the list of consecutive failed downloads.
                        number_of_files_not_found += 1
                        #if a large number of downloads failed, then the SIMAP server is probably not working. 
                        #Wait some time and try again later.
                        if number_of_files_not_found > 30:
                            utils.sleep_24_hours()
                        if number_of_files_not_found == 20:
                            utils.sleep_24_hours()
                        if number_of_files_not_found == 15:
                            utils.sleep_6_hours()
                        if number_of_files_not_found == 10:
                            utils.sleep_6_hours()
                        utils.sleep_120_seconds()
                    else:
                        #if download is successful or file exists, the SIMAP server must be working, 
                    #therefore reset the number_of_files_not_found
                        number_of_files_not_found = 0

                    #since we can't add files to the compressed tarfile, only when both the feature table 
                    #and xml file are downloaded should we pack and compress them
                    if feature_table_XML_exists and homologues_XML_exists:
                        with tarfile.open(df.loc[acc, 'SIMAP_tarfile'], mode='w:gz') as tar:
                            #add the files to the compressed tarfile
                            logging.info(
                                '%s XML files will be moved into the tarball, original XML files deleted' % protein_name)
                            tar.add(df.loc[acc, 'SIMAP_feature_table_XML_file_path'],
                                    arcname=df.loc[acc, 'SIMAP_feature_table_XML_file'])
                            tar.add(df.loc[acc, 'SIMAP_homologues_XML_file_path'],
                                    arcname=df.loc[acc, 'SIMAP_homologues_XML_file'])
                        #delete the original files
                        try:
                            os.remove(df.loc[acc, 'SIMAP_feature_table_XML_file_path'])
                            os.remove(df.loc[acc, 'SIMAP_homologues_XML_file_path'])
                        except FileNotFoundError:
                            pass

                            #now add the downloaded files to the tarball, and delete the original XML files

                        #            directory = r'E:/Stephis/Projects/Programming/Python/files/learning/tarfile'
                        #            newtar = directory + r'/newtar.tar.gz'
                        #            SIMAP_feature_table_XML_file_basename = 'P29274_AA2AR_HUMAN_feature_table.xml'
                        #            SIMAP_feature_table_XML_file = '%s/%s' % (directory, SIMAP_feature_table_XML_file_basename)
                        #            SIMAP_homologues_XML_file_basename = 'P29274_AA2AR_HUMAN_homologues.xml'
                        #            SIMAP_homologues_XML_file = '%s/%s' % (directory, SIMAP_homologues_XML_file_basename)
                        #
                        #            with tarfile.open(df.loc[acc, 'SIMAP_tarfile'], 'w:gz') as tar:
                        #                tar.add(SIMAP_feature_table_XML_file, arcname = SIMAP_feature_table_XML_file_basename)
                        #                tar.add(SIMAP_homologues_XML_file, arcname = SIMAP_homologues_XML_file_basename)
                        #            with tarfile.open(newtar, 'r:gz') as tar:
                        #                for tarinfo in tar:
                        #                    print(tarinfo.isreg())
                        #                    print(tarinfo.name)
                        #                    print(tarinfo.size)
                        #
                        #            df['SIMAP_tarfile'] = df.simap_filename_base + '_SIMAP.tar.gz'
                        #            df['SIMAP_feature_table_XML_file'] = df.A2_protein_name + '_feature_table.xml'
                        #            df['SIMAP_feature_table_XML_file_path'] = df.simap_filename_base + '_feature_table.xml'
                        #            df['SIMAP_homologues_XML_file'] = df.A2_protein_name + '_homologues.xml'
                        #            df['SIMAP_homologues_XML_file_path'] = df.A2_protein_name + '_homologues.xml'
                        #            df['SIMAP_csv_from_XML'] = df.simap_filename_base + '_homologues.csv'
                        #            df['SIMAP_csv_from_XML_path'] = df.simap_filename_base + '_homologues.csv'
                        #            df['csv_SIMAP_homologues_kept_for_statistical_analysis'] = df.simap_filename_base + '_homologues_for_stat.csv'

    logging.info('retrieve_simap_feature_table_and_homologues_from_list_in_csv is finished')

#A## variables are included only to help navigate the document in PyCharm
A07_parse_SIMAP_to_csv = settings["run_settings"]["parse_SIMAP_to_csv"]
if A07_parse_SIMAP_to_csv:
    counter_XML_to_CSV = 0
    logging.info('~~~~~~~~~~~~  starting parse_SIMAP_to_csv   ~~~~~~~~~~~~')
    #test if the dataframe has already been created, otherwise reopen from csv file
    if 'df' in globals():    
        if isinstance(df, pd.DataFrame):
            logging.info('first protein acc = %s, df already exists' % df.iloc[0][0])
    else:
        logging.info('df loaded from %s' % dfout04_uniprotcsv_incl_paths)
        df = pd.read_csv(dfout04_uniprotcsv_incl_paths, sep=",", quoting=csv.QUOTE_NONNUMERIC, index_col=0)

    #filter to remove sequences where no TMDs are found
    df = df.loc[df['list_of_TMDs'].notnull()]
    #filter to remove sequences where no TMDs are found (if string)
    df = df.loc[df['list_of_TMDs'] != 'nan']

    #iterate over the dataframe.
    for acc in df.index:
        #set up counters
        number_of_hits_missing_protein_node = 0
        num_hits_with_SW_align_node = 0
        number_of_hits_missing_smithWatermanAlignment_node = 0
        #number_of_hits_kept_for_statistical_analysis = 0  # number_of_hits
        organism_domain = df.loc[acc, 'organism_domain']
        protein_name = df.loc[acc, 'A2_protein_name']
        #try:
        logging.info('%s' % protein_name)
        #check which files exist
        feature_table_XML_exists, homologues_XML_exists, SIMAP_tarfile_exists = utils.check_tarfile(df, acc)
        #check if the feature table and homologue files actually exist
        if os.path.isfile(df.loc[acc, 'SIMAP_tarfile']):
            SIMAP_tarfile_exists = True
        else:
            SIMAP_tarfile_exists = False
            #at the moment, we'll only create the tarfile if both the feature table and
        #homologue XML files downloaded successfully, but this might change depending on preference
        if SIMAP_tarfile_exists:
            try:
                with tarfile.open(df.loc[acc, 'SIMAP_tarfile'], mode='r:gz') as tar:
                    if df.loc[acc, 'SIMAP_feature_table_XML_file'] in [tarinfo.name for tarinfo in tar]:
                        feature_table_in_tarfile = True
                    #else:
                    #    feature_table_in_tarfile = False
                    if df.loc[acc, 'SIMAP_homologues_XML_file'] in [tarinfo.name for tarinfo in tar]:
                        homologues_XML_in_tarfile = True
                        #else:
                        #    homologues_XML_in_tarfile = False
            except (EOFError, tarfile.ReadError):
                #file may be corrupted, if script stopped unexpectedly before compression was finished
                logging.info('%s seems to be corrupted.' 
                             'File will be deleted, and will need to be re-downloaded next time program is run' %
                             df.loc[acc, 'SIMAP_tarfile'])
                SIMAP_tarfile_exists = False
                feature_table_in_tarfile = False
                homologues_XML_in_tarfile = False
                os.remove(df.loc[acc, 'SIMAP_tarfile'])
        if not SIMAP_tarfile_exists:
            feature_table_in_tarfile = False
            homologues_XML_in_tarfile = False

        if all([feature_table_in_tarfile, homologues_XML_in_tarfile]):
            '''get the Phobius and TMHMM predictions from the feature table of the query sequence 
            NOT USED, PHOBIUS PRED OFTEN MISSING, in the future the TMD region taken from uniprot record
            '''
            #phobius_TMD_start, phobius_TMD_end, phobius_TMD_length = get_phobius_TMD_region(simap_feature_table_root)
            #TMHMM_TMD_start, TMHMM_TMD_end = get_TMHMM_TMD_region(simap_feature_table_root)
            #df.loc[acc,'phobius_TMD_start'] = phobius_TMD_start
            #df.loc[acc,'phobius_TMD_end'] = phobius_TMD_end
            #df.loc[acc,'phobius_TMD_length'] = phobius_TMD_length
            #df.loc[acc,'TMHMM_TMD_start'] = TMHMM_TMD_start
            #df.loc[acc,'TMHMM_TMD_end'] = TMHMM_TMD_end

            #create a new file to store all of the simap data, and write the csv header
            #SIMAP_csv_from_XML_path = r"E:/Databases/simap/%s/%s_homologues.csv" % (organism_domain, protein_name[:30])

            #if the setting is "False", and you don't want to overwrite the files, skip this section
            if settings["run_settings"]["calculate_AAIMON_ratios"]["overwrite_homologue_csv_files"]:
                create_homol_csv = True
            else:
                #check if output file already exists
                if os.path.isfile(df.loc[acc, 'SIMAP_csv_from_XML_tarfile']):
                    try:
                        with tarfile.open(df.loc[acc, 'SIMAP_csv_from_XML_tarfile'], 'r:gz') as tar:
                            #create a list of files
                            files_in_output_tarball = [t.name for t in tar]
                            #check that the analysed files are actually there
                            if df.loc[acc, 'SIMAP_csv_from_XML'] in files_in_output_tarball: 
                                #read output csv in tarfile
                                dfs = pd.read_csv(tar.extractfile(df.loc[acc, 'SIMAP_csv_from_XML']), index_col = 0)
                                description_of_first_hit = dfs.loc[1, 'A4_description']
                                logging.info('%s homologues already converted to csv. (%s)' % (acc, description_of_first_hit))
                                #filter to include only desired hits
                                '''OLD STUFF, from when XML to CSV was not saved separately
                                dfs_filt = dfs.query(
                                    'gapped_ident_above_cutoff == True and hit_contains_SW_node == True and disallowed_words_not_in_descr == True')
                                #avoid a divide by zero error in the unlikely case that there are no_identical_residues_in_alignment
                                dfs_filt_AAIMON = dfs_filt.loc[dfs_filt['nonTMD_perc_ident'] != 0]
                                list_of_TMDs = eval(df.loc[acc, 'list_of_TMDs'])
                                for TMD in list_of_TMDs:
                                    #following the general filters, filter to only analyse sequences with TMD identity above cutoff, 
                                    #and a nonTMD_perc_ident above zero ,to avoid a divide by zero error
                                    dfs_filt_AAIMON = dfs_filt_AAIMON.loc[
                                        dfs['%s_perc_ident' % TMD] >= settings['simap_match_filters'][
                                            'min_identity_of_TMD_initial_filter']]
                                    #add to original dataframe with the list of sequences
                                    df.loc[acc, '%s_AAIMON_ratio_mean' % TMD] = dfs_filt_AAIMON[
                                        '%s_AAIMON_ratio' % TMD].mean()
                                    df.loc[acc, '%s_AAIMON_ratio_std' % TMD] = dfs_filt_AAIMON[
                                        '%s_AAIMON_ratio' % TMD].std()
                                    logging.info('AIMON MEAN %s: %0.2f' % (TMD, df.loc[acc, '%s_AAIMON_ratio_mean' % TMD]))
                                '''
                            else:
                                logging.info('%s not in output file tarball, tarball will be deleted' % df.loc[acc, 'SIMAP_csv_from_XML'])
                                os.remove(df.loc[acc, 'SIMAP_csv_from_XML_tarfile'])
                                create_homol_csv = True
                        logging.info('%s already converted to csv, moving to next sequence' % 
                                    df.loc[acc, 'SIMAP_csv_from_XML'])
                        create_homol_csv = False
                    except (EOFError, KeyError, tarfile.ReadError):
                        #file may be corrupted, if script stopped unexpectedly before compression was finished
                        logging.info(
                            '%s seems to be corrupted. File will be deleted.' % df.loc[acc, 'SIMAP_csv_from_XML_tarfile'])
                        os.remove(df.loc[acc, 'SIMAP_csv_from_XML_tarfile'])
                        create_homol_csv = True
                else:
                    logging.info('%s not found, create_homol_csv = True' % df.loc[acc, 'SIMAP_csv_from_XML_tarfile'])
                    create_homol_csv = True
            #if the files don't exist, or you want to overwrite them
            if create_homol_csv:
                #extract the tarfile so that it can be read as xml
                with tarfile.open(df.loc[acc, 'SIMAP_tarfile'], 'r:gz') as tar:
                    SIMAP_homologues_XML_file_extracted = tar.extractfile(df.loc[acc, 'SIMAP_homologues_XML_file'])

                    #parse the XML file with elementtree, define the 'root' of the XML file
                    simap_homologue_tree = ET.parse(SIMAP_homologues_XML_file_extracted)
                    simap_homologue_root = simap_homologue_tree.getroot()

                    #print the simap search details (database, e-value cutoff etc)
                    #dict_with_query_data = print_query_details_from_homologue_XML(simap_homologue_root, dict_with_query_data)
                    for parameters in simap_homologue_root[0][0][0][0].iter('parameters'):
                        df.loc[acc, 'SIMAP_input_seq_details_dict'] = str(parameters[0][0].attrib)
                        for SIMAP_filter in parameters.iter('filter'):
                            SIMAP_filter_string = SIMAP_filter.text
                        df.loc[acc, 'SIMAP_filter_string'] = str(SIMAP_filter_string)
                        for resultSpecification in parameters.iter('resultSpecification'):
                            SIMAP_resultSpecification_dict = resultSpecification.attrib
                        df.loc[acc, 'SIMAP_resultSpecification_dict'] = '"%s"' % SIMAP_resultSpecification_dict
                        for databases in parameters.iter('databases'):
                            database_details_dict = databases[0].attrib
                        df.loc[acc, 'database_details_dict'] = '"%s"' % database_details_dict
                        df.loc[acc, 'simap_version'] = simap_homologue_root[0][0][0][0][0].attrib['version']
                        df.loc[acc, 'SIMAP_total_hits'] = int(simap_homologue_root[0][0][0][1][0].attrib['total'])
                    if df.loc[acc, 'simap_version'] != '4.0':
                        logging.warning('WARNING! Your XML file is simap version %s,'
                                        'however this SIMAP parser was developed for SIMAP version 4.0.' % 
                                         df.loc[acc, 'simap_version'])
                    counter_XML_to_CSV += 1
                    logging.info('%s homologous sequences parsed from SIMAP XML to csv' % 
                                 int(df.loc[acc, 'SIMAP_total_hits']))

                    query_sequence_node = simap_homologue_root[0][0][0][0][2][0][0]
                    ''' xxxx CURRENTLY THE df is filled with nan values, 
                        but that doesn't make sense as the script seems to work
                    '''
                    df.loc[acc, 'query_md5'] = query_sequence_node.attrib['md5']
                    df.loc[acc, 'query_length'] = int(query_sequence_node.attrib['length'])
                    df.loc[acc, 'query_selfscore'] = query_sequence_node.attrib['selfscore']
                    df.loc[acc, 'query_sequenceid'] = query_sequence_node.attrib['sequenceid']
                    df.loc[acc, 'total_number_of_simap_hits'] = query_sequence_node[0].attrib['number_hits']
                    df.loc[acc, 'query_sequence_from_homologue_XML_file'] = query_sequence_node[0][0].text
                    df.loc[acc, 'number_of_hits_in_homologue_XML_file'] = int(
                        simap_homologue_root[0][0][0][1][0].attrib['total'])
                    '''
                    Create an updated csv_file_with_uniprot_data to include the data from SIMAP regarding the query
                    '''
                    #save current df to csv
                    with open(dfout05_simapcsv, 'w') as f:
                        df.to_csv(f, sep=",", quoting=csv.QUOTE_NONNUMERIC)

                    #create new files to store the fasta sequences, and save the query sequence (the row here is the uniprot number in th df index)
                    for row in df.index:
                        #for a protein with TMDs, the list of TMD names should be ['TM01','TM02']
                        list_of_TMDs = eval(df.loc[row, 'list_of_TMDs'])
                        
                    #for each hit, save all the relevant data in the form of a dictionary,
                    # so it can be added to a csv file or used in other calculations
                    simap_homologue_hits = simap_homologue_root[0][0][0][1][0]

                    #see if there are any hits at all
                    try:
                        test2 = simap_homologue_root[0][0][0][1][0][0]
                        xml_contains_simap_homologue_hits = True
                    except IndexError:
                        xml_contains_simap_homologue_hits = False

                    if xml_contains_simap_homologue_hits:
                        #load the amino acid substitution matrices from the settings file
                        list_of_aa_sub_matrices = settings['aa_substitution_scoring']['matrices']

                        #import the amino acid substitution matrices
                        utils.import_amino_acid_substitution_matrices()

                        #add the similarity ratios to the csv_header_for_SIMAP_homologue_file.
                        # These will depend on the individual settings
                        #                    if settings['run_settings']['calculate_AAIMON_ratios']['calculate_TMD_conservation_with_aa_matrices']:
                        #                        for j in range(settings["aa_substitution_scoring"]["gap_open_penalty_min"],
                        #                                       settings["aa_substitution_scoring"]["gap_open_penalty_max"],
                        #                                       settings["aa_substitution_scoring"]["gap_open_penalty_increment"]):
                        #                            gap_open_penalty = j
                        #                            gap_extension_penalty = j
                        #print('gap_open_penalty:%s' % j)
                        #                            for matrix_name in list_of_aa_sub_matrices:
                        #                                column_name = 'sim_ratio_%s_gapo%i' % (matrix_name.replace("'", "")[0:-7], j)
                        #                                csv_header_for_SIMAP_homologue_file.append(column_name)

                        #import the necessary matrices
                        #for matrix_name in list_of_aa_sub_matrices:
                        #matrix = matrix_name[0:-7]
                        #print (matrix)
                        #from Bio.SubsMat.MatrixInfo import matrix as matrix_name

                        #print('SIMAP_csv_from_XML_path %s' % df.loc[acc, 'SIMAP_csv_from_XML_path'])

                        SIMAP_csv_from_XML_path = df.loc[acc, 'SIMAP_csv_from_XML_path']
                        #fasta_file_path = df.loc[acc, 'fasta_file_path']

                        #create an empty file
                        open(SIMAP_csv_from_XML_path, 'w').close()

                        #reopen to add match details iteratively from dictionary
                        with open(SIMAP_csv_from_XML_path, 'a') as csvfile:

                            #set up a bool to catch those files where not a single hit actually gives data
                            at_least_one_hit_contains_SW_node = False

                            for hit in simap_homologue_hits:
                                match_details_dict = {}

                                #add desired hit information to the dictionary for transfer to csv
                                A1_hit_number = int(hit.attrib['number'])
                                match_details_dict['A1_hit_number'] = A1_hit_number
                                match_details_dict['A3_md5'] = hit[1].attrib['md5']

                                #define the major nodes in the XML-file
                                try:
                                    protein_node = hit[1][1]
                                    hit_contains_protein_node = True
                                except IndexError:
                                    hit_contains_protein_node = False
                                    number_of_hits_missing_protein_node += 1
                                    logging.warning('%s hit %s contains no protein node' % (protein_name, 
                                                                                            match_details_dict['A3_md5'])
                                                                                            )
                                if hit_contains_protein_node:
                                    try:
                                        smithWatermanAlignment_node = hit[0][0][14]
                                        hit_contains_SW_node = True
                                        num_hits_with_SW_align_node += 1
                                    except IndexError:
                                        hit_contains_SW_node = False
                                    match_details_dict['hit_contains_SW_node'] = hit_contains_SW_node
                                    #add the description. Make an empty string if it is the first (query) hit, preventing the presence of np.nan in the later dataframe
                                    if A1_hit_number == 1:
                                        A4_description = '%s_SIMAP_query_sequence' % protein_name
                                    else:
                                        A4_description = protein_node.attrib['description']
                                    match_details_dict['A4_description'] = A4_description
                                    try:
                                        databaseId = int(protein_node[1].attrib['databaseId'])
                                        match_details_dict['databaseId'] = int(protein_node[1].attrib['databaseId'])
                                    except KeyError:
                                        databaseId = 0
                                        #match_details_dict['databaseId'] = int(0)
                                    #databaseId = int(protein_node[1].attrib['databaseId'])
                                    databasenode = protein_node[1]
                                    match_details_dict['database'] = databasenode.attrib['name']
                                    try:
                                        taxonomyNode = protein_node[2]
                                        match_details_dict['A2_organism'] = taxonomyNode.attrib['name']
                                        match_details_dict['taxonomy_node_id'] = taxonomyNode.attrib['node_id']
                                        match_details_dict['taxonomy_rank'] = taxonomyNode.attrib['rank']
                                    except IndexError:
                                        #sequence is probably synthetic, as it has no database node
                                        match_details_dict['A4_description'] += ', synthetic'
                                        match_details_dict['A2_organism'] = 'synthetic'
                                        match_details_dict['taxonomy_node_id'] = 'synthetic'
                                        match_details_dict['taxonomy_rank'] = 'synthetic'
                                    match_details_dict['len_full_match_seq'] = len(hit[1][0][0].text)
                                    #len_full_match_seq = len(full_match_seq)
                                    alignment_node = hit[0][0]
                                    #E-value for hit
                                    match_details_dict['FASTA_expectation'] = float(alignment_node[1].text)
                                    #convert identity from e.g. 80 (80%) to 0.8
                                    match_details_dict['FASTA_identity'] = float(alignment_node[3].text) / 100
                                    #strangely, I think gappedIdentity is the identity EXCLUDING gaps, which is a better value to base judgements on. convert identity from e.g. 80 (80%) to 0.8
                                    match_details_dict['FASTA_gapped_identity'] = float(alignment_node[4].text) / 100
                                    '''xxx notes on the gapped identity
                                    N.B The FASTA_gapped_identity data here is from the FASTA algorithm, that precedes the SW algorithm. 
                                    Occasionally they don’t match!!!
                                    I calculate the TMD identity manually from the SW alignment, BUT 
                                    currently for the calculation of membranous/nonmembranous I use the gappedIdentity from the FASTA output 
                                    (the SW output inly has identity including gaps)
                                    -    if I counted the gaps from the SW alignment, I COULD recalculate the gappedIdentity for the SW alignment
                                    -    OR: I could simply remove the data where the FASTA and SW don’t match.
                                    '''
                                    #FASTA overlap should be the length of the aligned region after running the FASTA algorithm (alignment is not shown by SIMAP)
                                    match_details_dict['FASTA_overlap'] = int(alignment_node[5].text)
                                    match_details_dict['FASTA_query_coverage'] = float(alignment_node[11].text)
                                    match_details_dict['FASTA_match_coverage'] = float(alignment_node[12].text)
                                    #find the start and the stop of the hsp
                                    querySeq = alignment_node[6]
                                    match_details_dict['FASTA_query_start'] = int(querySeq.attrib['start'])
                                    match_details_dict['FASTA_query_end'] = int(querySeq.attrib['end'])
                                    matchSeq = alignment_node[7]
                                    match_details_dict['FASTA_match_start'] = int(matchSeq.attrib['start'])
                                    match_details_dict['FASTA_match_end'] = int(matchSeq.attrib['end'])
                                    #some parameters that are needed for identity calculations later
                                    #FASTA_num_ident_res = FASTA_identity / 100.0 * FASTA_overlap
                                    #check if the TMD is in the smith waterman alignment. Note that start and stop in the alignment node is based on FASTA alignment, 
                                    #which is not shown. Occasionally, this will not match the SW alignment.
                                    #xxx it might be better to insert a function that determines of the TMD is after the start of the SW alignment
                                    #is_start_of_TMD_in_FASTA = True if FASTA_query_start <= TMDstart else False
                                    #is_end_of_TMD_in_FASTA = True if TMDend <= FASTA_query_end else False
                                    #is_TMD_in_FASTA_alignment = True if all(
                                    #    [is_start_of_TMD_in_FASTA, is_end_of_TMD_in_FASTA]) else False
                                    #mmmmm TEMP
                                    #logging.info('gapped_identity_too_low: %s' % gapped_identity_too_low)
                                    '''***********************if the TMD region is actually covered by the hsp, then conduct some further analyses of the match TMD region*************************'''
                                    if hit_contains_SW_node:
                                        #check that at least one hit gives data
                                        at_least_one_hit_contains_SW_node = True
                                        query_alignment_sequence = ''
                                        '''For the moment, there is no need to put the whole match hsp sequence into the csv file'''
                                        #for smithWatermanAlignment in alignment_node.iter('smithWatermanAlignment'):
                                        match_details_dict['SW_query_score_ratio'] = smithWatermanAlignment_node[0].text
                                        match_details_dict['SW_match_score_ratio'] = smithWatermanAlignment_node[1].text
                                        match_details_dict['SW_query_coverage'] = smithWatermanAlignment_node[2].text
                                        match_details_dict['SW_match_coverage'] = smithWatermanAlignment_node[3].text
                                        match_details_dict['SW_coverage_ratio'] = smithWatermanAlignment_node[4].text
                                        match_details_dict['A5_alignment_pretty'] = smithWatermanAlignment_node[8].text
                                        match_details_dict['SW_alignment_seq1offset'] = int(smithWatermanAlignment_node.attrib['alignment-seq1offset'])
                                        match_details_dict['SW_alignment_seq2offset'] = int(smithWatermanAlignment_node.attrib['alignment-seq2offset'])
                                        match_details_dict['SW_identity'] = float(smithWatermanAlignment_node.attrib['identity'])
                                        match_details_dict['SW_similarity'] = float(smithWatermanAlignment_node.attrib['similarity'])
                                        #Get the full sequences. Note that they greatly increase the size of the csv file.
                                        match_details_dict['query_alignment_sequence'] = smithWatermanAlignment_node[5].text
                                        match_details_dict['alignment_markup'] = smithWatermanAlignment_node[6].text
                                        match_details_dict['match_alignment_sequence'] = smithWatermanAlignment_node[7].text
                                        #create a list of TMD names to be used in the loops below (TM01, TM02 etc)
                                        list_of_TMDs = eval(df.loc[acc, 'list_of_TMDs'])
                                        #run the search using the regular expression that will find the TMD even if it contains gaps
                                        for TMD in list_of_TMDs:
                                            #if is_TMD_in_FASTA_alignment:
                                            query_TMD_sequence = df.loc[acc, '%s_seq' % TMD]
                                            query_TMD_length = len(query_TMD_sequence)
                                    else:
                                        number_of_hits_missing_smithWatermanAlignment_node += 1
                                    if A1_hit_number == 1:
                                        #sort
                                        csv_header_for_SIMAP_homologue_file = sorted(list(match_details_dict.keys()))
                                        #save the csv header to the csv file
                                        writer = csv.writer(csvfile, delimiter=',', quotechar='"', lineterminator='\n',
                                                            quoting=csv.QUOTE_NONNUMERIC, doublequote=True)
                                        writer.writerow(csv_header_for_SIMAP_homologue_file)
                                    #save the match_details_dict as a line in the csv file
                                    writer = csv.DictWriter(csvfile, fieldnames=csv_header_for_SIMAP_homologue_file,
                                                            extrasaction='ignore', delimiter=',', quotechar='"',
                                                            lineterminator='\n', quoting=csv.QUOTE_NONNUMERIC,
                                                            doublequote=True)
                                    writer.writerow(match_details_dict)

                        with tarfile.open(df.loc[acc, 'SIMAP_csv_from_XML_tarfile'], 'w:gz') as tar_SIMAP_out:
                            tar_SIMAP_out.add(SIMAP_csv_from_XML_path, arcname=df.loc[acc, 'SIMAP_csv_from_XML'])
                        os.remove(SIMAP_csv_from_XML_path)
    logging.info(
        'number_of_hits_missing_smithWatermanAlignment_node: %i' % number_of_hits_missing_smithWatermanAlignment_node)
    logging.info('number_of_hits_missing_protein_node: %i' % number_of_hits_missing_protein_node)
    logging.info('****parse_SIMAP_to_csv finished!!****\n%g files parsed from SIMAP XML to csv' % counter_XML_to_CSV)

#A## variables are included only to help navigate the document in PyCharm
A08_calculate_AAIMON_ratios = settings["run_settings"]["calculate_AAIMON_ratios"]["run"]
if A08_calculate_AAIMON_ratios:
    logging.info('~~~~~~~~~~~~starting calculate_AAIMON_ratios~~~~~~~~~~~~')
    overwrite_prev_calculated_AAIMON_ratios = settings["run_settings"]["calculate_AAIMON_ratios"][
                                              "overwrite_prev_calculated_AAIMON_ratios"]
    #test if the dataframe has already been created, otherwise re-open from csv file containing the simap data
    try:
        logging.info('first protein acc = %s, df already exists,'
                     'continuing with parse_list_proteins_to_csv_and_fasta' % df.iloc[0][0])
    except NameError:
        if os.path.isfile(dfout08_simap_AAIMON):
            #backup_original_file
            dfout08_simap_AAIMON_backup_before_adding_data = dfout08_simap_AAIMON[:-4] + 'backup.csv'
            import shutil
            #copy file, keeping original modification and access info
            shutil.copy2(dfout08_simap_AAIMON, dfout08_simap_AAIMON_backup_before_adding_data)
            #open file as dataframe
            df = pd.read_csv(dfout08_simap_AAIMON, sep=",", quoting=csv.QUOTE_NONNUMERIC, index_col=0)
            logging.info('df loaded from %s' % dfout08_simap_AAIMON)
        else:
            df = pd.read_csv(dfout05_simapcsv, sep=",", quoting=csv.QUOTE_NONNUMERIC, index_col=0)
            logging.info('df loaded from %s' % dfout05_simapcsv)
    #filter to remove sequences where no TMDs are found, 
    df = df.loc[df['list_of_TMDs'].notnull()]
    #filter to remove sequences where no TMDs are found (if string)
    df.loc[df['list_of_TMDs'] != 'nan']
    #filter to remove sequences where no TMDs are found (if string)
    df = df.loc[df['list_of_TMDs'] != 'nan']
    #determine if the dataframe already contains some previously analysed data
    dataframe_contains_prev_calc_AAIMON_ratios = True if 'TM01_AAIMON_ratio_mean' in df.columns else False
    logging.info('dataframe_contains_prev_calc_AAIMON_ratios = %s' % dataframe_contains_prev_calc_AAIMON_ratios)
    #iterate over the dataframe. Note that acc = uniprot accession here.
    for acc in df.index:
        #assume af first that there is no previous data, and that the calculations can be re-run
        prev_calc_AAIMON_ratio_for_this_protein_exists = False
        if overwrite_prev_calculated_AAIMON_ratios == False:
            if dataframe_contains_prev_calc_AAIMON_ratios == True:
                cell_is_empty = np.isnan(df.loc[acc, 'TM01_AAIMON_ratio_mean'])
                #if the cell is empty (np.nan), the AAIMON ratios haven't been calculated for that protein yet
                if cell_is_empty:
                    prev_calc_AAIMON_ratio_for_this_protein_exists = False
                else:
                    #if the data is there, it should be a float and not a np.nan, therefore there is no need to repeat all the calculations
                    prev_calc_AAIMON_ratio_for_this_protein_exists = True
                    logging.info(
                        'calculate_AAIMON_ratios skipped, prev_calc_AAIMON_ratio_for_this_protein_exists = True for %s' %
                        df.loc[acc, 'A1_uniprot_accession'])
        else:
            #if the settings says to overwrite the data, ignore its existence
            prev_calc_AAIMON_ratio_for_this_protein_exists = False
        if overwrite_prev_calculated_AAIMON_ratios == True or prev_calc_AAIMON_ratio_for_this_protein_exists == False:
            organism_domain = df.loc[acc, 'organism_domain']
            protein_name = df.loc[acc, 'A2_protein_name']
            logging.info('%s' % protein_name)
            SIMAP_csv_from_XML_tarfile = df.loc[acc, 'SIMAP_csv_from_XML_tarfile']
            if os.path.exists(SIMAP_csv_from_XML_tarfile):
                try:
                    with tarfile.open(df.loc[acc, 'SIMAP_csv_from_XML_tarfile'], mode='r:gz') as tar:
                        if df.loc[acc, 'SIMAP_csv_from_XML'] in [tarinfo.name for tarinfo in tar]:
                            SIMAP_csv_from_XML_exists = True
                        else:
                            SIMAP_csv_from_XML_exists = False
                except (EOFError, tarfile.ReadError, OSError):
                    #file may be corrupted, if script stopped unexpectedly before compression was finished
                    logging.info('%s seems to be corrupted. File will be deleted.'
                                 'Will need to be parsed again next time program is run' %
                                  df.loc[acc, 'SIMAP_csv_from_XML_tarfile'])
                    SIMAP_csv_from_XML_exists = False
                    os.remove(df.loc[acc, 'SIMAP_csv_from_XML_tarfile'])
            else:
                SIMAP_csv_from_XML_exists = False
            if SIMAP_csv_from_XML_exists:
                with tarfile.open(SIMAP_csv_from_XML_tarfile, 'r:gz') as tar_in:
                    SIMAP_csv_from_XML_extracted = tar_in.extractfile(df.loc[acc, 'SIMAP_csv_from_XML'])
                    #reopen the csv file containing all the homologue data for that particular protein as a pandas dataframe (labelled Data Frame SIMAP, or dfs)
                    dfs = pd.read_csv(SIMAP_csv_from_XML_extracted, sep=",", index_col=0, quoting=csv.QUOTE_NONNUMERIC)
                    #create a column with the length of the Smith Waterman alignment (same length for query, markup and match)
                    if True in dfs['hit_contains_SW_node']:
                        at_least_one_hit_contains_SW_node = True
                    else:
                        at_least_one_hit_contains_SW_node = False
                    if at_least_one_hit_contains_SW_node:
                        try:
                            dfs['len_query_alignment_sequence'] = dfs['query_alignment_sequence'].str.len()
                        except KeyError:
                            #dataframe does not contain query_alignment_sequence, 
                            #which means that the XML file is probably damaged somehow
                            logging.warning('SIMAP_csv_from_XML seems to be damaged for %s' % protein_name)
                            dfs['query_alignment_sequence'] = ''
                            dfs['len_query_alignment_sequence'] = 0
                    else:
                        dfs['query_alignment_sequence'] = ''
                        dfs['len_query_alignment_sequence'] = 0

                    #add the list of words to the globals, to be accessed by utils.find_disallowed_words
                    words_not_allowed_in_description = settings["simap_match_filters"]["words_not_allowed_in_description"]
                    #collect disallowed words in hit protein description (patent, synthetic, etc)
                    dfs['list_disallowed_words_in_descr'] = dfs['A4_description'].dropna().apply(
                        utils.find_disallowed_words, args=(words_not_allowed_in_description,))
                    #create a boolean column to select hits that do not contain these words in the description
                    dfs['disallowed_words_not_in_descr'] = dfs['list_disallowed_words_in_descr'] == '[]'

                    list_of_TMDs = eval(df.loc[acc, 'list_of_TMDs'])

                    '''slice the TMD regions out of the alignment markup and match sequences
                    the method first finds the indices of the TMD region in the query, and uses these indices to slice
                    the filters are not applied yet, so the sequences can be viewed in the csv file for analysis
                    '''
                    for TMD in list_of_TMDs:
                        #create regex string for each TMD
                        query_TMD_sequence = df.loc[acc, '%s_seq' % TMD]
                        TMD_for_regular_expression_search = utils.create_regex_string(query_TMD_sequence)
                        #select to only include data where the XML contained a SW node, and then apply function for a regex search
                        dfs['%s_start_end_tuple_in_SW_alignment' % TMD] = dfs.query('hit_contains_SW_node == True')[
                                'query_alignment_sequence'].apply(utils.get_start_and_end_of_TMD_in_query,
                                                                  args=(TMD_for_regular_expression_search,)
                                                                  )
                        '''the output of the regex search is a tuple with three components. 
                        1) a match (e.g. True), 
                        2) the start of the match (e.g. 124), 
                        3) the end of the match (e.g. 144)
                        '''
                        #for convenience, separate the components into separate columns
                        columns_from_regex_output = ['%s_in_SW_alignment' % TMD, '%s_start_in_SW_alignment' % TMD,
                                                     '%s_end_in_SW_alignment' % TMD]
                        #n = the index (1,2,3) in the tuple 
                        #col = column item in the list (e.g. 'TM09_in_SW_alignment')
                        for n, col in enumerate(columns_from_regex_output):
                            #first filter to analyse only columns that contain a SW node
                            df_match = dfs.query('hit_contains_SW_node == True')
                            #add a new column which is named TM01_start, etc, and insert the appropriate integer (start or stop) or bool from the tuple
                            dfs[col] = df_match['%s_start_end_tuple_in_SW_alignment' % TMD].apply(lambda x: x[n])

                        # Defining juxtaregions 
                    
                    for TMD in list_of_TMDs:
                        last_TMD_of_acc = list_of_TMDs[-1]
                                                   
                        if TMD == "TM01":
                            dfs['start_juxta_before_TM01'] = np.where(dfs['TM01_start_in_SW_alignment']>0,0,np.nan)
                            dfs['end_juxta_before_TM01'] = np.where(dfs['TM01_start_in_SW_alignment']==0,np.nan,dfs['TM01_start_in_SW_alignment'])
                            dfs['start_juxta_after_TM01'] = dfs['TM01_end_in_SW_alignment']
                            dfs['end_juxta_after_TM01'] = dfs["TM01_end_in_SW_alignment"]+((dfs["TM02_start_in_SW_alignment"]-dfs["TM01_end_in_SW_alignment"])/2).apply(lambda x :int(x) if not np.isnan(x) else np.nan)
                            
                            #dfs['seq_juxta_after_TM01_in_query'] = dfs[dfs['start_juxta_after_TM01'].notnull()].apply(r_utils.slice_juxta_after_TMD_in_query, args = (TMD,), axis=1)
                            #dfs['seq_juxta_after_TM01_in_match'] = dfs[dfs['end_juxta_after_TM01'].notnull()].apply(r_utils.slice_juxta_after_TMD_in_match, args = (TMD,), axis=1) 

                            
                        if not TMD == "TM01" and not TMD == last_TMD_of_acc:
                            dfs['start_juxta_after_%s'%TMD] = np.where(r_utils.isNaN(dfs['TM%.2d_start_in_SW_alignment'%(int(TMD[2:])+1)])==True,np.nan,dfs['%s_end_in_SW_alignment'%TMD])
                            dfs['end_juxta_before_%s'%TMD] = np.where(dfs["%s_start_in_SW_alignment"%TMD]!=0,dfs["%s_start_in_SW_alignment"%TMD],np.nan)
                            dfs['end_juxta_after_%s'%TMD] = dfs["%s_end_in_SW_alignment"%TMD]+((dfs["TM%.2d_start_in_SW_alignment"%(int(TMD[2:])+1)]-dfs["%s_end_in_SW_alignment"%TMD])/2).apply(lambda x :int(x) if not np.isnan(x) else np.nan)
                            dfs['start_juxta_before_%s'%TMD] = np.where(dfs["end_juxta_after_TM%.2d"%(int(TMD[2:])-1)] == dfs['end_juxta_before_%s'%TMD] ,dfs["end_juxta_after_TM%.2d"%(int(TMD[2:])-1)],dfs["end_juxta_after_TM%.2d"%(int(TMD[2:])-1)])
                            #dfs['seq_juxta_after_%s_in_query'%TMD] = dfs[dfs['start_juxta_after_%s' % TMD].notnull()].apply(r_utils.slice_juxta_after_TMD_in_query, args = (TMD,), axis=1)
                            #dfs['seq_juxta_after_%s_in_match'%TMD] = dfs[dfs['end_juxta_after_%s' % TMD].notnull()].apply(r_utils.slice_juxta_after_TMD_in_match, args = (TMD,), axis=1)  

                            
                        if TMD == last_TMD_of_acc:   
                            dfs['start_juxta_before_%s'%TMD] = dfs['end_juxta_after_TM%.2d'%(int(TMD[2:])-1)]
                            dfs['end_juxta_before_%s'%TMD] = dfs['%s_start_in_SW_alignment'%TMD]
                            
                            dfs['start_juxta_after_%s'%TMD] = np.where(dfs['%s_end_in_SW_alignment'%TMD] == dfs['len_query_alignment_sequence'],np.nan,dfs['%s_end_in_SW_alignment'%TMD])
                            dfs['end_juxta_after_%s'%TMD] = np.where(r_utils.isNaN(dfs['start_juxta_after_%s'%TMD]) == True,np.nan,dfs['len_query_alignment_sequence'])
                            
                            
                            #dfs['seq_juxta_after_%s_in_query'%TMD] = dfs[dfs['start_juxta_after_%s' % TMD].notnull()].apply(r_utils.slice_juxta_after_TMD_in_query, args = (TMD,), axis=1)
                            
                            #dfs['seq_juxta_after_%s_in_query'%TMD] = dfs.query_alignment_sequence[int(dfs['start_juxta_after_TM10']):int(dfs['end_juxta_after_TM10'])]
                            #dfs['seq_juxta_after_%s_in_match'%TMD] =
                            
                    for TMD in list_of_TMDs:
                        last_TMD_of_acc = list_of_TMDs[-1]                    
                        dfs['seq_juxta_before_%s_in_query'%TMD] = dfs[dfs['start_juxta_before_%s' % TMD].notnull()].apply(r_utils.slice_juxta_before_TMD_in_query, args = (TMD,), axis=1)
                        dfs['seq_juxta_before_%s_in_match'%TMD] = dfs[dfs['start_juxta_before_%s' % TMD].notnull()].apply(r_utils.slice_juxta_before_TMD_in_match, args = (TMD,), axis=1)
                        if not TMD == last_TMD_of_acc:
                            dfs['seq_juxta_after_%s_in_query'%TMD] = dfs[dfs['end_juxta_after_%s' % TMD].notnull()].apply(r_utils.slice_juxta_after_TMD_in_query, args = (TMD,), axis=1)
                            dfs['seq_juxta_after_%s_in_match'%TMD] = dfs[dfs['end_juxta_after_%s' % TMD].notnull()].apply(r_utils.slice_juxta_after_TMD_in_match, args = (TMD,), axis=1)
                        else:
                            dfs['seq_juxta_after_%s_in_query'%TMD] = np.nan
                            dfs['seq_juxta_after_%s_in_match'%TMD] = np.nan
                            for hit in dfs.index:
                                if not r_utils.isNaN(dfs['start_juxta_after_%s'%TMD])[hit]:
                                    dfs['seq_juxta_after_%s_in_match'%TMD][hit] = dfs.match_alignment_sequence[hit][int(dfs["start_juxta_after_%s" % TMD][hit]):int(dfs["end_juxta_after_%s" % TMD][hit])]
                                    dfs['seq_juxta_after_%s_in_query'%TMD][hit] = dfs.query_alignment_sequence[hit][int(dfs["start_juxta_after_%s" % TMD][hit]):int(dfs["end_juxta_after_%s" % TMD][hit])]
                           
                            
                        #redefine the number of amino acids before and after the TMD to be inserted into the FastA files
                        aa_before_tmd = settings["simap_match_filters"]["fasta_files"]["aa_before_tmd"]
                        aa_after_tmd = settings["simap_match_filters"]["fasta_files"]["aa_after_tmd"]

                        #define the start of theTMD + surrounding sequence
                        dfs['%s_start_in_SW_alignment_plus_surr' % TMD] = dfs['%s_start_in_SW_alignment' % TMD] - aa_before_tmd
                        #replace negative values with zero
                        dfs.loc[dfs['%s_start_in_SW_alignment_plus_surr' % TMD] < 0, '%s_start_in_SW_alignment_plus_surr' % TMD] = 0
                        #define the end of the TMD + surrounding sequence
                        dfs['%s_end_in_SW_alignment_plus_surr' % TMD] = dfs['%s_end_in_SW_alignment' % TMD] + aa_after_tmd
                        #replace values longer than the actual sequence with the length of the sequence
                        dfs.loc[dfs['%s_end_in_SW_alignment_plus_surr' % TMD] > dfs['len_query_alignment_sequence'], 
                                    '%s_end_in_SW_alignment_plus_surr' % TMD] = dfs['len_query_alignment_sequence']


                    number_of_TMDs_containing_some_homologue_data = 0
                    for TMD in list_of_TMDs:
                        #in some cases, there is no data to obtain as the hit_contains_SW_node = False for too many sequences, giving no start_in_SW_alignment
                        number_of_rows_containing_data = dfs[dfs['%s_start_in_SW_alignment' % TMD].notnull()].shape[0]
                        if number_of_rows_containing_data == 0:
                            logging.info('%s does not have any valid homologues for %s. '
                                         'Re-downloading simap homologue XML may be necessary.' % (protein_name, TMD))
                        if number_of_rows_containing_data != 0:
                            number_of_TMDs_containing_some_homologue_data += 1
                            len_query_TMD = len(df.loc[acc, '%s_seq' % TMD])
                            #use small throwaway functions to slice each TMD from the query, markup and match sequences
                            #notnull() removes missing data
                            #explanation: dfs['new_column_with_selected_seq'] = dfs[dfs['only_rows_containing_data]].apply(utils.slice_function_that_specifies_columns_with_start_and_stop)
                            dfs['%s_SW_query_seq' % TMD] = dfs[dfs['%s_start_in_SW_alignment' % TMD].notnull()].apply(utils.slice_SW_query_TMD_seq, args = (TMD,), axis=1)
                            dfs['%s_SW_markup_seq' % TMD] = dfs[dfs['%s_start_in_SW_alignment' % TMD].notnull()].apply(utils.slice_SW_markup_TMD, args = (TMD,), axis=1)
                            dfs['%s_SW_match_seq' % TMD] = dfs[dfs['%s_start_in_SW_alignment' % TMD].notnull()].apply(utils.slice_SW_match_TMD_seq, args = (TMD,), axis=1)
                            #and the same for the TMD + surrounding sequence, useful to examine the TMD interface
                            dfs['%s_SW_query_seq_plus_surr' % TMD] = dfs[dfs['%s_start_in_SW_alignment_plus_surr' % TMD].notnull()].apply(utils.slice_SW_query_TMD_seq_plus_surr, args = (TMD,), axis=1)
                            dfs['%s_SW_markup_seq_plus_surr' % TMD] = dfs[dfs['%s_start_in_SW_alignment_plus_surr' % TMD].notnull()].apply(utils.slice_SW_markup_TMD_plus_surr, args = (TMD,), axis=1)
                            dfs['%s_SW_match_seq_plus_surr' % TMD] = dfs[dfs['%s_start_in_SW_alignment_plus_surr' % TMD].notnull()].apply(utils.slice_SW_match_TMD_seq_plus_surr, args = (TMD,), axis=1)
                            #check if there are non-IUPAC amino acids in the sequence (frequently large gaps from NG sequencing data)
                            dfs['X_in_match_seq'] = 'X' in ['match_alignment_sequence']
                            #count the number of gaps in the query and match sequences
                            dfs['%s_SW_query_num_gaps' % TMD] = dfs['%s_SW_query_seq' % TMD].dropna().apply(lambda x: x.count('-'))
                            dfs['%s_SW_match_num_gaps' % TMD] = dfs['%s_SW_match_seq' % TMD].dropna().apply(lambda x: x.count('-'))
                            #calculate the length of the match TMD seq excluding gaps
                            dfs['%s_SW_m_seq_len' % TMD] = dfs['%s_SW_match_seq' % TMD].str.len()
                            #for the alignment length, take the smallest value from the length of query or match 
                            #this will exclude gaps from the length in the following calculations, preventing false "low conservation" where the query TMD is much longer than the match TMD)
                            #note that for most calculations this is somewhat redundant, because the max number of acceptable gaps in sequence is probable ~2
                            dfs['%s_SW_align_len' % TMD] = dfs['%s_SW_m_seq_len' % TMD].apply(lambda x: x if x < len_query_TMD else len_query_TMD)
                            #create a boolean column that allows filtering by the accepted number of gaps, according to the settings file
                            dfs['%s_SW_query_acceptable_num_gaps' % TMD] = dfs['%s_SW_query_num_gaps' % TMD] <= settings["simap_match_filters"]["number_of_gaps_allowed_in_query_TMD"]
                            dfs['%s_SW_match_acceptable_num_gaps' % TMD] = dfs['%s_SW_match_num_gaps' % TMD] <= settings["simap_match_filters"]["number_of_gaps_allowed_in_match_TMD"]
                            #count identical residues between query and match TMDs by counting the number of pipes in the markup string                     
                            dfs['%s_SW_num_ident_res' % TMD] = dfs['%s_SW_markup_seq' % TMD].dropna().apply(lambda x: x.count('|'))
                            dfs['%s_SW_num_sim_res' % TMD] = dfs['%s_SW_markup_seq' % TMD].dropna().apply(lambda x: x.count(':'))
                            #check that the TMD seq in match is not just 100% gaps!
                            dfs['%s_in_SW_align_match' % TMD] = dfs['%s_SW_num_ident_res' % TMD].dropna() != 0
                            dfs['%s_in_SW_align_match' % TMD].fillna(value=False)

                            #the percentage identity of that TMD is defined as the number of identical residues (pipes in markup) divided by the length of the the aligned residues (excluding gaps, based on the length of the shortest TMD, either match or query)                
                            #note that the nonTMD percentage identity is calculated the same way
                            dfs['%s_perc_ident' % TMD] = dfs['%s_SW_num_ident_res' % TMD] / dfs['%s_SW_align_len' % TMD]
                            #calculate percentage similar residues
                            dfs['%s_perc_sim' % TMD] = dfs['%s_SW_num_sim_res' % TMD] / dfs['%s_SW_align_len' % TMD]
                            #add together to obtain the percentage similar + identical residues
                            dfs['%s_perc_sim_plus_ident' % TMD] = dfs['%s_perc_ident' % TMD] + dfs['%s_perc_sim' % TMD]
                            #add to main dataframe
                            dfs['%s_perc_ident' % TMD] = dfs['%s_perc_ident' % TMD]
                            #calculate the average number of gaps per residue in the TMD alignment
                            #(number of gaps)/(length of sequence excluding gaps)
                            dfs['%s_SW_q_gaps_per_q_residue' % TMD] = dfs['%s_SW_query_num_gaps' % TMD].dropna() / len_query_TMD

                    if number_of_TMDs_containing_some_homologue_data == 0:
                        #there was no data obtained from the csv file, which probably means the original XML file was not properly downloaded
                        #there is no need to continue the script for this protein
                        logging.warning('%s does not contain any valid homologues at all. '
                                        'CSV and/or XML file is damaged. Re-downloading simap homologue XML may be necessary. '
                                        'AAIMON ratios will not be calculated for this protein.' % protein_name)

                    if number_of_TMDs_containing_some_homologue_data > 0:
                        #create a boolean column that describel whether the sequence is above the minimum gapped identity
                        minimum_identity_of_full_protein = settings["simap_match_filters"]["minimum_identity_of_full_protein"]
                        dfs['gapped_ident_above_cutoff'] = dfs['FASTA_gapped_identity'] > minimum_identity_of_full_protein

                        '''¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦             nonTMD calculations                        ¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦¦
                        '''
                        #filter to only analyse sequences with tho following:
                        #1) full protein identity above cutoff 
                        #2)containing smith waterman alignment in XML file 
                        #3) hit description lacks disallowedd words (patent, synthetic, etc)
                        dfs_filt = dfs.query('gapped_ident_above_cutoff == True and '
                                             'hit_contains_SW_node == True and '
                                             'disallowed_words_not_in_descr == True')

                        #check if all tmds are in SW alignment
                        #create a list of columns to reindex the DataFrame
                        list_columns_TMD_in_SW_alignment = []
                        for TMD in list_of_TMDs:
                            #TMD found by regex
                            list_columns_TMD_in_SW_alignment.append('%s_in_SW_alignment' % TMD)
                            #TMD matching useful sequence
                            list_columns_TMD_in_SW_alignment.append('%s_in_SW_align_match' % TMD)

                        #create a slice of the filtered dataframe that only contains the relevant columns (N.B. copy=False, this will provide a view, not a copy)
                        df2 = dfs_filt.reindex(index=dfs.index, columns=list_columns_TMD_in_SW_alignment, copy=False)
                        #create a new column in the original dataframe that shows that ALL TMDs have been found in the SW alignment
                        dfs['all_tmds_in_SW_alignment'] = df2.dropna().all(axis=1)
                        #filter to contain only hits with all tmds
                        #create a copy of the original dfs dataframe containing only hits where all tmds are found in the match
                        dfs_nonTMD = dfs.loc[dfs['all_tmds_in_SW_alignment'].notnull()].query('all_tmds_in_SW_alignment == True')
                        #filter to contain only hits that do not include disallowed words
                        dfs_nonTMD = dfs_nonTMD.query('disallowed_words_not_in_descr == True')
                        #print(dfs_nonTMD2.shape)
                        #filter to contain only hits where the index for the TMD is present
                        first_TMD_start_index = '%s_start_in_SW_alignment' % list_of_TMDs[0]

                        dfs_nonTMD = dfs_nonTMD.loc[dfs_nonTMD[first_TMD_start_index].notnull()]
                        #print(dfs_nonTMD2.shape)
                        #print(dfs_nonTMD2[first_TMD_end_index])
                        dfs_nonTMD['nonTMD_index_tuple_first'] = dfs_nonTMD[first_TMD_start_index].apply(lambda x: (0, int(x)))

                        #create start and stop indices for all sections between tmds
                        #the start of the last nonTMD section will be the end of the last TMD
                        dfs_nonTMD['nonTMD_index_tuple_last0'] = dfs_nonTMD['%s_end_in_SW_alignment' % list_of_TMDs[-1]].dropna().astype('int32')
                        #the end of the last nonTMD section will be the end of the full alignment sequence
                        dfs_nonTMD['nonTMD_index_tuple_last1'] = dfs_nonTMD['len_query_alignment_sequence'].dropna().astype('int32')
                        #join to make a tuple
                        #dfs_nonTMD['nonTMD_index_tuple_last'] = dfs_nonTMD[['nonTMD_index_tuple_last0', 'nonTMD_index_tuple_last1']].apply(tuple, axis=1)

                        #create the index tuple
                        dfs_nonTMD['nonTMD_index_tuple_last'] = dfs_nonTMD.apply(utils.create_indextuple_nonTMD_last, axis=1)

                        #for each TMD EXCEPT the last, which ends at the sequence end, create the indices for the nonTMD region (after the TMD)
                        for TM_Nr in range(len(list_of_TMDs) - 1):
                            #the TMD is the equivalent item in the list
                            TMD = list_of_TMDs[TM_Nr]
                            #the next TMD, which contains the end index, is the next item in the list
                            next_TMD = list_of_TMDs[TM_Nr + 1]
                            #select only the columns in the dataframe that are of interest, and change the data type to integer
                            index_columns = ['%s_end_in_SW_alignment' % TMD, '%s_start_in_SW_alignment' % next_TMD]
                            dfs_nonTMD[index_columns] = dfs_nonTMD[index_columns].astype('int64')
                            #create a tuple containing the indices for the nonTMD sequence regions in between each TMD (middle indices)                 
                            dfs_nonTMD['nonTMD_index_%s' % TMD] = tuple(zip(dfs_nonTMD['%s_end_in_SW_alignment' % TMD],
                                                                  dfs_nonTMD['%s_start_in_SW_alignment' % next_TMD]))

                        #now join all the indices together to make one tuple of tuples for the non-TMD region
                        #dfs_nonTMD = dfs.query('all_tmds_in_SW_alignment == True')
                        dfs_nonTMD['nested_tuple_indices_all_nonTMD_regions'] = dfs_nonTMD[['nonTMD_index_tuple_last']].apply(tuple, axis=1)

                        #create a view of the dataframe that contains only the desired columns
                        #first create a list of the desired columns
                        #start with the first tuple, from 0 to the start of the first TMD
                        list_of_nonTMD_index_columns = ['nonTMD_index_tuple_first']
                        #create a nonTMD region for each of the TMDs (except the last one)
                        list_from_TMs = ['nonTMD_index_%s' % TMD2 for TMD2 in list_of_TMDs[:-1]]
                        #join lists
                        list_of_nonTMD_index_columns = list_of_nonTMD_index_columns + list_from_TMs
                        list_of_nonTMD_index_columns += ['nonTMD_index_tuple_last']

                        #create the new view by reindexing the dataframe with the list of desired columns
                        dfs_tuples = dfs_nonTMD.reindex(index=dfs_nonTMD.index, columns=list_of_nonTMD_index_columns,
                                                        copy=False)  #.astype('int32')
                        #now for convenience, these tuples can be combined together to form one column, with a tuple of tuples
                        #first convert all values in each row to a list, excluding the index column
                        list_tuple_indices_all_nonTMD_regions = list(dfs_tuples.itertuples(index=False))
                        #convert to a series, and reindex with the original index from the dataframe
                        tuples_series = pd.Series(list_tuple_indices_all_nonTMD_regions, index=dfs_nonTMD.index)

                        dfs_nonTMD['nested_tuple_indices_all_nonTMD_regions'] = tuples_series
                        #change to a string, in case this solves the weird effect with only the last tuple shown
                        dfs_nonTMD['nested_tuple_indices_all_nonTMD_regions'] = dfs_nonTMD[
                            'nested_tuple_indices_all_nonTMD_regions'].astype(str)
                        #you can test that the original index is maintained as follows:
                        #add the series as a new column in the original dataframe. Missing data (when not all TMDs found) will be filled using np.nan
                        dfs['nested_tuple_indices_all_nonTMD_regions'] = dfs_nonTMD['nested_tuple_indices_all_nonTMD_regions']

                        #filter to remove incomplete sequences
                        dfs_nonTMD = dfs_nonTMD.query('all_tmds_in_SW_alignment == True')
                        #define the string for slicing as a numpy array
                        #use the numpy vectorize function, which effectively applies the function in a for loop (not optimized for speed)
                        #dfs_nonTMD['nonTMD_seq_query'] = np.vectorize(utils.slice_with_nested_tuple)(np.array(dfs_nonTMD['query_alignment_sequence']),np.array(dfs_nonTMD['nested_tuple_indices_all_nonTMD_regions']))
                        #dfs_nonTMD['nonTMD_markup'] = np.vectorize(utils.slice_with_nested_tuple)(np.array(dfs_nonTMD['alignment_markup']),np.array(dfs_nonTMD['nested_tuple_indices_all_nonTMD_regions']))
                        #dfs_nonTMD['nonTMD_seq_match'] = np.vectorize(utils.slice_with_nested_tuple)(np.array(dfs_nonTMD['match_alignment_sequence']),np.array(dfs_nonTMD['nested_tuple_indices_all_nonTMD_regions']))

                        #due to problems with np.vectorize and the pandas methods, slice the sequences one at a time with a simple 'for loop'
                        #create empty columns
                        dfs['nonTMD_seq_query'] = ''
                        dfs['nonTMD_markup'] = ''
                        dfs['nonTMD_seq_match'] = ''
                        #for each hit, perform the slice
                        for hit in dfs_nonTMD.index:
                            dfs_nonTMD.loc[hit, 'nonTMD_seq_query'] = utils.slice_with_nested_tuple(dfs_nonTMD.loc[hit, 'query_alignment_sequence'],
                                                                      dfs.loc[hit, 'nested_tuple_indices_all_nonTMD_regions']
                                                                      )
                            dfs_nonTMD.loc[hit, 'nonTMD_markup'] = utils.slice_with_nested_tuple(dfs_nonTMD.loc[hit, 'alignment_markup'],
                                                                   dfs_nonTMD.loc[hit, 'nested_tuple_indices_all_nonTMD_regions']
                                                                   )
                            dfs_nonTMD.loc[hit, 'nonTMD_seq_match'] = utils.slice_with_nested_tuple(dfs_nonTMD.loc[hit, 'match_alignment_sequence'],
                                                                      dfs_nonTMD.loc[hit, 'nested_tuple_indices_all_nonTMD_regions']
                                                                      )
                        #transfer to original dataframe (index should still match the original, partial seqs will be filled with np.nan)
                        dfs['nonTMD_seq_query'] = dfs_nonTMD['nonTMD_seq_query']
                        dfs['nonTMD_markup'] = dfs_nonTMD['nonTMD_markup']
                        dfs['nonTMD_seq_match'] = dfs_nonTMD['nonTMD_seq_match']
                        '''FILLNA REMOVED TO MAKE SORTING EASIER
                        '''
                        #fill in np.nan with empty strings
                        #dfs['nonTMD_seq_query'] = dfs['nonTMD_seq_query'].fillna(value = '')
                        #dfs['nonTMD_markup'] = dfs['nonTMD_markup'].fillna(value = '')
                        #dfs['nonTMD_seq_match'] = dfs['nonTMD_seq_match'].fillna(value = '')                   

                        #calculate identical residues in the nonTMD region (simply count the pipes '|' in the markup sequence)
                        dfs['nonTMD_num_ident_res'] = dfs['nonTMD_markup'].dropna().apply(lambda x: x.count('|'))
                        #calculate similar residues in the nonTMD region (simply count the colons ':' in the markup sequence)
                        dfs['nonTMD_num_sim_res'] = dfs['nonTMD_markup'].dropna().apply(lambda x: x.count(':'))
                        #add the identical and similar residues together to get the total number of similar + identical residues
                        dfs['nonTMD_num_sim_plus_ident_res'] = dfs['nonTMD_num_ident_res'] + dfs['nonTMD_num_sim_res']

                        #count the gaps in the nonTMD sequence of the query
                        dfs['nonTMD_q_num_gaps'] = dfs['nonTMD_seq_query'].dropna().apply(lambda x: x.count('-'))
                        #count the gaps in the nonTMD sequence of the match
                        dfs['nonTMD_m_num_gaps'] = dfs['nonTMD_seq_match'].dropna().apply(lambda x: x.count('-'))
                        #calculate the length of the nonTMD sequences, which may include gaps
                        dfs['len_nonTMD_seq_query'] = dfs['nonTMD_seq_query'].dropna().str.len()
                        dfs['len_nonTMD_seq_match'] = dfs['nonTMD_seq_match'].dropna().str.len()
                        #calculate the number aligned sequences, excluding gaps (length of query, or length of match, whichever is shorter)
                        dfs['len_nonTMD_align'] = dfs[['len_nonTMD_seq_query', 'len_nonTMD_seq_match']].dropna(
                            how='all').min(axis=1)

                        #calculate the length of the nonTMD sequence excluding gaps
                        dfs['len_nonTMD_q_excl_gaps'] = dfs['len_nonTMD_seq_query'] - dfs['nonTMD_q_num_gaps']
                        dfs['len_nonTMD_m_excl_gaps'] = dfs['len_nonTMD_seq_match'] - dfs['nonTMD_m_num_gaps']
                        #calculate the lenth of the alignment by finding which seq excl gaps is smaller
                        dfs['len_nonTMD_align'] = dfs[['len_nonTMD_q_excl_gaps', 'len_nonTMD_m_excl_gaps']].min(axis=1)

                        #calculate the percentage identity of the nonTMD region (number of identical residues divided by the length excluding gaps)
                        #used for the Amino Acid Identity : Membranous over Nonmembranous (AAIMON ratio)
                        #note that the length = length of the aligned residues excluding gaps
                        dfs['nonTMD_perc_ident'] = dfs['nonTMD_num_ident_res'] / dfs['len_nonTMD_align']
                        dfs['nonTMD_perc_sim'] = dfs['nonTMD_num_sim_res'] / dfs['len_nonTMD_align']
                        dfs['nonTMD_perc_sim_plus_ident'] = dfs['nonTMD_num_sim_plus_ident_res'] / dfs[
                            'len_nonTMD_align']
                        #calculate the average number of gaps per residue in the nonTMD alignment
                        #filter to analyse only sequences that are valid (length > 0)
                        dfs_filt_gaps = dfs.loc[dfs['len_nonTMD_q_excl_gaps'] != 0]
                        #calculate number of gaps in query AND match
                        dfs_filt_gaps['nonTMD_qm_num_gaps'] = dfs_filt_gaps['nonTMD_q_num_gaps'] + dfs_filt_gaps['nonTMD_m_num_gaps']
                        #add to simap dataframe
                        dfs['nonTMD_qm_num_gaps'] = dfs_filt_gaps['nonTMD_qm_num_gaps']
                        #gaps per query residue for both query and match = ((gaps in query + gaps in match)/2))/length of query excluding gaps
                        dfs['nonTMD_qm_gaps_per_q_residue'] = dfs_filt_gaps['nonTMD_qm_num_gaps'] / 2 / dfs_filt_gaps['len_nonTMD_q_excl_gaps']

                        '''  _________________________________________SSR ratio calculations____________________________________________________
                        '''

                        #                    def calc_score_ss_qTMD(dfs):
                        #                        score = sum(utils.score_pairwise(seq1=dfs['%s_SW_query_seq' % TMD], seq2=dfs['%s_SW_query_seq' % TMD],
                        #                                                     matrix=aa_sub_matrix,
                        #                                                     gap_open_penalty=gap_open_penalty,
                        #                                                     gap_extension_penalty=gap_extension_penalty))
                        #                        return(score)
                        #
                        #                    def calc_score_ss_mTMD(dfs):
                        #                        score = sum(utils.score_pairwise(seq1=dfs['%s_SW_match_seq' % TMD], seq2=dfs['%s_SW_match_seq' % TMD],
                        #                                                     matrix=aa_sub_matrix,
                        #                                                     gap_open_penalty=gap_open_penalty,
                        #                                                     gap_extension_penalty=gap_extension_penalty))
                        #                        return(score)
                        #
                        #                    def calc_score_qTMD_mTMD(dfs):
                        #                        score = sum(utils.score_pairwise(seq1=dfs['%s_SW_query_seq' % TMD], seq2=dfs['%s_SW_match_seq' % TMD],
                        #                                                     matrix=aa_sub_matrix,
                        #                                                     gap_open_penalty=gap_open_penalty,
                        #                                                     gap_extension_penalty=gap_extension_penalty))
                        #                        return(score)
                        #
                        #                    def calc_ss_q_nonTMD(dfs):
                        #                        score = sum(utils.score_pairwise(seq1=dfs['nonTMD_seq_query'], seq2=dfs['nonTMD_seq_query'],
                        #                                                     matrix=aa_sub_matrix,
                        #                                                     gap_open_penalty=gap_open_penalty,
                        #                                                     gap_extension_penalty=gap_extension_penalty))
                        #                        return(score)
                        #
                        #                    def calc_ss_m_nonTMD(dfs):
                        #                        score = sum(utils.score_pairwise(seq1=dfs['nonTMD_seq_match'], seq2=dfs['nonTMD_seq_match'],
                        #                                                     matrix=aa_sub_matrix,
                        #                                                     gap_open_penalty=gap_open_penalty,
                        #                                                     gap_extension_penalty=gap_extension_penalty))
                        #                        return(score)
                        #
                        #                    def calc_q_m_nonTMD(dfs):
                        #                        score = sum(utils.score_pairwise(seq1=dfs['nonTMD_seq_query'], seq2=dfs['nonTMD_seq_match'],
                        #                                                     matrix=aa_sub_matrix,
                        #                                                     gap_open_penalty=gap_open_penalty,
                        #                                                     gap_extension_penalty=gap_extension_penalty))
                        #                        return(score)
                        #
                        #                    for gap_open_penalty in range(settings["aa_substitution_scoring"]["gap_open_penalty_min"],settings["aa_substitution_scoring"]["gap_open_penalty_max"],
                        #                                   settings["aa_substitution_scoring"]["gap_open_penalty_increment"]):
                        #                        #print(gap_open_penalty)
                        #                        #for simplicity, give the gap open and gap extend the same value
                        #                        gap_extension_penalty = gap_open_penalty
                        #                        for matrix_name in list_of_aa_sub_matrices:
                        #                            #so long as the matrix is imported into python, eval will convert the matrix name to an object
                        #                            aa_sub_matrix = eval(matrix_name)
                        #                            #update the matrix (unsure what this deos! Taken directly from Stackoverflow)
                        #                            aa_sub_matrix.update(((b, a), val) for (a, b), val in list(aa_sub_matrix.items()))
                        #                            column_basename = 'sim_ratio_%s_gapo%i' % (matrix_name.replace("'", "")[0:-7], gap_open_penalty)
                        #                            #print(column_name)
                        #                            #scores of query and match. Jan used match-query and query-match because in some cases in his algorithm it gave different scores. In this case, the score is independent of the sequence order, and 2*score_qTMD_mTMD can be used instead of score_qTMD_mTMD + score_mTMD_qTMD
                        #                            #score_qTMD_mTMD = sum(utils.score_pairwise(seq1=SW_query_TMD_seq, seq2=SW_match_TMD_seq,
                        #                            #                         matrix=aa_sub_matrix,
                        #                            #                         gap_open_penalty=gap_open_penalty,
                        #                            #                         gap_extension_penalty=gap_extension_penalty))
                        #                            #print(score_qTMD_mTMD)
                        #                            dfs_nonTMD = dfs.query('"X" not in match_alignment_sequence')
                        #                            dfs_nonTMD = dfs_nonTMD.loc[dfs_nonTMD['match_alignment_sequence'].notnull()]
                        #                            dfs_nonTMD = dfs_nonTMD.loc[dfs_nonTMD['nonTMD_seq_query'].notnull()]
                        #                            dfs_nonTMD = dfs_nonTMD.loc[dfs_nonTMD['match_alignment_sequence'].apply(lambda x : 'X' not in x)]
                        #
                        #                            #score/self-score ratio of nonTMD query
                        #                            dfs_nonTMD[column_basename + '_ss_q_nonTMD'] = dfs_nonTMD.apply(calc_ss_q_nonTMD, axis = 1)
                        #                            #score/self-score ratio of match
                        #                            dfs_nonTMD[column_basename + '_ss_m_nonTMD'] = dfs_nonTMD.apply(calc_ss_m_nonTMD, axis = 1)
                        #                            #scores of query and match. Jan used match-query and query-match because in some cases in his algorithm it gave different scores. In this case, the score is independent of the sequence order, and 2*score_qTMD_mTMD can be used instead of score_qTMD_mTMD + score_mTMD_qTMD
                        #                            dfs_nonTMD[column_basename + '_q_m_nonTMD'] = dfs_nonTMD.apply(calc_q_m_nonTMD, axis = 1)
                        #                            #calculate the score/selfscore ratio
                        #                            dfs_nonTMD[column_basename + '_ssr_nonTMD'] = dfs_nonTMD[column_basename + '_q_m_nonTMD'] * 2 / (dfs_nonTMD[column_basename + '_ss_q_nonTMD'] + dfs_nonTMD[column_basename + '_ss_m_nonTMD'])
                        #                            #add to main dataframe
                        #                            dfs[column_basename + '_ssr_nonTMD'] = dfs_nonTMD[column_basename + '_ssr_nonTMD']
                        #
                        #                            for TMD in list_of_TMDs:
                        #                                column_name = TMD + '_' + column_basename
                        #                                dfs_nonTMD = dfs_nonTMD.loc[dfs['%s_SW_query_seq' % TMD].notnull()]
                        #                                #dfs_nonTMD = dfs_nonTMD.loc[dfs['X_in_match_seq'] == False]
                        #                                #score/self-score ratio of query
                        #                                dfs[column_name + '_ss_qTMD'] = dfs_nonTMD.apply(calc_score_ss_qTMD, axis = 1)
                        #                                #score/self-score ratio of match
                        #                                dfs[column_name + '_ss_mTMD'] = dfs_nonTMD.apply(calc_score_ss_mTMD, axis = 1)
                        #                                #scores of query and match. Jan used match-query and query-match because in some cases in his algorithm it gave different scores. In this case, the score is independent of the sequence order, and 2*score_qTMD_mTMD can be used instead of score_qTMD_mTMD + score_mTMD_qTMD
                        #                                dfs[column_name + '_qTMD_mTMD'] = dfs_nonTMD.apply(calc_score_qTMD_mTMD, axis = 1)
                        #                                #score/self-score ratio
                        #                                dfs[column_name + '_ssrTMD'] = dfs[column_name + '_qTMD_mTMD'] * 2 / (dfs[column_name + '_ss_qTMD'] + dfs[column_name + '_ss_mTMD'])
                        #                                #calculate the ssrTMD/ssr_nonTMD
                        #                                dfs_filt3 = dfs.loc[dfs[column_name + '_ssrTMD'].notnull()]
                        #                                dfs_filt3 = dfs.loc[dfs[column_name + '_ssrTMD'] > 0]
                        #                                dfs_filt3 = dfs.loc[dfs[column_basename + '_ssr_nonTMD'] > 0]
                        #                                dfs_filt3[column_name + '_ssrTMD_over_nonTMD'] = dfs[column_name + '_ssrTMD'] / dfs[column_basename + '_ssr_nonTMD']
                        #                                #add to main dataframe
                        #                                dfs[column_name + '_ssrTMD_over_nonTMD'] = dfs_filt3[column_name + '_ssrTMD_over_nonTMD']
                        '''  _________________________________________END SSR ratio calculations____________________________________________________
                        '''
                        #re-filter the original dataframe to create another copy with the desired sequences
                        dfs_filt = dfs.query('gapped_ident_above_cutoff == True and '
                                             'hit_contains_SW_node == True and '
                                             'disallowed_words_not_in_descr == True and '
                                             'X_in_match_seq == False')
                        
                        '''Calculate average values, add to original dataframe.
                           1) values associated with the FASTA output of SIMAP
                        '''
                        #fasta identity
                        df.loc[acc, 'FASTA_ident_mean'] = float('%0.2f' % dfs['FASTA_identity'].mean())
                        #number of identical residues in FASTA alignment
                        dfs['FASTA_num_ident_res'] = dfs_filt['FASTA_identity'] / 100 * dfs_filt['FASTA_overlap']
                        df.loc[acc, 'FASTA_num_ident_res'] = float('%0.2f' % dfs_filt['FASTA_identity'].mean())
                        
                        '''2) values associated with the nonTMD region
                        '''                        
                        #add the average values regarding the nonTMD region to the original file/dataframe with each protein
                        df.loc[acc, 'len_nonTMD_seq_match_mean'] = float('%0.2f' % dfs_filt['len_nonTMD_seq_match'].dropna().mean())
                        df.loc[acc, 'nonTMD_perc_ident_mean'] = float('%0.3f' % dfs_filt['nonTMD_perc_ident'].dropna().mean())
                        df.loc[acc, 'nonTMD_perc_sim_mean'] = float('%0.3f' % dfs_filt['nonTMD_perc_sim'].dropna().mean())
                        df.loc[acc, 'nonTMD_perc_sim_plus_ident_mean'] = float('%0.3f' % dfs_filt['nonTMD_perc_sim_plus_ident'].dropna().mean())
                        df.loc[acc, 'len_nonTMD_align_mean'] = float('%0.2f' % dfs_filt['len_nonTMD_align'].dropna().mean())
                        df.loc[acc, 'nonTMD_qm_gaps_per_q_residue_mean'] = float('%0.2f' % dfs_filt['nonTMD_qm_gaps_per_q_residue'].dropna().mean())
                        logging.info('nonTMD_qm_gaps_per_q_residue : %0.5f' % df.loc[acc, 'nonTMD_qm_gaps_per_q_residue_mean'])

                        '''3) values associated each TMD, such as average AAIMON ratio
                        '''  
                        #calculate AAISMON etc for each TMD
                        for TMD in list_of_TMDs:
                            len_query_TMD = len(df.loc[acc, '%s_seq' % TMD])
                            #following the general filters, filter to only analyse sequences with TMD identity above cutoff, and a nonTMD_perc_ident above zero ,to avoid a divide by zero error
                            min_identity_of_TMD_initial_filter = settings['simap_match_filters']['min_identity_of_TMD_initial_filter']
                            dfs_filt_AAIMON = dfs_filt.loc[dfs['%s_perc_ident' % TMD] >= min_identity_of_TMD_initial_filter]
                            #avoid a divide by zero error in the unlikely case that there are no_identical_residues_in_alignment
                            dfs_filt_AAIMON = dfs_filt_AAIMON.loc[dfs_filt_AAIMON['nonTMD_perc_ident'] != 0]
                            #calculate the Amino Acid Identity : Membranous Over Nonmembranous
                            dfs_filt_AAIMON['%s_AAIMON_ratio' % TMD] = dfs_filt_AAIMON['%s_perc_ident' % TMD] / dfs_filt_AAIMON['nonTMD_perc_ident']
                            #calculate the Amino Acid Similarity : Membranous Over Nonmembranous (AASMON) (includes similarity + identity based on the matrix used in the SW alignment of SIMAP)
                            dfs_filt_AAIMON['%s_AASMON_ratio' % TMD] = dfs_filt_AAIMON['%s_perc_sim_plus_ident' % TMD] / dfs_filt_AAIMON['nonTMD_perc_sim_plus_ident']
                            df.loc[acc, '%s_SW_q_gaps_per_q_residue_mean' % TMD] = dfs_filt_AAIMON['%s_SW_q_gaps_per_q_residue' % TMD].dropna().mean()
                            logging.info('%s_SW_q_gaps_per_q_residue Average: %0.3e' % 
                                        (TMD, df.loc[acc, '%s_SW_q_gaps_per_q_residue_mean' % TMD])
                                        )

                            #add to original dataframe with the list of uniprot sequences

                            df.loc[acc, '%s_perc_ident_mean' % TMD] = dfs_filt_AAIMON['%s_perc_ident' % TMD].mean()
                            df.loc[acc, '%s_perc_sim_mean' % TMD] = dfs_filt_AAIMON['%s_perc_sim' % TMD].mean()
                            df.loc[acc, '%s_perc_sim_plus_ident_mean' % TMD] = dfs_filt_AAIMON['%s_perc_sim_plus_ident' % TMD].mean()
                            df.loc[acc, '%s_AAIMON_ratio_mean' % TMD] = float(dfs_filt_AAIMON['%s_AAIMON_ratio' % TMD].mean())
                            df.loc[acc, '%s_AAIMON_ratio_std' % TMD] = dfs_filt_AAIMON['%s_AAIMON_ratio' % TMD].std()
                            df.loc[acc, '%s_AASMON_ratio_mean' % TMD] = dfs_filt_AAIMON['%s_AASMON_ratio' % TMD].mean()
                            df.loc[acc, '%s_AASMON_ratio_std' % TMD] = dfs_filt_AAIMON['%s_AASMON_ratio' % TMD].std()
                            logging.info('AAIMON MEAN %s: %0.2f' % (TMD, df.loc[acc, '%s_AAIMON_ratio_mean' % TMD]))
                            logging.info('AAISON MEAN %s: %0.2f' % (TMD, df.loc[acc, '%s_AASMON_ratio_mean' % TMD]))

                            #add to the dataframe with the SIMAP data for that particular protein
                            dfs['%s_AAIMON_ratio' % TMD] = dfs_filt_AAIMON['%s_AAIMON_ratio' % TMD]
                            dfs['%s_AASMON_ratio' % TMD] = dfs_filt_AAIMON['%s_AASMON_ratio' % TMD]

                            dfs['%s_ratio_length_of_TMD_to_rest_of_alignment' % TMD] = dfs['%s_SW_query_seq' % TMD].str.len() / dfs['FASTA_overlap']
                            dfs['%s_ratio_length_of_query_TMD_to_rest_of_match_protein' % TMD] = dfs['%s_SW_query_seq' % TMD].str.len() / dfs['len_full_match_seq']

                            df.loc[acc, '%s_ratio_length_of_TMD_to_rest_of_alignment_mean' % TMD] = float('%0.2f' % dfs['%s_ratio_length_of_TMD_to_rest_of_alignment' % TMD].dropna().mean())
                            df.loc[acc, '%s_ratio_length_of_query_TMD_to_rest_of_match_protein_mean' % TMD] = float('%0.2f' % dfs['%s_ratio_length_of_query_TMD_to_rest_of_match_protein' % TMD].dropna().mean())

                        '''re-filter the original dataframe to create another copy with the desired sequences
                        note that some values were added after filtering in the last round, 
                        but all were added to the dataframe dfs, not the copy dfs_filt
                        '''
                        dfs_filt = dfs.query('gapped_ident_above_cutoff == True and '
                                             'hit_contains_SW_node == True and '
                                             'disallowed_words_not_in_descr == True and '
                                             'X_in_match_seq == False')

                        #use linspace to get a fixid number of points between tha min and the max for the histogram
                        #set up evenly distributed bins between the chosen min and max
                        #if possible, 1.0 should be in the centre of a bin, to catch cases where a lot of homologues have a ratio that approximates 1
                        linspace_binlist = np.linspace(settings["hist_settings_single_protein"]["smallest_bin"],
                                                       settings["hist_settings_single_protein"]["largest_bin"],
                                                       settings["hist_settings_single_protein"]["number_of_bins"])
                        #add 30 as the last bin, to make sure 100% of the data is added to the histogram, including major outliers
                        binlist = np.append(linspace_binlist,
                                            settings["hist_settings_single_protein"]["final_highest_bin"])

                        #se default font size for text in the plot
                        fontsize = 4
                        #use a dictionary to organise the saving of multiple plots in multiple figures, with a certain number of plots per figure
                        n_plots_per_fig = 4
                        nrows_in_each_fig = 2
                        ncols_in_each_fig = 2
                        dict_organising_subplots = utils.create_dict_organising_subplots(
                                                   n_plots_per_fig=n_plots_per_fig,
                                                   n_rows=nrows_in_each_fig,
                                                   n_cols=ncols_in_each_fig
                                                   )
                        with tarfile.open(df.loc[acc, 'output_tarfile_path'], mode='w:gz') as tar_out:
                            #calculate ratio of Amino acid Identity Membranous Over Nonmembranous  (AAIMON ratio)                 
                            for TMD in list_of_TMDs:
                                len_query_TMD = len(df.loc[acc, '%s_seq' % TMD])
                                #following the general filters, filter to only analyse sequences with TMD identity above cutoff, and a nonTMD_perc_ident above zero ,to avoid a divide by zero error
                                min_identity_of_TMD_initial_filter = settings['simap_match_filters']['min_identity_of_TMD_initial_filter']
                                dfs_filt_AAIMON = dfs_filt.loc[dfs['%s_perc_ident' % TMD] >= min_identity_of_TMD_initial_filter]
                                #avoid a divide by zero error in the unlikely case that there are no_identical_residues_in_alignment
                                dfs_filt_AAIMON = dfs_filt_AAIMON.loc[dfs_filt_AAIMON['nonTMD_perc_ident'] != 0]
                                #find the TMD number (starting from 1)                                
                                TMD_Nr = list_of_TMDs.index(TMD) + 1
                                #use the dictionary to obtain the figure number, plot number in figure, plot indices, etc
                                newfig, savefig, Plot_Nr, plot_nr_in_fig, row_nr, col_nr = dict_organising_subplots[TMD_Nr]
                                #if the TMD is the last one, the figure should be saved
                                if TMD_Nr == len(list_of_TMDs):
                                    savefig = True
                                #if a new figure should be created (either because the orig is full, or the last TMD is analysed)
                                if newfig:
                                    #create a new figure
                                    fig, axarr = plt.subplots(nrows=nrows_in_each_fig,
                                                              ncols=ncols_in_each_fig)  # sharex=True
                                #create numpy array of membranous over nonmembranous conservation ratios (identity)
                                hist_data_I = np.array(dfs_filt_AAIMON['%s_AAIMON_ratio' % TMD].dropna())
                                #use numpy to create a histogram
                                freq_counts, bin_array = np.histogram(hist_data_I, bins=binlist)
                                #assuming all of the bins are exactly the same size, make the width of the column equal to 70% of each bin
                                col_width = float('%0.3f' % (0.7 * (bin_array[1] - bin_array[0])))
                                #when align='center', the central point of the bar in the x-axis is simply the middle of the bins ((bin_0-bin_1)/2, etc)
                                centre_of_bar_in_x_axis = (bin_array[:-2] + bin_array[1:-1]) / 2
                                #add the final bin, which is physically located just after the last regular bin but represents all higher values
                                centre_of_bar_in_x_axis = np.append(centre_of_bar_in_x_axis,
                                                                    centre_of_bar_in_x_axis[-1] +
                                                                    centre_of_bar_in_x_axis[0]
                                                                    )
                                barcontainer_I = axarr[row_nr, col_nr].bar(left=centre_of_bar_in_x_axis,
                                                                           height=freq_counts, align='center',
                                                                           width=col_width, color="#0489B1",
                                                                           alpha=0.5
                                                                           )  # edgecolor='black',
                                #create numpy array of membranous over nonmembranous conservation ratios (identity + similarity)
                                hist_data_S = np.array(dfs_filt_AAIMON['%s_AASMON_ratio' % TMD].dropna())
                                #use numpy to create a histogram
                                freq_counts, bin_array = np.histogram(hist_data_S, bins=binlist)
                                #create a line graph rather than a bar graph for the AAISON (ident + similarity)
                                linecontainer_S = axarr[row_nr, col_nr].plot(centre_of_bar_in_x_axis, freq_counts,
                                                                             color="#0101DF", alpha=0.5)
                                #other colours that are compatible with colourblind readers: #8A084B Dark red, #B45F04 deep orange, reddish purple #4B088A 
                                #http://html-color-codes.info/
                                #label the x-axis for each plot, based on the TMD
                                axarr[row_nr, col_nr].set_xlabel('%s conservation ratio (membranous over nonmembranous)' % TMD, 
                                                                 fontsize=fontsize)
                                if savefig:
                                    #take x-axis min from settings
                                    xlim_min = settings["hist_settings_single_protein"]["smallest_bin"]
                                    #take x-axis max from settings
                                    xlim_max = settings["hist_settings_single_protein"]["largest_bin"]
                                    #apply the following formatting changes to all plots in the figure
                                    for ax in axarr.flat:
                                        #set x-axis min
                                        ax.set_xlim(xlim_min, xlim_max)
                                        #set x-axis ticks
                                        #use the slide selection to select every second item in the list as an xtick(axis label)
                                        ax.set_xticks([float('%0.1f' % c) for c in centre_of_bar_in_x_axis[::3]])
                                        ax.set_ylabel('freq', rotation='vertical', fontsize=fontsize)
                                        #change axis font size
                                        ax.tick_params(labelsize=fontsize)
                                        #create legend?#http://stackoverflow.com/questions/9834452/how-do-i-make-a-single-legend-for-many-subplots-with-matplotlib
                                        ax.legend(['AASMON (identity + similarity)', 'AAIMON (identity)'],
                                                  loc='upper right', fontsize=fontsize)
                                        #add background grid
                                        ax.grid(True, color='0.75', alpha=0.5)
                                    #automatically tighten the layout of plots in the figure
                                    fig.tight_layout()
                                    #save files
                                    fig.savefig(df.loc[acc, 'AAIMON_hist_BASENAMEPATH'] + '_%01d.png' % Plot_Nr,
                                                format='png', dpi=200)
                                    fig.savefig(df.loc[acc, 'AAIMON_hist_BASENAMEPATH'] + '_%01d.pdf' % Plot_Nr,
                                                format='pdf')
                                    #close figure
                                    plt.close('all')
                                    #add to tarfile
                                    tar_out.add(df.loc[acc, 'AAIMON_hist_BASENAMEPATH'] + '_%01d.png' % Plot_Nr,
                                                arcname=df.loc[acc, 'AAIMON_hist_BASENAME'] + '_%01d.png' % Plot_Nr)
                                    tar_out.add(df.loc[acc, 'AAIMON_hist_BASENAMEPATH'] + '_%01d.pdf' % Plot_Nr,
                                                arcname=df.loc[acc, 'AAIMON_hist_BASENAME'] + '_%01d.pdf' % Plot_Nr)
                                    #delete files
                                    #os.remove(df.loc[acc,'AAIMON_hist_BASENAMEPATH'] + '_%01d.png' % Plot_Nr, format='png', dpi=200)
                                    os.remove(df.loc[acc, 'AAIMON_hist_BASENAMEPATH'] + '_%01d.pdf' % Plot_Nr)


                                #start with the same dataframe copy that has filtered for gapped identity, etc (the AIMAN ratio is not necessary for the FastA saving)
                                #remove hits lacking sequence and also remove hits with too many gaps in TMD from either query or match
                                dfs_filt_FastA = dfs_filt.loc[dfs['%s_SW_match_seq' % TMD].notnull()].loc[
                                    dfs['%s_SW_query_acceptable_num_gaps' % TMD]].loc[
                                    dfs['%s_SW_match_acceptable_num_gaps' % TMD]]
                                #setup the file names again. Note that the file should already exist, and the query sequence included.
                                fasta_file = df.loc[acc, 'fasta_file_BASENAME'] + '%s.fas' % TMD
                                fasta_file_path = df.loc[acc, 'fasta_file_BASENAMEPATH'] + '%s.fas' % TMD
                                #setup filenames for fastA plus surrounding sequence (interfacial region)
                                fasta_file_plus_surr = df.loc[acc, 'fasta_file_plus_surr_path_BASENAME'] + '%s.fas' % TMD
                                fasta_file_plus_surr_path = df.loc[acc, 'fasta_file_plus_surr_path_BASENAMEPATH'] + '%s.fas' % TMD

                                with open(fasta_file_path, 'w') as f:
                                    #add original query seq to fasta file. Note that the first SIMAP hit is excluded below.
                                    f.write('>0000_%s_%s_uniprot_query\n%s\n' % (
                                        df.loc[acc, 'A2_protein_name'], TMD, df.loc[acc, '%s_seq' % TMD]))
                                    #for each hit, add the sequence to the fastA file
                                    for hit in dfs_filt_FastA.loc[dfs['%s_SW_match_seq' % TMD].notnull()].index[1:]:
                                        #add the original query seq
                                        f.write('>%04d_%s_%s\n%s\n' % (hit,str(dfs_filt_FastA.loc[hit, 'A2_organism'])[:30],
                                                                       str(dfs_filt_FastA.loc[hit, 'A4_description'])[:30],
                                                                       dfs_filt_FastA.loc[hit, '%s_SW_match_seq' % TMD]))
                                        #logging.info('saved ' + fasta_file_path)
                                tar_out.add(fasta_file_path, arcname=fasta_file)
                                os.remove(fasta_file_path)

                                dfs_filt_FastA_plus_surr = dfs_filt.loc[dfs['%s_SW_match_seq_plus_surr' % TMD].notnull()]
                                with open(fasta_file_plus_surr_path, 'w') as f:
                                    #add the original query seq
                                    f.write('>00_%s_query_seq\n%s\n' % (df.loc[acc, 'A2_protein_name'], 
                                            df.loc[acc, '%s_with_surrounding_seq' % TMD]))
                                    for hit in dfs_filt_FastA.loc[dfs_filt_FastA['%s_SW_match_seq_plus_surr' % TMD].notnull()].index:
                                        f.write('>%04d_%s_%s\n%s\n' % (hit,str(dfs_filt_FastA.loc[hit, 'A2_organism'])[:30],
                                                                       str(dfs_filt_FastA.loc[hit, 'A4_description'])[:30],
                                                                       dfs_filt_FastA.loc[hit, '%s_SW_match_seq_plus_surr' % TMD]))
                                        #logging.info('saved ' + fasta_file_plus_surr_path)
                                tar_out.add(fasta_file_plus_surr_path, arcname=fasta_file_plus_surr)
                                os.remove(fasta_file_plus_surr_path)

                            #remove columns to make output csv smaller
                            if settings['run_settings']['calculate_AAIMON_ratios']['drop_columns_to_reduce_csv_filesize']:
                                dfs = dfs.drop(['match_alignment_sequence', 'query_alignment_sequence', 'alignment_markup',
                                     'nonTMD_seq_query', 'nonTMD_markup'], axis=1)
                            dfs.to_csv(df.loc[acc, 'SIMAP_csv_analysed_path'], sep=",", quoting=csv.QUOTE_NONNUMERIC)
                            tar_out.add(df.loc[acc, 'SIMAP_csv_analysed_path'], arcname=df.loc[acc, 'SIMAP_csv_analysed'])
                        #delete original uncompressed file
                        os.remove(df.loc[acc, 'SIMAP_csv_analysed_path'])
                        df.loc[acc, 'num_hits_with_SW_align_node'] = dfs['hit_contains_SW_node'].value_counts()[True]
                        logging.info('num_hits_with_SW_align_node: \n%s' % df.loc[acc, 'num_hits_with_SW_align_node'])
                        df.loc[acc, 'num_FastA_seqs_saved'] = int(dfs_filt_FastA.shape[0])
                        logging.info('num_FastA_seqs_saved; %i\n' % df.loc[acc, 'num_FastA_seqs_saved'])
                    #save to csv after each protein is analysed, incrementally adding the extra data
                    with open(dfout08_simap_AAIMON, 'w') as csv_out:
                        df.to_csv(csv_out, sep=",", quoting=csv.QUOTE_NONNUMERIC)
                        
    logging.info('calculate_AAIMON_ratios is finished.')

    
#A## variables are included only to help navigate the document in PyCharm
A08a_calculate_gap_densities = settings["run_settings"]["calculate_gap_densities"]["run"]

if A08a_calculate_gap_densities:
    logging.info('~~~~~~~~~~~~starting calculate_gap_densities~~~~~~~~~~~~')
    
    # If script previously has been run, continues with proteins not beeing processed yet, or overwrites previous gap analysis 
    overwrite_previous_gap_analysis = settings["run_settings"]["calculate_gap_densities"]["overwrite_previous_gap_analysis"]     
    
    # Maximum number of gaps for tmds to be considered 
    allowed_gaps_per_tmd = settings["run_settings"]["calculate_gap_densities"]["allowed_gaps_per_tmd"]

    # 24 for beta barrel proteins, can be altered if only several TMDs to consider 
    max_number_of_tmds = settings["run_settings"]["calculate_gap_densities"]["max_number_of_tmds"]

    #test if the dataframe has already been created, otherwise re-open from uniprot csv file 
    try:
        logging.info('first protein acc = %s, df already exists,'
                     'continuing with parse_list_proteins_to_csv_and_fasta' % df.iloc[0][0])
    except NameError:
        if os.path.isfile(dfout10_uniprot_gaps):
            df = pd.read_csv(dfout10_uniprot_gaps, sep=",", quoting=csv.QUOTE_NONNUMERIC, index_col=[0])
            logging.info('df loaded from %s' % dfout10_uniprot_gaps)        
        else:
            df = pd.read_csv(dfout08_simap_AAIMON, sep=",", quoting=csv.QUOTE_NONNUMERIC, index_col=[0])
            logging.info('no gap file found, df loaded from %s' % dfout08_simap_AAIMON)  

            # If no previous analysis had been done, uniprot csv is opened and additional columns are created 
            df["gaps_analysed"]= np.nan     # Column, containing true or false
            for n in range (1,max_number_of_tmds+1):
                df["TM%.2d_occuring_gaps"%n]=np.nan     # List of unique gappositions, which occur in the tmd                 
                df["TM%.2d_amount_possible_gappositions"%n]=np.nan      
                df['total_amount_of_TM%.2d'%n] = np.nan     # How often TMD is considered 
                df['juxta_TM%.2d_intracellular_possible_gappositions'%n] = np.nan   # List of unique gappositions, which occur in the intracellular loop-part/end
                df['juxta_TM%.2d_extracellular_possible_gappositions'%n] = np.nan   # List of unique gappositions, which occur in the intracellular loop-part/end
                df['juxta_TM%.2d_intracellular_num_gaps'%n] = np.nan
                df['juxta_TM%.2d_exracellular_num_gaps'%n] = np.nan
                df['len_juxta_TM%.2d_intracellular'%n] = np.nan
                df['len_juxta_TM%.2d_extracellular'%n] = np.nan

        for acc in df.index:

            # The next steps (the main analysis) is only executed, if previous analysis can be overwritten or no analysis has yet been done 
            if (overwrite_previous_gap_analysis == True) or (df.loc[acc,"gaps_analysed"] != True):            
                logging.info("%s"%acc) 
                list_of_TMDs = ast.literal_eval(df.loc[acc,"list_of_TMDs"])

                # Checks if outputfiles (tar) exist 
                if os.path.exists("/nas/teeselab/students/rimma/databases/simap/%s/%s_outputfiles.tar.gz"%(acc[0:2],acc)):

                # opens the analysed csv for each protein and loads it into a dataframe 
                    with tarfile.open("/nas/teeselab/students/rimma/databases/simap/%s/%s_outputfiles.tar.gz"%(acc[0:2],acc), mode= 'r:gz')as tar:
                    
                    # checks if the analysed file exists, otherwise prints that it does not exist
                        if "%s_analysed.csv"%acc in tar.getnames():
                    
                        # loads file into analysed csv 
                            analysed_csv = tar.extractfile('%s_analysed.csv'%acc)
                            analysed_df = pd.read_csv(analysed_csv,low_memory=False,index_col=[0])
           
                            # checks if first amino acid is located inside (or periplasmatic) or outside, returns a boolean, true or false
                            # if first residue is located inside, every even tmd (tmd2,tmd4,tmd6...) is reversed, otherwise every odd tmd is reversed 
                            # output is a boolean for each tmd, depending on the number and on the first amino acid
                            
                            if df.n_term_ec[acc] == False:
                                reverse_tmd = False
                            else:
                                reverse_tmd = True                            
                                       
           
                            # for each TMD in the proteins, creates new lists which will contain gappositions, lists are saved in a column and created newly for each tmd  
                            for tmd in list_of_TMDs:
                                print(tmd)
                                tmd_int = int(tmd[-2:]) # Integer of TMD number 
                                len_of_query = len(analysed_df["%s_SW_query_seq"%tmd][1]) # Length of first query sequence, which does (usually) not contain any gaps
                                len_of_query_reversed= ((1/len_of_query)+1) # Reversed length, important if TMD needs to be reversed afterwards 
                                list_of_gaps_in_tmd = []
                                list_of_gaps_intracellular = []
                                list_of_gaps_extracellular = []
        
 
                                for hit in analysed_df.index:
                                                               
        #'''
        #Start of the main gap analysis 
        #Code searches for "-" in the TMD sequence and returns the index!! (not the position)
        #'''            
                    
                                    # Following if conditions only refer to gaps in the query!
                                    # Query gaps are counted as "in between positions", for example: 4,5 refers to a gap between position 4 and 5;
                                    # if two gaps occur one after another: only one position (between two amino acids is considered)
                                                
                                    # Filter to make sure, that there are 1 or 2 gaps in the query sequence and up to the max allowed gaps in the match
                                    if (analysed_df["%s_SW_query_num_gaps"%tmd][hit] != 0.0) and (analysed_df["%s_SW_query_num_gaps"%tmd][hit] <= 2.0)\
                                        and (analysed_df["%s_SW_match_num_gaps"%tmd][hit] <= int("%s"%allowed_gaps_per_tmd)):
                                                
                                        # Stores the endpoints in a temp list; endpoints are used, to switch from python indices to numbers 
                                        list_of_gaps_per_hit_in_query = [m.start() for m in re.finditer("-",analysed_df.loc[hit,"%s_SW_query_seq"%tmd]) if m.start()]
                                        print (list_of_gaps_per_hit_in_query)
                                        
                                        # if there are two gaps in the query (and 2 allowed), code checks if they are side by side (difference of 1) 
                                        # and appends the gap position, else appends both gap positions
                                        # 0.5 is substracted in order to make them "in between" position; 
                                        #if two gaps are observed, 1.5 is substracted from the second one, since the residue positions are moved due to the first gap
                                        if len(list_of_gaps_per_hit_in_query) == 2 and allowed_gaps_per_tmd==2:
                                            
                                            if list_of_gaps_per_hit_in_query[1]- list_of_gaps_per_hit_in_query[0] ==1:
                                                list_of_gaps_in_tmd.append(list_of_gaps_per_hit_in_query[0]-0.5)
                                                
                                            else: 
                                                list_of_gaps_in_tmd.append(list_of_gaps_per_hit_in_query[0]-0.5)
                                                list_of_gaps_in_tmd.append(list_of_gaps_per_hit_in_query[1]-1.5)
                                                
                                        # if there is only one gap in query or only one gap is allowed, it appends the first (and only) gap from the list_of_gaps_per_hit_in_query to the list_of_TMDs       
                                        else:
                                            if len(list_of_gaps_per_hit_in_query) == 1:
                                                list_of_gaps_in_tmd.append (list_of_gaps_per_hit_in_query[0]-0.5)

                                  
                                    # Following if conditions only refer to gaps in the match!
                                    # Query gaps are counted as deletions of positions; for example: 4 refers to a gap on position 4;
                                    # if two gaps occur one after another, both are considered since two actual amino acids from the original query are deleted
                                    # Since the gap positions are dependend on the query sequence, query-gap positions in the same alignment have to be considered as well 
                                      
                                    
                                    # Filter to make sure, that there are 1 or 2 gaps in the match sequence and up to the max allowed gaps in the query
                                    if (analysed_df["%s_SW_query_num_gaps"%tmd][hit] <=2.0) and (analysed_df["%s_SW_match_num_gaps"%tmd][hit] <= 2.0)\
                                        and (analysed_df["%s_SW_match_num_gaps"%tmd][hit] != 0.0):
                                        
                                        # It's not sure that the list of hits in query was already determined, maybe there were no gaps, anyway here it is important how many 
                                        list_of_gaps_per_hit_in_query = [m.start() for m in re.finditer("-",analysed_df.loc[hit,"%s_SW_query_seq"%tmd]) if m.start()]
                                        list_of_gaps_per_hit_in_match = [m.start() for m in re.finditer("-",analysed_df.loc[hit,"%s_SW_match_seq"%tmd])if m.start()]
                                        #print (list_of_gaps_per_hit_in_query)
                                        #print (list_of_gaps_per_hit_in_match)
                                        
                                        if len(list_of_gaps_per_hit_in_match)>0:
                                            for n in list(reversed(list_of_gaps_per_hit_in_match)):
                                                substracted_value = len([m<n for m in list_of_gaps_per_hit_in_query])
                                                list_of_gaps_in_tmd.append(abs(n-substracted_value))                            

        #######
        # Start of the Juxta Consideration
        # In the case of n_term being located intracellular: 
        # there are 4 groups: 1. juxta_before_odd_TMDs + 2.juxta_after_even_TMDs 3. Juxta_before_even_TMDs + 4. Juxta_after_odd_TMDs
        # 1 + 2 --> Intracellular 
        # 3 + 4 --> Extracellular           
        # If the n_term is extracellular, that it's the other way round. 1+2 --> Extracellular 3+4 --> Intracellular                    
        ### The data will already be flipped in order to align extracellular and intracellular parts, extracellular: + , intracellular: -                                                       

                                    # juxta before_odd_TMDs:
                                    
                                    if r_utils.isOdd(tmd_int)==True:  # also für 1 , 3 ...
                                      
                                        
                                        # makes sure that the search is done in a string
                                        if type(analysed_df.loc[hit,"seq_juxta_before_TM%.2d_in_query"%tmd_int])== str:
                               
                                            # list of gap indices 
                                            list_of_gaps_in_query_before_odd = [m.start()+0.5 for m in re.finditer("-",analysed_df.loc[hit,"seq_juxta_before_TM%.2d_in_query"%tmd_int][::-1])if m.start()+0.5 < 31]
                                            
                                            # if one gap is found, code checks location and appends it
                                            if len (list_of_gaps_in_query_before_odd)==1:
                                                if reverse_tmd == False:                                      
                                                    list_of_gaps_intracellular.append(list_of_gaps_in_query_before_odd[0])   
                                                else: 
                                                    list_of_gaps_extracellular.append(list_of_gaps_in_query_before_odd[0])
                                            # if more than one gap is found, code checks if the gaps are one after another in the query (therefore the following_gap is used)!
                                            if len (list_of_gaps_in_query_before_odd)>1.0:
                                                following_gap = 0
                                                rev_value = list_of_gaps_in_query_before_odd[0]
                                                for n in list_of_gaps_in_query_before_odd:
                                                    if n-following_gap == rev_value:
                                                        if reverse_tmd == False: 
                                                            list_of_gaps_intracellular.append(n-following_gap)  
                                                            
                                                            following_gap = following_gap+1
                                                        else: 
                                                            list_of_gaps_extracellular.append(n-following_gap)                                                   
                                                            following_gap = following_gap+1 
                                                    else: 
                                                        if reverse_tmd == False:
                                                            list_of_gaps_intracellular.append(n-following_gap)
                                                            following_gap = following_gap+1
                                                            rev_value = n 
                                                        else: 
                                                            list_of_gaps_extracellular.append(n-following_gap)
                                                            following_gap = following_gap+1
                                                            rev_value = n   
                                            
                                            
                                            if type(analysed_df.loc[hit,"seq_juxta_before_TM%.2d_in_match"%tmd_int])== str:
                                            # Makes a list of gaps of the match of the odd juxta before the TMD
                                                list_of_gaps_in_query_before_odd = [m.start()+1 for m in re.finditer("-",analysed_df.loc[hit,"seq_juxta_before_TM%.2d_in_query"%tmd_int][::-1])if m.start() < 32]

                                                list_of_gaps_in_match_before_odd = [m.start()+1 for m in re.finditer("-",analysed_df.loc[hit,"seq_juxta_before_TM%.2d_in_match"%tmd_int][::-1])if m.start() < 32]
                                               # This step is essential, to control, if there is a gap before, in the query region 
                                                for n in list(reversed(list_of_gaps_in_match_before_odd)):
                                                    greater_values = sum(i< n for i in list_of_gaps_in_query_before_odd)    
                                                    if reverse_tmd== False:
                                                        list_of_gaps_intracellular.append(n-greater_values)
                                                                                                    
                                                    else:
                                                        list_of_gaps_extracellular.append(n-greater_values)
                                                                                                
# juxta after odd TMDs:
                                        
                                            if type(analysed_df.loc[hit,"seq_juxta_after_TM%.2d_in_query"%tmd_int])== str:

                                                list_of_gaps_in_query_after_odd = [m.start()+0.5 for m in re.finditer("-",analysed_df.loc[hit,"seq_juxta_after_TM%.2d_in_query"%tmd_int])if m.start()+0.5 < 31]

                                                # if one gap is found, code checks location and appends it
                                                if len (list_of_gaps_in_query_after_odd)==1:
                                                    if reverse_tmd == False:                                      
                                                        list_of_gaps_extracellular.append(list_of_gaps_in_query_after_odd[0])                                            
                                                    else: 
                                                        list_of_gaps_intracellular.append(list_of_gaps_in_query_after_odd[0])
                                                                                                    
                                                # if more than one gap is found, code checks if the gaps are one after another in the query!
                                                if len (list_of_gaps_in_query_after_odd)>1.0:
                                                    following_gap = 0
                                                    rev_value = list_of_gaps_in_query_after_odd[0]
                                                    for n in list_of_gaps_in_query_after_odd:
                                                        if n+following_gap == rev_value:
                                                            if reverse_tmd == False: 
                                                                list_of_gaps_extracellular.append(n-following_gap)                                                   
                                                                following_gap = following_gap+1
                                                            else: 
                                                                list_of_gaps_intracellular.append(n-following_gap)                                                   
                                                                following_gap = following_gap+1                                                         
                                                        else: 
                                                            if reverse_tmd == False:
                                                                list_of_gaps_extracellular.append(n-following_gap)
                                                                following_gap = following_gap+1
                                                                rev_value = n 
                                                            else: 
                                                                list_of_gaps_intracellular.append(n-following_gap)
                                                                following_gap = following_gap+1
                                                                rev_value = n         
       
                                                                
                                        else:  # for 2,4
                                    
                                        # juxta before even TMDs:
                                        
                                        # makes sure that the search is done in a string
                                            if type(analysed_df.loc[hit,"seq_juxta_before_TM%.2d_in_query"%tmd_int])== str:
                                                
                                                # list of gap indices 
                                                list_of_gaps_in_query_before_even = [m.start()+0.5 for m in re.finditer("-",analysed_df.loc[hit,"seq_juxta_before_TM%.2d_in_query"%tmd_int])if m.start()+0.5 < 31]
                                                
                                                # if one gap is found, code checks location and appends it
                                                if len (list_of_gaps_in_query_before_even)==1:
                                                    if reverse_tmd == False:                                      
                                                        list_of_gaps_extracellular.append(list_of_gaps_in_query_before_even[0])                                            
                                                    else: 
                                                        list_of_gaps_intracellular.append(list_of_gaps_in_query_before_even[0])
                                                                                                    
                                                # if more than one gap is found, code checks if the gapy are one after another in the query!
                                                if len (list_of_gaps_in_query_before_even)>1.0:
                                                    following_gap = 0
                                                    rev_value = list_of_gaps_in_query_before_even[0]
                                                    for n in list_of_gaps_in_query_before_even:
                                                        if n+following_gap == rev_value:
                                                            if reverse_tmd == False: 
                                                                list_of_gaps_extracellular.append(n-following_gap)                                                   
                                                                following_gap = following_gap+1
                                                            else: 
                                                                list_of_gaps_intracellular.append(n-following_gap)                                                   
                                                                following_gap = following_gap+1                                                         
                                                        else: 
                                                            if reverse_tmd == False:
                                                                list_of_gaps_extracellular.append(n-following_gap)
                                                                following_gap = following_gap+1
                                                                rev_value = n 
                                                            else: 
                                                                list_of_gaps_intracellular.append(n-following_gap)
                                                                following_gap = following_gap+1
                                                                rev_value = n                         
                                                                    
                                            if type(analysed_df.loc[hit,"seq_juxta_before_TM%.2d_in_match"%tmd_int])== str:                
                                            # Makes a list of gaps of the match of the odd juxta before the TMD
                                                list_of_gaps_in_query_before_even = [m.start()+0.5 for m in re.finditer("-",analysed_df.loc[hit,"seq_juxta_before_TM%.2d_in_query"%tmd_int])if m.start()+0.5 < 31]

                                                list_of_gaps_in_match_before_even = [m.start()+1 for m in re.finditer("-",analysed_df.loc[hit,"seq_juxta_before_TM%.2d_in_match"%tmd_int])if m.start() < 31]

                                                for n in list(reversed(list_of_gaps_in_match_before_even)):
                                                    greater_values = sum(i< n for i in list_of_gaps_in_query_before_even)    
                                                    if reverse_tmd== False:
                                                        list_of_gaps_extracellular.append(n-greater_values)
                                                    else:
                                                        list_of_gaps_intracellular.append(n-greater_values)
                                            
                                           
                                           
                                            # juxta after even TMDs:
                                            
                                            if type(analysed_df.loc[hit,"seq_juxta_after_TM%.2d_in_query"%tmd_int])== str:

                                                list_of_gaps_in_query_after_even = [m.start()+0.5 for m in re.finditer("-",analysed_df.loc[hit,"seq_juxta_after_TM%.2d_in_query"%tmd_int][::-1]) if m.start()+0.5 < 31]
                                                
                                                # if one gap is found, code checks location and appends it
                                                if len (list_of_gaps_in_query_after_even)==1:
                                                    if reverse_tmd == False:                                      
                                                        list_of_gaps_intracellular.append(list_of_gaps_in_query_after_even[0])                                            
                                                    else: 
                                                        list_of_gaps_extracellular.append(list_of_gaps_in_query_after_even[0])
                                                                                                    
                                                # if more than one gap is found, code checks if the gaps are one after another in the query!
                                                if len (list_of_gaps_in_query_after_even)>1.0:
                                                    following_gap = 0
                                                    rev_value = list_of_gaps_in_query_after_even[0]
                                                    for n in list_of_gaps_in_query_after_even:
                                                        if n-following_gap == rev_value:
                                                            if reverse_tmd == False: 
                                                                list_of_gaps_intracellular.append(n-following_gap)                                                   
                                                                following_gap = following_gap+1
                                                            else: 
                                                                list_of_gaps_extracellular.append(n-following_gap)                                                   
                                                                following_gap = following_gap+1                                                         
                                                        else: 
                                                            if reverse_tmd == False:
                                                                list_of_gaps_intracellular.append(n-following_gap)
                                                                following_gap = following_gap+1
                                                                rev_value = n 
                                                            else: 
                                                                list_of_gaps_extracellular.append(n-following_gap)
                                                                following_gap = following_gap+1
                                                                rev_value = n                              
                                                                    
                                                           # Makes a list of gaps of the match of the odd juxta before the TMD
                                            
                                            if type(analysed_df.loc[hit,"seq_juxta_after_TM%.2d_in_match"%tmd_int])== str and type(analysed_df.loc[hit,"seq_juxta_after_TM%.2d_in_query"%tmd_int])==str:
                                            
                                                list_of_gaps_in_query_after_even = [m.start()+0.5 for m in re.finditer("-",analysed_df.loc[hit,"seq_juxta_after_TM%.2d_in_query"%tmd_int][::-1])if m.start()+0.5 < 31]
                                                list_of_gaps_in_match_after_even = [m.start()+1 for m in re.finditer("-",analysed_df.loc[hit,"seq_juxta_after_TM%.2d_in_match"%tmd_int][::-1])if m.start() < 31]
                                         
                                                for n in list(reversed(list_of_gaps_in_match_after_even)):
                                                    greater_values = sum(i< n for i in list_of_gaps_in_query_after_even)    
                                                    if reverse_tmd== False:
                                                        list_of_gaps_intracellular.append(n-greater_values)
                                                        
                                                    else:
                                                        list_of_gaps_extracellular.append(n-greater_values)
                                    
                                    # sets of lists are created, to assure, that each gapposition contributes only once to the possible gap positions 
                                    unique_list_of_gaps_in_tmd = list(set(list_of_gaps_in_tmd))
                                    unique_list_of_gaps_intracellular = list(set(list_of_gaps_intracellular))
                                    unique_list_of_gaps_extracellular = list(set(list_of_gaps_extracellular))
                                
                                    # Saves the calculated lists into cells in the columns 
                                    df["%s_occuring_gaps"%tmd][acc]=str(unique_list_of_gaps_in_tmd)
                                    df["%s_amount_possible_gappositions"%tmd][acc]=len(unique_list_of_gaps_in_tmd)                             
                                    df['juxta_%s_intracellular_possible_gappositions'%tmd][acc] = str(unique_list_of_gaps_intracellular)
                                    df['juxta_%s_extracellular_possible_gappositions'%tmd][acc] = str(unique_list_of_gaps_extracellular)
                                    df['juxta_%s_intracellular_num_gaps'%tmd][acc] = len(unique_list_of_gaps_intracellular)
                                    df['juxta_%s_exracellular_num_gaps'%tmd][acc] = len(unique_list_of_gaps_extracellular)
                                            
                            # At the end, sets analysed to true, this is important to not overwrite                                     
                            df.gaps_analysed[acc] = "True"
                            logging.info("--Analysed")
                            with open(dfout10_uniprot_gaps, 'w') as csv_out:
                                df.to_csv(csv_out, sep=",", quoting=csv.QUOTE_NONNUMERIC)
                                                       
                        else:
                            logging.info("Analysed csv for %s does not exist" %acc)
                else: 
                    logging.info("Output file for %s does not exist" %acc)

            else:
                logging.info("Gap analysis for %s already done" %acc)                    
                                                                    

#A## variables are included only to help navigate the document in PyCharm
A08b_calculate_gap_densities = settings["run_settings"]["create_graph_of_gap_density"]["run"]

if A08b_calculate_gap_densities:
    logging.info('~~~~~~~~~~~~starting creating graphs of gap density~~~~~~~~~~~~')

    # Tests if values to process do exist, otherwise prints out an error, that analysis has not been done yet

    if os.path.isfile(dfout10_uniprot_gaps):
        df = pd.read_csv(dfout10_uniprot_gaps, sep=",", quoting=csv.QUOTE_NONNUMERIC, index_col=[0])
        logging.info('df loaded from %s' % dfout10_uniprot_gaps)        
    else:
        logging.info('no analysis has been done yet, please set calculate gap densities to true')  
        
    # Creates  variables, which are defined in the settings file
    amount_of_bins_in_tmd_region = settings["run_settings"]["create_graph_of_gap_density"]["amount_of_bins_in_tmd_region"]
    amount_of_bins_in_juxta_regions = settings["run_settings"]["create_graph_of_gap_density"]["amount_of_bins_in_juxta_regions"]
    labeling_of_tmd_region = settings["run_settings"]["create_graph_of_gap_density"]["labeling_of_tmd_region"]
    labeling_of_intracellular_region = settings["run_settings"]["create_graph_of_gap_density"]["labeling_of_intracellular_region"]
    labeling_of_outer_region = settings["run_settings"]["create_graph_of_gap_density"]["labeling_of_outer_region"]
    fontsize = settings["run_settings"]["create_graph_of_gap_density"]["fontsize_of_labels"]
    
    # Creates to empty lists, which contain the gap positions of tmd regions, depending on the orientation in the membrane 
    to_reverse = []
    not_to_reverse = [] 
        
    for acc in df.index:      
        print(acc)   
        
        # Checks for gap positions in TMD 1,3,5,7... and appends it to the right lists (depending on the orientation of the tmd), if the tmd needs to be reversed
        # the "reversed value" is appended
        for num_TMD in range(1,int(df.loc[acc,"number_of_TMDs_in_uniprot_feature_list"])+1,2):

            if not r_utils.isNaN(df.loc[acc,"TM%.2d_occuring_gaps"%num_TMD]):
                for n in ast.literal_eval(df.loc[acc,"TM%.2d_occuring_gaps"%num_TMD]):                                
                    if df.loc[acc,"n_term_ec"] == False:
                        not_to_reverse.append ((n/(len(df.loc[acc,"TM%.2d_seq"%num_TMD])-1))*amount_of_bins_in_tmd_region)
                    if df.loc[acc,"n_term_ec"] == True:
                        to_reverse.append(amount_of_bins_in_tmd_region-((n/(len(df.loc[acc,"TM%.2d_seq"%num_TMD])-1))*amount_of_bins_in_tmd_region))
            

        # Does the same for gap positions in TMD 2,4,6,8...                
        for num_TMD in range(2,int(df.loc[acc,"number_of_TMDs_in_uniprot_feature_list"])+1,2):   
            if not r_utils.isNaN(df.loc[acc,"TM%.2d_occuring_gaps"%num_TMD]):
                for n in ast.literal_eval(df.loc[acc,"TM%.2d_occuring_gaps"%num_TMD]):                  
                    if df.loc[acc,"n_term_ec"] == False:                        
                        to_reverse.append (amount_of_bins_in_tmd_region-((n/(len(df.loc[acc,"TM%.2d_seq"%num_TMD])-1))*amount_of_bins_in_tmd_region))                      
                    if df.loc[acc,"n_term_ec"] == True:                         
                        not_to_reverse.append((n/(len(df.loc[acc,"TM%.2d_seq"%num_TMD])-1))*amount_of_bins_in_tmd_region)

        # Merging the two lists together
        tmd_gaps = to_reverse+not_to_reverse
        
        # Creating a plot
        fig, axarr = plt.subplots()
        
        # Creating the tmd region, which is represented as a histogram; the amount of bins can be selected in the settings file
        n, bins, patches = plt.hist(tmd_gaps,bins=amount_of_bins_in_tmd_region,range=(0,amount_of_bins_in_tmd_region),normed=True,histtype='bar',color="darkblue",zorder=3, edgecolor = "white")

        # Labeling and adjusting graph 
        axarr.set_ylabel("Propability for gapposition",fontsize=fontsize,color="black")
        axarr.set_xlabel("Depth of Residue",fontsize=fontsize,color="black")
        axarr.set_xlim(0,amount_of_bins_in_tmd_region)                                   
        axarr.set_xticks(np.linspace(0.8,19.2,2))
        axarr.set_xticklabels((labeling_of_intracellular_region,labeling_of_outer_region),fontsize=fontsize,color="black")
        axarr.tick_params(reset=True, labelsize=fontsize,color="black")
        for tic in axarr.xaxis.get_major_ticks():
            tic.tick1On = tic.tick2On = False
        axarr.grid(False, 'minor', color='0.99', linestyle='-', linewidth=0.7)
        axarr.spines['top'].set_visible(False)
        axarr.spines['right'].set_visible(False)
        axarr.yaxis.tick_left()
        axarr.patch.set_facecolor('white')
        axarr.yaxis.grid(True,zorder = 0,color = "grey",linestyle=":")      
        
        # Saving the gap propensity 
        fig.savefig("/nas/teeselab/students/rimma/omp/summaries/List%.2d_gaps_TMD.png"%list_number, format='png', dpi=200)
                                            
                                                            
                                                            
          
        # data for intracellular part --> left 
        # Creates a nested list, of all lists, which contain gap positions, afterwards flattens it and convertes it to an array
        nested_list_of_gaps_intracellular = [ast.literal_eval(m) for n in range (1,25) for m in df['juxta_TM%.2d_intracellular_possible_gappositions'%n].dropna().tolist()]
        hist_data_juxta_intracellular = np.array (list(itertools.chain(*nested_list_of_gaps_intracellular)))

        # Checks for the smallest value (important for y-axis)
        min_value = int(abs(hist_data_juxta_intracellular.min()))
        

        # data for extracellular part --> right
        nested_list_of_gaps_extracellular = [ast.literal_eval(m) for n in range (1,25) for m in df['juxta_TM%.2d_extracellular_possible_gappositions'%n].dropna().tolist()]
        hist_data_juxta_extracellular = np.array (list(itertools.chain(*nested_list_of_gaps_extracellular)))

        # data for tmd --> middle
        # Creates an array out of the list
        hist_data_tmds = np.array (tmd_gaps)
        # times 2, because TMDs in gap and in query are considered! --> double amount
        total_amount_of_TMDs_in_protein = df.loc[df.gaps_analysed==True,"number_of_TMDs_in_uniprot_feature_list"].sum()*2

     
        # Creating final graph
        fig, axarr = plt.subplots()

      
        ### Intracellular part

        freq_counts_I, bin_array_I = np.histogram(hist_data_juxta_intracellular, bins = hist_data_juxta_intracellular.max())
        centre_of_bar_in_x_axis_I = -((bin_array_I[:-2] + bin_array_I[1:-1]) / 2)

        bar_width_I = centre_of_bar_in_x_axis_I[3] - centre_of_bar_in_x_axis_I[2]

        centre_of_bar_in_x_axis_I = np.append(centre_of_bar_in_x_axis_I, centre_of_bar_in_x_axis_I[-1] + bar_width_I)

        axarr.bar(left=centre_of_bar_in_x_axis_I, height= [((freq_counts_I.tolist()[n])/(r_utils.positionfreq_in_list_intra(n,hist_data_juxta_intracellular))) for n in range (0,len(freq_counts_I))],  width=0.4, color="mediumblue", linewidth = 0,zorder = 3)  # edgecolor='black',

        ######### TMD

        hist_data_tmds = np.array (tmd_gaps)
        freq_counts_II, bin_array_II = np.histogram(hist_data_tmds,bins = 20,range = (0,20))

        centre_of_bar_in_x_axis_II = (bin_array_II[:-2] + bin_array_II[1:-1]) / 2

        bar_width_II =  centre_of_bar_in_x_axis_II[3] - centre_of_bar_in_x_axis_II[2]

        centre_of_bar_in_x_axis_II = np.append(centre_of_bar_in_x_axis_II, centre_of_bar_in_x_axis_II[-1] + bar_width_II)

        axarr.bar(left=centre_of_bar_in_x_axis_II, height=[n/total_amount_of_TMDs_in_protein for n in freq_counts_II.tolist()], align='center', width=0.5, color="blue",linewidth =0,zorder=3)  # edgecolor='black',

        #####Extracellular part

        freq_counts_III, bin_array_III = np.histogram(hist_data_juxta_extracellular,bins = hist_data_juxta_extracellular.max())
        centre_of_bar_in_x_axis_III = ((bin_array_III[:-2] + bin_array_III[1:-1]) / 2)

        bar_width_III = centre_of_bar_in_x_axis_III[3] - centre_of_bar_in_x_axis_III[2]

        centre_of_bar_in_x_axis_III = np.append(centre_of_bar_in_x_axis_III, centre_of_bar_in_x_axis_III[-1] + bar_width_III)

        axarr.bar(left=centre_of_bar_in_x_axis_III+19.5, height= [((freq_counts_III.tolist()[n])/(r_utils.positionfreq_in_list_extra(n,hist_data_juxta_extracellular))) for n in range (0,len(freq_counts_III))],  width=0.4, color="mediumblue", linewidth = 0,zorder=3)  # edgecolor='black',



        ##### Style and labels
        
        axarr.set_xlim(-amount_of_bins_in_juxta_regions,amount_of_bins_in_juxta_regions+10)
        axarr.set_xlabel("residue position relative to TM domain",fontsize=fontsize)
        axarr.set_ylabel("gap propensity",color="mediumblue",fontsize=fontsize )
      
        # Setting the x-axis; Labels depend on selected size of juxta regions
        axarr.set_xticks([n for n in range(-amount_of_bins_in_juxta_regions,0,5)]+[-1,10,21]+ [n for n in range(25,amount_of_bins_in_juxta_regions+25,5)])
        labels = ["%s"%n for n in range(-amount_of_bins_in_juxta_regions,0,5)]+["-1","TM_domain","1"]+ ["%s"%n for n in range(5,amount_of_bins_in_juxta_regions+1,5)]
        axarr.set_xticklabels(labels,color="mediumblue")

        axarr.tick_params(reset=True, labelsize=fontsize,color="mediumblue")

        for tic in axarr.xaxis.get_major_ticks():
            tic.tick1On = tic.tick2On = False
            
        axarr.grid(False)
        axarr.yaxis.tick_left()

        axarr.patch.set_facecolor('white')
        axarr.yaxis.grid(True,zorder =0,linestyle = ":", color = "g")
        
        # for some reason it needs to be imported again 
        import matplotlib.patches as patches

        # Creating light background for TMD region and Annotations for labels
        axarr.add_patch(patches.Rectangle((0, 0),amount_of_bins_in_tmd_region,10,alpha=0.1,linewidth=None,facecolor="lightseagreen"))  
        axarr.annotate('Intracellular', xy = (0,0),xytext=(-20,1.7),alpha=0.5, fontsize=fontsize )
        axarr.annotate('Extracellular', xy = (0,0),xytext=(30,1.7),alpha=0.5,fontsize=fontsize )
        axarr.annotate('TM', xy = (10,0),xytext =(8.6,1.3),alpha=0.5,fontsize=fontsize )
        axarr.annotate('helix', xy = (10,0),xytext =(8.2,1.2),alpha=0.5,fontsize=fontsize )

        axarr.spines['left'].set_color('mediumblue')


        # Creating second y-axis 
        axarr2 = axarr.twinx()
        axarr2.spines['right'].set_color("black")
        axarr2.grid(False)
        axarr2.plot([n+20 for n in range(0,30)],[r_utils.positionfreq_in_list_extra(n,hist_data_juxta_extracellular) for n in range(0,30)],"black",linewidth=1.2)
        axarr2.plot([-n for n in range(0,30)],[r_utils.positionfreq_in_list_intra(n,hist_data_juxta_intracellular) for n in range (0,30)],"black",linewidth=1.2)
        axarr2.set_ylabel("Frequency of considered position in dataset",color="black",fontsize=fontsize)
        axarr2.yaxis.label.set_size(fontsize)
        axarr2.tick_params(axis='y', colors='black',labelsize=fontsize)
        axarr2.spines["left"].set_color("mediumblue")


        fig.savefig("/nas/teeselab/students/rimma/omp/summaries/List%.2d_gap_propensity.png"%list_number, format='png', dpi=200)


                                                            
'''+++++++++++++++ Summary figures describing the conservation ratios of proteins in the list ++++++++++++++++++'''
A09_save_figures_describing_proteins_in_list = settings["run_settings"]["save_figures_describing_proteins_in_list"]
if A09_save_figures_describing_proteins_in_list:

    
    backgroundcolour = 'white'
    #plt.style.use('ggplot')
    #add non-functional object to aid document navigation in some IDEs (e.g. Spyder)
    fig_title = ''
    
    '''
    Prepare subplots and default fontsizes etc
    '''
    #set default font size for plot
    fontsize = 4
    datapointsize = 2
    alpha = 0.1
    #use a dictionary to organise the saving of multiple plots in multiple figures, with a certain number of plots per figure
    n_plots_per_fig = 4
    nrows_in_each_fig = 2
    ncols_in_each_fig = 2
    dict_organising_subplots = utils.create_dict_organising_subplots(n_plots_per_fig=n_plots_per_fig,
                                                                     n_rows=nrows_in_each_fig,
                                                                     n_cols=ncols_in_each_fig)

    '''
    Prepare data for following figures
    '''
    
    #test if the dataframe has already been created, otherwise re-open from csv file containing the simap data
    try:
        logging.info('first protein acc = %s, df already exists, '
                     'continuing with save_figures_describing_proteins_in_list' % df.iloc[0][0])
    except NameError:
        logging.info('df loaded from %s' % dfout01_uniprotcsv)
        df = pd.read_csv(dfout08_simap_AAIMON, sep=",", quoting=csv.QUOTE_NONNUMERIC,index_col=0)
    #filter to remove sequences where no TMDs are found (will contain either np.nan, or 'nan') 
    df = df.loc[df['list_of_TMDs'].notnull()]
    df = df.loc[df['list_of_TMDs'] != 'nan']
    #iterate over the dataframe. Note that acc = uniprot accession here.    
    linspace_binlist = np.linspace(settings["hist_settings_mult_proteins"]["smallest_bin"],
                                   settings["hist_settings_mult_proteins"]["largest_bin"],
                                   settings["hist_settings_mult_proteins"]["number_of_bins"])

    #add 30 as the last bin, to make sure 100% of the data is added to the histogram, including major outliers
    binlist = np.append(linspace_binlist, settings["hist_settings_mult_proteins"]["final_highest_bin"])   
    
    for acc in df.index:
        dict_AAIMON_ratio_mean = {}
        dict_AAIMON_ratio_std = {}
        dict_AASMON_ratio_mean = {}
        dict_AASMON_ratio_std = {}
        for TMD in eval(df.loc[acc, 'list_of_TMDs']):
            dict_AAIMON_ratio_mean[TMD] = df.loc[acc, '%s_AAIMON_ratio_mean' % TMD]
            dict_AAIMON_ratio_std[TMD] = df.loc[acc, '%s_AAIMON_ratio_std' % TMD]
            dict_AASMON_ratio_mean[TMD] = df.loc[acc, '%s_AASMON_ratio_mean' % TMD]
            dict_AASMON_ratio_std[TMD] = df.loc[acc, '%s_AASMON_ratio_std' % TMD]
        df.loc[acc, 'AAIMON_ratio_mean_all_TMDs'] = np.mean(list(dict_AAIMON_ratio_mean.values()))
        df.loc[acc, 'AAIMON_ratio_std_all_TMDs'] = np.mean(list(dict_AAIMON_ratio_std.values()))
        df.loc[acc, 'AASMON_ratio_mean_all_TMDs'] = np.mean(list(dict_AASMON_ratio_mean.values()))
        df.loc[acc, 'AASMON_ratio_std_all_TMDs'] = np.mean(list(dict_AASMON_ratio_std.values()))
    '''
    Fig01: Histogram of mean AAIMON and AASMON ratios, SP vs MP
    '''
    Plot_Nr = 1
    #add non-functional "if" line to aid document navigation in some IDEs (e.g. Spyder)
    if fig_title == "Fig01: Histogram of mean AAIMON and AASMON ratios, SP vs MP":
        pass
    sys.stdout.write('Figures Processed' + str(Plot_Nr) + ', ')
    title = 'Mean ratios'
    newfig, savefig, Plot_Nr, plot_nr_in_fig, row_nr, col_nr = dict_organising_subplots[Plot_Nr]
    #for the first figure, create empty objects that will be overwritten into matplotlib objects. 
    fig, axarr = 'empty', 'objects'
    #create a new figure
    fig, axarr = utils.create_new_fig_if_necessary(newfig, fig, axarr, nrows_in_each_fig, ncols_in_each_fig, dpi = 300)
    #create numpy array of membranous over nonmembranous conservation ratios (identity)
    hist_data_AAIMON_mean = np.array(df['AAIMON_ratio_mean_all_TMDs'].dropna())
    #use numpy to create a histogram
    freq_counts_I, bin_array_I = np.histogram(hist_data_AAIMON_mean, bins=binlist)
    #assuming all of the bins are exactly the same size, make the width of the column equal to 70% of each bin
    col_width = float('%0.3f' % (0.95 * (bin_array_I[1] - bin_array_I[0])))
    #when align='center', the central point of the bar in the x-axis is simply the middle of the bins ((bin_0-bin_1)/2, etc)
    centre_of_bar_in_x_axis = (bin_array_I[:-2] + bin_array_I[1:-1]) / 2
    #add the final bin, which is physically located just after the last regular bin but represents all higher values
    bar_width = centre_of_bar_in_x_axis[3] - centre_of_bar_in_x_axis[2]
    centre_of_bar_in_x_axis = np.append(centre_of_bar_in_x_axis, centre_of_bar_in_x_axis[-1] + bar_width)
    barcontainer_AAIMON_mean = axarr[row_nr, col_nr].bar(left=centre_of_bar_in_x_axis, height=freq_counts_I,
                                                         align='center', width=col_width, color="mediumblue")
                                                         #alpha=0.5)  # edgecolor='black',
    #create numpy array of membranous over nonmembranous conservation ratios (identity + similarity)
    hist_data_AAISON_mean = np.array(df['AASMON_ratio_mean_all_TMDs'].dropna())
    #use numpy to create a histogram
    freq_counts_S, bin_array_S = np.histogram(hist_data_AAISON_mean, bins=binlist)
    #barcontainer_S = axarr[row_nr,col_nr].bar(left=centre_of_bar_in_x_axis, height=freq_counts_S, align='center', width=col_width, color="#0101DF", edgecolor="#0101DF", alpha = 0.5)
    #create a line graph rather than a bar graph for the AAISON (ident + similarity)
    linecontainer_AAISON_mean = axarr[row_nr, col_nr].plot(centre_of_bar_in_x_axis, freq_counts_S, color="black")
                                                           #alpha=0.5)
    #other colours that are compatible with colourblind readers: #8A084B Dark red, #B45F04 deep orange, reddish purple #4B088A 
    #http://html-color-codes.info/
    #label the x-axis for each plot, based on the TMD
    axarr[row_nr, col_nr].set_xlabel('average conservation ratio (membranous over nonmembranous)', fontsize=fontsize)
    #move the x-axis label closer to the x-axis
    axarr[row_nr, col_nr].xaxis.set_label_coords(0.45, -0.085)
    #pylab.rcParams['figure.figsize'] = (50.0, 40.0)
    #pylab.rcParams['figure.figsize'] = (20.0, 16.0)
    #plt.show()
    xlim_min = settings["hist_settings_mult_proteins"]["xlim_min01"]
    #take x-axis max from settings
    xlim_max = settings["hist_settings_mult_proteins"]["xlim_max01"]
    #set x-axis min
    axarr[row_nr, col_nr].set_xlim(xlim_min, xlim_max)
    #set x-axis ticks
    #use the slide selection to select every second item in the list as an xtick(axis label)
    axarr[row_nr, col_nr].set_xticks([float('%0.1f' % c) for c in centre_of_bar_in_x_axis[::3]])
    axarr[row_nr, col_nr].set_ylabel('freq', rotation='vertical', fontsize=fontsize)
    #change axis font size
    axarr[row_nr, col_nr].tick_params(labelsize=fontsize)
    #create legend?#http://stackoverflow.com/questions/9834452/how-do-i-make-a-single-legend-for-many-subplots-with-matplotlib
    #legend_obj = axarr[row_nr, col_nr].legend(['AASMON (identity + similarity)', 'AAIMON (identity)'], loc='upper right',
                                 #fontsize=fontsize)
    #add figure number to top left of subplot
    #axarr[row_nr, col_nr].annotate(s=str(Plot_Nr) + '.', xy = (0.04,0.9), fontsize=5, xytext=None, xycoords='axes fraction', alpha=0.75)
    #add figure title to top left of subplot
    #axarr[row_nr, col_nr].annotate(s=title, xy = (0.1,0.9), fontsize=5, xytext=None, xycoords='axes fraction', alpha=0.75)    
    #improve ggplot style for a canvas (fig) with 4 figures (plots)
    
    #axarr.set_axis_bgcolor("white")
    axarr[row_nr, col_nr].yaxis.grid(True,zorder =0,linestyle = ":", color = "grey")
    axarr[row_nr, col_nr]
    
    for tic in axarr[row_nr, col_nr].xaxis.get_major_ticks():
        tic.tick1On  = False
    for tic in axarr[row_nr, col_nr].yaxis.get_major_ticks():
        tic.tick1On  = False
    

    axarr[row_nr, col_nr].spines['top'].set_visible(False)
    axarr[row_nr, col_nr].spines['right'].set_visible(False)
    
    #utils.improve_ggplot_for_4_plots(axarr,row_nr,col_nr,backgroundcolour, legend_obj)


    '''
    Fig02: Histogram of standard deviations for AAIMON and AASMON (std among homologues for each protein)
    '''
    Plot_Nr = 2
    #add non-functional "if" line to aid document navigation in some IDEs (e.g. Spyder)
    if fig_title == "Fig02: Histogram of standard deviations for AAIMON and AASMON":
        pass    
    sys.stdout.write(str(Plot_Nr) + ', ')
    title = 'Standard Deviaton, SP vs MP'
    newfig, savefig, Plot_Nr, plot_nr_in_fig, row_nr, col_nr = dict_organising_subplots[Plot_Nr]
    #if a new figure should be created (either because the orig is full, or the last TMD is analysed)
    fig, axarr = utils.create_new_fig_if_necessary(newfig, fig, axarr, nrows_in_each_fig, ncols_in_each_fig, dpi = 300)
    #create numpy array of membranous over nonmembranous conservation ratios (identity)
    hist_data_AAIMON_std = np.array(df['AAIMON_ratio_std_all_TMDs'].dropna())
    #use numpy to create a histogram
    number_of_bins = 50
    freq_counts_S, bin_array_S = np.histogram(hist_data_AAIMON_std, bins=number_of_bins)
    #assuming all of the bins are exactly the same size, make the width of the column equal to 70% of each bin
    col_width = float('%0.3f' % (0.95 * (bin_array_S[1] - bin_array_S[0])))
    #when align='center', the central point of the bar in the x-axis is simply the middle of the bins ((bin_0-bin_1)/2, etc)
    centre_of_bar_in_x_axis = (bin_array_S[:-1] + bin_array_S[1:]) / 2
    barcontainer_AAIMON_std = axarr[row_nr, col_nr].bar(left=centre_of_bar_in_x_axis, height=freq_counts_S,
                                                        align='center', width=col_width, color="#0489B1",
                                                        alpha=0.5)  # edgecolor='black',
    #create numpy array of membranous over nonmembranous conservation ratios (identity + similarity)
    hist_data_AAISON_std = np.array(df['AASMON_ratio_std_all_TMDs'].dropna())
    #use numpy to create a histogram
    #N.B. use the bins from the previous plot
    freq_counts, bin_array = np.histogram(hist_data_AAISON_std, bins=bin_array_S)
    #barcontainer_S = axarr[row_nr,col_nr].bar(left=centre_of_bar_in_x_axis, height=freq_counts, align='center', width=col_width, color="#0101DF", edgecolor="#0101DF", alpha = 0.5)
    #create a line graph rather than a bar graph for the AAISON (ident + similarity)
    linecontainer_AAISON_std = axarr[row_nr, col_nr].plot(centre_of_bar_in_x_axis, freq_counts, color="#0101DF",
                                                          alpha=0.5)
    #other colours that are compatible with colourblind readers: #8A084B Dark red, #B45F04 deep orange, reddish purple #4B088A 
    #http://html-color-codes.info/
    #label the x-axis for each plot, based on the TMD
    axarr[row_nr, col_nr].set_xlabel('average standard deviation', fontsize=fontsize)
    #move the x-axis label closer to the x-axis
    axarr[row_nr, col_nr].xaxis.set_label_coords(0.45, -0.085)
    #xlim_min = settings["hist_settings_mult_proteins"]["xlim_min01"]
    #take x-axis max from settings
    #xlim_max = settings["hist_settings_mult_proteins"]["xlim_max01"]
    #set x-axis min
    #axarr[row_nr,col_nr].set_xlim(xlim_min,xlim_max)
    #set x-axis ticks
    #use the slide selection to select every second item in the list as an xtick(axis label)
    axarr[row_nr, col_nr].set_xticks([float('%0.1f' % c) for c in centre_of_bar_in_x_axis[::3]])
    axarr[row_nr, col_nr].set_ylabel('freq', rotation='vertical', fontsize=fontsize)
    #change axis font size
    axarr[row_nr, col_nr].tick_params(labelsize=fontsize)
    #create legend?#http://stackoverflow.com/questions/9834452/how-do-i-make-a-single-legend-for-many-subplots-with-matplotlib
    legend_obj = axarr[row_nr, col_nr].legend(['AASMON (identity + similarity)', 'AAIMON (identity)'], loc='upper right',
                                 fontsize=fontsize)
    #add figure number to top left of subplot
    axarr[row_nr, col_nr].annotate(s=str(Plot_Nr) + '.', xy = (0.04,0.9), fontsize=5, xytext=None, xycoords='axes fraction', alpha=0.75)
    #add figure title to top left of subplot
    axarr[row_nr, col_nr].annotate(s=title, xy = (0.1,0.9), fontsize=5, xytext=None, xycoords='axes fraction', alpha=0.75)    
    #improve ggplot style for a canvas (fig) with 4 figures (plots)
    utils.improve_ggplot_for_4_plots(axarr,row_nr,col_nr,backgroundcolour, legend_obj)
 

    '''
    Fig03: Scattergram comparing mean AAIMON and AASMON
    '''
    Plot_Nr = 3
    #add non-functional "if" line to aid document navigation in some IDEs (e.g. Spyder)
    if fig_title == "Fig03: Scattergram comparing mean AAIMON and AASMON":
        pass  
    sys.stdout.write(str(Plot_Nr) + ', ')
    title = 'AAIMON vs AASMON'
    newfig, savefig, Plot_Nr, plot_nr_in_fig, row_nr, col_nr = dict_organising_subplots[Plot_Nr]
    #if a new figure should be created (either because the orig is full, or the last TMD is analysed)
    fig, axarr = utils.create_new_fig_if_necessary(newfig, fig, axarr, nrows_in_each_fig, ncols_in_each_fig, dpi = 300)
    #pylab.rcParams['figure.figsize'] = (50.0, 40.0)
    x = np.array(df['AAIMON_ratio_mean_all_TMDs'])
    y = np.array(df['AASMON_ratio_mean_all_TMDs'])
    scattercontainer_AAIMON_AASMON_mean = axarr[row_nr, col_nr].scatter(x=x, y=y, color="#8A084B", alpha=alpha, s=datapointsize)
    #label the x-axis for each plot, based on the TMD
    axarr[row_nr, col_nr].set_xlabel('AAIMON_ratio (aa identity)', fontsize=fontsize)
    #move the x-axis label closer to the x-axis
    axarr[row_nr, col_nr].xaxis.set_label_coords(0.45, -0.085)
    axarr[row_nr, col_nr].set_ylabel('AASMON_ratio (aa identity + similarity)', fontsize=fontsize)
    #change axis font size
    axarr[row_nr, col_nr].tick_params(labelsize=fontsize)
    #create legend?#http://stackoverflow.com/questions/9834452/how-do-i-make-a-single-legend-for-many-subplots-with-matplotlib
    legend_obj = axarr[row_nr, col_nr].legend(['mean'], loc='upper right', fontsize=fontsize)
    #add figure number to top left of subplot
    axarr[row_nr, col_nr].annotate(s=str(Plot_Nr) + '.', xy = (0.04,0.9), fontsize=5, xytext=None, xycoords='axes fraction', alpha=0.75)
    #add figure title to top left of subplot
    axarr[row_nr, col_nr].annotate(s=title, xy = (0.1,0.9), fontsize=5, xytext=None, xycoords='axes fraction', alpha=0.75)    
    #improve ggplot style for a canvas (fig) with 4 figures (plots)
    utils.improve_ggplot_for_4_plots(axarr,row_nr,col_nr,backgroundcolour, legend_obj)
    
    '''
    Fig04: Scattergram comparing standard deviation AAIMON and AASMON
    '''
    Plot_Nr = 4
    #add non-functional "if" line to aid document navigation in some IDEs (e.g. Spyder)
    if fig_title == "Fig04: Scattergram comparing standard deviation AAIMON and AASMON":
        pass  
    title = 'standard deviation AAIMON vs AASMON'
    newfig, savefig, Plot_Nr, plot_nr_in_fig, row_nr, col_nr = dict_organising_subplots[Plot_Nr]
    #if a new figure should be created (either because the orig is full, or the last TMD is analysed)
    fig, axarr = utils.create_new_fig_if_necessary(newfig, fig, axarr, nrows_in_each_fig, ncols_in_each_fig, dpi = 300)
    x = np.array(df['AAIMON_ratio_std_all_TMDs'])
    y = np.array(df['AASMON_ratio_std_all_TMDs'])
    scattercontainer_AAIMON_AASMON_std = axarr[row_nr, col_nr].scatter(x=x, y=y, color="#B45F04", alpha=alpha, s=datapointsize)
    #other colours that are compatible with colourblind readers: #8A084B Dark red, #B45F04 deep orange, reddish purple #4B088A 
    #http://html-color-codes.info/
    #label the x-axis for each plot, based on the TMD
    axarr[row_nr, col_nr].set_xlabel('AAIMON_ratio', fontsize=fontsize)
    #move the x-axis label closer to the x-axis
    axarr[row_nr, col_nr].xaxis.set_label_coords(0.45, -0.085)
    axarr[row_nr, col_nr].set_ylabel('AASMON_ratio', fontsize=fontsize)
    #set x-axis ticks
    #use the slide selection to select every second item in the list as an xtick(axis label)
    #axarr[row_nr,col_nr].set_xticks([float('%0.1f' % c) for c in centre_of_bar_in_x_axis[::3]])
    #axarr[row_nr,col_nr].set_ylabel('freq',rotation = 'vertical', fontsize = fontsize)                                    
    #change axis font size
    axarr[row_nr, col_nr].tick_params(labelsize=fontsize)
    #create legend?#http://stackoverflow.com/questions/9834452/how-do-i-make-a-single-legend-for-many-subplots-with-matplotlib
    legend_obj = axarr[row_nr, col_nr].legend(['standard deviation'], loc='upper right', fontsize=fontsize)
    #add figure number to top left of subplot
    axarr[row_nr, col_nr].annotate(s=str(Plot_Nr) + '.', xy = (0.04,0.9), fontsize=5, xytext=None, xycoords='axes fraction', alpha=0.75)
    #add figure title to top left of subplot
    axarr[row_nr, col_nr].annotate(s=title, xy = (0.1,0.9), fontsize=5, xytext=None, xycoords='axes fraction', alpha=0.75)    
    #improve ggplot style for a canvas (fig) with 4 figures (plots)
    utils.improve_ggplot_for_4_plots(axarr,row_nr,col_nr,backgroundcolour, legend_obj)
    
    #save the figure if necessary (i.e., if the maximum number of plots per figure has been obtained)
    utils.savefig_if_necessary(savefig, fig, Plot_Nr, base_filepath = base_filename_summaries)

    '''
    Fig05: Scattergram comparing number_of_TMDs_in_uniprot_feature_list with mean AAIMON
    '''
    Plot_Nr = 5
    sys.stdout.write(str(Plot_Nr) + ', ')
    title = 'num_TMDs vs AAIMON'
    newfig, savefig, Plot_Nr, plot_nr_in_fig, row_nr, col_nr = dict_organising_subplots[Plot_Nr]
    #if a new figure should be created (either because the orig is full, or the last TMD is analysed)
    fig, axarr = utils.create_new_fig_if_necessary(newfig, fig, axarr, nrows_in_each_fig, ncols_in_each_fig, dpi = 300)
    #for backwards compatibility, check for old name.
    if 'number_of_TMDs_in_uniprot_feature_list' in df.columns:
        x = np.array(df['number_of_TMDs_in_uniprot_feature_list'])
    else:
        x = np.array(df['number_of_TMDs_in_uniprot_feature_list'])
    y = np.array(df['AAIMON_ratio_mean_all_TMDs'])
    scattercontainer_AAIMON_AASMON_std = axarr[row_nr, col_nr].scatter(x=x, y=y, color="#0489B1", alpha=alpha, s=datapointsize)
    #other colours that are compatible with colourblind readers: #8A084B Dark red, #B45F04 deep orange, reddish purple #4B088A 
    #http://html-color-codes.info/
    #label the x-axis for each plot, based on the TMD
    axarr[row_nr, col_nr].set_xlabel('number of TMDs in protein', fontsize=fontsize)
    #move the x-axis label closer to the x-axis
    axarr[row_nr, col_nr].xaxis.set_label_coords(0.45, -0.085)
    axarr[row_nr, col_nr].set_ylabel('Average AAIMON ratio for all TMDs', fontsize=fontsize)
    #change axis font size
    axarr[row_nr, col_nr].tick_params(labelsize=fontsize)
    #create legend?#http://stackoverflow.com/questions/9834452/how-do-i-make-a-single-legend-for-many-subplots-with-matplotlib
    legend_obj = axarr[row_nr, col_nr].legend(['AAIMON_ratio_mean'], loc='upper right', fontsize=fontsize)
    #add background grid
    axarr[row_nr, col_nr].grid(True, color='0.75', alpha=0.3)
    #add figure number to top left of subplot
    axarr[row_nr, col_nr].annotate(s=str(Plot_Nr) + '.', xy = (0.04,0.9), fontsize=5, xytext=None, xycoords='axes fraction', alpha=0.75)
    #add figure title to top left of subplot
    axarr[row_nr, col_nr].annotate(s=title, xy = (0.1,0.9), fontsize=5, xytext=None, xycoords='axes fraction', alpha=0.75)    
    #improve ggplot style for a canvas (fig) with 4 figures (plots)
    utils.improve_ggplot_for_4_plots(axarr,row_nr,col_nr,backgroundcolour, legend_obj)
    
    '''
    Fig06: Scattergram comparing query_length with mean AAIMON
    '''
    Plot_Nr = 6
    #add non-functional "if" line to aid document navigation in some IDEs (e.g. Spyder)
    if fig_title == "Fig06: Scattergram comparing query_length with mean AAIMON":
        pass  
    sys.stdout.write(str(Plot_Nr) + ', ')
    title = 'query_length vs AAIMON'
    newfig, savefig, Plot_Nr, plot_nr_in_fig, row_nr, col_nr = dict_organising_subplots[Plot_Nr]
    #if a new figure should be created (either because the orig is full, or the last TMD is analysed)
    fig, axarr = utils.create_new_fig_if_necessary(newfig, fig, axarr, nrows_in_each_fig, ncols_in_each_fig, dpi = 300)
    #pylab.rcParams['figure.figsize'] = (100.0, 80.0)
    x = np.array(df['query_length'])
    y = np.array(df['AAIMON_ratio_mean_all_TMDs'])
    scattercontainer_AAIMON_AASMON_std = axarr[row_nr, col_nr].scatter(x=x, y=y, color="#0489B1", alpha=alpha, s=datapointsize)
    #label the x-axis for each plot, based on the TMD
    axarr[row_nr, col_nr].set_xlabel('Length of protein', fontsize=fontsize)
    #move the x-axis label closer to the x-axis
    axarr[row_nr, col_nr].xaxis.set_label_coords(0.45, -0.085)
    axarr[row_nr, col_nr].set_ylabel('AAIMON_ratio', fontsize=fontsize)
    #change axis font size
    axarr[row_nr, col_nr].tick_params(labelsize=fontsize)
    #add figure number to top left of subplot
    axarr[row_nr, col_nr].annotate(s=str(Plot_Nr) + '.', xy = (0.04,0.9), fontsize=5, xytext=None, xycoords='axes fraction', alpha=0.75)
    #add figure title to top left of subplot
    axarr[row_nr, col_nr].annotate(s=title, xy = (0.1,0.9), fontsize=5, xytext=None, xycoords='axes fraction', alpha=0.75)    
    #improve ggplot style for a canvas (fig) with 4 figures (plots)
    utils.improve_ggplot_for_4_plots(axarr,row_nr,col_nr,backgroundcolour, legend_obj)

    '''
    Fig07: Scattergram comparing len_nonTMD_align_mean with mean AAIMON
    '''
    Plot_Nr = 7
    sys.stdout.write(str(Plot_Nr) + ', ')
    title = 'length nonTMD region'
    newfig, savefig, Plot_Nr, plot_nr_in_fig, row_nr, col_nr = dict_organising_subplots[Plot_Nr]
    #if a new figure should be created (either because the orig is full, or the last TMD is analysed)
    fig, axarr = utils.create_new_fig_if_necessary(newfig, fig, axarr, nrows_in_each_fig, ncols_in_each_fig, dpi = 300)
    #pylab.rcParams['figure.figsize'] = (100.0, 80.0)
    x = np.array(df['len_nonTMD_align_mean'])
    y = np.array(df['AAIMON_ratio_mean_all_TMDs'])
    scattercontainer_AAIMON_AASMON_std = axarr[row_nr, col_nr].scatter(x=x, y=y, color="#0489B1", alpha=alpha, s=datapointsize)
    #label the x-axis for each plot, based on the TMD
    axarr[row_nr, col_nr].set_xlabel('Average length of nonTMD region in homologues', fontsize=fontsize)
    #move the x-axis label closer to the x-axis
    axarr[row_nr, col_nr].xaxis.set_label_coords(0.45, -0.085)
    axarr[row_nr, col_nr].set_ylabel('AAIMON_ratio', fontsize=fontsize)
    #change axis font size
    axarr[row_nr, col_nr].tick_params(labelsize=fontsize)
    #add figure number to top left of subplot
    axarr[row_nr, col_nr].annotate(s=str(Plot_Nr) + '.', xy = (0.04,0.9), fontsize=5, xytext=None, xycoords='axes fraction', alpha=0.75)
    #add figure title to top left of subplot
    axarr[row_nr, col_nr].annotate(s=title, xy = (0.1,0.9), fontsize=5, xytext=None, xycoords='axes fraction', alpha=0.75)    
    #improve ggplot style for a canvas (fig) with 4 figures (plots)
    utils.improve_ggplot_for_4_plots(axarr,row_nr,col_nr,backgroundcolour, legend_obj)
    
    
    '''
    Fig08: Scattergram comparing total_number_of_simap_hits with mean AAIMON
    note that the total hits comes from SIMAP, so this doesn't say anything about how much data is available for each protein
    '''
    Plot_Nr = 8
    #add non-functional "if" line to aid document navigation in some IDEs (e.g. Spyder)
    if fig_title == "Fig08: Scattergram comparing total_number_of_simap_hits with mean AAIMON":
        pass  
    sys.stdout.write(str(Plot_Nr) + ', ')
    title = 'number SIMAP hits'
    newfig, savefig, Plot_Nr, plot_nr_in_fig, row_nr, col_nr = dict_organising_subplots[Plot_Nr]
    #if a new figure should be created (either because the orig is full, or the last TMD is analysed)
    fig, axarr = utils.create_new_fig_if_necessary(newfig, fig, axarr, nrows_in_each_fig, ncols_in_each_fig, dpi = 300)
    #pylab.rcParams['figure.figsize'] = (100.0, 80.0)
    x = np.array(df['total_number_of_simap_hits'])
    y = np.array(df['AAIMON_ratio_mean_all_TMDs'])
    scattercontainer_AAIMON_AASMON_std = axarr[row_nr, col_nr].scatter(x=x, y=y, color="#0489B1", alpha=alpha, s=datapointsize)
    #label the x-axis for each plot, based on the TMD
    axarr[row_nr, col_nr].set_xlabel('total number of homologues', fontsize=fontsize)
    #move the x-axis label closer to the x-axis
    axarr[row_nr, col_nr].xaxis.set_label_coords(0.45, -0.085)
    axarr[row_nr, col_nr].set_ylabel('AAIMON_ratio', fontsize=fontsize)
    #change axis font size
    axarr[row_nr, col_nr].tick_params(labelsize=fontsize)
    #add figure number to top left of subplot
    axarr[row_nr, col_nr].annotate(s=str(Plot_Nr) + '.', xy = (0.04,0.9), fontsize=5, xytext=None, xycoords='axes fraction', alpha=0.75)
    #add figure title to top left of subplot
    axarr[row_nr, col_nr].annotate(s=title, xy = (0.1,0.9), fontsize=5, xytext=None, xycoords='axes fraction', alpha=0.75)    
    #improve ggplot style for a canvas (fig) with 4 figures (plots)
    utils.improve_ggplot_for_4_plots(axarr,row_nr,col_nr,backgroundcolour, legend_obj)
    #save the figure if necessary (i.e., if the maximum number of plots per figure has been obtained)
    utils.savefig_if_necessary(savefig, fig, Plot_Nr, base_filepath = base_filename_summaries)


    '''
    SHOW THE DATA FOR EACH TMD IN THE DATASET
    '''
    #create list of colours to use in figures
    colour_lists = utils.create_colour_lists()
    tableau20 = colour_lists['tableau20']
    

    '''
    Fig09: Histogram of mean AAIMON ratios for each TMD separately
    '''
    Plot_Nr = 9
    sys.stdout.write(str(Plot_Nr) + ', ')
    title = 'Histogram of mean AAIMON ratios'    
    
   
    
    df_mean_AAIMON_each_TM, max_num_TMDs, legend = utils.create_df_with_mean_AAIMON_each_TM(df) 
    
    newfig, savefig, Plot_Nr, plot_nr_in_fig, row_nr, col_nr = dict_organising_subplots[Plot_Nr]
    #if a new figure should be created (either because the orig is full, or the last TMD is analysed)
    fig, axarr = utils.create_new_fig_if_necessary(newfig, fig, axarr, nrows_in_each_fig, ncols_in_each_fig, dpi = 300)
    
    title = 'AAIMON each TMD separately'
    num_bins = 30
    #"#0489B1"
    alpha = 0.25
    col_width_value = 0.95
    ylabel = 'freq'
    xlabel = 'average conservation ratio (membranous over nonmembranous)'
            
    for n, TM in enumerate(df_mean_AAIMON_each_TM.columns):
        #define the colour for that TMD
        #if there are more TMDs than colours, simply start from the beginning of the list again
        if n < len(tableau20):
            color_num = n
        else:
            color_num = n - len(tableau20)
        color = tableau20[color_num]
        
        hist_data = np.array(df_mean_AAIMON_each_TM[TM].dropna())
        '''
        Calculated the bins for a histogram, even for highly non-normal data
        '''
        #calculate 5th percentile
        percentile_5 = np.percentile(hist_data, 5)#hist_data.min()
        #calculate 9th percentile
        percentile_95 = np.percentile(hist_data, 95)#hist_data.max()
        #calculate difference
        percentile_95_minus_5 = percentile_95 - percentile_5
        #create buffer for bins
        extra_xaxis_range = percentile_95_minus_5 / 4
        #lowest bin is the 5th percentile minus the buffer, except where that is below zero
        data_min = percentile_5 - extra_xaxis_range#hist_data.min()
        #ata_min = 0 if data_max < 0 else data_max
        #highest bin is the 95th percentile
        data_max = percentile_95 + extra_xaxis_range#hist_data.max()
        #create bins using the min and max
        binlist = np.linspace(data_min,data_max,num_bins)
        #print(binlist)
        #use numpy to create a histogram
        freq_counts_I, bin_array_I = np.histogram(hist_data, bins = binlist)
        #print(freq_counts_I)
        #print(bin_array_I)
        #assuming all of the bins are exactly the same size, make the width of the column equal to 70% of each bin
        col_width = float('%0.3f' % (col_width_value * (bin_array_I[1] - bin_array_I[0])))
        #when align='center', the central point of the bar in the x-axis is simply the middle of the bins ((bin_0-bin_1)/2, etc)
        centre_of_bar_in_x_axis = (bin_array_I[:-2] + bin_array_I[1:-1]) / 2
        #add the final bin, which is physically located just after the last regular bin but represents all higher values
        bar_width = centre_of_bar_in_x_axis[3] - centre_of_bar_in_x_axis[2]
        centre_of_bar_in_x_axis = np.append(centre_of_bar_in_x_axis, centre_of_bar_in_x_axis[-1] + bar_width)
        barcontainer = axarr[row_nr, col_nr].bar(left=centre_of_bar_in_x_axis, height=freq_counts_I,
                                                             align='center', width=col_width, facecolor=color,
                                                             alpha=alpha, edgecolor='black', linewidth=0.2)  # edgecolor='black',
        #barcontainer = axarr[row_nr, col_nr].bar(left=centre_of_bar_in_x_axis, height=freq_counts_I,
        #                                                     align='center', width=col_width, facecolor=color,
        #                                                     alpha=alpha, edgecolor='black', linewidth=0.1)  # edgecolor='black',    #label the x-axis for each plot, based on the TMD
    axarr[row_nr, col_nr].set_xlabel(xlabel, fontsize=fontsize)
    #move the x-axis label closer to the x-axis
    axarr[row_nr, col_nr].xaxis.set_label_coords(0.45, -0.085)
    axarr[row_nr, col_nr].set_ylabel(ylabel, rotation='vertical', fontsize=fontsize)
    #change axis font size
    axarr[row_nr, col_nr].tick_params(labelsize=fontsize)
    #create legend?#http://stackoverflow.com/questions/9834452/how-do-i-make-a-single-legend-for-many-subplots-with-matplotlib
    legend_obj = axarr[row_nr, col_nr].legend(legend, loc='upper right',fontsize=fontsize)
    #add title
    #axarr[row_nr, col_nr].set_title(title,fontsize=fontsize)
    #add figure number to top left of subplot
    axarr[row_nr, col_nr].annotate(s=str(Plot_Nr) + '.', xy = (0.04,0.9), fontsize=5, xytext=None, xycoords='axes fraction', alpha=0.75)
    #add figure title to top left of subplot
    axarr[row_nr, col_nr].annotate(s=title, xy = (0.1,0.9), fontsize=5, xytext=None, xycoords='axes fraction', alpha=0.75)    
    #improve ggplot style for a canvas (fig) with 4 figures (plots)
    utils.improve_ggplot_for_4_plots(axarr,row_nr,col_nr,backgroundcolour, legend_obj)
    
    
    '''
    Fig10: Line histogram of mean AAIMON ratios for each TMD separately
    '''        
    Plot_Nr = 10
    #add non-functional "if" line to aid document navigation in some IDEs (e.g. Spyder)
    if fig_title == "Fig10: Line histogram of mean AAIMON ratios for each TMD separately":
        pass  
    sys.stdout.write(str(Plot_Nr) + ', ')
    title = 'Line histogram each TMD'    
    newfig, savefig, Plot_Nr, plot_nr_in_fig, row_nr, col_nr = dict_organising_subplots[Plot_Nr]
    #if a new figure should be created (either because the orig is full, or the last TMD is analysed)
    fig, axarr = utils.create_new_fig_if_necessary(newfig, fig, axarr, nrows_in_each_fig, ncols_in_each_fig, dpi = 300)
    num_bins = 50
    #"#0489B1"
    alpha = 0.9
    col_width_value = 0.95
    ylabel = 'freq'
    xlabel = 'average conservation ratio (membranous over nonmembranous)'
    
    for n, TM in enumerate(df_mean_AAIMON_each_TM.columns):
        #define the colour for that TMD
        #if there are more TMDs than colours, simply start from the beginning of the list again
        if n < len(tableau20):
            color_num = n
        else:
            color_num = n - len(tableau20)
        color = tableau20[color_num]
        
        hist_data = np.array(df_mean_AAIMON_each_TM[TM].dropna())
        '''
        Calculated the bins for a histogram, even for highly non-normal data
        '''
        #calculate 5th percentile
        percentile_5 = np.percentile(hist_data, 5)#hist_data.min()
        #calculate 9th percentile
        percentile_95 = np.percentile(hist_data, 95)#hist_data.max()
        #calculate difference
        percentile_95_minus_5 = percentile_95 - percentile_5
        #create buffer for bins
        extra_xaxis_range = percentile_95_minus_5 / 4
        #lowest bin is the 5th percentile minus the buffer, except where that is below zero
        data_min = percentile_5 - extra_xaxis_range#hist_data.min()
        #ata_min = 0 if data_max < 0 else data_max
        #highest bin is the 95th percentile
        data_max = percentile_95 + extra_xaxis_range#hist_data.max()
        #create bins using the min and max
        binlist = np.linspace(data_min,data_max,num_bins)
        #print(binlist)
        #use numpy to create a histogram
        freq_counts_I, bin_array_I = np.histogram(hist_data, bins = binlist)
        #print(freq_counts_I)
        #print(bin_array_I)
        #assuming all of the bins are exactly the same size, make the width of the column equal to 70% of each bin
        col_width = float('%0.3f' % (col_width_value * (bin_array_I[1] - bin_array_I[0])))
        #when align='center', the central point of the bar in the x-axis is simply the middle of the bins ((bin_0-bin_1)/2, etc)
        centre_of_bar_in_x_axis = (bin_array_I[:-2] + bin_array_I[1:-1]) / 2
        #add the final bin, which is physically located just after the last regular bin but represents all higher values
        bar_width = centre_of_bar_in_x_axis[3] - centre_of_bar_in_x_axis[2]
        centre_of_bar_in_x_axis = np.append(centre_of_bar_in_x_axis, centre_of_bar_in_x_axis[-1] + bar_width)
        linecontainer = axarr[row_nr, col_nr].plot(centre_of_bar_in_x_axis, freq_counts_I, alpha=alpha, linewidth=0.3)  # edgecolor='black',
    #label the x-axis for each plot, based on the TMD
    axarr[row_nr, col_nr].set_xlabel(xlabel, fontsize=fontsize)
    #move the x-axis label closer to the x-axis
    axarr[row_nr, col_nr].xaxis.set_label_coords(0.45, -0.085)
    axarr[row_nr, col_nr].set_ylabel(ylabel, rotation='vertical', fontsize=fontsize)
    #change axis font size
    axarr[row_nr, col_nr].tick_params(labelsize=fontsize)
    #create legend?#http://stackoverflow.com/questions/9834452/how-do-i-make-a-single-legend-for-many-subplots-with-matplotlib
    legend_obj = axarr[row_nr, col_nr].legend(legend, loc='upper right',fontsize=fontsize)
    #add title
    #axarr[row_nr, col_nr].set_title(title,fontsize=fontsize)
    #add figure number to top left of subplot
    axarr[row_nr, col_nr].annotate(s=str(Plot_Nr) + '.', xy = (0.04,0.9), fontsize=5, xytext=None, xycoords='axes fraction', alpha=0.75)
    #add figure title to top left of subplot
    axarr[row_nr, col_nr].annotate(s=title, xy = (0.1,0.9), fontsize=5, xytext=None, xycoords='axes fraction', alpha=0.75)    
    #improve ggplot style for a canvas (fig) with 4 figures (plots)
    utils.improve_ggplot_for_4_plots(axarr,row_nr,col_nr,backgroundcolour, legend_obj)
    


    #these graphs are only applicable for multi-pass proteins. Use where at least 2 proteins have a 7th TMD
    if df_mean_AAIMON_each_TM['TM07'].dropna().shape[0] >= 2:
        dataset_contains_multipass_prots = True
    else:
        dataset_contains_multipass_prots = False
    if dataset_contains_multipass_prots:
        '''
        Fig11: Line histogram of mean AAIMON ratios for selected TMDs, highlighting difference for TM07
        '''         
        Plot_Nr = 11
        sys.stdout.write(str(Plot_Nr) + ', ')
        title = 'Select TMDs, all data'    
        cols_for_analysis = ['TM01','TM07','TM08','last_TM_AAIMON_ratio_mean']
        newfig, savefig, Plot_Nr, plot_nr_in_fig, row_nr, col_nr = dict_organising_subplots[Plot_Nr]
        #if a new figure should be created (either because the orig is full, or the last TMD is analysed)
        fig, axarr = utils.create_new_fig_if_necessary(newfig, fig, axarr, nrows_in_each_fig, ncols_in_each_fig, dpi = 300)
        num_bins = 50
        #"#0489B1"
        alpha = 0.7
        col_width_value = 0.95
        ylabel = 'freq'
        xlabel = 'average conservation ratio (membranous over nonmembranous)'
     
        for n, TM in enumerate(cols_for_analysis):
            #define the colour for that TMD
            #if there are more TMDs than colours, simply start from the beginning of the list again
            if n < len(tableau20):
                color_num = n
            else:
                color_num = n - len(tableau20)
            color = tableau20[color_num]
            
            hist_data = np.array(df_mean_AAIMON_each_TM[TM].dropna())
            '''
            Calculated the bins for a histogram, even for highly non-normal data
            '''
            #calculate 5th percentile
            percentile_5 = np.percentile(hist_data, 5)#hist_data.min()
            #calculate 9th percentile
            percentile_95 = np.percentile(hist_data, 95)#hist_data.max()
            #calculate difference
            percentile_95_minus_5 = percentile_95 - percentile_5
            #create buffer for bins
            extra_xaxis_range = percentile_95_minus_5 / 4
            #lowest bin is the 5th percentile minus the buffer, except where that is below zero
            data_min = percentile_5 - extra_xaxis_range#hist_data.min()
            #ata_min = 0 if data_max < 0 else data_max
            #highest bin is the 95th percentile
            data_max = percentile_95 + extra_xaxis_range#hist_data.max()
            #create bins using the min and max
            binlist = np.linspace(data_min,data_max,num_bins)
            #print(binlist)
            #use numpy to create a histogram
            freq_counts_I, bin_array_I = np.histogram(hist_data, bins = binlist)
            #print(freq_counts_I)
            #print(bin_array_I)
            #assuming all of the bins are exactly the same size, make the width of the column equal to 70% of each bin
            col_width = float('%0.3f' % (col_width_value * (bin_array_I[1] - bin_array_I[0])))
            #when align='center', the central point of the bar in the x-axis is simply the middle of the bins ((bin_0-bin_1)/2, etc)
            centre_of_bar_in_x_axis = (bin_array_I[:-2] + bin_array_I[1:-1]) / 2
            #add the final bin, which is physically located just after the last regular bin but represents all higher values
            bar_width = centre_of_bar_in_x_axis[3] - centre_of_bar_in_x_axis[2]
            centre_of_bar_in_x_axis = np.append(centre_of_bar_in_x_axis, centre_of_bar_in_x_axis[-1] + bar_width)
            linecontainer = axarr[row_nr, col_nr].plot(centre_of_bar_in_x_axis, freq_counts_I, alpha=alpha, linewidth=0.3)  # edgecolor='black',
        #label the x-axis for each plot, based on the TMD
        axarr[row_nr, col_nr].set_xlabel(xlabel, fontsize=fontsize)
        #move the x-axis label closer to the x-axis
        axarr[row_nr, col_nr].xaxis.set_label_coords(0.45, -0.085)
        axarr[row_nr, col_nr].set_ylabel(ylabel, rotation='vertical', fontsize=fontsize)
        #change axis font size
        axarr[row_nr, col_nr].tick_params(labelsize=fontsize)
        #create legend?#http://stackoverflow.com/questions/9834452/how-do-i-make-a-single-legend-for-many-subplots-with-matplotlib
        legend_obj = axarr[row_nr, col_nr].legend(cols_for_analysis, loc='upper right',fontsize=fontsize)
        #add title
        #axarr[row_nr, col_nr].set_title(title,fontsize=fontsize)
        #add figure number to top left of subplot
        axarr[row_nr, col_nr].annotate(s=str(Plot_Nr) + '.', xy = (0.04,0.9), fontsize=5, xytext=None, xycoords='axes fraction', alpha=0.75)
        #add figure title to top left of subplot
        axarr[row_nr, col_nr].annotate(s=title, xy = (0.1,0.9), fontsize=5, xytext=None, xycoords='axes fraction', alpha=0.75)    
        #improve ggplot style for a canvas (fig) with 4 figures (plots)
        utils.improve_ggplot_for_4_plots(axarr,row_nr,col_nr,backgroundcolour, legend_obj)
        
    
        '''
        Fig12: TMD 1-5 only
        '''
        Plot_Nr = 12
        #add non-functional "if" line to aid document navigation in some IDEs (e.g. Spyder)
        if fig_title == "Fig12: TMD 1-5 only":
            pass  
        sys.stdout.write(str(Plot_Nr) + ', ')
        col_start = 0
        col_end = 5
        cols_to_analyse = df_mean_AAIMON_each_TM.columns[col_start:col_end]
        title = 'TMD 1 to 5, all data'    
        newfig, savefig, Plot_Nr, plot_nr_in_fig, row_nr, col_nr = dict_organising_subplots[Plot_Nr]
        #if a new figure should be created (either because the orig is full, or the last TMD is analysed)
        fig, axarr = utils.create_new_fig_if_necessary(newfig, fig, axarr, nrows_in_each_fig, ncols_in_each_fig, dpi = 300)
        num_bins = 30
        #"#0489B1"
        alpha = 0.9
        col_width_value = 0.95
        ylabel = 'freq'
        xlabel = 'average conservation ratio (membranous over nonmembranous)'
        
        for n, TM in enumerate(cols_to_analyse):
            #define the colour for that TMD
            #if there are more TMDs than colours, simply start from the beginning of the list again
            if n < len(tableau20):
                color_num = n
            else:
                color_num = n - len(tableau20)
            color = tableau20[color_num]
            
            hist_data = np.array(df_mean_AAIMON_each_TM[TM].dropna())
            '''
            Calculated the bins for a histogram, even for highly non-normal data
            '''
            #calculate 5th percentile
            percentile_5 = np.percentile(hist_data, 5)#hist_data.min()
            #calculate 9th percentile
            percentile_95 = np.percentile(hist_data, 95)#hist_data.max()
            #calculate difference
            percentile_95_minus_5 = percentile_95 - percentile_5
            #create buffer for bins
            extra_xaxis_range = percentile_95_minus_5 / 4
            #lowest bin is the 5th percentile minus the buffer, except where that is below zero
            data_min = percentile_5 - extra_xaxis_range#hist_data.min()
            #ata_min = 0 if data_max < 0 else data_max
            #highest bin is the 95th percentile
            data_max = percentile_95 + extra_xaxis_range#hist_data.max()
            #create bins using the min and max
            binlist = np.linspace(data_min,data_max,num_bins)
            #print(binlist)
            #use numpy to create a histogram
            freq_counts_I, bin_array_I = np.histogram(hist_data, bins = binlist)
            #print(freq_counts_I)
            #print(bin_array_I)
            #assuming all of the bins are exactly the same size, make the width of the column equal to 70% of each bin
            col_width = float('%0.3f' % (col_width_value * (bin_array_I[1] - bin_array_I[0])))
            #when align='center', the central point of the bar in the x-axis is simply the middle of the bins ((bin_0-bin_1)/2, etc)
            centre_of_bar_in_x_axis = (bin_array_I[:-2] + bin_array_I[1:-1]) / 2
            #add the final bin, which is physically located just after the last regular bin but represents all higher values
            bar_width = centre_of_bar_in_x_axis[3] - centre_of_bar_in_x_axis[2]
            centre_of_bar_in_x_axis = np.append(centre_of_bar_in_x_axis, centre_of_bar_in_x_axis[-1] + bar_width)
            linecontainer = axarr[row_nr, col_nr].plot(centre_of_bar_in_x_axis, freq_counts_I, alpha=alpha, linewidth=0.3)  # edgecolor='black',
        #label the x-axis for each plot, based on the TMD
        axarr[row_nr, col_nr].set_xlabel(xlabel, fontsize=fontsize)
        #move the x-axis label closer to the x-axis
        axarr[row_nr, col_nr].xaxis.set_label_coords(0.45, -0.085)
        axarr[row_nr, col_nr].set_ylabel(ylabel, rotation='vertical', fontsize=fontsize)
        #change axis font size
        axarr[row_nr, col_nr].tick_params(labelsize=fontsize)
        #create legend?#http://stackoverflow.com/questions/9834452/how-do-i-make-a-single-legend-for-many-subplots-with-matplotlib
        legend_obj = axarr[row_nr, col_nr].legend(cols_to_analyse, loc='upper right',fontsize=fontsize)
        #add title
        #axarr[row_nr, col_nr].set_title(title,fontsize=fontsize)
        #add figure number to top left of subplot
        axarr[row_nr, col_nr].annotate(s=str(Plot_Nr) + '.', xy = (0.04,0.9), fontsize=5, xytext=None, xycoords='axes fraction', alpha=0.75)
        #add figure title to top left of subplot
        axarr[row_nr, col_nr].annotate(s=title, xy = (0.1,0.9), fontsize=5, xytext=None, xycoords='axes fraction', alpha=0.75)    
        #improve ggplot style for a canvas (fig) with 4 figures (plots)
        utils.improve_ggplot_for_4_plots(axarr,row_nr,col_nr,backgroundcolour, legend_obj)
        
        #save the figure as it is
        savefig = True
        if savefig:
            #automatically tighten the layout of plots in the figure
            #fig.tight_layout()
            #save files
            fig.savefig(base_filename_summaries + '_%01d.png' % Plot_Nr, format='png', dpi=400)
            fig.savefig(base_filename_summaries + '_%01d.pdf' % Plot_Nr, format='pdf')    
    
    
    
        '''
        Fig13: TMD 5-10 only
        '''
        Plot_Nr = 13
        title = 'TMD 5-10, all data' 
        sys.stdout.write(str(Plot_Nr) + ', ')
        col_start = 5
        col_end = 10
        cols_to_analyse = ['TM01'] + list(df_mean_AAIMON_each_TM.columns[col_start:col_end])
       
        newfig, savefig, Plot_Nr, plot_nr_in_fig, row_nr, col_nr = dict_organising_subplots[Plot_Nr]
        #if a new figure should be created (either because the orig is full, or the last TMD is analysed)
        fig, axarr = utils.create_new_fig_if_necessary(newfig, fig, axarr, nrows_in_each_fig, ncols_in_each_fig, dpi = 300)
        num_bins = 30
        #"#0489B1"
        alpha = 0.9
        col_width_value = 0.95
        ylabel = 'freq'
        xlabel = 'average conservation ratio (membranous over nonmembranous)'
        
        
        for n, TM in enumerate(cols_to_analyse):
            #define the colour for that TMD
            #if there are more TMDs than colours, simply start from the beginning of the list again
            if n < len(tableau20):
                color_num = n
            else:
                color_num = n - len(tableau20)
            color = tableau20[color_num]
            
            hist_data = np.array(df_mean_AAIMON_each_TM[TM].dropna())
            '''
            Calculated the bins for a histogram, even for highly non-normal data
            '''
            #calculate 5th percentile
            percentile_5 = np.percentile(hist_data, 5)#hist_data.min()
            #calculate 9th percentile
            percentile_95 = np.percentile(hist_data, 95)#hist_data.max()
            #calculate difference
            percentile_95_minus_5 = percentile_95 - percentile_5
            #create buffer for bins
            extra_xaxis_range = percentile_95_minus_5 / 4
            #lowest bin is the 5th percentile minus the buffer, except where that is below zero
            data_min = percentile_5 - extra_xaxis_range#hist_data.min()
            #ata_min = 0 if data_max < 0 else data_max
            #highest bin is the 95th percentile
            data_max = percentile_95 + extra_xaxis_range#hist_data.max()
            #create bins using the min and max
            binlist = np.linspace(data_min,data_max,num_bins)
            #print(binlist)
            #use numpy to create a histogram
            freq_counts_I, bin_array_I = np.histogram(hist_data, bins = binlist)
            #print(freq_counts_I)
            #print(bin_array_I)
            #assuming all of the bins are exactly the same size, make the width of the column equal to 70% of each bin
            col_width = float('%0.3f' % (col_width_value * (bin_array_I[1] - bin_array_I[0])))
            #when align='center', the central point of the bar in the x-axis is simply the middle of the bins ((bin_0-bin_1)/2, etc)
            centre_of_bar_in_x_axis = (bin_array_I[:-2] + bin_array_I[1:-1]) / 2
            #add the final bin, which is physically located just after the last regular bin but represents all higher values
            bar_width = centre_of_bar_in_x_axis[3] - centre_of_bar_in_x_axis[2]
            centre_of_bar_in_x_axis = np.append(centre_of_bar_in_x_axis, centre_of_bar_in_x_axis[-1] + bar_width)
            linecontainer = axarr[row_nr, col_nr].plot(centre_of_bar_in_x_axis, freq_counts_I, alpha=alpha, linewidth=0.3)  # edgecolor='black',
        #label the x-axis for each plot, based on the TMD
        axarr[row_nr, col_nr].set_xlabel(xlabel, fontsize=fontsize)
        #move the x-axis label closer to the x-axis
        axarr[row_nr, col_nr].xaxis.set_label_coords(0.45, -0.085)
        axarr[row_nr, col_nr].set_ylabel(ylabel, rotation='vertical', fontsize=fontsize)
        #change axis font size
        axarr[row_nr, col_nr].tick_params(labelsize=fontsize)
        #create legend?#http://stackoverflow.com/questions/9834452/how-do-i-make-a-single-legend-for-many-subplots-with-matplotlib
        legend_obj = axarr[row_nr, col_nr].legend(cols_to_analyse, loc='upper right',fontsize=fontsize)
        #add title
        #axarr[row_nr, col_nr].set_title(title,fontsize=fontsize)
        #add figure number to top left of subplot
        axarr[row_nr, col_nr].annotate(s=str(Plot_Nr) + '.', xy = (0.04,0.9), fontsize=5, xytext=None, xycoords='axes fraction', alpha=0.75)
        #add figure title to top left of subplot
        axarr[row_nr, col_nr].annotate(s=title, xy = (0.1,0.9), fontsize=5, xytext=None, xycoords='axes fraction', alpha=0.75)    
        #improve ggplot style for a canvas (fig) with 4 figures (plots)
        utils.improve_ggplot_for_4_plots(axarr,row_nr,col_nr,backgroundcolour, legend_obj)
        
    
        '''
        Fig14: TMD 10-15 only
        '''    
        Plot_Nr = 14
        #add non-functional "if" line to aid document navigation in some IDEs (e.g. Spyder)
        if fig_title == "Fig14: TMD 10-15 only":
            pass  
        title = 'TMD 10-15, all data'
        sys.stdout.write(str(Plot_Nr) + ', ')
        col_start = 10
        col_end = 15
        cols_to_analyse = ['TM01'] + list(df_mean_AAIMON_each_TM.columns[col_start:col_end])   
        newfig, savefig, Plot_Nr, plot_nr_in_fig, row_nr, col_nr = dict_organising_subplots[Plot_Nr]
        #if a new figure should be created (either because the orig is full, or the last TMD is analysed)
        fig, axarr = utils.create_new_fig_if_necessary(newfig, fig, axarr, nrows_in_each_fig, ncols_in_each_fig, dpi = 300)
        num_bins = 10
        #"#0489B1"
        alpha = 0.9
        col_width_value = 0.95
        ylabel = 'freq'
        xlabel = 'average conservation ratio (membranous over nonmembranous)'
        
        for n, TM in enumerate(cols_to_analyse):
            #define the colour for that TMD
            #if there are more TMDs than colours, simply start from the beginning of the list again
            if n < len(tableau20):
                color_num = n
            else:
                color_num = n - len(tableau20)
            color = tableau20[color_num]
            
            hist_data = np.array(df_mean_AAIMON_each_TM[TM].dropna())
            '''
            Calculated the bins for a histogram, even for highly non-normal data
            '''
            #calculate 5th percentile
            percentile_5 = np.percentile(hist_data, 5)#hist_data.min()
            #calculate 9th percentile
            percentile_95 = np.percentile(hist_data, 95)#hist_data.max()
            #calculate difference
            percentile_95_minus_5 = percentile_95 - percentile_5
            #create buffer for bins
            extra_xaxis_range = percentile_95_minus_5 / 4
            #lowest bin is the 5th percentile minus the buffer, except where that is below zero
            data_min = percentile_5 - extra_xaxis_range#hist_data.min()
            #ata_min = 0 if data_max < 0 else data_max
            #highest bin is the 95th percentile
            data_max = percentile_95 + extra_xaxis_range#hist_data.max()
            #create bins using the min and max
            binlist = np.linspace(data_min,data_max,num_bins)
            #print(binlist)
            #use numpy to create a histogram
            freq_counts_I, bin_array_I = np.histogram(hist_data, bins = binlist)
            #normalize the frequency counts
            freq_counts_normalised = freq_counts_I/freq_counts_I.max()
    
            #assuming all of the bins are exactly the same size, make the width of the column equal to 70% of each bin
            col_width = float('%0.3f' % (col_width_value * (bin_array_I[1] - bin_array_I[0])))
            #when align='center', the central point of the bar in the x-axis is simply the middle of the bins ((bin_0-bin_1)/2, etc)
            centre_of_bar_in_x_axis = (bin_array_I[:-2] + bin_array_I[1:-1]) / 2
            #add the final bin, which is physically located just after the last regular bin but represents all higher values
            bar_width = centre_of_bar_in_x_axis[3] - centre_of_bar_in_x_axis[2]
            centre_of_bar_in_x_axis = np.append(centre_of_bar_in_x_axis, centre_of_bar_in_x_axis[-1] + bar_width)
            linecontainer = axarr[row_nr, col_nr].plot(centre_of_bar_in_x_axis, freq_counts_normalised, alpha=alpha, linewidth=0.3)  # edgecolor='black',
        #label the x-axis for each plot, based on the TMD
        axarr[row_nr, col_nr].set_xlabel(xlabel, fontsize=fontsize)
        #move the x-axis label closer to the x-axis
        axarr[row_nr, col_nr].xaxis.set_label_coords(0.45, -0.085)
        axarr[row_nr, col_nr].set_ylabel(ylabel, rotation='vertical', fontsize=fontsize)
        #change axis font size
        axarr[row_nr, col_nr].tick_params(labelsize=fontsize)
        #create legend?#http://stackoverflow.com/questions/9834452/how-do-i-make-a-single-legend-for-many-subplots-with-matplotlib
        legend_obj = axarr[row_nr, col_nr].legend(cols_to_analyse, loc='upper right',fontsize=fontsize)
        #add title
        #axarr[row_nr, col_nr].set_title(title,fontsize=fontsize)
        #add figure number to top left of subplot
        axarr[row_nr, col_nr].annotate(s=str(Plot_Nr) + '.', xy = (0.04,0.9), fontsize=5, xytext=None, xycoords='axes fraction', alpha=0.75)
        #add figure title to top left of subplot
        axarr[row_nr, col_nr].annotate(s=title, xy = (0.1,0.9), fontsize=5, xytext=None, xycoords='axes fraction', alpha=0.75)    
        #improve ggplot style for a canvas (fig) with 4 figures (plots)
        utils.improve_ggplot_for_4_plots(axarr,row_nr,col_nr,backgroundcolour, legend_obj)
        
        '''
        Fig15: TMD 15-20 only
        '''    
        Plot_Nr = 15
        title = 'TMD 15-20, all data' 
        sys.stdout.write(str(Plot_Nr) + ', ')
        col_start = 15
        col_end = 20
        cols_to_analyse = ['TM01'] + list(df_mean_AAIMON_each_TM.columns[col_start:col_end])
        newfig, savefig, Plot_Nr, plot_nr_in_fig, row_nr, col_nr = dict_organising_subplots[Plot_Nr]
        #if a new figure should be created (either because the orig is full, or the last TMD is analysed)
        fig, axarr = utils.create_new_fig_if_necessary(newfig, fig, axarr, nrows_in_each_fig, ncols_in_each_fig, dpi = 300)
        num_bins = 10
        #"#0489B1"
        alpha = 0.9
        col_width_value = 0.95
        ylabel = 'freq'
        xlabel = 'average conservation ratio (membranous over nonmembranous)'
        
        for n, TM in enumerate(cols_to_analyse):
            #define the colour for that TMD
            #if there are more TMDs than colours, simply start from the beginning of the list again
            if n < len(tableau20):
                color_num = n
            else:
                color_num = n - len(tableau20)
            color = tableau20[color_num]
            
            hist_data = np.array(df_mean_AAIMON_each_TM[TM].dropna())
            '''
            Calculated the bins for a histogram, even for highly non-normal data
            '''
            #calculate 5th percentile
            percentile_5 = np.percentile(hist_data, 5)#hist_data.min()
            #calculate 9th percentile
            percentile_95 = np.percentile(hist_data, 95)#hist_data.max()
            #calculate difference
            percentile_95_minus_5 = percentile_95 - percentile_5
            #create buffer for bins
            extra_xaxis_range = percentile_95_minus_5 / 4
            #lowest bin is the 5th percentile minus the buffer, except where that is below zero
            data_min = percentile_5 - extra_xaxis_range#hist_data.min()
            #ata_min = 0 if data_max < 0 else data_max
            #highest bin is the 95th percentile
            data_max = percentile_95 + extra_xaxis_range#hist_data.max()
            #create bins using the min and max
            binlist = np.linspace(data_min,data_max,num_bins)
            #print(binlist)
            #use numpy to create a histogram
            freq_counts_I, bin_array_I = np.histogram(hist_data, bins = binlist)
            #normalize the frequency counts
            freq_counts_normalised = freq_counts_I/freq_counts_I.max()
            #assuming all of the bins are exactly the same size, make the width of the column equal to 70% of each bin
            col_width = float('%0.3f' % (col_width_value * (bin_array_I[1] - bin_array_I[0])))
            #when align='center', the central point of the bar in the x-axis is simply the middle of the bins ((bin_0-bin_1)/2, etc)
            centre_of_bar_in_x_axis = (bin_array_I[:-2] + bin_array_I[1:-1]) / 2
            #add the final bin, which is physically located just after the last regular bin but represents all higher values
            bar_width = centre_of_bar_in_x_axis[3] - centre_of_bar_in_x_axis[2]
            centre_of_bar_in_x_axis = np.append(centre_of_bar_in_x_axis, centre_of_bar_in_x_axis[-1] + bar_width)
            linecontainer = axarr[row_nr, col_nr].plot(centre_of_bar_in_x_axis, freq_counts_normalised, alpha=alpha, linewidth=0.3)  # edgecolor='black',
        #label the x-axis for each plot, based on the TMD
        axarr[row_nr, col_nr].set_xlabel(xlabel, fontsize=fontsize)
        #move the x-axis label closer to the x-axis
        axarr[row_nr, col_nr].xaxis.set_label_coords(0.45, -0.085)
        axarr[row_nr, col_nr].set_ylabel(ylabel, rotation='vertical', fontsize=fontsize)
        #change axis font size
        axarr[row_nr, col_nr].tick_params(labelsize=fontsize)
        #create legend?#http://stackoverflow.com/questions/9834452/how-do-i-make-a-single-legend-for-many-subplots-with-matplotlib
        legend_obj = axarr[row_nr, col_nr].legend(cols_to_analyse, loc='upper right',fontsize=fontsize)
        #add title
        #axarr[row_nr, col_nr].set_title(title,fontsize=fontsize)
        #add figure number to top left of subplot
        axarr[row_nr, col_nr].annotate(s=str(Plot_Nr) + '.', xy = (0.04,0.9), fontsize=5, xytext=None, xycoords='axes fraction', alpha=0.75)
        #add figure title to top left of subplot
        axarr[row_nr, col_nr].annotate(s=title, xy = (0.1,0.9), fontsize=5, xytext=None, xycoords='axes fraction', alpha=0.75)    
        #improve ggplot style for a canvas (fig) with 4 figures (plots)
        utils.improve_ggplot_for_4_plots(axarr,row_nr,col_nr,backgroundcolour, legend_obj)
            
       
        '''
        full dataset
        '''
        #df_nonGPCR = df.loc[df['Gprotein'] == False]
        
        #create numpy array of membranous over nonmembranous conservation ratios (identity)
        hist_data_AAIMON_mean = np.array(df['AAIMON_ratio_mean_all_TMDs'].dropna())
        #use numpy to create a histogram
        freq_counts_I, bin_array_I = np.histogram(hist_data_AAIMON_mean, bins=binlist)
        #normalize the frequency counts
        freq_counts_normalised = freq_counts_I/freq_counts_I.max()
        #assuming all of the bins are exactly the same size, make the width of the column equal to 70% of each bin
        col_width = float('%0.3f' % (0.95 * (bin_array_I[1] - bin_array_I[0])))
        #when align='center', the central point of the bar in the x-axis is simply the middle of the bins ((bin_0-bin_1)/2, etc)
        centre_of_bar_in_x_axis = (bin_array_I[:-2] + bin_array_I[1:-1]) / 2
        #add the final bin, which is physically located just after the last regular bin but represents all higher values
        bar_width = centre_of_bar_in_x_axis[3] - centre_of_bar_in_x_axis[2]
        centre_of_bar_in_x_axis = np.append(centre_of_bar_in_x_axis, centre_of_bar_in_x_axis[-1] + bar_width)
        barcontainer_AAIMON_mean = axarr[row_nr, col_nr].bar(left=centre_of_bar_in_x_axis, height=freq_counts_normalised,
                                                             align='center', width=col_width, color='#B45F04',
                                                             alpha=0.5, linewidth = 0.1)  # edgecolor='black',
        #other colours that are compatible with colourblind readers: #8A084B Dark red, #B45F04 deep orange, reddish purple #4B088A 
        #http://html-color-codes.info/
        #label the x-axis for each plot, based on the TMD
        axarr[row_nr, col_nr].set_xlabel('average conservation ratio (membranous over nonmembranous)', fontsize=fontsize)
        #move the x-axis label closer to the x-axis
        axarr[row_nr, col_nr].xaxis.set_label_coords(0.45, -0.085)
        xlim_min = 0.8
        xlim_max = 1.5
        axarr[row_nr, col_nr].set_xlim(xlim_min, xlim_max)
        #set x-axis ticks
        #use the slide selection to select every second item in the list as an xtick(axis label)
        axarr[row_nr, col_nr].set_xticks([float('%0.1f' % c) for c in centre_of_bar_in_x_axis[::3]])
        axarr[row_nr, col_nr].set_ylabel('freq', rotation='vertical', fontsize=fontsize)
        #change axis font size
        axarr[row_nr, col_nr].tick_params(labelsize=fontsize)
        #create legend?#http://stackoverflow.com/questions/9834452/how-do-i-make-a-single-legend-for-many-subplots-with-matplotlib
        legend_obj = axarr[row_nr, col_nr].legend(['AAIMON GPCR', 'AAIMON ALL'], loc='upper right',
                                     fontsize=fontsize)
        #add figure number to top left of subplot
        axarr[row_nr, col_nr].annotate(s=str(Plot_Nr) + '.', xy = (0.04,0.9), fontsize=5, xytext=None, xycoords='axes fraction', alpha=0.75)
        #add figure title to top left of subplot
        axarr[row_nr, col_nr].annotate(s=title, xy = (0.1,0.9), fontsize=5, xytext=None, xycoords='axes fraction', alpha=0.75)    
        #improve ggplot style for a canvas (fig) with 4 figures (plots)
        utils.improve_ggplot_for_4_plots(axarr,row_nr,col_nr,backgroundcolour, legend_obj)             
            
        '''
        Fig20: Boxplot of all TMDs
        '''
        Plot_Nr = 20
        title = 'Boxplot of all TMDs'
        #add non-functional "if" line to aid document navigation in some IDEs (e.g. Spyder)
        if fig_title == "Fig20: Boxplot of all TMDs":
            pass 
        sys.stdout.write(str(Plot_Nr) + ', ')
        newfig, savefig, Plot_Nr, plot_nr_in_fig, row_nr, col_nr = dict_organising_subplots[Plot_Nr]
        #if a new figure should be created (either because the orig is full, or the last TMD is analysed)
        fig, axarr = utils.create_new_fig_if_necessary(newfig, fig, axarr, nrows_in_each_fig, ncols_in_each_fig, dpi = 300)
        
        num_bins = 30
        #"#0489B1"
        alpha = 0.25
        col_width_value = 0.95
        ylabel = 'freq'
        xlabel = 'average conservation ratio (membranous over nonmembranous)'
        #legend = 
        
        max_num_TMDs = df.number_of_TMDs_in_uniprot_feature_list.max()
        legend = []
        data_to_plot = []
        for i in range(1, max_num_TMDs):
            TM = 'TM%02d' % i
            hist_data_AAIMON_each_TM = df['TM%02d_AAIMON_ratio_mean' % i].dropna()
            if len(hist_data_AAIMON_each_TM) > 0:
                data_to_plot.append(hist_data_AAIMON_each_TM)
                legend.append(TM)
        
        meanpointprops = dict(marker='o', markerfacecolor='black', markersize=2)#markeredgecolor='0.75',
        
        flierprops = dict(marker='o', markerfacecolor='green', markersize=12,
                      linestyle='none')
        boxplotcontainer = axarr[row_nr, col_nr].boxplot(data_to_plot, sym='+', whis = 1.5, showmeans=True, meanprops = meanpointprops)
        axarr[row_nr, col_nr].tick_params(labelsize=fontsize)
        for box in boxplotcontainer['boxes']:
            # change outline color
            box.set( color='black', linewidth=0.4)#'7570b3'
            # change fill color
            #box.set( facecolor = '#1b9e77' )
            box.set_linewidth(0.4)
            
        ## change color and linewidth of the whiskers
        for whisker in boxplotcontainer['whiskers']:
            whisker.set(color='black', linewidth=0.4, dashes = (1,1))
        
        ## change color and linewidth of the caps
        for cap in boxplotcontainer['caps']:
            cap.set(color='black', linewidth=0.4)
        
        ## change color and linewidth of the medians
        for median in boxplotcontainer['medians']:
            median.set(color='black', linewidth=0.4)
        
        # change the style of fliers and their fill
        for flier in boxplotcontainer['fliers']:
            flier.set(marker='o', color = '0.8', alpha=0.1, markerfacecolor='0.3', markersize=3)
    
        ## Remove top axes and right axes ticks
        axarr[row_nr, col_nr].get_xaxis().tick_bottom()
        axarr[row_nr, col_nr].get_yaxis().tick_left()
        ## Custom x-axis labels
        axarr[row_nr, col_nr].set_xticklabels(legend, rotation = 45)
        #add figure number to top left of subplot
        axarr[row_nr, col_nr].annotate(s=str(Plot_Nr) + '.', xy = (0.04,0.9), fontsize=5, xytext=None, xycoords='axes fraction', alpha=0.75)
        #add figure title to top left of subplot
        axarr[row_nr, col_nr].annotate(s=title, xy = (0.1,0.9), fontsize=5, xytext=None, xycoords='axes fraction', alpha=0.75)    
        #improve ggplot style for a canvas (fig) with 4 figures (plots)
        utils.improve_ggplot_for_4_plots(axarr,row_nr,col_nr,backgroundcolour, legend_obj)
        
        if savefig:
            #automatically tighten the layout of plots in the figure
            #fig.tight_layout()
            #save files
            fig.savefig(base_filename_summaries + '_%01d.png' % Plot_Nr, format='png', dpi=400)
            fig.savefig(base_filename_summaries + '_%01d.pdf' % Plot_Nr, format='pdf')
  
    '''
    FINAL - MAKE SURE LAST GRAPHS ARE SAVED
    '''        
    #save the figure as it is
    savefig = True
    #save the figure if necessary (i.e., if the maximum number of plots per figure has been obtained)
    utils.savefig_if_necessary(savefig, fig, Plot_Nr, base_filepath = base_filename_summaries)
    '''
    save the updated dataframe, containing the various extra columns used for the figure
    '''
    with open(dfout09_simap_AAIMON_02, 'w') as csv_out:
        df.to_csv(csv_out, sep = ",", quoting = csv.QUOTE_NONNUMERIC)
    logging.info('A07b_save_figures_describing_proteins_in_list is finished')
    

'''+++++++++++++++ Summary figures describing the conservation ratios of proteins in the list ++++++++++++++++++'''
A10_compare_lists = settings["run_settings"]["compare_lists"]["run"]
if A10_compare_lists:

#    if 'df1' in globals():
#        if isinstance(df1, pd.DataFrame):
#            reload_data_from_summary_files = False
#        else:
#            reload_data_from_summary_files = True
#    else:
#        reload_data_from_summary_files = True
 
     
    #if reload_data_from_summary_files == True:
    
    protein_lists = settings["run_settings"]["compare_lists"]["protein_lists"]
    protein_list_names = settings["run_settings"]["compare_lists"]["protein_list_names"]
    
    protein_lists_joined = '_'.join(['%02d' % n for n in protein_lists])
    base_path_summ_two_lists = os.path.join(main_folder, 'summaries', 'compare_lists')
    if os.path.isdir(base_path_summ_two_lists) == False:
        os.mkdir(base_path_summ_two_lists)
    base_filename_summ_two_lists = os.path.join(base_path_summ_two_lists,'Lists_%s' % protein_lists_joined)
        
    df_list = []
    for index, list_num in enumerate(protein_lists):
        base_filename_summ = os.path.join(main_folder, 'summaries', 'List%02d' % list_num)
        dfout09_simap_AAIMON_02 = '%s_simap_AAIMON_02.csv' % base_filename_summ
        if index == 0:
            df1 = pd.read_csv(dfout09_simap_AAIMON_02, sep=",", quoting=csv.QUOTE_NONNUMERIC,index_col=0)
            df_list.append(df1)
        if index == 1:
            df2 = pd.read_csv(dfout09_simap_AAIMON_02, sep=",", quoting=csv.QUOTE_NONNUMERIC,index_col=0)
            df_list.append(df2)
        if index == 2:
            df3 = pd.read_csv(dfout09_simap_AAIMON_02, sep=",", quoting=csv.QUOTE_NONNUMERIC,index_col=0)
            df_list.append(df3)
        if index == 3:
            df4 = pd.read_csv(dfout09_simap_AAIMON_02, sep=",", quoting=csv.QUOTE_NONNUMERIC,index_col=0)
            df_list.append(df4)
        if index == 4:
            df5 = pd.read_csv(dfout09_simap_AAIMON_02, sep=",", quoting=csv.QUOTE_NONNUMERIC,index_col=0)
            df_list.append(df5)
        if index == 5:
            df6 = pd.read_csv(dfout09_simap_AAIMON_02, sep=",", quoting=csv.QUOTE_NONNUMERIC,index_col=0)
            df_list.append(df6)
        if index == 6:
            logging.warning('ERROR! Too many lists included for analysis in A10_compare_lists')

    '''
    Prepare fonts, colours etc for following figures
    '''
    #set default font size for plot
    fontsize = 4
    datapointsize = 2
    alpha = 0.1
    #use a dictionary to organise the saving of multiple plots in multiple figures, with a certain number of plots per figure
    n_plots_per_fig = 4
    nrows_in_each_fig = 2
    ncols_in_each_fig = 2
    dict_organising_subplots = utils.create_dict_organising_subplots(n_plots_per_fig=n_plots_per_fig,
                                                                     n_rows=nrows_in_each_fig,
                                                                     n_cols=ncols_in_each_fig)

    colour_lists = utils.create_colour_lists()
    TUM_colours_list_with_greys = colour_lists['TUM_colours_list_with_greys']

    '''
    Fig01: Histogram of mean AAIMON ratios
    '''
    Plot_Nr = 1
    title = 'Histogram of mean AAIMON ratios'    
    newfig, savefig, Plot_Nr, plot_nr_in_fig, row_nr, col_nr = dict_organising_subplots[Plot_Nr]
    #if the plot is the last one, the figure should be saved
    #if a new figure should be created (either because the orig is full, or the last TMD is analysed)
    if newfig:
        #clear any previous plots
        plt.close('all')
        #create a new figure
        fig, axarr = plt.subplots(nrows=nrows_in_each_fig, ncols=ncols_in_each_fig, dpi=300)  # sharex=True
    
    for n, df in enumerate(df_list):
        data_column = 'AAIMON_ratio_mean_all_TMDs'  
        color = TUM_colours_list_with_greys[n+2]#"#0489B1"
        alpha = 0.7
        col_width_value = 0.95
        ylabel = 'freq'
        xlabel = 'average conservation ratio (membranous over nonmembranous)'
        utils.create_hist_from_df_col(df=df,
                                title=title,row_nr=row_nr,col_nr=col_nr,axarr=axarr,
                                settings=settings,
                                data_column = data_column,
                                color=color,
                                alpha=alpha,
                                col_width_value=col_width_value,
                                fontsize=fontsize,
                                xlabel=xlabel,
                                ylabel=ylabel,legend=protein_list_names
                                )
    
    Plot_Nr = 2
    title = 'Histogram of mean AASMON ratios'  
    newfig, savefig, Plot_Nr, plot_nr_in_fig, row_nr, col_nr = dict_organising_subplots[Plot_Nr]
    #if the plot is the last one, the figure should be saved
    #if a new figure should be created (either because the orig is full, or the last TMD is analysed)
    if newfig:
        #clear any previous plots
        plt.close('all')
        #create a new figure
        fig, axarr = plt.subplots(nrows=nrows_in_each_fig, ncols=ncols_in_each_fig, dpi=300)  # sharex=True
    
    for n, df in enumerate(df_list):
        data_column = 'AASMON_ratio_mean_all_TMDs'  
        color = TUM_colours_list_with_greys[n+2]#"#0489B1"
        alpha = 0.7
        col_width_value = 0.95
        ylabel = 'freq'
        xlabel = 'average conservation ratio (membranous over nonmembranous)'
        utils.create_hist_from_df_col(df=df,
                                    title=title,row_nr=row_nr,col_nr=col_nr,axarr=axarr,
                                    settings=settings,
                                    data_column = data_column,
                                    color=color,
                                    alpha=alpha,
                                    col_width_value=col_width_value,
                                    fontsize=fontsize,
                                    xlabel=xlabel,
                                    ylabel=ylabel,legend=protein_list_names
                                    )
    Plot_Nr = 3
    title = 'len_nonTMD_seq_match_mean' 
    newfig, savefig, Plot_Nr, plot_nr_in_fig, row_nr, col_nr = dict_organising_subplots[Plot_Nr]
    #if the plot is the last one, the figure should be saved
    #if a new figure should be created (either because the orig is full, or the last TMD is analysed)
    if newfig:
        #clear any previous plots
        plt.close('all')
        #create a new figure
        fig, axarr = plt.subplots(nrows=nrows_in_each_fig, ncols=ncols_in_each_fig, dpi=300)  # sharex=True
    
    for n, df in enumerate(df_list):
        data_column = 'len_nonTMD_seq_match_mean'
        color = TUM_colours_list_with_greys[n+2]#"#0489B1"
        alpha = 0.7
        col_width_value = 0.95
        ylabel = 'freq'
        xlabel = data_column
    
        utils.create_hist_from_df_col_with_auto_binlist(df=df,
                                                title=title,row_nr=row_nr,col_nr=col_nr,axarr=axarr,num_bins=60,
                                                settings=settings,
                                                data_column = data_column,
                                                color=color,
                                                alpha=alpha,
                                                col_width_value=col_width_value,
                                                fontsize=fontsize,
                                                xlabel=xlabel,
                                                ylabel=ylabel,legend=protein_list_names
                                                )
    Plot_Nr = 4
    title = 'len_nonTMD_align_mean' 
    newfig, savefig, Plot_Nr, plot_nr_in_fig, row_nr, col_nr = dict_organising_subplots[Plot_Nr]
    #if the plot is the last one, the figure should be saved
    #if a new figure should be created (either because the orig is full, or the last TMD is analysed)
    if newfig:
        #clear any previous plots
        plt.close('all')
        #create a new figure
        fig, axarr = plt.subplots(nrows=nrows_in_each_fig, ncols=ncols_in_each_fig, dpi=300)  # sharex=True
    
    for n, df in enumerate(df_list):
        data_column = 'len_nonTMD_align_mean'
        color = TUM_colours_list_with_greys[n+2]#"#0489B1"
        alpha = 0.7
        col_width_value = 0.95
        ylabel = 'freq'
        xlabel = data_column
    
        utils.create_hist_from_df_col_with_auto_binlist(df=df,
                                                title=title,row_nr=row_nr,col_nr=col_nr,axarr=axarr,num_bins=100,
                                                settings=settings,
                                                data_column = data_column,
                                                color=color,
                                                alpha=alpha,
                                                col_width_value=col_width_value,
                                                fontsize=fontsize,
                                                xlabel=xlabel,
                                                ylabel=ylabel,legend=protein_list_names
                                                )        
    #save the figure as it is
    if savefig:
        utils.save_fig_with_subplots(fig=fig, axarr = axarr, base_filename = base_filename_summ_two_lists,
                                     fig_nr = Plot_Nr, fontsize=fontsize)


    Plot_Nr = 5
    title = 'nonTMD_perc_ident_mean' 
    newfig, savefig, Plot_Nr, plot_nr_in_fig, row_nr, col_nr = dict_organising_subplots[Plot_Nr]
    #if the plot is the last one, the figure should be saved
    #if a new figure should be created (either because the orig is full, or the last TMD is analysed)
    if newfig:
        #close any old figures
        plt.close('all')
        #create a new figure
        fig, axarr = plt.subplots(nrows=nrows_in_each_fig, ncols=ncols_in_each_fig, dpi=300)  # sharex=True
    
    for n, df in enumerate(df_list):
        data_column = 'nonTMD_perc_ident_mean'
        color = TUM_colours_list_with_greys[n+2]#"#0489B1"
        alpha = 0.7
        col_width_value = 0.95
        ylabel = 'freq'
        xlabel = data_column
    
        utils.create_hist_from_df_col_with_auto_binlist(df=df,
                                                title=title,row_nr=row_nr,col_nr=col_nr,axarr=axarr,num_bins=60,
                                                settings=settings,
                                                data_column = data_column,
                                                color=color,
                                                alpha=alpha,
                                                col_width_value=col_width_value,
                                                fontsize=fontsize,
                                                xlabel=xlabel,
                                                ylabel=ylabel,legend=protein_list_names
                                                )


    Plot_Nr = 6
    title = 'nonTMD_perc_sim_mean' 
    newfig, savefig, Plot_Nr, plot_nr_in_fig, row_nr, col_nr = dict_organising_subplots[Plot_Nr]
    #if the plot is the last one, the figure should be saved
    #if a new figure should be created (either because the orig is full, or the last TMD is analysed)
    if newfig:
        #clear any previous plots
        plt.close('all')
        #create a new figure
        fig, axarr = plt.subplots(nrows=nrows_in_each_fig, ncols=ncols_in_each_fig, dpi=300)  # sharex=True
    
    for n, df in enumerate(df_list):
        data_column = 'nonTMD_perc_sim_mean'
        color = TUM_colours_list_with_greys[n+2]#"#0489B1"
        alpha = 0.7
        col_width_value = 0.95
        ylabel = 'freq'
        xlabel = data_column
    
        utils.create_hist_from_df_col_with_auto_binlist(df=df,
                                                title=title,row_nr=row_nr,col_nr=col_nr,axarr=axarr,num_bins=60,
                                                settings=settings,
                                                data_column = data_column,
                                                color=color,
                                                alpha=alpha,
                                                col_width_value=col_width_value,
                                                fontsize=fontsize,
                                                xlabel=xlabel,
                                                ylabel=ylabel,legend=protein_list_names
                                                )        

    Plot_Nr = 7
    title = 'nonTMD_perc_sim_plus_ident_mean' 
    newfig, savefig, Plot_Nr, plot_nr_in_fig, row_nr, col_nr = dict_organising_subplots[Plot_Nr]
    #if the plot is the last one, the figure should be saved
    #if a new figure should be created (either because the orig is full, or the last TMD is analysed)
    if newfig:
        #clear any previous plots
        plt.close('all')
        #create a new figure
        fig, axarr = plt.subplots(nrows=nrows_in_each_fig, ncols=ncols_in_each_fig, dpi=300)  # sharex=True
    
    for n, df in enumerate(df_list):
        data_column = 'nonTMD_perc_sim_plus_ident_mean'
        color = TUM_colours_list_with_greys[n+2]#"#0489B1"
        alpha = 0.7
        col_width_value = 0.95
        ylabel = 'freq'
        xlabel = data_column
    
        utils.create_hist_from_df_col_with_auto_binlist(df=df,
                                                title=title,row_nr=row_nr,col_nr=col_nr,axarr=axarr,num_bins=60,
                                                settings=settings,
                                                data_column = data_column,
                                                color=color,
                                                alpha=alpha,
                                                col_width_value=col_width_value,
                                                fontsize=fontsize,
                                                xlabel=xlabel,
                                                ylabel=ylabel,legend=protein_list_names
                                                )          

                                     
    Plot_Nr = 8
    title = 'nonTMD_qm_gaps_per_q_residue_mean' 
    newfig, savefig, Plot_Nr, plot_nr_in_fig, row_nr, col_nr = dict_organising_subplots[Plot_Nr]
    #if the plot is the last one, the figure should be saved
    #if a new figure should be created (either because the orig is full, or the last TMD is analysed)
    if newfig:
        #clear any previous plots
        plt.close('all')
        #create a new figure
        fig, axarr = plt.subplots(nrows=nrows_in_each_fig, ncols=ncols_in_each_fig, dpi=300)  # sharex=True
    
    for n, df in enumerate(df_list):
        data_column = 'nonTMD_qm_gaps_per_q_residue_mean'
        color = TUM_colours_list_with_greys[n+2]#"#0489B1"
        alpha = 0.7
        col_width_value = 0.95
        ylabel = 'freq'
        xlabel = data_column
    
        utils.create_hist_from_df_col_with_auto_binlist(df=df,
                                                title=title,row_nr=row_nr,col_nr=col_nr,axarr=axarr,num_bins=60,
                                                settings=settings,
                                                data_column = data_column,
                                                color=color,
                                                alpha=alpha,
                                                col_width_value=col_width_value,
                                                fontsize=fontsize,
                                                xlabel=xlabel,
                                                ylabel=ylabel,legend=protein_list_names
                                                )   
                                                

    #save the figure as it is
    if savefig:
        utils.save_fig_with_subplots(fig=fig, axarr = axarr, base_filename = base_filename_summ_two_lists,
                                     fig_nr = Plot_Nr, fontsize=fontsize)



    Plot_Nr = 9
    title = 'TM01_AAIMON_ratio_mean' 
    newfig, savefig, Plot_Nr, plot_nr_in_fig, row_nr, col_nr = dict_organising_subplots[Plot_Nr]
    #if the plot is the last one, the figure should be saved
    #if a new figure should be created (either because the orig is full, or the last TMD is analysed)
    if newfig:
        #clear any previous plots
        plt.close('all')
        #create a new figure
        fig, axarr = plt.subplots(nrows=nrows_in_each_fig, ncols=ncols_in_each_fig, dpi=300)  # sharex=True
    
    for n, df in enumerate(df_list):
        data_column = 'TM01_AAIMON_ratio_mean'
        color = TUM_colours_list_with_greys[n+2]#"#0489B1"
        alpha = 0.7
        col_width_value = 0.95
        ylabel = 'freq'
        xlabel = data_column
    
        utils.create_hist_from_df_col_with_auto_binlist(df=df,
                                                title=title,row_nr=row_nr,col_nr=col_nr,axarr=axarr,num_bins=60,
                                                settings=settings,
                                                data_column = data_column,
                                                color=color,
                                                alpha=alpha,
                                                col_width_value=col_width_value,
                                                fontsize=fontsize,
                                                xlabel=xlabel,
                                                ylabel=ylabel,legend=protein_list_names
                                                )                                                   
                                                
 
#
#

                      
    #save the figure as it is
    savefig = True
    if savefig:
        utils.save_fig_with_subplots(fig=fig, axarr = axarr, base_filename = base_filename_summ_two_lists,
                                     fig_nr = Plot_Nr, fontsize=fontsize)
    logging.info('A10_compare_lists is finished')
    logging.info('df1 AAIMON_ratio_mean_all_TMDs : %0.5f' % df1['AAIMON_ratio_mean_all_TMDs'].mean())
    logging.info('df2 AAIMON_ratio_mean_all_TMDs : %0.5f' % df2['AAIMON_ratio_mean_all_TMDs'].mean())
    logging.info('df1 AASMON_ratio_mean_all_TMDs : %0.5f' % df1['AASMON_ratio_mean_all_TMDs'].mean())
    logging.info('df2 AASMON_ratio_mean_all_TMDs : %0.5f' % df2['AASMON_ratio_mean_all_TMDs'].mean())
#        for ax in axarr.flat:                 
#            #change axis font size
#            ax.tick_params(labelsize = fontsize)
#            #hide spines
#            ax.spines["top"].set_visible(False)  
#            ax.spines["right"].set_visible(False) 
#            ax.tick_params(
#                axis='x',          # changes apply to the x-axis
#                which='both',      # both major and minor ticks are affected
#                bottom='off',      # ticks along the bottom edge are off
#                top='off',         # ticks along the top edge are off
#                labelbottom='on') # labels along the bottom edge are off    
#        #automatically tighten the layout of plots in the figure
#        fig.subplots_adjust(bottom = 0)
#        fig.subplots_adjust(top = 1)
#        fig.subplots_adjust(right = 1)
#        fig.subplots_adjust(left = 0)
#        #save files
#        fig.savefig(base_filename_summ_two_lists + '_%01d.png' % Plot_Nr, format='png', dpi=400)
#        fig.savefig(base_filename_summ_two_lists + '_%01d.pdf' % Plot_Nr, format='pdf')
    
'''+++++++++++++++TMD CONSERVATION (OLD)++++++++++++++++++'''

#A## variables are included only to help navigate the document in PyCharm
B01_OLD_calculate_TMD_conservation = settings["run_settings"]["calculate_TMD_conservation"]["calculate_TMD_conservation"]
if B01_OLD_calculate_TMD_conservation:
    #convert file to dataframe (eg List70_simap.csv)
    df_dfout05_simapcsv = pd.read_csv(dfout05_simapcsv, sep=",", index_col=0,
                                      quoting=csv.QUOTE_NONNUMERIC)  # index_col = 0, only wher reopening a saved pandas file
    #add desired columns (only necessary if for loop method is kept)
    original_columns = list(df_dfout05_simapcsv.columns)
    new_columns = ['md5_checksums_in_homologue_list_that_are_also_in_query_list']
    #columns_added_after_SIMAP_analysis as above
    new_unique_column_list = set(original_columns + new_columns)
    #add extra columns
    df_dfout05_simapcsv = df_dfout05_simapcsv.reindex(index=df.index, columns=new_unique_column_list)
    #sort columns
    df_dfout05_simapcsv = df_dfout05_simapcsv.sort_index(axis=1)

    #specify some data types:
    df_dfout05_simapcsv['kept_after_redundancy_check'] = \
        df_dfout05_simapcsv['kept_after_redundancy_check'].astype(np.bool)
    df_dfout05_simapcsv['protein_kept_for_statistical_analysis'] = \
        df_dfout05_simapcsv['protein_kept_for_statistical_analysis'].astype(np.bool)
    df_dfout05_simapcsv['query_md5'] = df_dfout05_simapcsv['query_md5'].astype('<U32')
    df_dfout05_simapcsv['md5_checksums_in_homologue_list_that_are_also_in_query_list'] = \
        df_dfout05_simapcsv[
            'md5_checksums_in_homologue_list_that_are_also_in_query_list'].astype('<U3000')

    #create a list of the files to ba analysed from the original csv file
    #list_files_with_features, list_files_with_homologues, list_of_protein_names, list_of_organism_domains = create_list_of_files_from_csv_with_uniprot_data(csv_file_with_uniprot_data, list_of_keys_to_keep)

    '''OLD BIN CREATION'''
    #create the bins used for the histogram data. I add 3 empty field 0,0,0 to account for the three columns
    #The first bin is 0. The bins range UPWARDS, so 0 = =-0.4,
    #    list_of_bins_for_histogram = []
    #    for i in range(20,250,10):
    #        i2 = i/100.0
    #        list_of_bins_for_histogram.append(i2)
    #    list_of_bins_for_histogram.append(30)
    #logging.info(list_of_bins_for_histogram)

    #    list_of_bins_for_histogram2 = np.linspace(0.2,2.5,24)
    #    logging.info(list_of_bins_for_histogram2)

    #use linspace to get a fixid number of points between tha min and the max for the histogram
    list_of_bins_for_histogram = np.linspace(settings["hist_settings_single_protein"]["smallest_bin"],
                                             settings["hist_settings_single_protein"]["largest_bin"],
                                             settings["hist_settings_single_protein"]["number_of_bins"])
    #add 30 as the last bin, to make sure 100% of the data is added to the histogram, including major outliers
    list_of_bins_for_histogram = np.append(list_of_bins_for_histogram,
                                           settings["hist_settings_single_protein"]["final_highest_bin"])

    list_of_keys_for_hit_identity_analysis = ['match_TMD_added_to_FastA_alignment', 'ratio_ident_mem_to_nonmem',
                                              'match_TMD_kept_for_statistical_analysis']

    #minimum_number_hits_for_data_analysis = 20

    number_query_seq_processed = 0
    number_query_seq_added_to_histogram = 0
    list_of_proteins_kept_for_statistical_analysis = []
    list_of_proteins_whose_XML_file_doesnt_exist = []
    list_of_proteins_with_not_enough_valid_hits = []
    list_of_csv_files_with_homologues_that_are_old = []
    dict_of_arrays_mem_nonmem_ratios = {}

    for i in df_dfout05_simapcsv.index:
        #take the organism domain (Eukaryota, Bacteria, Archaea) from the full organism classification list
        first_two_letters_of_uniprot_acc = df_dfout05_simapcsv.loc[i, 'first_two_letters_of_uniprot_acc']
        protein_name = df_dfout05_simapcsv.loc[i, 'A2_protein_name']
        SIMAP_feature_table_XML_file_path = os.path.join(df.loc[i, 'simap_filename_base'], first_two_letters_of_uniprot_acc,
                                                         '%s_feature_table.xml' % protein_name)
        SIMAP_homologues_XML_file_path = os.path.join(df.loc[i, 'simap_filename_base'], first_two_letters_of_uniprot_acc,
                                                      '%s_homologues.xml' % protein_name)
        #for i in range(len(list_of_protein_names)):
        #print dots to show that program is running
        sys.stdout.write(" .")
        sys.stdout.flush()
        csv_file_with_hit_details = os.path.join(df.loc[i, 'simap_filename_base'], first_two_letters_of_uniprot_acc,
                                                 '%s_homologues.csv' % protein_name)

        #check if the feature table and homologue files actually exist
        if os.path.isfile(csv_file_with_hit_details):
            csv_file_with_hit_details_exists = True
            #get the size of the file (actual size, not windows "size on disk")
            statinfo = os.stat(csv_file_with_hit_details)
            csv_file_with_hit_details_filesize = statinfo.st_size
            if csv_file_with_hit_details_filesize > 2400:
                csv_file_with_hit_details_contains_hits = True
            else:
                csv_file_with_hit_details_contains_hits = False
        else:
            csv_file_with_hit_details_exists = False
            list_of_proteins_whose_XML_file_doesnt_exist.append(protein_name)
            #logging.info('%s_homologues.csv file is missing' % protein_name)
        if csv_file_with_hit_details_exists and csv_file_with_hit_details_contains_hits:
            #check if any of the files are old: in the old files, the headers had spaces in between each word, rather than an underscore. Later there might be other woys to search for older csv files with homologue data.
            csv_file_with_homologues_is_old = False
            with open(csv_file_with_hit_details, 'r') as f:
                reader = csv.reader(f)
                rownumber = 0
                for row in reader:
                    if rownumber == 0:
                        if row[0] == 'hit number':
                            csv_file_with_homologues_is_old = True
                            list_of_csv_files_with_homologues_that_are_old.append(protein_name)
            if csv_file_with_homologues_is_old:
                logging.warning("%s csv file is old!! Need to recreate csv from XML file of homologues (%s,%s)" % (
                    protein_name[:6], df_dfout05_simapcsv.loc[i, 'SIMAP_homologues_XML_file_path'],
                    df_dfout05_simapcsv.loc[i, 'SIMAP_csv_analysed_path']))
                #The uniprot TMD start and stop is necessary. Take this from the csv with uniprot data. Thi file and list_of_keys_to_keep should be global objects.
                "since parse_single_homologue_XML_to_csv_and_fasta is no longer  function, "
                "the old file cannot simply be updated here"
            #                nested_dict_with_uniprot_seq = create_nested_dict_from_csv(csv_file_with_uniprot_data, list_of_keys_to_keep)
            #                #since this is a nested dictionary, I have to scroll through it until I obtain the correct protein:
            #                for i in range(1,len(nested_dict_with_uniprot_seq)+1):
            #                    protein_in_csv_cell = df_dfout05_simapcsv.loc[i,'A1_uniprot_accession']
            #                    print(protein_in_csv_cell)
            #                    print(protein_name[:6])
            #
            #                    if protein_in_csv_cell == protein_name[:6]:
            #                        row_in_csv_with_protein = i
            #now update the csv with homologues, just for this sequence
            #                parse_single_homologue_XML_to_csv_and_fasta(df_dfout05_simapcsv.loc[i,'SIMAP_homologues_XML_file_path'],df_dfout05_simapcsv.loc[i,'SIMAP_csv_analysed_path'], nested_dict_with_uniprot_seq[row_in_csv_with_protein])

            #if not csv_file_with_homologues_is_old: #xxx hopefully they auto-update now!!!
            '''Use Pandas for the analysis of the hits, rather than the nested dictionary approach used previously
            '''
            df_SIMAP_hits = pd.read_csv(csv_file_with_hit_details, sep=",", na_values=['nan'],
                                        quoting=csv.QUOTE_NONNUMERIC)  # index_col = 0, only wher reopening a saved pandas file
            #convert all np.nan (empty values) to empty strings
            #for col in df_SIMAP_hits.columns[pd.isnull(df_SIMAP_hits).all()]:
            #    df_SIMAP_hits[col] = df_SIMAP_hits[col].astype(object).fillna("UNKNOWN")
            #df_SIMAP_hits.fillna(value='empty')
            #            original_columns_df_SIMAP_hits = list(df_SIMAP_hits.columns)
            #            added_columns_df_SIMAP_hits = ['match_kept', 'free_column']
            #            new_column_list = list(original_columns_df_SIMAP_hits + added_columns_df_SIMAP_hits)
            #            df_SIMAP_hits = df_SIMAP_hits.reindex(index = df_SIMAP_hits.index, columns = new_column_list)
            #            df_SIMAP_hits.loc[:,'match_kept'] = False
            #            df_SIMAP_hits['match_kept'].astype = np.bool_
            list_of_columns_to_fill_empty_values = ['disallowed_words_in_description']
            for col in list_of_columns_to_fill_empty_values:
                df_SIMAP_hits[col] = df_SIMAP_hits[col].astype(object).fillna("none")


            #filter each homologue from SIMAP
            number_of_gaps_allowed_in_match_TMD = settings["simap_match_filters"]["number_of_gaps_allowed_in_match_TMD"]
            number_of_gaps_allowed_in_query_TMD = settings["simap_match_filters"]["number_of_gaps_allowed_in_query_TMD"]
            minimum_identity_of_full_protein = settings["simap_match_filters"]["minimum_identity_of_full_protein"]
            min_identity_of_TMD_final_filter = settings["simap_match_filters"]["min_identity_of_TMD_final_filter"]
            for j in range(len(df_SIMAP_hits)):
                df_SIMAP_hits.loc[j, 'match_TMD_kept_for_statistical_analysis'] = True if all([
                    df_SIMAP_hits.loc[j, 'FASTA_expectation'] <= settings['simap_match_filters']['e_value_filter'],
                    df_SIMAP_hits.loc[j, 'disallowed_words_in_description'] == 'none',
                    settings['simap_match_filters']['database'] == 'all',
                    df_SIMAP_hits.loc[j, 'number_of_gaps_in_match_TMD'] <= number_of_gaps_allowed_in_match_TMD,
                    df_SIMAP_hits.loc[j, 'number_of_gaps_in_query_TMD'] <= number_of_gaps_allowed_in_query_TMD,
                    df_SIMAP_hits.loc[j, 'SW_identity'] >= minimum_identity_of_full_protein,
                    df_SIMAP_hits.loc[j, 'percentage_identity_of_TMD'] >= min_identity_of_TMD_final_filter,
                    #df_SIMAP_hits.loc[j, 'number_of_X_in_TMD'] <= settings["simap_match_filters"]["number_of_X_allowed_in_TMD"],
                    df_SIMAP_hits.loc[j, 'SW_match_TMD_seq'] != 'TMD_not_in_SW_alignment'
                ]) else False
                #len(df_SIMAP_hits.loc[j,'query_aln_seq_excl_TMD']) >= settings["simap_match_filters"]["min_len_query_aln_seq_excl_TMD"]
                #df_SIMAP_hits.loc[j,'number_of_X_in_match_seq'] <= settings["simap_match_filters"]["number_of_X_allowed_in_seq"],

            #simply select the true data in a smaller dataframe
            #use the value_counts method to count the number of hits kept for statistical analysis
            df_SIMAP_homologues_kept_for_statistical_analysis_stat_only = df_SIMAP_hits[df_SIMAP_hits.match_TMD_kept_for_statistical_analysis == True]
            #both stat and fasta filters. NEEDS TO BE UPDATED LATER TO ONLY REFLECT WORKING STAT FILTERS
            df_SIMAP_homologues_kept_for_statistical_analysis = \
                df_SIMAP_homologues_kept_for_statistical_analysis_stat_only[
                    df_SIMAP_homologues_kept_for_statistical_analysis_stat_only.match_TMD_added_to_FastA_alignment == True]

            with open(df.loc[i, 'csv_SIMAP_homologues_kept_for_statistical_analysis'], 'w') as f:
                df_SIMAP_homologues_kept_for_statistical_analysis.to_csv(f, sep=",", quoting=csv.QUOTE_NONNUMERIC)

                #NOT working! gives error whel len(df) = 0
                #number_of_hits_kept_for_statistical_analysis = df_SIMAP_hits['match_TMD_kept_for_statistical_analysis'].value_counts()[True]
                #now check that the DataFrame really only holds the selected data (match_TMD_kept_for_statistical_analysis).
                #if number_of_hits_kept_for_statistical_analysis != len(df_SIMAP_homologues_kept_for_statistical_analysis):
                #logging.warning('number_of_hits_kept_for_statistical_analysis != len(df_SIMAP_homologues_kept_for_statistical_analysis) for %s)' % protein_name)

            #check if the csv_file_av_cons_ratios_hits file already exists
            csv_mem_nonmem_ratios_is_old = False
            if os.path.isfile(df.loc[i, 'csv_file_av_cons_ratios_hits']):
                #load into a pandas dataframe
                df_mem_nonmem_ratios = pd.read_csv(df.loc[i, 'csv_file_av_cons_ratios_hits'],
                                                   sep=",", index_col=0, quoting=csv.QUOTE_NONNUMERIC)
                gap_penalties_in_file = list(df_mem_nonmem_ratios.index)
                matrices_in_file = list(df_mem_nonmem_ratios.columns)
                if gap_penalties_in_file == list(range(settings["aa_substitution_scoring"]["gap_open_penalty_min"],
                                                       settings["aa_substitution_scoring"]["gap_open_penalty_max"],
                                                       settings["aa_substitution_scoring"][
                                                           "gap_open_penalty_increment"])) and \
                                matrices_in_file == settings['aa_substitution_scoring']['matrices']:
                    csv_mem_nonmem_ratios_is_old = False
                else:
                    csv_mem_nonmem_ratios_is_old = True
                    logging.info('%s is old, repeating calculation' % df.loc[
                        i, 'csv_file_av_cons_ratios_hits'])

            if settings["run_settings"]["calculate_TMD_conservation"]["overwrite_csv_file_av_cons_ratios_hits"] or csv_mem_nonmem_ratios_is_old:
                #if not os.path.isfile(df.loc[i,'csv_file_av_cons_ratios_hits']):
                #create a new dataframe to hold the array of mem/nonmem ratios for the varying matrices and  gap penalties
                df_mem_nonmem_ratios = pd.DataFrame(0.0,
                                                    index=range(
                                                        settings["aa_substitution_scoring"]["gap_open_penalty_min"],
                                                        settings["aa_substitution_scoring"]["gap_open_penalty_max"],
                                                        settings["aa_substitution_scoring"][
                                                            "gap_open_penalty_increment"]),
                                                    columns=settings['aa_substitution_scoring']['matrices'])

                #load the amino acid substitution matrices from the settings file
                list_of_aa_sub_matrices = settings['aa_substitution_scoring']['matrices']
                dict_of_aa_matrices = {key: eval(key) for key in list_of_aa_sub_matrices}

                #for each gap penalty
                for k in range(settings["aa_substitution_scoring"]["gap_open_penalty_min"],
                               settings["aa_substitution_scoring"]["gap_open_penalty_max"],
                               settings["aa_substitution_scoring"]["gap_open_penalty_increment"]):
                    gap_open_penalty = k
                    #for each aa sub matrix, represented as a key in the dictionary
                    for key in dict_of_aa_matrices:
                        matrix_name = repr(key).replace("'", "")
                        column_name = 'sim_ratio_%s_gapo%i' % (matrix_name.replace("'", "")[0:-7], k)
                        mean_sim_ratio = df_SIMAP_homologues_kept_for_statistical_analysis[column_name].mean()
                        #if I wanted to, I could save a histogram for each protein!!
                        #hist_sim_ratio = df_SIMAP_homologues_kept_for_statistical_analysis[column_name].hist()

                        df_mem_nonmem_ratios.loc[gap_open_penalty, matrix_name] = mean_sim_ratio

                #save as a separate file?
                with open(df.loc[i, 'csv_file_av_cons_ratios_hits'], 'w') as f:
                    df_mem_nonmem_ratios.to_csv(f, sep=",", quoting=csv.QUOTE_NONNUMERIC)
                    print('%s created' % df.loc[i, 'csv_file_av_cons_ratios_hits'])
            else:
                #csv_file_av_cons_ratios_hits already created, load from file
                logging.info('%s already exists, moving to next' % df.loc[
                    i, 'csv_file_av_cons_ratios_hits'])
                df_mem_nonmem_ratios = pd.read_csv(df.loc[i, 'csv_file_av_cons_ratios_hits'],
                                                   sep=",", quoting=csv.QUOTE_NONNUMERIC)

            df_dfout05_simapcsv.loc[i, 'number_of_valid_hits'] = len(
                df_SIMAP_homologues_kept_for_statistical_analysis)
            df_dfout05_simapcsv.loc[
                i, 'df_mem_nonmem_ratios'] = df_mem_nonmem_ratios.to_string()
            #df_dfout05_simapcsv.loc[i,'df_mem_nonmem_ratios'] = df_mem_nonmem_ratios
            #add to the uniprot df/csv the mem/nonmem and stdev calculated from the SW alignment and counting the identical residues
            df_dfout05_simapcsv.loc[i, 'mean_ratio_ident_mem_to_nonmem'] = \
                df_SIMAP_homologues_kept_for_statistical_analysis['ratio_ident_mem_to_nonmem'].mean()
            df_dfout05_simapcsv.loc[i, 'stdev_ratio_ident_mem_to_nonmem'] = \
                df_SIMAP_homologues_kept_for_statistical_analysis['ratio_ident_mem_to_nonmem'].std()

            #nested_dict_for_hit_identity_analysis = create_nested_dict_from_csv(csv_file_with_hit_details, 'all')#xxx had problems with KeyError: 'match_TMD_added_to_FastA_alignment',
            nested_dict_for_hit_identity_analysis = utils.create_nested_dict_from_csv(csv_file_with_hit_details,
                                                                                      list_of_keys_for_hit_identity_analysis)

            #create a list of the ratio_ident_mem_to_nonmem for each hit, to be used in the histogram for each protein
            list_ratio_ident_mem_to_nonmem = []
            number_of_valid_hits = 0
            for j in range(1, len(nested_dict_for_hit_identity_analysis)):
                #if the TMD from this hit was good enough to add to the alignment
                match_TMD_added_to_FastA_alignment = utils.convert_string_to_boolean_value(
                    nested_dict_for_hit_identity_analysis[j]['match_TMD_added_to_FastA_alignment'])
                match_TMD_kept_for_statistical_analysis = utils.convert_string_to_boolean_value(
                    nested_dict_for_hit_identity_analysis[j]['match_TMD_kept_for_statistical_analysis'])

                #logging.info(match_added)
                #start with the second hit, as the first is the query sequence
                if j >= 2:
                    if df_SIMAP_hits.loc[j, 'match_TMD_kept_for_statistical_analysis']:
                        if nested_dict_for_hit_identity_analysis[j]['ratio_ident_mem_to_nonmem'] != '':
                            number_of_valid_hits += 1
                            ratio_ident_mem_to_nonmem = float(
                                nested_dict_for_hit_identity_analysis[j]['ratio_ident_mem_to_nonmem'])
                            list_ratio_ident_mem_to_nonmem.append(ratio_ident_mem_to_nonmem)
            #logging.info(list_ratio_ident_mem_to_nonmem)

            #determine if protein is kept for statistical analysis. Currently the only filter is the number of valid hits.
            if df_dfout05_simapcsv.loc[i, 'number_of_valid_hits'] >= \
                    settings["hist_settings_single_protein"]["minimum_number_hits_for_data_analysis"]:
                df_dfout05_simapcsv.loc[i, 'protein_kept_for_statistical_analysis'] = True
                #if there is enough valid hits, add the protein name to the list, so the data can be used later
                list_of_proteins_kept_for_statistical_analysis.append(protein_name)
            else:
                df_dfout05_simapcsv.loc[i, 'protein_kept_for_statistical_analysis'] = False
                list_of_proteins_with_not_enough_valid_hits.append(
                    df_dfout05_simapcsv.loc[i, 'A1_uniprot_accession'])

            logging.info('%s\tnumber_of_valid_hits: %i, protein_kept_for_statistical_analysis: %s' % (
                protein_name, df_dfout05_simapcsv.loc[i, 'number_of_valid_hits'],
                df_dfout05_simapcsv.loc[i, 'protein_kept_for_statistical_analysis']))
            '''Previously, a histogram was created for each protein
            '''
            if df_dfout05_simapcsv.loc[i, 'protein_kept_for_statistical_analysis']:
                #numpy will automatically create a histogram
                hist_ratio_ident_mem_to_nonmem = np.histogram(list_ratio_ident_mem_to_nonmem,
                                                              bins=list_of_bins_for_histogram)
                # If True, the result is the value of the probability density function at the bin, normalized such that the integral over the range is 1.
                hist_ratio_ident_mem_to_nonmem_normalised = np.histogram(list_ratio_ident_mem_to_nonmem, density=True,
                                                                         bins=list_of_bins_for_histogram)
                #logging.info(protein_name, hist)

                #prepare csv header
                #this bin list should be the same as the list_of_bins_for_histogram
                bin_list = np.array(hist_ratio_ident_mem_to_nonmem)[1].tolist()
                list_of_bins_as_strings = []
                for i8 in range(len(bin_list) - 1):
                    bin_string = '%.1f-%.1f' % (bin_list[i8], bin_list[i8 + 1])
                    list_of_bins_as_strings.append(bin_string)

                csv_added_header = ['protein_name', 'organism_domain', 'number_of_valid_hits',
                                    'mean_ratio_ident_mem_to_nonmem', 'stdev_ratio_ident_mem_to_nonmem']
                header_with_bins_for_csv = csv_added_header + list_of_bins_as_strings

                #for the first file, create a new file and save the csv header
                if number_query_seq_added_to_histogram == 0:
                    utils.save_list_as_row_in_csv(header_with_bins_for_csv, csv_file_with_histogram_data, 'w')
                    utils.save_list_as_row_in_csv(header_with_bins_for_csv, csv_file_with_histogram_data_normalised,'w')

                #rpitrh = ratio_ident_mem_to_nonmem
                array_rpitrh = np.array(hist_ratio_ident_mem_to_nonmem)
                hist_ratio_ident_mem_to_nonmem_list = np.array(hist_ratio_ident_mem_to_nonmem)[0].tolist()
                hist_ratio_ident_mem_to_nonmem_normalised_list = np.array(hist_ratio_ident_mem_to_nonmem_normalised)[
                    0].tolist()

                #calculate mean and stdev
                mean_ratio_ident_mem_to_nonmem = np.mean(list_ratio_ident_mem_to_nonmem)
                stdev_ratio_ident_mem_to_nonmem = np.std(list_ratio_ident_mem_to_nonmem)

                #prepare data as list for csv
                list_protein_name = [protein_name, organism_domain, number_of_valid_hits,
                                     mean_ratio_ident_mem_to_nonmem, stdev_ratio_ident_mem_to_nonmem]
                hist_ratio_ident_mem_to_nonmem_csv_row = list_protein_name + hist_ratio_ident_mem_to_nonmem_list
                hist_ratio_ident_mem_to_nonmem_normalised_csv_row = list_protein_name + hist_ratio_ident_mem_to_nonmem_normalised_list

                #save data to csv
                utils.save_list_as_row_in_csv(hist_ratio_ident_mem_to_nonmem_csv_row, csv_file_with_histogram_data, 'a')
                utils.save_list_as_row_in_csv(hist_ratio_ident_mem_to_nonmem_normalised_csv_row,
                                              csv_file_with_histogram_data_normalised, 'a')

                number_query_seq_added_to_histogram += 1
        number_query_seq_processed += 1
    logging.info('\n\nSUMMARY OF TMD CONSERVATION ANALYSIS:\n%i query sequences processed' % number_query_seq_processed)
    logging.info(
        '%i sequences had enough valid hits to be added added to histogram (rpitrh)' % number_query_seq_added_to_histogram)
    logging.info('%i files excluded:' % (number_query_seq_processed - number_query_seq_added_to_histogram))
    logging.info('\t%i XML files missing:\n\t%s' % (
        len(list_of_proteins_whose_XML_file_doesnt_exist), list_of_proteins_whose_XML_file_doesnt_exist))
    logging.info('\t%i proteins had less than %s valid homologues and were excluded from analysis:\n\t%s' % (
        len(list_of_proteins_with_not_enough_valid_hits),
        settings["hist_settings_single_protein"]["minimum_number_hits_for_data_analysis"],
        list_of_proteins_with_not_enough_valid_hits))
    logging.info('list_of_csv_files_with_homologues_that_are_old: %s' % list_of_csv_files_with_homologues_that_are_old)

    '''
     REMOVE REDUNDANT/DUPLICATE RECORDS
    '''

    '''new pandas method: take the list of md5 directly from the pandas dataframe
    '''

    #make a list of the md5 checksums, including nan values
    list_of_query_md5_checksums = []
    for ii2 in df_dfout05_simapcsv.index:
        list_of_query_md5_checksums.append(df_dfout05_simapcsv.loc[ii2, 'query_md5'])

    #print(list_of_query_md5_checksums)
    #convert to pandas data series
    series_of_query_md5_checksums = pd.Series(list_of_query_md5_checksums)
    #remove nan values
    series_of_query_md5_checksums_without_nan_values = series_of_query_md5_checksums.dropna()
    #convert back to a list
    list_of_query_md5_checksums = list(series_of_query_md5_checksums_without_nan_values)
    if 'nan' in list_of_query_md5_checksums:
        list_of_query_md5_checksums.remove('nan')
    #remove any empty values (note that the filter function tries to link to elementtree, and not the basic python function)
    list_of_query_md5_checksums = [x for x in list_of_query_md5_checksums if x]
    print('list_of_query_md5_checksums: %s' % list_of_query_md5_checksums)
    '''
    For each csv file containing all the homologues, create a list of the md5s (md5_list_from_simap_hits).
    Then use the 'set' function to identify which md5 checksums are shared by both the list of query proteins,
    and the list of homologue hits for that particular protein.
    '''
    #create an empty array to hold the data. For the data type thecolumn containing the list of md5s will be large, so use unicode strings with at least 600 characters (>U300)
    #    number_of_columns = 6
    #    dtype_array_with_query_md5_data = [('number', '<i4'), ('query_name', '<U32'),
    #                                       ('query_md5', '<U33'),
    #                                       ('md5_checksums_in_homologue_list_that_are_also_in_query_list', '<U3000'),
    #                                       ('number_of_valid_hits', '<U32'),
    #                                       ('kept_after_redundancy_check', np.bool)]
    #    array_with_query_md5_data = np.ones(len(df_dfout05_simapcsv),
    #                                        dtype = dtype_array_with_query_md5_data)
    #logging.info(array_with_query_md5_data)
    '''Create a list of all the md5 checksums of all SIMAP hits
    '''
    list_of_proteins_with_damaged_file_0_valid_hits = []
    list_of_proteins_with_no_csv_file_with_hit_details = []

    for i in df_dfout05_simapcsv.index:
        csv_file_with_hit_details = df_dfout05_simapcsv.loc[
            i, 'SIMAP_csv_analysed_path']
        #if the file exists
        if os.path.isfile(csv_file_with_hit_details):
            #create a dataframe from the csv file, replace empty values with the tewt 'nan'
            df_SIMAP_hits = pd.read_csv(csv_file_with_hit_details, sep=",", na_values=['nan'],
                                        quoting=csv.QUOTE_NONNUMERIC)
            #if there are actually hits, and the SIMAP download was sucessful
            if len(df_SIMAP_hits) != 0:
                #create a list of all the md5 values from the simap homolague file
                md5_list_from_simap_hits = list(df_SIMAP_hits['A3_md5'])

                query_md5 = df_dfout05_simapcsv.loc[i, 'query_md5']
                #use the .intersection function to identify other query sequences that are found as hits in the homologues for this query
                #create a set, exclude the first sequence, as that is the query
                set_md5_list_from_simap_hits = set(md5_list_from_simap_hits[1:])
                #find elemnts common to the set and the list
                md5_checksums_in_homologue_list_that_are_also_in_query_list = list(
                    set_md5_list_from_simap_hits.intersection(set(list_of_query_md5_checksums)))
                #add data to the array, so it can be saved as a csv.
                #            df_dfout05_simapcsv.loc[i,'number'] = i
                #            df_dfout05_simapcsv.loc[i,'query_name'] = protein_name
                #            df_dfout05_simapcsv.loc[i,'query_md5'] = query_md5
                #for the list of redundant md5s, replace commas with | to improve compatibility with csv files. Where there are no redundant md5 accessions, leave blank
                if md5_checksums_in_homologue_list_that_are_also_in_query_list != []:
                    df_dfout05_simapcsv.loc[
                        i, 'md5_checksums_in_homologue_list_that_are_also_in_query_list'] = ' | '.join(
                        md5_checksums_in_homologue_list_that_are_also_in_query_list)
                else:
                    df_dfout05_simapcsv.loc[
                        i, 'md5_checksums_in_homologue_list_that_are_also_in_query_list'] = ''
                    #logging.info(protein_name)
                    #logging.info(md5_list_from_simap_hits)
                    #logging.info('md5_checksums_found_in_both_lists: %s' % md5_checksums_found_in_both_lists)
            else:
                #if there are no valid hits, replace the np.nan with 0, add the uniprot accession to a list so that the SIMAP retrieval can be repeated
                df_dfout05_simapcsv.loc[
                    i, 'md5_checksums_in_homologue_list_that_are_also_in_query_list'] = ''
                list_of_proteins_with_damaged_file_0_valid_hits.append(
                    df_dfout05_simapcsv.loc[i, 'A1_uniprot_accession'])
                logging.info('file damaged or no hits (len(df_SIMAP_hits) = 0 : %s)' % csv_file_with_hit_details)
        else:
            #the file is not there, replace the np.nan with an empty field
            df_dfout05_simapcsv.loc[i, 'query_md5'] = ''
            list_of_proteins_with_no_csv_file_with_hit_details.append(
                df_dfout05_simapcsv.loc[i, 'A1_uniprot_accession'])
            logging.info('file does not exist: %s' % csv_file_with_hit_details)
    logging.info(
        'list_of_proteins_with_damaged_file_0_valid_hits (repeating SIMAP download is necessary) : %s' % list_of_proteins_with_damaged_file_0_valid_hits)
    logging.info('list_of_proteins_with_no_csv_file_with_hit_details (if len>%i, try repeating SIMAP download) : %s' % (
        settings["simap_homologue_download"]["max_query_sequence_length"],
        list_of_proteins_with_no_csv_file_with_hit_details))

    #specify data type, repeat of earlier specification for df
    df_dfout05_simapcsv['query_md5'] = df_dfout05_simapcsv[
        'query_md5'].astype('<U32')
    df_dfout05_simapcsv['md5_checksums_in_homologue_list_that_are_also_in_query_list'] = \
        df_dfout05_simapcsv[
            'md5_checksums_in_homologue_list_that_are_also_in_query_list'].astype('<U3000')

    #convert the md5s from a string (separated by |) to a list.
    for i in df_dfout05_simapcsv.index:
        if df_dfout05_simapcsv.loc[
            i, 'md5_checksums_in_homologue_list_that_are_also_in_query_list'] != '':
            md5_checksums_in_homologue_list_that_are_also_in_query_list = \
                df_dfout05_simapcsv.loc[
                    i, 'md5_checksums_in_homologue_list_that_are_also_in_query_list'].split(' | ')
    '''
    Save the updated df with SIMAP data later??
    '''
    #innocent until proven guilty! all seqs (at least those with SIMAP data) are nonredundant by default
    df_dfout05_simapcsv['kept_after_redundancy_check'].astype(np.bool)
    #df_dfout05_simapcsv.loc[:,'kept_after_redundancy_check'] = True
    #label all with some valid hits nonredundant, all those with valid hits are innocent until proven guilty
    df_dfout05_simapcsv['kept_after_redundancy_check'] = \
        df_dfout05_simapcsv['SIMAP_total_hits'] > 0

    #%%
    #do the redundancy check by loading all of the md5s from SIMAP homologues for each protein, and checking against the md5s of all of the query sequences in that list [note.. I had many troubles clearing these of np.nan!]
    for i in df_dfout05_simapcsv.index:
        #only examine the files with redundant sequences to see which has the most valid hits
        if df_dfout05_simapcsv.loc[
            i, 'md5_checksums_in_homologue_list_that_are_also_in_query_list'] != '':
            #split the md5s into a list
            md5_checksums_in_homologue_list_that_are_also_in_query_list = \
                df_dfout05_simapcsv.loc[
                    i, 'md5_checksums_in_homologue_list_that_are_also_in_query_list'].split(' | ')
            #add the original query md5 to the list
            md5_checksums_in_homologue_list_that_are_also_in_query_list.append(
                df_dfout05_simapcsv.loc[i, 'query_md5'])
            #logging.info(md5_checksums_in_homologue_list_that_are_also_in_query_list)

            #create a new subarray to hold the data from the redundant md5s
            #subarray = np.zeros(len(md5_checksums_in_homologue_list_that_are_also_in_query_list))

            #array_of_redundant_md5s = np.array(md5_checksums_in_homologue_list_that_are_also_in_query_list)
            #array_with_redundant_md5_data = np.ones(len(array_of_redundant_md5s)*4, dtype = [('redundant_md5', '<U33'), ('row_in_df_dfout05_simapcsv', '<i4'), ('number_of_valid_hits', '<i4'), ('kept_after_redundancy_check', 'b')]).reshape(len(array_of_redundant_md5s),4)
            #create an array, just for those redundant md5s
            array_with_redundant_md5_data = np.ones(len(md5_checksums_in_homologue_list_that_are_also_in_query_list),
                                                    dtype=[('query_name', '<U32'), ('redundant_md5', '<U33'), (
                                                        'row_in_df_dfout05_simapcsv', '<i4'),
                                                           ('number_of_valid_hits', '<i4'),
                                                           ('kept_after_redundancy_check', '<U5')])
            #for readability, convert to a pandas dataframe
            df_redundant_md5 = pd.DataFrame(array_with_redundant_md5_data)
            '''New Pandas dataframe: transfer directly to df with uniprot data
            '''
            #if the kept_after_redundancy_check is currently labelled "True", due to some valid hits
            if df_dfout05_simapcsv.loc[i, 'kept_after_redundancy_check'] == True:
                #simply add the mds to the dataframe
                df_redundant_md5['redundant_md5'] = md5_checksums_in_homologue_list_that_are_also_in_query_list
                #set up the subarray
                for h in df_redundant_md5.index:
                    redundant_md5 = df_redundant_md5.loc[h, 'redundant_md5']
                    #find the location of the redundant md5 in the dataframe of query sequences
                    row_in_df_dfout05_simapcsv = int(
                        df_dfout05_simapcsv[
                            df_dfout05_simapcsv['query_md5'] == redundant_md5].index)
                    df_redundant_md5.loc[
                        h, 'row_in_df_dfout05_simapcsv'] = row_in_df_dfout05_simapcsv
                    #add the number of valid hits to the subarray, for readability
                    df_redundant_md5.loc[h, 'number_of_valid_hits'] = \
                        df_dfout05_simapcsv.loc[
                            row_in_df_dfout05_simapcsv, 'number_of_valid_hits']
                    df_redundant_md5.loc[h, 'query_name'] = df_dfout05_simapcsv.loc[
                        row_in_df_dfout05_simapcsv, 'A2_protein_name']
                #find the protein seq with the most hits
                index_maximum_number_of_valid_hits = df_redundant_md5.loc[:, 'number_of_valid_hits'].argmax()
                #label them as 'not kept' as default
                df_redundant_md5.loc[:, 'kept_after_redundancy_check'] = False
                #label the best as 'kept'
                df_redundant_md5.loc[index_maximum_number_of_valid_hits, 'kept_after_redundancy_check'] = True
                #now transfer this data to the original array (df_dfout05_simapcsv)
                for m in df_redundant_md5.index:
                    row = df_redundant_md5.loc[m, 'row_in_df_dfout05_simapcsv']
                    kept_after_redundancy_check = df_redundant_md5.loc[m, 'kept_after_redundancy_check']
                    df_dfout05_simapcsv.loc[
                        row, 'kept_after_redundancy_check'] = kept_after_redundancy_check


    ##if this particular query has not already been rejected as a redundant sequence with a lower number of valid hits, continue and analyse the redundant md5s
    #if df_dfout05_simapcsv.loc[i,'kept_after_redundancy_check'] != 'False':
    #
    #    for i in range(len(md5_checksums_in_homologue_list_that_are_also_in_query_list)):
    #        #logging.info(array_with_redundant_md5_data[i][0])
    #        #add the redundant md5
    #        redundant_md5 = md5_checksums_in_homologue_list_that_are_also_in_query_list[i]
    #        array_with_redundant_md5_data[i]['redundant_md5'] = redundant_md5
    #        #find the location of the redundant md5 in the array of query sequences
    #        row_in_df_dfout05_simapcsv = list_of_query_md5_checksums.index(array_with_redundant_md5_data[i]['redundant_md5'])
    #        array_with_redundant_md5_data[i]['row_in_df_dfout05_simapcsv'] = row_in_df_dfout05_simapcsv
    #        #use this row_in_df_dfout05_simapcsv to find the number of valid hits
    #        array_with_redundant_md5_data[i]['number_of_valid_hits'] = array_with_query_md5_data[row_in_df_dfout05_simapcsv]['number_of_valid_hits']
    #        array_with_redundant_md5_data[i]['query_name'] = array_with_query_md5_data[row_in_df_dfout05_simapcsv]['query_name']
    #
    #
    #    #find the protein with the maximum number of hits
    #    index_maximum_number_of_valid_hits = array_with_redundant_md5_data[:]['number_of_valid_hits'].argmax()
    #    #label them as 'not kept' as default
    #    array_with_redundant_md5_data[:]['kept_after_redundancy_check'] = False
    #    #label the best as 'kept'
    #    array_with_redundant_md5_data[index_maximum_number_of_valid_hits]['kept_after_redundancy_check'] = True
    #    #now tranfer this data to the original array (array_with_query_md5_data)
    #    for i in range(len(array_with_redundant_md5_data)):
    #        row_in_df_dfout05_simapcsv = array_with_redundant_md5_data[i]['row_in_df_dfout05_simapcsv']
    #        array_with_query_md5_data[row_in_df_dfout05_simapcsv]['kept_after_redundancy_check'] = array_with_redundant_md5_data[i]['kept_after_redundancy_check']
    #logging.info(array_with_redundant_md5_data)
    #
    #find row in dataframe that contains the redundant md5
    #    index_row_in_df_containing_redundant_md5 = int(df_dfout05_simapcsv[df_dfout05_simapcsv['query_md5']==redundant_md5].index)
    #    #now update the df to show that this protein should be regarded as redundant
    #df_dfout05_simapcsv.loc[index_row_in_df_containing_redundant_md5,'kept_after_redundancy_check'] = array_with_redundant_md5_data[i]['kept_after_redundancy_check']
    #utils.save_structured_array_to_csv(array_with_query_md5_data,file_with_redundant_md5_checksums)
    #array6 = load_structured_array_from_csv(file_with_redundant_md5_checksums,dtype_array_with_query_md5_data)
    #array6 == array_with_query_md5_data
    #np.savetxt(filexxx3, array_with_query_md5_data, fmt='%s', delimiter=',', newline='\n', header='', footer='', comments='#')
    '''create list of protein names that are kept_after_redundancy_check
    this is necennary because of a corrupted csv file, where the order of the
    data rows doesn't match the array with the md5 redundancy data
    '''
    #    for datarow in array_with_query_md5_data:
    #        if datarow['kept_after_redundancy_check'] == True:
    #            list_of_proteins_kept_after_redundancy_check.append(datarow['query_name'])
    try:
        df_dfout07_simapnonred = \
            df_dfout05_simapcsv[
                df_dfout05_simapcsv.kept_after_redundancy_check == True]
    except KeyError:
        #for small lists of protein sequences, there are no redundancies and an error is received as all kept_after_redundancy_check values = True
        df_dfout07_simapnonred = df_dfout05_simapcsv
    list_of_proteins_kept_after_redundancy_check = list(
        df_dfout07_simapnonred.loc[:, 'A2_protein_name'])

    #sort the columns alphabetically (capitals first)
    df_dfout05_simapcsv = df_dfout05_simapcsv.sort_index(
        axis=1)
    #save the updated DataFrame of all proteins to file
    with open(dfout05_simapcsv, 'w') as f:
        df_dfout05_simapcsv.to_csv(f, sep=",", quoting=csv.QUOTE_NONNUMERIC)
    #save the updated DataFrame of nonredundant proteins to file
    with open(dfout07_simapnonred, 'w') as f:
        df_dfout07_simapnonred.to_csv(f, sep=",", quoting=csv.QUOTE_NONNUMERIC)

    '''
    Remove the data of 'redundant' sequences
    (sequences that fall into the SIMAP data of another query)
    The query with the MOST VALID HITS will be kept for further analysis.
    '''

    with open(csv_file_with_histogram_data_normalised, mode='r') as infile:
        reader = csv.reader(infile)
        with open(csv_file_with_histogram_data_normalised_redundant_removed, 'w') as outfile:
            pass
        with open(csv_file_with_histogram_data_normalised_redundant_removed, 'a') as outfile:
            writer = csv.writer(outfile, delimiter=',', quotechar='"', lineterminator='\n',
                                quoting=csv.QUOTE_NONNUMERIC, doublequote=True)
            rownumber = 0
            for row in reader:
                #logging.info(row)
                if rownumber == 0:
                    writer.writerow(row)
                else:
                    logging.info('redundancy check: %s' % row[0])
                    #logging.info("array_with_query_md5_data[rownumber-1]['query_name'] %s" % array_with_query_md5_data[rownumber-1]['query_name'])
                    #logging.info("array_with_query_md5_data[rownumber]['query_name'] %s\n" % array_with_query_md5_data[rownumber]['query_name'])
                    if row[0] in list_of_proteins_kept_after_redundancy_check:
                        #logging.info('match')
                        #logging.info(array5[rownumber-1]['kept_after_redundancy_check'])
                        #if convert_string_to_boolean_value(array_with_query_md5_data[rownumber-1]['kept_after_redundancy_check']):
                        #logging.info(array5[rownumber-1]['kept_after_redundancy_check'])
                        writer.writerow(row)
                rownumber += 1
    logging.info('calculate_TMD_conservation is finished.')

'''The gapped identity analysis gives the average TMD mem/nonmem conservation for homologues at different levels of similarity (eg. close homologues vs far homologues)
This is a funny system. The dictionary keys are the 5000 hits. This would be much more efficient and simple to re-write in pandas style.
'''
#A## variables are included only to help navigate the document in PyCharm
B02_OLD_calculate_TMD_conservation_by_gappedIdentity = settings["run_settings"][
    "calculate_TMD_conservation_by_gappedIdentity"]
if B02_OLD_calculate_TMD_conservation_by_gappedIdentity:
    #import matplotlib.pyplot as plt
    #csv_file_with_uniprot_data = r'E:\Stephis\Projects\Programming\Python\files\20131115_bacterial_TMD_conservation\List04-parser-test_uniprot_data.csv'

    csv_file_with_histogram_data_normalised_redundant_removed = r'E:\Stephis\Projects\Programming\Python\files\20131115_bacterial_TMD_conservation\List%02d_histogram_normalised_redundant_removed.csv' % list_number
    fieldnames_list_nonredundant_protein_names = ['A2_protein_name']
    nested_list_from_csv_nonred = utils.create_nested_dict_from_csv(
        csv_file_with_histogram_data_normalised_redundant_removed, fieldnames_list_nonredundant_protein_names)
    list_nonredundant_protein_names = []
    for i in range(1, len(nested_list_from_csv_nonred)):
        list_nonredundant_protein_names.append(nested_list_from_csv_nonred[i]['A2_protein_name'])
    #csv_file_TMD_cons_ratio_by_identity = r'E:\Stephis\Projects\Programming\Python\files\20131115_bacterial_TMD_conservation\List%s_TMD_cons_ratio_by_identitym.csv' % list_number

    #create a list of the files to ba analysed from the original csv file
    #list_files_with_features, list_files_with_homologues, list_of_protein_names, list_of_organism_domains = create_list_of_files_from_csv_with_uniprot_data(csv_file_with_uniprot_data, list_of_keys_to_keep)

    '''
    open csv file with homologues
    '''
    headers_to_retrieve_from_csv01 = ['A1_hit_number', 'match_TMD_added_to_FastA_alignment',
                                      'ratio_ident_mem_to_nonmem',
                                      'A3_md5', 'percentage_identity_of_TMD', 'FASTA_identity', 'FASTA_gapped_identity',
                                      'match_TMD_kept_for_statistical_analysis']

    #create new lists for each gapped identity cutoff
    list_of_gappedIdentity_cutoffs = settings["analysis_by_gappedIdentity"]["gappedIdentity_cutoffs"]

    for i in range(len(list_of_gappedIdentity_cutoffs) - 1):
        min_ = list_of_gappedIdentity_cutoffs[i]
        max_ = list_of_gappedIdentity_cutoffs[i + 1]
        csv_file_with_average_rpitrh_for_GI_cutoff = r'E:\Stephis\Projects\Programming\Python\files\20131115_bacterial_TMD_conservation\List%02d_TMD_cons_ratio_identity_%s-%s.csv' % (
            list_number, min_, max_)

        #create a new ewpty file
        headers = ['protein_name', 'organism_domain', 'number_of_valid_hits_within_cutoff',
                   'mean_ratio_ident_mem_to_nonmem_within_cutoff', 'stdev_ratio_ident_mem_to_nonmem_within_cutoff',
                   'min_', 'max_']
        utils.save_list_as_row_in_csv(headers, csv_file_with_average_rpitrh_for_GI_cutoff, 'w')

        #for i in range(len(list_of_protein_names)):
        for i in range(len(df_dfout05_simapcsv)):
            #take the organism domain (Eukaryota, Bacteria, Archaea) from the full organism classification list
            organism_domain = utils.convert_stringlist_to_list(
                df_dfout05_simapcsv.loc[i, 'uniprot_orgclass'])[0]
            protein_name = '%s_%s' % (df_dfout05_simapcsv.loc[i, 'A1_uniprot_accession'],
                                      df_dfout05_simapcsv.loc[
                                          i, 'uniprot_entry_name'])
            #open the csv file that contains a summary of all 'hits' for that protein
            csv_file_with_hit_details = 'E:\Databases\simap\%s\%s_homologues.csv' % (organism_domain, protein_name)

            #create a new list to hold the filtered TMD conservation values
            list_rpitrh_within_cutoffs = []

            if protein_name in list_nonredundant_protein_names:
                if os.path.isfile(csv_file_with_hit_details):
                    logging.info('conducting gappedIdentity cutoff analysis for %s...' % protein_name)

                    dict_from_homologue_csv_for_GI_analysis = utils.create_nested_dict_from_csv(
                        csv_file_with_hit_details, headers_to_retrieve_from_csv01)

                    for i in range(1, len(dict_from_homologue_csv_for_GI_analysis)):
                        #the old system used the match_TMD_added_to_FastA_alignment, but this needs to be updated like the match_TMD_kept_for_statistical_analysis in a pandas fashion
                        #match_TMD_added_to_FastA_alignment = utils.convert_string_to_boolean_value(dict_from_homologue_csv_for_GI_analysis[1]['match_TMD_added_to_FastA_alignment'])
                        match_TMD_kept_for_statistical_analysis = utils.convert_string_to_boolean_value(
                            dict_from_homologue_csv_for_GI_analysis[1]['match_TMD_kept_for_statistical_analysis'])
                        #print('%s, %s' % (dict_from_homologue_csv_for_GI_analysis[i], match_TMD_added_to_FastA_alignment)
                        if match_TMD_kept_for_statistical_analysis:
                            #I don't really understand why some ratio_ident_mem_to_nonmem = 0, even though they were match_TMD_kept_for_statistical_analysis
                            if dict_from_homologue_csv_for_GI_analysis[i]['ratio_ident_mem_to_nonmem'] != '':
                                FASTA_gapped_identity = float(
                                    dict_from_homologue_csv_for_GI_analysis[i]['FASTA_gapped_identity'])
                                #print(FASTA_gapped_identity)
                                if dict_from_homologue_csv_for_GI_analysis[i]['ratio_ident_mem_to_nonmem'] != '':
                                    ratio_ident_mem_to_nonmem = float(
                                        dict_from_homologue_csv_for_GI_analysis[i]['ratio_ident_mem_to_nonmem'])
                                else:
                                    ratio_ident_mem_to_nonmem = ''
                                if min_ < FASTA_gapped_identity <= max_:
                                    list_rpitrh_within_cutoffs.append(ratio_ident_mem_to_nonmem)
                number_of_valid_hits_within_cutoff = len(list_rpitrh_within_cutoffs)
                if number_of_valid_hits_within_cutoff >= settings["analysis_by_gappedIdentity"][
                    "minimum_number_hits_within_GI_range"]:
                    mean_ratio_ident_mem_to_nonmem_within_cutoff = np.mean(list_rpitrh_within_cutoffs)
                    stdev_ratio_ident_mem_to_nonmem_within_cutoff = np.std(list_rpitrh_within_cutoffs)
                else:
                    #there ARE not enough hits in this range, leave the cell empty
                    mean_ratio_ident_mem_to_nonmem_within_cutoff, stdev_ratio_ident_mem_to_nonmem_within_cutoff = '', ''

                #save the data for that protein as a line in the csv file for that cutoff
                list_to_save_in_csv = [protein_name, organism_domain, number_of_valid_hits_within_cutoff,
                                       mean_ratio_ident_mem_to_nonmem_within_cutoff,
                                       stdev_ratio_ident_mem_to_nonmem_within_cutoff, min_, max_]
                utils.save_list_as_row_in_csv(list_to_save_in_csv, csv_file_with_average_rpitrh_for_GI_cutoff, 'a')

#A## variables are included only to help navigate the document in PyCharm
A10_conduct_statistical_analysis_of_sim_ratios_saved_in_dfout05_simapcsv = settings["run_settings"][
    "conduct_statistical_analysis_of_sim_ratios_saved_in_dfout05_simapcsv"]
if A10_conduct_statistical_analysis_of_sim_ratios_saved_in_dfout05_simapcsv:
    #load the csv containing only nonredundant sequences as a pandas dataframe
    if settings["run_settings"]["conduct_stat_analysis_with_all_seqs_or_nonredundant_seqs"] == "all":
        df_dfout05_simapcsv_stat_analysis = pd.read_csv(dfout05_simapcsv, sep=",", index_col=0,
                                                        quoting=csv.QUOTE_NONNUMERIC)
    else:
        df_dfout05_simapcsv_stat_analysis = pd.read_csv(dfout07_simapnonred, sep=",", index_col=0,
                                                        quoting=csv.QUOTE_NONNUMERIC)

    #create a dictionary that contains the average conservation ratios for each protein, for each gap penalty and matrix
    dict_of_all_df_csv_file_av_cons_ratios_hits = {}
    for row in list(df_dfout05_simapcsv_stat_analysis.index):
        protein_name = df_dfout05_simapcsv_stat_analysis.loc[row, 'A2_protein_name']
        csv_file_av_cons_ratios_hits = df_dfout05_simapcsv_stat_analysis.loc[row, 'csv_file_av_cons_ratios_hits']
        if os.path.isfile(csv_file_av_cons_ratios_hits):
            df_csv_file_av_cons_ratios_hits = pd.read_csv(csv_file_av_cons_ratios_hits, sep=",", index_col=0,
                                                          quoting=csv.QUOTE_NONNUMERIC)
            #print(df_csv_file_av_cons_ratios_hits)
            #print(list(df_csv_file_av_cons_ratios_hits.columns))
            dict_of_all_df_csv_file_av_cons_ratios_hits[protein_name] = df_csv_file_av_cons_ratios_hits

    #convert the nested dictionary inta a 3D array (pandas panel) to hold the data of average sim_ratios for each protein
    wp_av_cons_ratios = pd.Panel.from_dict(dict_of_all_df_csv_file_av_cons_ratios_hits, orient='minor')

    #default = maor axis, each dataframe has the data for all matrixes for one protein (as originally in the dictionary)
    wp_av_cons_ratios = pd.Panel.from_dict(dict_of_all_df_csv_file_av_cons_ratios_hits)
    #create a pandas panel that contains all the average similarity ratios for each protein
    #orient='minor': each dataframe in the panel contains the data from one matrix, with all proteins as columns and all gap penalties as rows
    '''Example of axes:
    Items axis: ident_matrix to levin_matrix
    Major_axis axis: 0 to -10
    Minor_axis axis: P15509_CSF2R_HUMAN to Q6PJG9_LRFN4_HUMAN
    [each items axis is a separate dataframe]
    '''
    wp_av_cons_ratios = pd.Panel.from_dict(dict_of_all_df_csv_file_av_cons_ratios_hits, orient='minor')
    list_of_matrices = list(wp_av_cons_ratios.items)

    #calculate the average conservation ratio
    #create empty dataframe
    df_av_cons_ratio_all_proteins = pd.DataFrame(columns=list_of_matrices, index=list(wp_av_cons_ratios.major_axis))
    for matrix in list_of_matrices:
        #calculate, creating a new series
        mean_for_each_gap_penalty = wp_av_cons_ratios[matrix].mean(axis=1)
        #add the series to the final dataframe summary
        df_av_cons_ratio_all_proteins[matrix] = pd.Series(mean_for_each_gap_penalty)
    #save the averages to file
    with open(csv_av_cons_ratio_all_proteins, 'w') as f:
        df_av_cons_ratio_all_proteins.to_csv(f, sep=",", quoting=csv.QUOTE_NONNUMERIC)
        logging.info('%s   saved' % csv_av_cons_ratio_all_proteins)

    #calculate the standard deviation of the conservation ratios
    #create empty dataframe
    df_std_cons_ratio_all_proteins = pd.DataFrame(columns=list_of_matrices, index=list(wp_av_cons_ratios.major_axis))
    #calculate from panel
    for matrix in list_of_matrices:
        #print(matrix)
        #calculate, creating a new series
        std_for_each_gap_penalty = wp_av_cons_ratios[matrix].std(axis=1)
        #print('%s, %s' % (matrix,std_for_each_gap_penalty))
        #add the series to the final dataframe summary
        df_std_cons_ratio_all_proteins[matrix] = pd.Series(std_for_each_gap_penalty)
    #save the averages to file
    with open(csv_std_cons_ratio_all_proteins, 'w') as f:
        df_std_cons_ratio_all_proteins.to_csv(f, sep=",", quoting=csv.QUOTE_NONNUMERIC)
    logging.info('%s   saved' % csv_std_cons_ratio_all_proteins)

#A## variables are included only to help navigate the document in PyCharm
B03_OLD_fix_dfout05_simapcsv_by_adding_query_md5 = settings["run_settings"]["fix_dfout05_simapcsv_by_adding_query_md5"]
if B03_OLD_fix_dfout05_simapcsv_by_adding_query_md5:
    #convert the csv file to a pandas dataframe
    df_dfout05_simapcsv = pd.read_csv(dfout05_simapcsv,
                                      sep=",", na_values=['nan'],
                                      quoting=csv.QUOTE_NONNUMERIC)
    #create backup
    dfout05_simapcsv_backup = dfout05_simapcsv[
                              0:-4] + '_backup_after_query_md5_addition.csv'
    with open(dfout05_simapcsv_backup, 'w') as f:
        df_dfout05_simapcsv.to_csv(f, sep=",", quoting=csv.QUOTE_NONNUMERIC)

    #check that the query md5s are not there at first
    print(df_dfout05_simapcsv.query_md5)

    for i in range(len(df_dfout05_simapcsv)):
        organism_domain = df_dfout05_simapcsv.loc[i, 'organism_domain']
        protein_name = df_dfout05_simapcsv.loc[i, 'A2_protein_name']
        #try:
        logging.info('%s %s' % (i, protein_name))
        if os.path.isfile(df_dfout05_simapcsv.loc[i, 'SIMAP_feature_table_XML_file_path']):
            feature_table_exists = True
        else:
            feature_table_exists = False
            logging.info('feature table and/or homologue file is missing\n')
        if os.path.isfile(df_dfout05_simapcsv.loc[i, 'SIMAP_homologues_XML_file_path']):
            homologue_file_exists = True
        else:
            homologue_file_exists = False
            logging.info('feature table and/or homologue file is missing\n')
        if all([feature_table_exists, homologue_file_exists]):
            #parse the XML file with elementtree, define the 'root' of the XML file
            simap_homologue_tree = ET.parse(
                df_dfout05_simapcsv.loc[i, 'SIMAP_homologues_XML_file_path'])
            simap_homologue_root = simap_homologue_tree.getroot()
            query_sequence_node = simap_homologue_root[0][0][0][0][2][0][0]
            ''' xxxx CURRENTLY THE df_dfout05_simapcsv is filled with nan values, but that doesn't make sense as the script seems to work
            '''
            df_dfout05_simapcsv.loc[i, 'query_md5'] = query_sequence_node.attrib['md5']

            #check that the query md5s are now added
    print(df_dfout05_simapcsv.query_md5)

    #now save the improved dataframe to the csv
    with open(dfout05_simapcsv, 'w') as f:
        df_dfout05_simapcsv.to_csv(f, sep=",", quoting=csv.QUOTE_NONNUMERIC)

