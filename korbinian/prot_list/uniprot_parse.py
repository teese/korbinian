from Bio import SwissProt
import ast
import csv
import itertools
import korbinian
import korbinian.utils as utils
import logging
import numpy as np
import os
import pandas as pd
import re
import sys
# import debugging tools
from korbinian.utils import pr, pc, pn, aaa

def parse_flatfile_to_csv(selected_uniprot_records_flatfile, n_aa_before_tmd, n_aa_after_tmd, analyse_sp, logging, list_parsed_csv, slice=True):
    """ Parses a flatfile of UniProt records to csv.

    Parameters
    ----------
    selected_uniprot_records_flatfile : str
        Path to UniProt flatfile containing selected records for analysis.
    n_aa_before_tmd : int
        Number of amino acids before the TMD to be included when slicing the "TMD_plus_surr".
    n_aa_after_tmd : int
        Number of amino acids before the TMD to be included when slicing the "TMD_plus_surr".
    analyse_sp : bool
        Whether to analyse the signal peptides.
    logging : logging.Logger
        Logger for printing to console and logfile.
    list_parsed_csv : str
        Path to output csv file containing the list of proteins for analysis.

    Dataframes
    ----------
    dfu
        Dataframe for Uniprot
        index = acc for each protein
        columns = 'uniprot_acc', 'prot_descr', 'full_seq', etc

    Saved Files and Figures
    -----------------------
    list_summary_csv : csv
        CSV from dfu, with info for a protein on each row.

    """
    logging.info('~~~~~~~~~~~~                 starting parse_flatfile_to_csv              ~~~~~~~~~~~~')
    if not os.path.isfile(selected_uniprot_records_flatfile):
        return "parse_flatfile_to_csv could not be run. Uniprot flatfile not found. ({})".format(selected_uniprot_records_flatfile)
    uniprot_dict_all_proteins = {}
    with open(selected_uniprot_records_flatfile, "r") as f:
        records = SwissProt.parse(f)
        count_of_uniprot_records_processed = 0
        for m, record in enumerate(records):
            # create an empty output dictionary to hold the uniprot data for each record
            output_dict = {}
            # extract the subcellular location detail from the (poorly organized and unsorted) uniprot comments section
            comments_dict = {}

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
                # there are no comments in this uniprot file!
                logging.info('no comments in Uniprot file')
                output_dict['comments_subcellular_location_uniprot'] = ''

            # use regex to search for text describing subcellular locations
            # [ -]? accepts either space, hyphen, or no dividing character
            regex_word_dict = {'multipass': ['multi[ -]?(pass|span)', 'poly[ -]?topic'],
                               'singlepass': ['single[ -]?(pass|span)', 'bi[ -]?topic'],
                               'membrane': ['membran', 'lipid[ -](anchor|bound)'],
                               'typeI': ['type[ -](one|1|I)[ -]membran'],
                               'typeII': ['type[ -](two|2|II)[ -]membran']}
            # comments_subcellular_location_uniprot = 'Membrane; Bitopictype I membrane protein.'
            regex_subcell_loc_dict = {}
            for search_word in regex_word_dict:
                regex_subcell_loc_dict[search_word] = False
                regex_search_list = regex_word_dict[search_word]
                for regex_search_string in regex_search_list:
                    # search for the regex string, ignoring any mismatches in upper or lower case
                    comment_match = re.search(regex_search_string,
                                              output_dict['comments_subcellular_location_uniprot'],
                                              re.IGNORECASE)
                    if comment_match:
                        regex_subcell_loc_dict[search_word] = True
            # add all of the fields to the dictionary
            output_dict.update(regex_subcell_loc_dict)

            # print accession number
            sys.stdout.write("{}, ".format(record.accessions[0]))
            if count_of_uniprot_records_processed % 20 == 0:
                sys.stdout.write("\n".format(record.accessions[0]))
            sys.stdout.flush()

            # add data to dictionary
            output_dict['uniprot_acc'] = record.accessions[0]
            output_dict['organism'] = record.organism
            output_dict['uniprot_entry_name'] = record.entry_name
            output_dict['gene_name'] = record.gene_name
            output_dict['prot_descr'] = record.description
            output_dict['full_seq'] = record.sequence
            output_dict['uniprot_orgclass'] = record.organism_classification
            output_dict['uniprot_all_accessions'] = record.accessions
            output_dict['uniprot_KW'] = record.keywords
            output_dict['uniprot_features'] = record.features
            output_dict['seqlen'] = record.sequence_length

            # create a list of all the feature types (signal, transmem, etc)
            list_of_feature_types_in_uniprot_record = []
            for sublist in record.features:
                list_of_feature_types_in_uniprot_record.append(sublist[0])

            # list of the features that we want in the final csv
            desired_features_in_uniprot = ['TRANSMEM', 'VARIANT', 'CONFLICT', 'VAR_SEQ', 'VARSPLIC', 'TOPO_DOM']
            if analyse_sp == True:
                desired_features_in_uniprot.append('SIGNAL')
            desired_features_in_uniprot_dict = {}
            location_of_sp_in_feature_list = []
            location_of_tmds_in_feature_list = []
            location_of_non_tmds_in_feature_list = []

            # add bool if uniprot thinks that protein contains signal peptides
            if 'SIGNAL' in list_of_feature_types_in_uniprot_record:
                output_dict['uniprot_SiPe'] = True
            else:
                output_dict['uniprot_SiPe'] = False

            for feature in desired_features_in_uniprot:
                if feature in list_of_feature_types_in_uniprot_record:
                    # find the features in the feature list.
                    # For polytopic membrane proteins, there will be more than one tmd (labelled "TRANSMEM".
                    location_of_features_in_feature_list = [i for i, x in
                                                            enumerate(list_of_feature_types_in_uniprot_record) if
                                                            x == feature]
                    desired_features_in_uniprot_dict[feature] = location_of_features_in_feature_list
                    if feature == 'SIGNAL':
                        location_of_sp_in_feature_list = location_of_features_in_feature_list
                        # location_of_sp_in_feature_list.sort()
                    if feature == 'TRANSMEM':
                        location_of_tmds_in_feature_list = location_of_features_in_feature_list
                        # sort list to be sure that the "transmem" notation is definitely ordered correctly,
                        # as this order determines the TMD name
                        location_of_tmds_in_feature_list.sort()
                    if feature == 'TOPO_DOM':
                        location_of_non_tmds_in_feature_list = location_of_features_in_feature_list
                        # sort list to be sure that the "transmem" notation is definitely ordered correctly,
                        # as this order determines the TMD name
                        location_of_non_tmds_in_feature_list.sort()

            # count the number of SP
            output_dict['number_of_SP'] = len(location_of_sp_in_feature_list)
            # count the number of "TRANSMEM" TMDs listed in the feature-list
            output_dict['number_of_TMDs'] = len(location_of_tmds_in_feature_list)

            # information about location of first non-tmd (extracellular or periplasmic/cytoplasmic)
            if len(location_of_non_tmds_in_feature_list) > 0:
                output_dict['loc_start'] = record.features[location_of_non_tmds_in_feature_list[0]][3]
                output_dict['n_term_ec'] = "Extracellular" in output_dict["loc_start"]
            else:
                output_dict['loc_start'] = np.nan
                output_dict['n_term_ec'] = np.nan

            # number of TMDs excluding signal peptides (which might be added later)
            number_of_TMDs_excl_SP = len(location_of_tmds_in_feature_list)
            output_dict['number_of_TMDs_excl_SP'] = number_of_TMDs_excl_SP

            list_of_TMDs = ["TM{:02d}".format(n) for n in range(1, number_of_TMDs_excl_SP + 1)]
            output_dict['list_of_TMDs'] = list_of_TMDs
            output_dict['list_of_TMDs_excl_SP'] = list_of_TMDs

            if number_of_TMDs_excl_SP > 0:

                #list_of_TMDs = []
                TM_indices = []

                for n, TMD_location in enumerate(location_of_tmds_in_feature_list):
                    # consecutively number the TMDs based on the "TRANSMEM" location in the feature list
                    #TMD = 'TM{:02d}'.format(n+1)
                    #list_of_TMDs.append(TMD)

                    TM_start = record.features[TMD_location][1]
                    TM_end = record.features[TMD_location][2]

                    # remove any strings or floats
                    for TM_index in [TM_start, TM_end]:
                        TM_index = TM_index if isinstance(TM_index, int) else np.nan

                    # add to nested tuple
                    TM_indices.append((TM_start, TM_end))

                    # DEPRECATED
                    # # add the start and stop of each TMD, and the comments
                    # output_dict['%s_start'%TMD] = record.features[TMD_location][1]
                    # output_dict['%s_end'%TMD] = record.features[TMD_location][2]
                    # output_dict['%s_description'%TMD] = record.features[TMD_location][3]
                    # if isinstance(output_dict['%s_start'%TMD], str) or isinstance(output_dict['%s_end'%TMD], str):
                    #     logging.info("{} strings found in feature indices: {},{}".format(output_dict['uniprot_acc'], output_dict['%s_start'%TMD], output_dict['%s_end'%TMD]))
                    #     output_dict['%s_start' % TMD], output_dict['%s_end'%TMD] = np.nan, np.nan

                # information about SP location
                if output_dict['number_of_SP'] != 0:
                    for SP_location in location_of_sp_in_feature_list:
                        SP = 'SP01'
                        list_of_TMDs.append(SP)
                    for SP_location in location_of_sp_in_feature_list:
                        output_dict['SP01_start'] = record.features[SP_location][1]
                        output_dict['SP01_end'] = record.features[SP_location][2]
                        output_dict['SP01_description'] = record.features[SP_location][3]
                    if isinstance(output_dict['SP01_start'], str) or isinstance(output_dict['SP01_end'], str):
                        logging.info("{} strings found in feature indices: {},{}".format(output_dict['uniprot_acc'], output_dict['SP01_start'], output_dict['SP01_end']))
                        output_dict['SP01_start'], output_dict['SP01_end'] = np.nan, np.nan

                # add the list of TMD names to the dictionary and dataframe
                #output_dict['list_of_TMDs'] = list_of_TMDs
                output_dict["TM_indices"] = TM_indices

                # create a numpy array of any sequence variants are in the TMD (and SP) region
                list_of_variant_types_in_uniprot = ['VARIANT', 'CONFLICT', 'VARSPLIC', 'VAR_SEQ']
                for n, TMD in enumerate(list_of_TMDs):
                    TM_start = TM_indices[n][0]
                    TM_start = TM_indices[n][1]

                    # array_of_all_variants_in_tmd = np.zeros(4)
                    array_of_all_variants_in_tmd = np.array([])
                    for variant_type in list_of_variant_types_in_uniprot:
                        if variant_type in desired_features_in_uniprot_dict.keys():
                            # if that variant is in the uniprot data for that protein, create a list of the indices showing where that variant is found
                            list_of_variant_locations = list(desired_features_in_uniprot_dict[variant_type])
                            # get the specific start, end and details of that variant
                            for v in range(len(list_of_variant_locations)):
                                # get start
                                start_of_variant_in_seq = record.features[list_of_variant_locations[v]][1]
                                # get end
                                end_of_variant_in_seq = record.features[list_of_variant_locations[v]][2]
                                # get description
                                variant_description = record.features[list_of_variant_locations[v]][3]
                                variant_feature_identifier = record.features[list_of_variant_locations[v]][4]
                                # check if the variant is in the tmd
                                start_of_variant_is_after_start_of_tmd = True if start_of_variant_in_seq > TM_start else False
                                end_of_variant_is_before_end_of_tmd = True if end_of_variant_in_seq < TM_start else False
                                variant_is_in_tmd = True if all([start_of_variant_is_after_start_of_tmd, end_of_variant_is_before_end_of_tmd]) else False
                                # if the variants are the tmd region, add to numpy array
                                if variant_is_in_tmd:
                                    # create array of the variant data
                                    variant_array = np.array(
                                        [variant_type, start_of_variant_in_seq, end_of_variant_in_seq,
                                         variant_description, variant_feature_identifier])
                                    if array_of_all_variants_in_tmd.size != 0:
                                        # add array with the data for this variant to the array/list for all variants
                                        array_of_all_variants_in_tmd = np.row_stack(
                                            (array_of_all_variants_in_tmd, variant_array))
                                    else:
                                        # if the array is empty, replace the array for all variants with the array for the first variant
                                        array_of_all_variants_in_tmd = variant_array
                    # if there were variants added (array is not empty), convert to string and add them to the output dictionary
                    if array_of_all_variants_in_tmd.size:
                        output_dict['%s_seq_variants'%TMD] = str(array_of_all_variants_in_tmd)

            count_of_uniprot_records_processed += 1
            # nest each dictionary containing the data for each protein into a large dictionary that contains all data from all proteins
            uniprot_dict_all_proteins[output_dict['uniprot_acc']] = output_dict

        # convert that nested dict into a pandas dataframe, transverse
        dfu = pd.DataFrame(uniprot_dict_all_proteins).sort_index().T
        # count records in dataframe
        count_of_initial_uniprot_records = dfu.shape[0]
        # make a unique list of all TMD combinations in list([TM01], [TM01, TM03], etc)
        unique_TMD_combinations_orig = list(dfu.list_of_TMDs.astype(str).unique())
        # convert to python list
        unique_TMD_combinations_lists = [ast.literal_eval(s) for s in unique_TMD_combinations_orig if "nan" not in s]
        # grab all unique values into a large list(e.g. TM01, TM02, TM03 until last TM of protein with most TMs)
        unique_TMD_combinations_single_list = [i for i in itertools.chain.from_iterable(unique_TMD_combinations_lists)]
        # sort
        list_all_TMDs_in_dataset = sorted(list(set(unique_TMD_combinations_single_list)))

        # extract the organism domain (e.g. Eukaryota)
        dfu['uniprot_orgclass'] = dfu['uniprot_orgclass'].astype(str)
        dfu['organism_domain'] = dfu.uniprot_orgclass.apply(lambda x: x.strip("'[]").split("', '")[0])
        # convert python datatypes to strings, as these currently give a TypeError when saving to excel
        dfu['uniprot_all_accessions'] = dfu['uniprot_all_accessions'].astype(str)
        dfu['uniprot_KW'] = dfu['uniprot_KW'].astype(str)
        dfu['uniprot_features'] = dfu['uniprot_features'].astype(str)
        dfu['list_of_TMDs'] = dfu['list_of_TMDs'].astype(str)
        dfu['topol_source'] = "UniProt"         #Hotfix: "UniProt" instead of "uniprot" or
                                                #else all proteins will be filtered out later

        # save to a csv
        utils.make_sure_path_exists(list_parsed_csv, isfile=True)
        # count records in dataframe
        count_of_uniprot_records_added_to_csv = dfu.shape[0]

        dfu.to_csv(list_parsed_csv, sep=",", quoting=csv.QUOTE_NONNUMERIC)

    return '\n%i valid UniProt records parsed to csv (from %i initial)\n~~~~~~~~~~~~                 finished parse_flatfile_to_csv              ~~~~~~~~~~~~' % (count_of_uniprot_records_added_to_csv, count_of_initial_uniprot_records)


def slice_TMD_seqs(selected_uniprot_records_flatfile, n_aa_before_tmd, n_aa_after_tmd, analyse_sp, logging, list_parsed_csv, slice=True):
    list_all_TMDs_in_dataset = "tobeadded"

    for TMD in list_all_TMDs_in_dataset:
        dfu = korbinian.prot_list.prot_list.get_indices_TMD_plus_surr_for_summary_file(dfu, TMD, n_aa_before_tmd, n_aa_after_tmd)

    ''' ~~   SLICE TMDS FROM SEQ    ~~ '''
    for TMD in list_all_TMDs_in_dataset:
        df_proteins_with_start_and_end = dfu.loc[dfu['%s_start' % TMD].notnull()].loc[dfu['%s_end' % TMD].notnull()]
        if df_proteins_with_start_and_end.empty:
            # For this TMD, none of the proteins have an appropriate start and stop. skip to next TMD
            logging.info("For {}, none of the proteins had valid indices! This may cause errors later.".format(TMD))
            continue
        # slice TMD
        dfu['%s_seq' % TMD] = df_proteins_with_start_and_end.apply(utils.slice_uniprot_TMD_seq, args=(TMD,), axis=1)
        # slice TMD plus surrounding seq
        dfu['%s_seq_plus_surr' % TMD] = df_proteins_with_start_and_end.apply(utils.slice_uniprot_TMD_plus_surr_seq, args=(TMD,), axis=1)

        ''' ~~   SLICE nonTMD sequence FROM UNIPROT SEQ    ~~ '''
        if slice == True:
            sys.stdout.write ('\nslicing nonTMD sequences:')
            valid_acc_list = dfu.loc[dfu['list_of_TMDs'].notnull()].loc[dfu['list_of_TMDs'] != "nan"].index
            for n, acc in enumerate(valid_acc_list):
                if n % 10 == 0 and n is not 0:
                    sys.stdout.write('.')
                    sys.stdout.flush()
                    if n % 200 == 0 and n is not 0:
                        sys.stdout.write('\n')
                        sys.stdout.flush()
                list_of_TMDs = ast.literal_eval(dfu.loc[acc, 'list_of_TMDs'])
                if 'SP01' in list_of_TMDs:
                    list_of_TMDs.remove('SP01')
                # sequence from N-term. to first TMD
                nonTMD_first = dfu.loc[acc, 'full_seq'][0: (dfu.loc[acc, 'TM01_start']-1).astype('int64')]
                sequence = nonTMD_first
                # only for multipass proteins, generate sequences between TMDs
                if len(list_of_TMDs) > 1:
                    for TM_Nr in range(len(list_of_TMDs) - 1):
                        # the TMD is the equivalent item in the list
                        TMD = list_of_TMDs[TM_Nr]
                        # the next TMD, which contains the end index, is the next item in the list
                        next_TMD = list_of_TMDs[TM_Nr + 1]
                        between_TM_and_TMplus1 = dfu.loc[acc, 'full_seq'][dfu.loc[acc, '%s_end' %TMD].astype('int64'): dfu.loc[acc, '%s_start' %next_TMD].astype('int64')-1]
                        sequence += between_TM_and_TMplus1
                last_TMD = list_of_TMDs[-1]
                # sequence from last TMD to C-term.
                nonTMD_last = dfu.loc[acc, 'full_seq'][dfu.loc[acc, '%s_end' %last_TMD].astype('int64'):dfu.loc[acc, 'seqlen']]
                sequence += nonTMD_last
                dfu.loc[acc, 'nonTMD_seq'] = sequence
                dfu.loc[acc, 'len_nonTMD'] = len(sequence)


def create_dictionary_of_comments(uniprot_record_handle):
    """ Create a dictionary of the comments from a UniProt record.

    Parameters
    ----------
    uniprot_record_handle : Bio.SwissProt.Record
        UniProt record handle, via BioPython.

    Returns
    -------
    output_dictionary : dict
        Output dictionary containing the comments.
    """
    output_dictionary = {} #create a new empty dictionary
    try:
        for comment in uniprot_record_handle.comments:
            # splits comments based on first ":" symbol, creates a list called split_comment
            split_comment = comment.strip().split(': ', 1)
             # several comments have the same name. need to check if it is already in the dictionary
            if split_comment[0] in output_dictionary:
                # list the different comments, one after another
                output_dictionary[split_comment[0]] += ", %s" % split_comment[1]
            else:
                output_dictionary[split_comment[0]] = split_comment[1]
    except AttributeError:
        #there are no comments in this uniprot file!
        logging.info('no comments in Uniprot file')
        output_dictionary = {}
    return output_dictionary

def create_dict_of_data_from_uniprot_record_DEPRECATED(record):
    """ DEPRECATED.

    use parse_flatfile_to_csv instead.
    """
    pass
    # DEPRECATED FUNCTION
    # convert the comments in uniprot to a dictionary.
    #the comments in uniprot are unordered and non-hierarchical and need to be processed with string manipulation
    # comments_dict = create_dictionary_of_comments(record)
    #
    # #create an empty output dictionary to holnd the uniprot data for each record
    # output_dict = {}
    #
    # #print first accession number
    # logging.info(record.accessions[0])
    #
    # #by default this is zero
    # output_dict['uniprot_TMD_start'] = 0
    # output_dict['uniprot_TMD_end'] = 0
    # output_dict['uniprot_TMD_description'] = 0
    # output_dict['uniprot_record_included_in_csv'] = True
    #
    # #add data to dictionary
    # output_dict['accession_uniprot'] = record.accessions[0]
    # output_dict['full_list_of_accessions_uniprot'] = record.accessions
    # output_dict['record_entry_name_uniprot'] = record.entry_name
    # output_dict['record_gene_name_uniprot'] = record.gene_name
    # output_dict['record_description_uniprot'] = record.description
    # output_dict['sequence_uniprot'] = record.sequence
    # output_dict['organism_classification'] = record.organism_classification
    # output_dict['organism'] = record.organism
    # output_dict['record_keywords_uniprot'] = record.keywords
    # output_dict['record_features_all_uniprot'] = record.features
    # output_dict['record_sequence_length_uniprot'] = record.sequence_length
    # output_dict['comments_subcellular_location_uniprot'] = comments_dict['SUBCELLULAR LOCATION']
    #
    # #create a list of all the feature types (signal, transmem, etc)
    # list_of_feature_types_in_uniprot_record = []
    # for sublist in record.features:
    #     list_of_feature_types_in_uniprot_record.append(sublist[0])
    #
    # #list of the features that we want in the final csv
    # desired_features_in_uniprot = ['TRANSMEM', 'VARIANT', 'CONFLICT', 'VAR_SEQ','VARSPLIC'] # SIGNAL ?
    # desired_features_in_uniprot_dict = {}
    #
    # location_of_tmds_in_feature_list = []
    #
    # for feature in desired_features_in_uniprot:
    #     if feature in list_of_feature_types_in_uniprot_record:
    #         #find the features in the feature list. For polytopic membrane protoins, there will be more than one tmd.
    #         location_of_features_in_feature_list = [i for i,x in enumerate(list_of_feature_types_in_uniprot_record) if x == feature]
    #         desired_features_in_uniprot_dict[feature] = location_of_features_in_feature_list
    #         if feature == 'TRANSMEM':
    #             location_of_tmds_in_feature_list = location_of_features_in_feature_list
    #
    # #determine if the TMD is actually annotated in the uniprot record
    # if 'TRANSMEM' not in desired_features_in_uniprot_dict.keys():
    #     TRANSMEM_missing_from_uniprot_features = True
    # else:
    #     TRANSMEM_missing_from_uniprot_features = False
    # output_dict['TRANSMEM_missing_from_uniprot_features'] = TRANSMEM_missing_from_uniprot_features
    # output_dict['more_than_one_tmd_in_uniprot_annotation'] = False
    #
    # if not TRANSMEM_missing_from_uniprot_features:
    #
    #     # number of TMDs excluding signal peptides (which might be added later)
    #     number_of_TMDs_excl_SP = len(location_of_tmds_in_feature_list)
    #     output_dict['number_of_TMDs_excl_SP'] = number_of_TMDs_excl_SP
    #     list_of_TMDs_excl_SP = ["TM{}".format(n) for n in range(1, number_of_TMDs_excl_SP + 1)]
    #     #add the TMD to the dictionary for single-pass membrane proteins
    #     if number_of_TMDs_excl_SP == 1:
    #         location = location_of_tmds_in_feature_list[0]
    #         #output_dict['number_of_TMDs_excl_SP'] = 1
    #         output_dict['uniprot_TMD_start'] = record.features[location][1]
    #         output_dict['uniprot_TMD_end'] = record.features[location][2]
    #         output_dict['uniprot_TMD_description'] = record.features[location][3]
    #
    #     #add the TMD to the dictionary for multi-pass membrane proteins
    #     if number_of_TMDs_excl_SP > 1:
    #         logging.info('more than one "TRANSMEM" feature in uniprot for %s' % output_dict['accession_uniprot'])
    #         #output_dict['number_of_TMDs_excl_SP'] = len(location_of_tmds_in_feature_list)
    #         output_dict['more_than_one_tmd_in_uniprot_annotation'] = True
    #         output_dict['uniprot_record_included_in_csv'] = False
    #
    #     if number_of_TMDs_excl_SP == 1:
    #         #create a numpy array of any sequence variants are in the TMD region
    #         list_of_variant_types_in_uniprot = ['VARIANT', 'CONFLICT','VARSPLIC', 'VAR_SEQ']
    #         #array_of_all_variants_in_tmd = np.zeros(4)
    #         array_of_all_variants_in_tmd = np.array([])
    #         for variant_type in list_of_variant_types_in_uniprot:
    #             if variant_type in desired_features_in_uniprot_dict.keys():
    #                 list_of_variant_locations = list(desired_features_in_uniprot_dict[variant_type])
    #                 for i in range(len(list_of_variant_locations)):
    #                     start_of_variant_in_seq = record.features[list_of_variant_locations[i]][1]
    #                     end_of_variant_in_seq = record.features[list_of_variant_locations[i]][2]
    #                     variant_description = record.features[list_of_variant_locations[i]][3]
    #                     variant_feature_identifier = record.features[list_of_variant_locations[i]][4]
    #                    #check if the variant is in the tmd
    #                     start_of_variant_is_after_start_of_tmd = True if start_of_variant_in_seq > output_dict['uniprot_TMD_start'] else False
    #                     end_of_variant_is_before_end_of_tmd = True if end_of_variant_in_seq < output_dict['uniprot_TMD_end'] else False
    #                     variant_is_in_tmd = True if all([start_of_variant_is_after_start_of_tmd and end_of_variant_is_before_end_of_tmd]) else False
    #
    #                     #add to numpy array that contains all the variants in the tmd region
    #                     if variant_is_in_tmd:
    #                         variant_array = np.array([variant_type, start_of_variant_in_seq, end_of_variant_in_seq, variant_description,variant_feature_identifier])
    #                         if array_of_all_variants_in_tmd.size == 0:
    #                             array_of_all_variants_in_tmd = np.row_stack((array_of_all_variants_in_tmd, variant_array))
    #                         else:
    #                             array_of_all_variants_in_tmd = variant_array
    #         #if there were variants added, add them to the output dictionary
    #         output_dict['array_of_all_variants_in_tmd'] = array_of_all_variants_in_tmd
    #
    # #if the tmd region is annotated, get the sequence
    # if not TRANSMEM_missing_from_uniprot_features:
    #     if number_of_TMDs_excl_SP == 1:
    #         output_dict['uniprot_TMD_sequence'] = record.sequence[output_dict['uniprot_TMD_start'] - 1:output_dict['uniprot_TMD_end']]
    #     else:
    #         output_dict['uniprot_record_included_in_csv'] = False
    #         output_dict['uniprot_TMD_start'] = 0
    #         output_dict['uniprot_TMD_end'] = 0
    #         output_dict['uniprot_TMD_description'] = 0
    # else:
    #     output_dict['uniprot_record_included_in_csv'] = False
    #     logging.info('%s: "TRANSMEM" not found in uniprot features, therefore not incuded in csv file for further analysis' % output_dict['accession_uniprot'])
    #     output_dict['uniprot_TMD_start'] = 0
    #     output_dict['uniprot_TMD_end'] = 0
    #     output_dict['uniprot_TMD_description'] = 0
    #
    # return output_dict
