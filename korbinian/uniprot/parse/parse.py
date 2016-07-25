from Bio import SwissProt
import re
import numpy as np
import pandas as pd
import csv
import korbinian.mtutils as utils

def create_csv_from_uniprot_flatfile(uniprot_flatfile_of_selected_records, settingsdict, logging, pathdict):
    # create_csv_from_uniprot_flatfile(input_file=uniprot_flatfile_of_selected_records, output_file=csv_file_with_uniprot_data)
    ## open uniprot flatfile
    # def create_csv_from_uniprot_flatfile(input_file, output_file):
    # global uniprot_record, uni_dict, record
    # input_file=uniprot_flatfile_of_selected_records
    # output_file=csv_file_with_uniprot_data
    logging.info('~~~~~~~~~~~~  starting A03_create_csv_from_uniprot_flatfile   ~~~~~~~~~~~~')
    uniprot_dict_all_proteins = {}
    with open(uniprot_flatfile_of_selected_records, "r")as f:
        records = SwissProt.parse(f)
        count_of_uniprot_records_processed = 0
        count_of_uniprot_records_added_to_csv = 0
        for record in records:
            # uni_dict = utils.create_dict_of_data_from_uniprot_record(record)
            # create an empty output dictionary to hold the uniprot data for each record
            output_dict = {}
            # extract the subcellular location detail from the (poorly organized and unsorted) uniprot comments section
            comments_dict = {}
            # utils.create_dictionary_of_comments(record, comments_dict)
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
                                              re.IGNORECASE
                                              )
                    if comment_match:
                        regex_subcell_loc_dict[search_word] = True
            # the dictionary could also be nested within one column
            # output_dict['regex_subcell_loc_dict'] = regex_subcell_loc_dict
            # add all of the fields to the dictionary
            output_dict.update(regex_subcell_loc_dict)

            # print accession number
            logging.info(record.accessions[0])

            # add data to dictionary
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

            # create a list of all the feature types (signal, transmem, etc)
            list_of_feature_types_in_uniprot_record = []
            for sublist in record.features:
                list_of_feature_types_in_uniprot_record.append(sublist[0])
                # logging.info(sublist)

            # list of the features that we want in the final csv
            desired_features_in_uniprot = ['TRANSMEM', 'VARIANT', 'CONFLICT', 'VAR_SEQ', 'VARSPLIC', 'TOPO_DOM']
            desired_features_in_uniprot_dict = {}
            location_of_tmds_in_feature_list = []
            location_of_non_tmds_in_feature_list = []

            for feature in desired_features_in_uniprot:
                if feature in list_of_feature_types_in_uniprot_record:
                    # find the features in the feature list.
                    # For polytopic membrane proteins, there will be more than one tmd (labelled "TRANSMEM".
                    location_of_features_in_feature_list = [i for i, x in
                                                            enumerate(list_of_feature_types_in_uniprot_record) if
                                                            x == feature]
                    desired_features_in_uniprot_dict[feature] = location_of_features_in_feature_list
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

            # count the number of "TRANSMEM" TMDs listed in the feature-list
            output_dict['number_of_TMDs'] = len(location_of_tmds_in_feature_list)

            # information about location of first non-tmd (extracellular or perplasmic/cytoplasmic)
            if len(location_of_non_tmds_in_feature_list) > 0:
                output_dict['loc_start'] = record.features[location_of_non_tmds_in_feature_list[0]][3]
                output_dict['n_term_ec'] = "Extracellular" in output_dict["loc_start"]
            else:
                output_dict['loc_start'] = np.nan
                output_dict['n_term_ec'] = np.nan

            if output_dict['number_of_TMDs'] > 0:
                list_of_TMDs = []
                for TMD_location in location_of_tmds_in_feature_list:
                    # consequtively number the TMDs based on the "TRANSMEM" location in the feature list
                    TMD = 'TM%02d' % (location_of_tmds_in_feature_list.index(TMD_location) + 1)
                    list_of_TMDs.append(TMD)
                    # add the start and stop of each TMD, and the comments
                    output_dict['%s_start' % TMD] = record.features[TMD_location][1]
                    output_dict['%s_end' % TMD] = record.features[TMD_location][2]
                    output_dict['%s_description' % TMD] = record.features[TMD_location][3]

                # add the list of TMD names to the dictionary and dataframe
                output_dict['list_of_TMDs'] = list_of_TMDs

                # create a numpy array of any sequence variants are in the TMD region
                list_of_variant_types_in_uniprot = ['VARIANT', 'CONFLICT', 'VARSPLIC', 'VAR_SEQ']
                for TMD in list_of_TMDs:
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
                                start_of_variant_is_after_start_of_tmd = True if start_of_variant_in_seq > output_dict[
                                    '%s_start' % TMD] else False
                                end_of_variant_is_before_end_of_tmd = True if end_of_variant_in_seq < output_dict[
                                    '%s_end' % TMD] else False
                                variant_is_in_tmd = True if all([start_of_variant_is_after_start_of_tmd,
                                                                 end_of_variant_is_before_end_of_tmd]) else False
                                # if the variants are the tmd region, add to numpy array
                                if variant_is_in_tmd:
                                    # create array of the variant data
                                    variant_array = np.array(
                                        [variant_type, start_of_variant_in_seq, end_of_variant_in_seq,
                                         variant_description, variant_feature_identifier])
                                    if array_of_all_variants_in_tmd != ([]):
                                        # add array with the data for this variant to the array/list for all variants
                                        array_of_all_variants_in_tmd = np.row_stack(
                                            (array_of_all_variants_in_tmd, variant_array))
                                    else:
                                        # if the array is empty, replace the array for all variants with the array for the first variant
                                        array_of_all_variants_in_tmd = variant_array
                    # if there were variants added (array is not empty), convert to string and add them to the output dictionary
                    if array_of_all_variants_in_tmd.size:
                        output_dict['%s_seq_variants' % TMD] = str(array_of_all_variants_in_tmd)

            count_of_uniprot_records_processed += 1
            # nest each dictionary containing the data for each protein into a large dictionary that contains all data from all proteins
            uniprot_dict_all_proteins[output_dict['A1_uniprot_accession']] = output_dict

        # convert that nested dict into a pandas dataframe
        df = pd.DataFrame(uniprot_dict_all_proteins).sort_index()
        # count records in dataframe
        count_of_uniprot_records_added_to_csv = len(df.columns)
        # flip rows and columns (transverse)
        df = df.T.copy()

        ''' ~~ DETERMINE START AND STOP INDICES FOR TMD PLUS SURROUNDING SEQ ~~ '''
        aa_before_tmd = settingsdict["variables"]["analyse.simap_match_filters.aa_before_tmd"]
        aa_after_tmd = settingsdict["variables"]["analyse.simap_match_filters.aa_after_tmd"]
        # determine max number of TMD columns that need to be created
        max_num_TMDs = df['number_of_TMDs'].max()
        # currently the loop is run for each TMD, based on the sequence with the most TMDs
        for i in range(1, max_num_TMDs + 1):
            TMD = 'TM%02d' % i
            TMD_seq_name = '%s_seq' % TMD
            # instead of integers showing the start or end of the TMD, some people write strings into the
            # UniProt database, such as "<5" or "?"
            # to avoid the bugs that this introduces, it is necessary to convert all strings to np.nan (as floats),
            # using the convert objects function. The numbers can then be converted back from floats to integers.
            df['%s_start' % TMD] = pd.to_numeric(df['%s_start' % TMD]).dropna().astype('int64')
            df['%s_end' % TMD] = pd.to_numeric(df['%s_end' % TMD]).dropna().astype('int64')
            # determine the position of the start of the surrounding sequence
            df['start_surrounding_seq_in_query_%s' % TMD] = df['%s_start' % TMD] - aa_before_tmd
            # replace negative values with zero. (slicing method was replaced with lambda function to avoid CopyWithSetting warning)
            df['start_surrounding_seq_in_query_%s' % TMD] = df['start_surrounding_seq_in_query_%s' % TMD].apply(lambda x: x if x > 0 else 0)
            df['end_surrounding_seq_in_query_%s' % TMD] = df['%s_end' % TMD] + aa_after_tmd
            # create a boolean series, describing whether the end_surrounding_seq_in_query is longer than the protein seq
            series_indices_longer_than_prot_seq = df.apply(utils.find_indices_longer_than_prot_seq, args=(TMD,), axis=1)
            # obtain the indices of proteins in the series
            uniprot_acc_indices_longer_than_prot_seq = series_indices_longer_than_prot_seq[
                series_indices_longer_than_prot_seq].index
            # use indices to select the main dataframe, and convert these end_surrounding_seq_in_query values to the uniprot_seqlen value
            df.loc[uniprot_acc_indices_longer_than_prot_seq, 'end_surrounding_seq_in_query_%s' % TMD] = df.loc[
                uniprot_acc_indices_longer_than_prot_seq, 'uniprot_seqlen']
            # df = df['end_surrounding_seq_in_query_%s' % TMD].apply(lambda x: x['end_surrounding_seq_in_query_%s' % TMD] if x < df['uniprot_seqlen'] else df['uniprot_seqlen'])
            # df['end_surrounding_seq_in_query_%s' % TMD][df['end_surrounding_seq_in_query_%s' % TMD] > df['uniprot_seqlen']] = df['uniprot_seqlen']

        ''' ~~   SLICE TMDS FROM UNIPROT SEQ    ~~ '''
        # iterate through each TMD, slicing out the relevant sequence.
        # If there is no TMD, the cells will contain np.nan
        for i in range(1, max_num_TMDs + 1):
            TMD = 'TM%02d' % i
            # slice TMD
            df['%s_seq' % TMD] = df[df['%s_start' % TMD].notnull()].apply(utils.slice_uniprot_TMD_seq, args=(TMD,), axis=1)
            # slice TMD plus surrounding seq
            df['%s_with_surrounding_seq' % TMD] = df[df['%s_start' % TMD].notnull()].apply(utils.slice_uniprot_TMD_plus_surr_seq, args=(TMD,), axis=1)
        # extract the organism domain (e.g. Eukaryota)
        df['uniprot_orgclass'] = df['uniprot_orgclass'].astype(str)
        df['organism_domain'] = df.uniprot_orgclass.apply(lambda x: x.strip("'[]").split("', '")[0])
        # convert python datatypes to strings, as these currently give a TypeError when saving to excel
        df['uniprot_all_accessions'] = df['uniprot_all_accessions'].astype(str)
        df['uniprot_KW'] = df['uniprot_KW'].astype(str)
        df['uniprot_features'] = df['uniprot_features'].astype(str)
        df['list_of_TMDs'] = df['list_of_TMDs'].astype(str)
        # save to a csv
        df.to_csv(pathdict["dfout01_uniprotcsv"], sep=",", quoting=csv.QUOTE_NONNUMERIC)
        # save to Excel
        writer = pd.ExcelWriter(pathdict["dfout03_uniprotxlsx"])
        df.to_excel(writer, sheet_name='dfout03')
        writer.save()
        writer.close()

    logging.info('A03_create_csv_from_uniprot_flatfile was successful:'
                 '\n\t%i uniprot records processed\n\t%i uniprot records parsed to csv' % (
                     count_of_uniprot_records_processed, count_of_uniprot_records_added_to_csv))