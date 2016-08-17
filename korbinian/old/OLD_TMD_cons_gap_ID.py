import korbinian.mtutils as utils
import numpy as np
import pandas as pd
import csv
import os
# FOR SOME REASON, SCRIPT HAD BOTH df and df. All df references were then converted to df.

def OLD_calculate_TMD_conservation_by_gappedIdentity(pathdict, set_, logging, list_number):
     
    # import matplotlib.pyplot as plt
    # csv_file_with_uniprot_data = r'E:\Stephis\Projects\Programming\Python\files\20131115_bacterial_TMD_conservation\List04-parser-test_uniprot_data.csv'
    df = pd.read_csv(pathdict["dfout07_simapnonred"], sep=",", index_col=0, quoting=csv.QUOTE_NONNUMERIC)
    pathdict[ "csv_file_with_histogram_data_normalised_redundant_removed"] = r'E:\Stephis\Projects\Programming\Python\files\20131115_bacterial_TMD_conservation\List%02d_histogram_normalised_redundant_removed.csv' % list_number
    fieldnames_list_nonredundant_protein_names = ['protein_name']
    nested_list_from_csv_nonred = utils.create_nested_dict_from_csv(pathdict["csv_file_with_histogram_data_normalised_redundant_removed"], fieldnames_list_nonredundant_protein_names)
    list_nonredundant_protein_names = []
    for i in range(1, len(nested_list_from_csv_nonred)):
        list_nonredundant_protein_names.append(nested_list_from_csv_nonred[i]['protein_name'])
    # csv_file_TMD_cons_ratio_by_identity = r'E:\Stephis\Projects\Programming\Python\files\20131115_bacterial_TMD_conservation\List%s_TMD_cons_ratio_by_identitym.csv' % list_number
    
    # create a list of the files to ba analysed from the original csv file
    # list_files_with_features, list_files_with_homologues, list_of_protein_names, list_of_organism_domains = create_list_of_files_from_csv_with_uniprot_data(csv_file_with_uniprot_data, list_of_keys_to_keep)
    
    '''
    open csv file with homologues
    '''
    headers_to_retrieve_from_csv01 = ['hit_num', 'match_TMD_added_to_FastA_alignment',
                                      'ratio_ident_mem_to_nonmem',
                                      'A3_md5', 'percentage_identity_of_TMD', 'FASTA_identity', 'FASTA_gapped_identity',
                                      'match_TMD_kept_for_statistical_analysis']
    
    # create new lists for each gapped identity cutoff
    list_of_gappedIdentity_cutoffs = set_["gappedIdentity_cutoffs"]

    for i in range(len(list_of_gappedIdentity_cutoffs) - 1):
        min_ = list_of_gappedIdentity_cutoffs[i]
        max_ = list_of_gappedIdentity_cutoffs[i + 1]
        csv_file_with_average_rpitrh_for_GI_cutoff = r'D:\Schweris\Projects\Programming\Python\files\20131115_bacterial_TMD_conservation\List%02d_TMD_cons_ratio_identity_%s-%s.csv' % (
            list_number, min_, max_)
    
        # create a new ewpty file
        headers = ['protein_name', 'organism_domain', 'number_of_valid_hits_within_cutoff',
                   'mean_ratio_ident_mem_to_nonmem_within_cutoff', 'stdev_ratio_ident_mem_to_nonmem_within_cutoff',
                   'min_', 'max_']
        utils.save_list_as_row_in_csv(headers, csv_file_with_average_rpitrh_for_GI_cutoff, 'w')
    
        # for i in range(len(list_of_protein_names)):
        for i in range(len(df)):
            # take the organism domain (Eukaryota, Bacteria, Archaea) from the full organism classification list
            organism_domain = utils.convert_stringlist_to_list(
                df.loc[i, 'uniprot_orgclass'])[0]
            protein_name = '%s_%s' % (df.loc[i, 'uniprot_acc'],
                                      df.loc[
                                          i, 'uniprot_entry_name'])
            # open the csv file that contains a summary of all 'hits' for that protein
            csv_file_with_hit_details = 'E:\Databases\simap\%s\%s_homologues.csv' % (organism_domain, protein_name)
    
            # create a new list to hold the filtered TMD conservation values
            list_rpitrh_within_cutoffs = []
    
            if protein_name in list_nonredundant_protein_names:
                if os.path.isfile(csv_file_with_hit_details):
                    logging.info('conducting gappedIdentity cutoff analysis for %s...' % protein_name)
    
                    dict_from_homologue_csv_for_GI_analysis = utils.create_nested_dict_from_csv(
                        csv_file_with_hit_details, headers_to_retrieve_from_csv01)
    
                    for i in range(1, len(dict_from_homologue_csv_for_GI_analysis)):
                        # the old system used the match_TMD_added_to_FastA_alignment, but this needs to be updated like the match_TMD_kept_for_statistical_analysis in a pandas fashion
                        # match_TMD_added_to_FastA_alignment = utils.convert_string_to_boolean_value(dict_from_homologue_csv_for_GI_analysis[1]['match_TMD_added_to_FastA_alignment'])
                        match_TMD_kept_for_statistical_analysis = utils.convert_string_to_boolean_value(
                            dict_from_homologue_csv_for_GI_analysis[1]['match_TMD_kept_for_statistical_analysis'])
                        # print('%s, %s' % (dict_from_homologue_csv_for_GI_analysis[i], match_TMD_added_to_FastA_alignment)
                        if match_TMD_kept_for_statistical_analysis:
                            # I don't really understand why some ratio_ident_mem_to_nonmem = 0, even though they were match_TMD_kept_for_statistical_analysis
                            if dict_from_homologue_csv_for_GI_analysis[i]['ratio_ident_mem_to_nonmem'] != '':
                                FASTA_gapped_identity = float(
                                    dict_from_homologue_csv_for_GI_analysis[i]['FASTA_gapped_identity'])
                                # print(FASTA_gapped_identity)
                                if dict_from_homologue_csv_for_GI_analysis[i]['ratio_ident_mem_to_nonmem'] != '':
                                    ratio_ident_mem_to_nonmem = float(
                                        dict_from_homologue_csv_for_GI_analysis[i]['ratio_ident_mem_to_nonmem'])
                                else:
                                    ratio_ident_mem_to_nonmem = ''
                                if min_ < FASTA_gapped_identity <= max_:
                                    list_rpitrh_within_cutoffs.append(ratio_ident_mem_to_nonmem)
                number_of_valid_hits_within_cutoff = len(list_rpitrh_within_cutoffs)
                if number_of_valid_hits_within_cutoff >= set_["minimum_number_hits_within_GI_range"]:
                    mean_ratio_ident_mem_to_nonmem_within_cutoff = np.mean(list_rpitrh_within_cutoffs)
                    stdev_ratio_ident_mem_to_nonmem_within_cutoff = np.std(list_rpitrh_within_cutoffs)
                else:
                    # there ARE not enough hits in this range, leave the cell empty
                    mean_ratio_ident_mem_to_nonmem_within_cutoff, stdev_ratio_ident_mem_to_nonmem_within_cutoff = '', ''
    
                # save the data for that protein as a line in the csv file for that cutoff
                list_to_save_in_csv = [protein_name, organism_domain, number_of_valid_hits_within_cutoff,
                                       mean_ratio_ident_mem_to_nonmem_within_cutoff,
                                       stdev_ratio_ident_mem_to_nonmem_within_cutoff, min_, max_]
                utils.save_list_as_row_in_csv(list_to_save_in_csv, csv_file_with_average_rpitrh_for_GI_cutoff, 'a')