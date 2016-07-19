import korbinian.mtutils as utils
import numpy as np
import pandas as pd
import csv
import os
import xml.etree.ElementTree as ET



def OLD_conduct_statistical_analysis_of_sim_ratios_saved_in_dfout05_simapcsv(pathdict, settingsdict, logging):
    #load the csv containing only nonredundant sequences as a pandas dataframe
    if settingsdict["variables"]["conduct_statistical_analysis_of_sim_ratios_saved_in_dfout05_simapcsv.conduct_stat_analysis_with_all_seqs_or_nonredundant_seqs"] == "all":
        df_dfout05_simapcsv_stat_analysis = pd.read_csv(pathdict["dfout05_simapcsv"], sep=",", index_col=0,
                                                        quoting=csv.QUOTE_NONNUMERIC)
    else:
        df_dfout05_simapcsv_stat_analysis = pd.read_csv(pathdict["dfout07_simapnonred"], sep=",", index_col=0,
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
    with open(pathdict["csv_av_cons_ratio_all_proteins"], 'w') as f:
        df_av_cons_ratio_all_proteins.to_csv(f, sep=",", quoting=csv.QUOTE_NONNUMERIC)
        logging.info('%s   saved' % pathdict["csv_av_cons_ratio_all_proteins"])

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
    with open(pathdict["csv_std_cons_ratio_all_proteins"], 'w') as f:
        df_std_cons_ratio_all_proteins.to_csv(f, sep=",", quoting=csv.QUOTE_NONNUMERIC)
    logging.info('%s   saved' % pathdict["csv_std_cons_ratio_all_proteins"])


def OLD_fix_dfout05_simapcsv_by_adding_query_md5(pathdict, logging):
    #convert the csv file to a pandas dataframe
    df_dfout05_simapcsv = pd.read_csv(pathdict["dfout05_simapcsv"],
                                      sep=",", na_values=['nan'],
                                      quoting=csv.QUOTE_NONNUMERIC)
    #create backup
    pathdict["dfout05_simapcsv_backup"] = pathdict["dfout05_simapcsv"][0:-4] + '_backup_after_query_md5_addition.csv'
    with open(pathdict["dfout05_simapcsv_backup"], 'w') as f:
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
    with open(pathdict["dfout05_simapcsv"], 'w') as f:
        df_dfout05_simapcsv.to_csv(f, sep=",", quoting=csv.QUOTE_NONNUMERIC)
