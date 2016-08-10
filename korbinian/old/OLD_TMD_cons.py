import ast
import csv
import korbinian.mtutils as utils
import numpy as np
import pandas as pd
import sys
import os
# FOR SOME REASON, SCRIPT HAD BOTH df_old and df_old. All df_old references were then converted to df_old.

def OLD_calculate_TMD_conservation(pathdict, set_, logging):
    # convert file to dataframe (eg List70_simap.csv)
    df_old = pd.read_csv(pathdict["df_oldout05_simapcsv"], sep=",", index_col=0,quoting=csv.QUOTE_NONNUMERIC)  # index_col = 0, only wher reopening a saved pandas file
    # # add desired columns (only necessary if for loop method is kept)
    # original_columns = list(df_old.columns)
    # new_columns = ['md5_checksums_in_homologue_list_that_are_also_in_query_list']
    # # columns_added_after_SIMAP_analysis as above
    # new_unique_column_list = set(original_columns + new_columns)
    # # add extra columns
    # df_old = df_old.reindex(index=df_old.index, columns=new_unique_column_list)
    # # sort columns
    # df_old = df_old.sort_index(axis=1)
    
    # specify some data types:
    df_old['kept_after_redundancy_check'] = \
        df_old['kept_after_redundancy_check'].astype(np.bool)
    df_old['protein_kept_for_statistical_analysis'] = \
        df_old['protein_kept_for_statistical_analysis'].astype(np.bool)
    df_old['query_md5'] = df_old['query_md5'].astype('<U32')
    df_old['md5_checksums_in_homologue_list_that_are_also_in_query_list'] = \
        df_old[
            'md5_checksums_in_homologue_list_that_are_also_in_query_list'].astype('<U3000')
    
    # create a list of the files to ba analysed from the original csv file
    # list_files_with_features, list_files_with_homologues, list_of_protein_names, list_of_organism_domains = create_list_of_files_from_csv_with_uniprot_data(csv_file_with_uniprot_data, list_of_keys_to_keep)
    
    '''OLD BIN CREATION'''
    # create the bins used for the histogram data. I add 3 empty field 0,0,0 to account for the three columns
    # The first bin is 0. The bins range UPWARDS, so 0 = =-0.4,
    #    list_of_bins_for_histogram = []
    #    for i in range(20,250,10):
    #        i2 = i/100.0
    #        list_of_bins_for_histogram.append(i2)
    #    list_of_bins_for_histogram.append(30)
    # logging.info(list_of_bins_for_histogram)
    
    #    list_of_bins_for_histogram2 = np.linspace(0.2,2.5,24)
    #    logging.info(list_of_bins_for_histogram2)
    list_of_bins_for_histogram = np.linspace(set_["1p_smallest_bin"],
                                             set_["1p_largest_bin"],
                                             set_["1p_number_of_bins"])
    # add 30 as the last bin, to make sure 100% of the data is added to the histogram, including major outliers
    list_of_bins_for_histogram = np.append(list_of_bins_for_histogram,
                                           set_["1p_final_highest_bin"])
    
    list_of_keys_for_hit_identity_analysis = ['match_TMD_added_to_FastA_alignment', 'ratio_ident_mem_to_nonmem',
                                              'match_TMD_kept_for_statistical_analysis']
    
    # minimum_number_hits_for_data_analysis = 20
    number_query_seq_processed = 0
    number_query_seq_added_to_histogram = 0
    list_of_proteins_kept_for_statistical_analysis = []
    list_of_proteins_whose_XML_file_doesnt_exist = []
    list_of_proteins_with_not_enough_valid_hits = []
    list_of_csv_files_with_homologues_that_are_old = []
    dict_of_arrays_mem_nonmem_ratios = {}
    
    for i in df_old.index:
        # take the organism domain (Eukaryota, Bacteria, Archaea) from the full organism classification list
        first_two_letters_of_uniprot_acc = df_old.loc[i, 'first_two_letters_of_uniprot_acc']
        protein_name = df_old.loc[i, 'protein_name']
        organism_domain = df_old.loc[i, 'organism_domain']
        SIMAP_feature_table_XML_file_path = os.path.join(df_old.loc[i, 'simap_filename_base'], first_two_letters_of_uniprot_acc,
                                                         '%s_feature_table.xml' % protein_name)
        SIMAP_homologues_XML_file_path = os.path.join(df_old.loc[i, 'simap_filename_base'], first_two_letters_of_uniprot_acc,
                                                      '%s_homologues.xml' % protein_name)
        # for i in range(len(list_of_protein_names)):
        # print dots to show that program is running
        sys.stdout.write(" .")
        sys.stdout.flush()
        csv_file_with_hit_details = os.path.join(df_old.loc[i, 'simap_filename_base'], first_two_letters_of_uniprot_acc,
                                                 '%s_homologues.csv' % protein_name)
    
        # check if the feature table and homologue files actually exist
        if os.path.isfile(csv_file_with_hit_details):
            csv_file_with_hit_details_exists = True
            # get the size of the file (actual size, not windows "size on disk")
            statinfo = os.stat(csv_file_with_hit_details)
            csv_file_with_hit_details_filesize = statinfo.st_size
            if csv_file_with_hit_details_filesize > 2400:
                csv_file_with_hit_details_contains_hits = True
            else:
                csv_file_with_hit_details_contains_hits = False
        else:
            csv_file_with_hit_details_exists = False
            list_of_proteins_whose_XML_file_doesnt_exist.append(protein_name)
            # logging.info('%s_homologues.csv file is missing' % protein_name)
        if csv_file_with_hit_details_exists and csv_file_with_hit_details_contains_hits:
            # check if any of the files are old: in the old files, the headers had spaces in between each word, rather than an underscore. Later there might be other woys to search for older csv files with homologue data.
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
                    protein_name[:6], df_old.loc[i, 'SIMAP_homologues_XML_file_path'],
                    df_old.loc[i, 'SIMAP_csv_analysed_path']))
                # The uniprot TMD start and stop is necessary. Take this from the csv with uniprot data. Thi file and list_of_keys_to_keep should be global objects.
                "since parse_single_homologue_XML_to_csv_and_fasta is no longer  function, "
                "the old file cannot simply be updated here"
            # nested_dict_with_uniprot_seq = create_nested_dict_from_csv(csv_file_with_uniprot_data, list_of_keys_to_keep)
            #                #since this is a nested dictionary, I have to scroll through it until I obtain the correct protein:
            #                for i in range(1,len(nested_dict_with_uniprot_seq)+1):
            #                    protein_in_csv_cell = df_old.loc[i,'uniprot_acc']
            #                    print(protein_in_csv_cell)
            #                    print(protein_name[:6])
            #
            #                    if protein_in_csv_cell == protein_name[:6]:
            #                        row_in_csv_with_protein = i
            # now update the csv with homologues, just for this sequence
            #                parse_single_homologue_XML_to_csv_and_fasta(df_old.loc[i,'SIMAP_homologues_XML_file_path'],df_old.loc[i,'SIMAP_csv_analysed_path'], nested_dict_with_uniprot_seq[row_in_csv_with_protein])
    
            # if not csv_file_with_homologues_is_old: #xxx hopefully they auto-update now!!!
            '''Use Pandas for the analysis of the hits, rather than the nested dictionary approach used previously
            '''
            df_old_SIMAP_hits = pd.read_csv(csv_file_with_hit_details, sep=",", na_values=['nan'],
                                        quoting=csv.QUOTE_NONNUMERIC)  # index_col = 0, only wher reopening a saved pandas file
            # convert all np.nan (empty values) to empty strings
            # for col in df_old_SIMAP_hits.columns[pd.isnull(df_old_SIMAP_hits).all()]:
            #    df_old_SIMAP_hits[col] = df_old_SIMAP_hits[col].astype(object).fillna("UNKNOWN")
            # df_old_SIMAP_hits.fillna(value='empty')
            #            original_columns_df_old_SIMAP_hits = list(df_old_SIMAP_hits.columns)
            #            added_columns_df_old_SIMAP_hits = ['match_kept', 'free_column']
            #            new_column_list = list(original_columns_df_old_SIMAP_hits + added_columns_df_old_SIMAP_hits)
            #            df_old_SIMAP_hits = df_old_SIMAP_hits.reindex(index = df_old_SIMAP_hits.index, columns = new_column_list)
            #            df_old_SIMAP_hits.loc[:,'match_kept'] = False
            #            df_old_SIMAP_hits['match_kept'].astype = np.bool_
            list_of_columns_to_fill_empty_values = ['disallowed_words_in_description']
            for col in list_of_columns_to_fill_empty_values:
                df_old_SIMAP_hits[col] = df_old_SIMAP_hits[col].astype(object).fillna("none")
    
            # filter each homologue from SIMAP
            for j in range(len(df_old_SIMAP_hits)):
                df_old_SIMAP_hits.loc[j, 'match_TMD_kept_for_statistical_analysis'] = True if all([
                    df_old_SIMAP_hits.loc[j, 'FASTA_expectation'] <= set_['e_value_filter'],
                    df_old_SIMAP_hits.loc[j, 'disallowed_words_in_description'] == 'none',
                    set_['database'] == 'all',
                    df_old_SIMAP_hits.loc[j, 'number_of_gaps_in_match_TMD'] <= set_["cr_max_n_gaps_in_match_TMD"],
                    df_old_SIMAP_hits.loc[j, 'number_of_gaps_in_query_TMD'] <= set_["cr_max_n_gaps_in_query_TMD"],
                    df_old_SIMAP_hits.loc[j, 'SW_identity'] >= set_["cr_minimum_identity_of_full_protein"],
                    df_old_SIMAP_hits.loc[j, 'percentage_identity_of_TMD'] >= set_["cr_min_identity_of_TMD_final_filter"],
                    # df_old_SIMAP_hits.loc[j, 'number_of_X_in_TMD'] <= set_['number_of_X_allowed_in_TMD'],
                    df_old_SIMAP_hits.loc[j, 'SW_match_TMD_seq'] != 'TMD_not_in_SW_alignment'
                ]) else False
                # len(df_old_SIMAP_hits.loc[j,'query_aln_seq_excl_TMD']) >= set_['min_len_query_aln_seq_excl_TMD']
                # df_old_SIMAP_hits.loc[j,'number_of_X_in_match_seq'] <= set_['number_of_X_allowed_in_seq'],
    
            # simply select the true data in a smaller dataframe
            # use the value_counts method to count the number of hits kept for statistical analysis
            df_old_SIMAP_homologues_kept_for_statistical_analysis_stat_only = df_old_SIMAP_hits[
                df_old_SIMAP_hits.match_TMD_kept_for_statistical_analysis == True]
            # both stat and fasta filters. NEEDS TO BE UPDATED LATER TO ONLY REFLECT WORKING STAT FILTERS
            df_old_SIMAP_homologues_kept_for_statistical_analysis = \
                df_old_SIMAP_homologues_kept_for_statistical_analysis_stat_only[
                    df_old_SIMAP_homologues_kept_for_statistical_analysis_stat_only.match_TMD_added_to_FastA_alignment == True]
    
            with open(df_old.loc[i, 'csv_SIMAP_homologues_kept_for_statistical_analysis'], 'w') as f:
                df_old_SIMAP_homologues_kept_for_statistical_analysis.to_csv(f, sep=",", quoting=csv.QUOTE_NONNUMERIC)
    
                # NOT working! gives error whel len(df_old) = 0
                # number_of_hits_kept_for_statistical_analysis = df_old_SIMAP_hits['match_TMD_kept_for_statistical_analysis'].value_counts()[True]
                # now check that the DataFrame really only holds the selected data (match_TMD_kept_for_statistical_analysis).
                # if number_of_hits_kept_for_statistical_analysis != len(df_old_SIMAP_homologues_kept_for_statistical_analysis):
                # logging.warning('number_of_hits_kept_for_statistical_analysis != len(df_old_SIMAP_homologues_kept_for_statistical_analysis) for %s)' % protein_name)
    
            # check if the csv_file_av_cons_ratios_hits file already exists
            csv_mem_nonmem_ratios_is_old = False
            if os.path.isfile(df_old.loc[i, 'csv_file_av_cons_ratios_hits']):
                # load into a pandas dataframe
                df_old_mem_nonmem_ratios = pd.read_csv(df_old.loc[i, 'csv_file_av_cons_ratios_hits'],
                                                   sep=",", index_col=0, quoting=csv.QUOTE_NONNUMERIC)
                gap_penalties_in_file = list(df_old_mem_nonmem_ratios.index)
                matrices_in_file = list(df_old_mem_nonmem_ratios.columns)
                if gap_penalties_in_file == list(
                        range(set_["gap_open_penalty_min"],
                              set_["gap_open_penalty_max"],
                              set_["gap_open_penalty_increment"])) and \
                                matrices_in_file == set_["matrices"]:
                    csv_mem_nonmem_ratios_is_old = False
                else:
                    csv_mem_nonmem_ratios_is_old = True
                    logging.info('%s is old, repeating calculation' % df_old.loc[
                        i, 'csv_file_av_cons_ratios_hits'])
    
            if set_["overwrite_csv_file_av_cons_ratios_hits"] or csv_mem_nonmem_ratios_is_old:
                # if not os.path.isfile(df_old.loc[i,'csv_file_av_cons_ratios_hits']):
                # create a new dataframe to hold the array of mem/nonmem ratios for the varying matrices and  gap penalties
                df_old_mem_nonmem_ratios = pd.DataFrame(0.0,
                                                    index=range(
                                                        set_["gap_open_penalty_min"],
                                                        set_["gap_open_penalty_max"],
                                                        set_["gap_open_penalty_increment"]),
                                                    columns=set_['["mp_matrices'])
    
                # load the amino acid substitution matrices from the settings file
                list_of_aa_sub_matrices = set_['["mp_matrices']
                dict_of_aa_matrices = {key: ast.literal_eval(key) for key in list_of_aa_sub_matrices}
    
                # for each gap penalty
                for k in range(set_["gap_open_penalty_min"],
                               set_["gap_open_penalty_max"],
                               set_["gap_open_penalty_increment"]):
                    gap_open_penalty = k
                    # for each aa sub matrix, represented as a key in the dictionary
                    for key in dict_of_aa_matrices:
                        matrix_name = repr(key).replace("'", "")
                        column_name = 'sim_ratio_%s_gapo%i' % (matrix_name.replace("'", "")[0:-7], k)
                        mean_sim_ratio = df_old_SIMAP_homologues_kept_for_statistical_analysis[column_name].mean()
                        # if I wanted to, I could save a histogram for each protein!!
                        # hist_sim_ratio = df_old_SIMAP_homologues_kept_for_statistical_analysis[column_name].hist()
    
                        df_old_mem_nonmem_ratios.loc[gap_open_penalty, matrix_name] = mean_sim_ratio
    
                # save as a separate file?
                with open(df_old.loc[i, 'csv_file_av_cons_ratios_hits'], 'w') as f:
                    df_old_mem_nonmem_ratios.to_csv(f, sep=",", quoting=csv.QUOTE_NONNUMERIC)
                    print('%s created' % df_old.loc[i, 'csv_file_av_cons_ratios_hits'])
            else:
                # csv_file_av_cons_ratios_hits already created, load from file
                logging.info('%s already exists, moving to next' % df_old.loc[
                    i, 'csv_file_av_cons_ratios_hits'])
                df_old_mem_nonmem_ratios = pd.read_csv(df_old.loc[i, 'csv_file_av_cons_ratios_hits'],
                                                   sep=",", quoting=csv.QUOTE_NONNUMERIC)
    
            df_old.loc[i, 'number_of_valid_hits'] = len(
                df_old_SIMAP_homologues_kept_for_statistical_analysis)
            df_old.loc[
                i, 'df_old_mem_nonmem_ratios'] = df_old_mem_nonmem_ratios.to_string()
            # df_old.loc[i,'df_old_mem_nonmem_ratios'] = df_old_mem_nonmem_ratios
            # add to the uniprot df_old/csv the mem/nonmem and stdev calculated from the SW alignment and counting the identical residues
            df_old.loc[i, 'mean_ratio_ident_mem_to_nonmem'] = \
                df_old_SIMAP_homologues_kept_for_statistical_analysis['ratio_ident_mem_to_nonmem'].mean()
            df_old.loc[i, 'stdev_ratio_ident_mem_to_nonmem'] = \
                df_old_SIMAP_homologues_kept_for_statistical_analysis['ratio_ident_mem_to_nonmem'].std()
    
            # nested_dict_for_hit_identity_analysis = create_nested_dict_from_csv(csv_file_with_hit_details, 'all')#xxx had problems with KeyError: 'match_TMD_added_to_FastA_alignment',
            nested_dict_for_hit_identity_analysis = utils.create_nested_dict_from_csv(csv_file_with_hit_details,
                                                                                      list_of_keys_for_hit_identity_analysis)
    
            # create a list of the ratio_ident_mem_to_nonmem for each hit, to be used in the histogram for each protein
            list_ratio_ident_mem_to_nonmem = []
            number_of_valid_hits = 0
            for j in range(1, len(nested_dict_for_hit_identity_analysis)):
                # if the TMD from this hit was good enough to add to the alignment
                match_TMD_added_to_FastA_alignment = utils.convert_string_to_boolean_value(
                    nested_dict_for_hit_identity_analysis[j]['match_TMD_added_to_FastA_alignment'])
                match_TMD_kept_for_statistical_analysis = utils.convert_string_to_boolean_value(
                    nested_dict_for_hit_identity_analysis[j]['match_TMD_kept_for_statistical_analysis'])
    
                # logging.info(match_added)
                # start with the second hit, as the first is the query sequence
                if j >= 2:
                    if df_old_SIMAP_hits.loc[j, 'match_TMD_kept_for_statistical_analysis']:
                        if nested_dict_for_hit_identity_analysis[j]['ratio_ident_mem_to_nonmem'] != '':
                            number_of_valid_hits += 1
                            ratio_ident_mem_to_nonmem = float(
                                nested_dict_for_hit_identity_analysis[j]['ratio_ident_mem_to_nonmem'])
                            list_ratio_ident_mem_to_nonmem.append(ratio_ident_mem_to_nonmem)
            # logging.info(list_ratio_ident_mem_to_nonmem)
    
            # determine if protein is kept for statistical analysis. Currently the only filter is the number of valid hits.
            if df_old.loc[i, 'number_of_valid_hits'] >= \
                    set_["1p_min_n_hits_for_data_analysis"]:
                df_old.loc[i, 'protein_kept_for_statistical_analysis'] = True
                # if there is enough valid hits, add the protein name to the list, so the data can be used later
                list_of_proteins_kept_for_statistical_analysis.append(protein_name)
            else:
                df_old.loc[i, 'protein_kept_for_statistical_analysis'] = False
                list_of_proteins_with_not_enough_valid_hits.append(
                    df_old.loc[i, 'uniprot_acc'])
    
            logging.info('%s\tnumber_of_valid_hits: %i, protein_kept_for_statistical_analysis: %s' % (
                protein_name, df_old.loc[i, 'number_of_valid_hits'],
                df_old.loc[i, 'protein_kept_for_statistical_analysis']))
            '''Previously, a histogram was created for each protein
            '''
            if df_old.loc[i, 'protein_kept_for_statistical_analysis']:
                # numpy will automatically create a histogram
                hist_ratio_ident_mem_to_nonmem = np.histogram(list_ratio_ident_mem_to_nonmem,
                                                              bins=list_of_bins_for_histogram)
                # If True, the result is the value of the probability density function at the bin, normalized such that the integral over the range is 1.
                hist_ratio_ident_mem_to_nonmem_normalised = np.histogram(list_ratio_ident_mem_to_nonmem, density=True,
                                                                         bins=list_of_bins_for_histogram)
                # logging.info(protein_name, hist)
    
                # prepare csv header
                # this bin list should be the same as the list_of_bins_for_histogram
                bin_list = np.array(hist_ratio_ident_mem_to_nonmem)[1].tolist()
                list_of_bins_as_strings = []
                for i8 in range(len(bin_list) - 1):
                    bin_string = '%.1f-%.1f' % (bin_list[i8], bin_list[i8 + 1])
                    list_of_bins_as_strings.append(bin_string)
    
                csv_added_header = ['protein_name', 'organism_domain', 'number_of_valid_hits',
                                    'mean_ratio_ident_mem_to_nonmem', 'stdev_ratio_ident_mem_to_nonmem']
                header_with_bins_for_csv = csv_added_header + list_of_bins_as_strings
    
                # for the first file, create a new file and save the csv header
                if number_query_seq_added_to_histogram == 0:
                    utils.save_list_as_row_in_csv(header_with_bins_for_csv, pathdict["csv_file_with_histogram_data"], 'w')
                    utils.save_list_as_row_in_csv(header_with_bins_for_csv,
                                                  pathdict["csv_file_with_histogram_data_normalised"], 'w')
    
                # rpitrh = ratio_ident_mem_to_nonmem
                array_rpitrh = np.array(hist_ratio_ident_mem_to_nonmem)
                hist_ratio_ident_mem_to_nonmem_list = np.array(hist_ratio_ident_mem_to_nonmem)[0].tolist()
                hist_ratio_ident_mem_to_nonmem_normalised_list = np.array(hist_ratio_ident_mem_to_nonmem_normalised)[
                    0].tolist()
    
                # calculate mean and stdev
                mean_ratio_ident_mem_to_nonmem = np.mean(list_ratio_ident_mem_to_nonmem)
                stdev_ratio_ident_mem_to_nonmem = np.std(list_ratio_ident_mem_to_nonmem)
    
                # prepare data as list for csv
                list_protein_name = [protein_name, organism_domain, number_of_valid_hits,
                                     mean_ratio_ident_mem_to_nonmem, stdev_ratio_ident_mem_to_nonmem]
                hist_ratio_ident_mem_to_nonmem_csv_row = list_protein_name + hist_ratio_ident_mem_to_nonmem_list
                hist_ratio_ident_mem_to_nonmem_normalised_csv_row = list_protein_name + hist_ratio_ident_mem_to_nonmem_normalised_list
    
                # save data to csv
                utils.save_list_as_row_in_csv(hist_ratio_ident_mem_to_nonmem_csv_row,
                                              pathdict["csv_file_with_histogram_data"], 'a')
                utils.save_list_as_row_in_csv(hist_ratio_ident_mem_to_nonmem_normalised_csv_row,
                                              pathdict["csv_file_with_histogram_data_normalised"], 'a')
    
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
        set_["1p_min_n_hits_for_data_analysis"],
        list_of_proteins_with_not_enough_valid_hits))
    logging.info('list_of_csv_files_with_homologues_that_are_old: %s' % list_of_csv_files_with_homologues_that_are_old)
    
    '''
     REMOVE REDUNDANT/DUPLICATE RECORDS
    '''
    
    '''new pandas method: take the list of md5 directly from the pandas dataframe
    '''
    
    # make a list of the md5 checksums, including nan values
    list_of_query_md5_checksums = []
    for ii2 in df_old.index:
        list_of_query_md5_checksums.append(df_old.loc[ii2, 'query_md5'])
    
    # print(list_of_query_md5_checksums)
    # convert to pandas data series
    series_of_query_md5_checksums = pd.Series(list_of_query_md5_checksums)
    # remove nan values
    series_of_query_md5_checksums_without_nan_values = series_of_query_md5_checksums.dropna()
    # convert back to a list
    list_of_query_md5_checksums = list(series_of_query_md5_checksums_without_nan_values)
    if 'nan' in list_of_query_md5_checksums:
        list_of_query_md5_checksums.remove('nan')
    # remove any empty values (note that the filter function tries to link to elementtree, and not the basic python function)
    list_of_query_md5_checksums = [x for x in list_of_query_md5_checksums if x]
    print('list_of_query_md5_checksums: %s' % list_of_query_md5_checksums)
    '''
    For each csv file containing all the homologues, create a list of the md5s (md5_list_from_simap_hits).
    Then use the 'set' function to identify which md5 checksums are shared by both the list of query proteins,
    and the list of homologue hits for that particular protein.
    '''
    # create an empty array to hold the data. For the data type thecolumn containing the list of md5s will be large, so use unicode strings with at least 600 characters (>U300)
    #    number_of_columns = 6
    #    dtype_array_with_query_md5_data = [('number', '<i4'), ('query_name', '<U32'),
    #                                       ('query_md5', '<U33'),
    #                                       ('md5_checksums_in_homologue_list_that_are_also_in_query_list', '<U3000'),
    #                                       ('number_of_valid_hits', '<U32'),
    #                                       ('kept_after_redundancy_check', np.bool)]
    #    array_with_query_md5_data = np.ones(len(df_old),
    #                                        dtype = dtype_array_with_query_md5_data)
    # logging.info(array_with_query_md5_data)
    '''Create a list of all the md5 checksums of all SIMAP hits
    '''
    list_of_proteins_with_damaged_file_0_valid_hits = []
    list_of_proteins_with_no_csv_file_with_hit_details = []
    
    for i in df_old.index:
        csv_file_with_hit_details = df_old.loc[
            i, 'SIMAP_csv_analysed_path']
        # if the file exists
        if os.path.isfile(csv_file_with_hit_details):
            # create a dataframe from the csv file, replace empty values with the tewt 'nan'
            df_old_SIMAP_hits = pd.read_csv(csv_file_with_hit_details, sep=",", na_values=['nan'],
                                        quoting=csv.QUOTE_NONNUMERIC)
            # if there are actually hits, and the SIMAP download was sucessful
            if len(df_old_SIMAP_hits) != 0:
                # create a list of all the md5 values from the simap homolague file
                md5_list_from_simap_hits = list(df_old_SIMAP_hits['A3_md5'])
    
                query_md5 = df_old.loc[i, 'query_md5']
                # use the .intersection function to identify other query sequences that are found as hits in the homologues for this query
                # create a set, exclude the first sequence, as that is the query
                set_md5_list_from_simap_hits = set(md5_list_from_simap_hits[1:])
                # find elemnts common to the set and the list
                md5_checksums_in_homologue_list_that_are_also_in_query_list = list(
                    set_md5_list_from_simap_hits.intersection(set(list_of_query_md5_checksums)))
                # add data to the array, so it can be saved as a csv.
                #            df_old.loc[i,'number'] = i
                #            df_old.loc[i,'query_name'] = protein_name
                #            df_old.loc[i,'query_md5'] = query_md5
                # for the list of redundant md5s, replace commas with | to improve compatibility with csv files. Where there are no redundant md5 accessions, leave blank
                if md5_checksums_in_homologue_list_that_are_also_in_query_list != []:
                    df_old.loc[
                        i, 'md5_checksums_in_homologue_list_that_are_also_in_query_list'] = ' | '.join(
                        md5_checksums_in_homologue_list_that_are_also_in_query_list)
                else:
                    df_old.loc[
                        i, 'md5_checksums_in_homologue_list_that_are_also_in_query_list'] = ''
                    # logging.info(protein_name)
                    # logging.info(md5_list_from_simap_hits)
                    # logging.info('md5_checksums_found_in_both_lists: %s' % md5_checksums_found_in_both_lists)
            else:
                # if there are no valid hits, replace the np.nan with 0, add the uniprot accession to a list so that the SIMAP retrieval can be repeated
                df_old.loc[
                    i, 'md5_checksums_in_homologue_list_that_are_also_in_query_list'] = ''
                list_of_proteins_with_damaged_file_0_valid_hits.append(
                    df_old.loc[i, 'uniprot_acc'])
                logging.info('file damaged or no hits (len(df_old_SIMAP_hits) = 0 : %s)' % csv_file_with_hit_details)
        else:
            # the file is not there, replace the np.nan with an empty field
            df_old.loc[i, 'query_md5'] = ''
            list_of_proteins_with_no_csv_file_with_hit_details.append(
                df_old.loc[i, 'uniprot_acc'])
            logging.info('file does not exist: %s' % csv_file_with_hit_details)
    logging.info(
        'list_of_proteins_with_damaged_file_0_valid_hits (repeating SIMAP download is necessary) : %s' % list_of_proteins_with_damaged_file_0_valid_hits)
    logging.info('list_of_proteins_with_no_csv_file_with_hit_details (if len>%i, try repeating SIMAP download) : %s' % (
        set_["max_query_sequence_length"],
        list_of_proteins_with_no_csv_file_with_hit_details))
    
    # specify data type, repeat of earlier specification for df_old
    df_old['query_md5'] = df_old[
        'query_md5'].astype('<U32')
    df_old['md5_checksums_in_homologue_list_that_are_also_in_query_list'] = \
        df_old[
            'md5_checksums_in_homologue_list_that_are_also_in_query_list'].astype('<U3000')
    
    # convert the md5s from a string (separated by |) to a list.
    for i in df_old.index:
        if df_old.loc[
            i, 'md5_checksums_in_homologue_list_that_are_also_in_query_list'] != '':
            md5_checksums_in_homologue_list_that_are_also_in_query_list = \
                df_old.loc[
                    i, 'md5_checksums_in_homologue_list_that_are_also_in_query_list'].split(' | ')
    '''
    Save the updated df_old with SIMAP data later??
    '''
    # innocent until proven guilty! all seqs (at least those with SIMAP data) are nonredundant by default
    df_old['kept_after_redundancy_check'].astype(np.bool)
    # df_old.loc[:,'kept_after_redundancy_check'] = True
    # label all with some valid hits nonredundant, all those with valid hits are innocent until proven guilty
    df_old['kept_after_redundancy_check'] = \
        df_old['SIMAP_total_hits'] > 0
    
    # %%
    # do the redundancy check by loading all of the md5s from SIMAP homologues for each protein, and checking against the md5s of all of the query sequences in that list [note.. I had many troubles clearing these of np.nan!]
    for i in df_old.index:
        # only examine the files with redundant sequences to see which has the most valid hits
        if df_old.loc[
            i, 'md5_checksums_in_homologue_list_that_are_also_in_query_list'] != '':
            # split the md5s into a list
            md5_checksums_in_homologue_list_that_are_also_in_query_list = \
                df_old.loc[
                    i, 'md5_checksums_in_homologue_list_that_are_also_in_query_list'].split(' | ')
            # add the original query md5 to the list
            md5_checksums_in_homologue_list_that_are_also_in_query_list.append(
                df_old.loc[i, 'query_md5'])
            # logging.info(md5_checksums_in_homologue_list_that_are_also_in_query_list)
    
            # create a new subarray to hold the data from the redundant md5s
            # subarray = np.zeros(len(md5_checksums_in_homologue_list_that_are_also_in_query_list))
    
            # array_of_redundant_md5s = np.array(md5_checksums_in_homologue_list_that_are_also_in_query_list)
            # array_with_redundant_md5_data = np.ones(len(array_of_redundant_md5s)*4, dtype = [('redundant_md5', '<U33'), ('row_in_df_old', '<i4'), ('number_of_valid_hits', '<i4'), ('kept_after_redundancy_check', 'b')]).reshape(len(array_of_redundant_md5s),4)
            # create an array, just for those redundant md5s
            array_with_redundant_md5_data = np.ones(len(md5_checksums_in_homologue_list_that_are_also_in_query_list),
                                                    dtype=[('query_name', '<U32'), ('redundant_md5', '<U33'), (
                                                        'row_in_df_old', '<i4'),
                                                           ('number_of_valid_hits', '<i4'),
                                                           ('kept_after_redundancy_check', '<U5')])
            # for readability, convert to a pandas dataframe
            df_old_redundant_md5 = pd.DataFrame(array_with_redundant_md5_data)
            '''New Pandas dataframe: transfer directly to df_old with uniprot data
            '''
            # if the kept_after_redundancy_check is currently labelled "True", due to some valid hits
            if df_old.loc[i, 'kept_after_redundancy_check'] == True:
                # simply add the mds to the dataframe
                df_old_redundant_md5['redundant_md5'] = md5_checksums_in_homologue_list_that_are_also_in_query_list
                # set up the subarray
                for h in df_old_redundant_md5.index:
                    redundant_md5 = df_old_redundant_md5.loc[h, 'redundant_md5']
                    # find the location of the redundant md5 in the dataframe of query sequences
                    row_in_df_old = int(
                        df_old[
                            df_old['query_md5'] == redundant_md5].index)
                    df_old_redundant_md5.loc[
                        h, 'row_in_df_old'] = row_in_df_old
                    # add the number of valid hits to the subarray, for readability
                    df_old_redundant_md5.loc[h, 'number_of_valid_hits'] = \
                        df_old.loc[
                            row_in_df_old, 'number_of_valid_hits']
                    df_old_redundant_md5.loc[h, 'query_name'] = df_old.loc[
                        row_in_df_old, 'protein_name']
                # find the protein seq with the most hits
                index_maximum_number_of_valid_hits = df_old_redundant_md5.loc[:, 'number_of_valid_hits'].argmax()
                # label them as 'not kept' as default
                df_old_redundant_md5.loc[:, 'kept_after_redundancy_check'] = False
                # label the best as 'kept'
                df_old_redundant_md5.loc[index_maximum_number_of_valid_hits, 'kept_after_redundancy_check'] = True
                # now transfer this data to the original array (df_old)
                for m in df_old_redundant_md5.index:
                    row = df_old_redundant_md5.loc[m, 'row_in_df_old']
                    kept_after_redundancy_check = df_old_redundant_md5.loc[m, 'kept_after_redundancy_check']
                    df_old.loc[
                        row, 'kept_after_redundancy_check'] = kept_after_redundancy_check
    
    ##if this particular query has not already been rejected as a redundant sequence with a lower number of valid hits, continue and analyse the redundant md5s
    # if df_old.loc[i,'kept_after_redundancy_check'] != 'False':
    #
    #    for i in range(len(md5_checksums_in_homologue_list_that_are_also_in_query_list)):
    #        #logging.info(array_with_redundant_md5_data[i][0])
    #        #add the redundant md5
    #        redundant_md5 = md5_checksums_in_homologue_list_that_are_also_in_query_list[i]
    #        array_with_redundant_md5_data[i]['redundant_md5'] = redundant_md5
    #        #find the location of the redundant md5 in the array of query sequences
    #        row_in_df_old = list_of_query_md5_checksums.index(array_with_redundant_md5_data[i]['redundant_md5'])
    #        array_with_redundant_md5_data[i]['row_in_df_old'] = row_in_df_old
    #        #use this row_in_df_old to find the number of valid hits
    #        array_with_redundant_md5_data[i]['number_of_valid_hits'] = array_with_query_md5_data[row_in_df_old]['number_of_valid_hits']
    #        array_with_redundant_md5_data[i]['query_name'] = array_with_query_md5_data[row_in_df_old]['query_name']
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
    #        row_in_df_old = array_with_redundant_md5_data[i]['row_in_df_old']
    #        array_with_query_md5_data[row_in_df_old]['kept_after_redundancy_check'] = array_with_redundant_md5_data[i]['kept_after_redundancy_check']
    # logging.info(array_with_redundant_md5_data)
    #
    # find row in dataframe that contains the redundant md5
    #    index_row_in_df_old_containing_redundant_md5 = int(df_old[df_old['query_md5']==redundant_md5].index)
    #    #now update the df_old to show that this protein should be regarded as redundant
    # df_old.loc[index_row_in_df_old_containing_redundant_md5,'kept_after_redundancy_check'] = array_with_redundant_md5_data[i]['kept_after_redundancy_check']
    # utils.save_structured_array_to_csv(array_with_query_md5_data,file_with_redundant_md5_checksums)
    # array6 = load_structured_array_from_csv(file_with_redundant_md5_checksums,dtype_array_with_query_md5_data)
    # array6 == array_with_query_md5_data
    # np.savetxt(filexxx3, array_with_query_md5_data, fmt='%s', delimiter=',', newline='\n', header='', footer='', comments='#')
    '''create list of protein names that are kept_after_redundancy_check
    this is necennary because of a corrupted csv file, where the order of the
    data rows doesn't match the array with the md5 redundancy data
    '''
    #    for datarow in array_with_query_md5_data:
    #        if datarow['kept_after_redundancy_check'] == True:
    #            list_of_proteins_kept_after_redundancy_check.append(datarow['query_name'])
    try:
        df_oldout07_simapnonred = \
            df_old[
                df_old.kept_after_redundancy_check == True]
    except KeyError:
        # for small lists of protein sequences, there are no redundancies and an error is received as all kept_after_redundancy_check values = True
        df_oldout07_simapnonred = df_old
    list_of_proteins_kept_after_redundancy_check = list(
        df_oldout07_simapnonred.loc[:, 'protein_name'])
    
    # sort the columns alphabetically (capitals first)
    df_old = df_old.sort_index(
        axis=1)
    # save the updated DataFrame of all proteins to file
    with open(pathdict["df_oldout05_simapcsv"], 'w') as f:
        df_old.to_csv(f, sep=",", quoting=csv.QUOTE_NONNUMERIC)
    # save the updated DataFrame of nonredundant proteins to file
    with open(pathdict["df_oldout07_simapnonred"], 'w') as f:
        df_oldout07_simapnonred.to_csv(f, sep=",", quoting=csv.QUOTE_NONNUMERIC)
    
    '''
    Remove the data of 'redundant' sequences
    (sequences that fall into the SIMAP data of another query)
    The query with the MOST VALID HITS will be kept for further analysis.
    '''
    
    with open(pathdict["csv_file_with_histogram_data_normalised"], mode='r') as infile:
        reader = csv.reader(infile)
        with open(pathdict["csv_file_with_histogram_data_normalised_redundant_removed"], 'w') as outfile:
            pass
        with open(pathdict["csv_file_with_histogram_data_normalised_redundant_removed"], 'a') as outfile:
            writer = csv.writer(outfile, delimiter=',', quotechar='"', lineterminator='\n',
                                quoting=csv.QUOTE_NONNUMERIC, doublequote=True)
            rownumber = 0
            for row in reader:
                # logging.info(row)
                if rownumber == 0:
                    writer.writerow(row)
                else:
                    logging.info('redundancy check: %s' % row[0])
                    # logging.info("array_with_query_md5_data[rownumber-1]['query_name'] %s" % array_with_query_md5_data[rownumber-1]['query_name'])
                    # logging.info("array_with_query_md5_data[rownumber]['query_name'] %s\n" % array_with_query_md5_data[rownumber]['query_name'])
                    if row[0] in list_of_proteins_kept_after_redundancy_check:
                        # logging.info('match')
                        # logging.info(array5[rownumber-1]['kept_after_redundancy_check'])
                        # if convert_string_to_boolean_value(array_with_query_md5_data[rownumber-1]['kept_after_redundancy_check']):
                        # logging.info(array5[rownumber-1]['kept_after_redundancy_check'])
                        writer.writerow(row)
                rownumber += 1
    logging.info('calculate_TMD_conservation is finished.')