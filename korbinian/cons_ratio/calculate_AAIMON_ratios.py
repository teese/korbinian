import ast
import numpy as np
import csv
import os
import pickle
import korbinian
import korbinian.mtutils as utils
import zipfile
import pandas as pd

def run_calculate_AAIMON_ratios(pathdict, set_, logging):
    logging.info('~~~~~~~~~~~~starting run_calculate_AAIMON_ratios~~~~~~~~~~~~')
    df = pd.read_csv(pathdict["list_summary_csv"], sep=",", quoting=csv.QUOTE_NONNUMERIC, index_col=0)
    # iterate over the dataframe for proteins with an existing list_of_TMDs. acc = uniprot accession.
    for acc in df.loc[df['list_of_TMDs'].notnull()].loc[df['list_of_TMDs'] != 'nan'].index:
        if os.path.exists(df.loc[acc, 'homol_df_orig_zip']):
            dfh = utils.open_df_from_pickle_zip(df.loc[acc, 'homol_df_orig_zip'],filename=df.loc[acc, 'homol_df_orig_pickle'], delete_corrupt=True)
            list_of_TMDs = df.loc[acc, 'list_of_TMDs']

        if set_["run_calculate_AAIMON_ratios"]:
            ########################################################################################
            #                                                                                      #
            #                        Calculate AAIMON conservation ratios                          #
            #                                                                                      #
            ########################################################################################
            # determine if the dataframe already contains some previously analysed data
            dataframe_contains_prev_calc_AAIMON_ratios = True if 'TM01_AAIMON_ratio_mean' in df.columns else False
            logging.info('dataframe_contains_prev_calc_AAIMON_ratios = %s' % dataframe_contains_prev_calc_AAIMON_ratios)

            #assume af first that there is no previous data, and that the calculations can be re-run
            prev_calc_AAIMON_ratio_for_this_protein_exists = False
            if set_["overwrite_prev_calculated_AAIMON_ratios"] == False:
                if dataframe_contains_prev_calc_AAIMON_ratios == True:
                    cell_is_empty = np.isnan(df.loc[acc, 'TM01_AAIMON_ratio_mean'])
                    #if the cell is empty (np.nan), the AAIMON ratios haven't been calculated for that protein yet
                    if cell_is_empty:
                        prev_calc_AAIMON_ratio_for_this_protein_exists = False
                    else:
                        #if the data is there, it should be a float and not a np.nan, therefore there is no need to repeat all the calculations
                        prev_calc_AAIMON_ratio_for_this_protein_exists = True
                        logging.info('calculate_AAIMON_ratios skipped, prev_calc_AAIMON_ratio_for_this_protein_exists = True for %s' % df.loc[acc, 'uniprot_acc'])
            else:
                #if the settings says to overwrite the data, ignore its existence
                prev_calc_AAIMON_ratio_for_this_protein_exists = False

            if set_["overwrite_prev_calculated_AAIMON_ratios"] == True or prev_calc_AAIMON_ratio_for_this_protein_exists == False:
                """Find disallowed words (e.g. synthetic, patent) [AAIMON] """
                # add the list of words to the globals, to be accessed by utils.find_disallowed_words
                cr_words_not_allowed_in_description = ast.literal_eval(set_["cr_words_not_allowed_in_description"])
                # collect disallowed words in hit protein description (patent, synthetic, etc)
                dfs['cr_list_disallowed_words_in_descr'] = dfs['description'].dropna().apply(utils.find_disallowed_words, args=(cr_words_not_allowed_in_description,))
                # create a boolean column to select hits that do not contain these words in the description
                dfs['cr_disallowed_words_not_in_descr'] = dfs['cr_list_disallowed_words_in_descr'] == '[]'

                # check if there are non-IUPAC amino acids in the sequence (frequently large gaps from NG sequencing data)
                dfs['X_in_match_seq'] = 'X' in dfs['match_alignment_sequence']

                for TMD in list_of_TMDs:
                    ########################################################################################
                    #                                                                                      #
                    #        Define juxtamembrane regions associated with each TMD  [AAIMON]               #
                    #                                                                                      #
                    ########################################################################################
                    # convert the tuple of (True, 32, 53) into separate dataframes.
                    # http://stackoverflow.com/questions/29550414/how-to-split-column-of-tuples-in-pandas-dataframe

                    # dfs_with_match = dfs['%s_start_end_list_in_SW_alignment' % TMD].dropna()
                    # df_match_start_end = pd.DataFrame(dfs_with_match.values.tolist())
                    # df_match_start_end.index = dfs_with_match.index
                    # dfs["%s_start_in_SW_alignment" % TMD] = df_match_start_end[1]
                    # dfs["%s_end_in_SW_alignment" % TMD] = df_match_start_end[2]

                    last_TMD_of_acc = list_of_TMDs[-1]

                    if set_["slice_juxtamembrane_regions"] == True:
                        if TMD == "TM01":
                            dfs['start_juxta_before_TM01'] = np.where(dfs['TM01_start_in_SW_alignment'] > 0, 0, np.nan)
                            dfs['end_juxta_before_TM01'] = np.where(dfs['TM01_start_in_SW_alignment'] == 0, np.nan,
                                                                    dfs['TM01_start_in_SW_alignment'])
                            dfs['start_juxta_after_TM01'] = dfs['TM01_end_in_SW_alignment']
                            if len(list_of_TMDs) == 1:
                                # if there is only one TMD, TM01 == last_TMD_of_acc
                                dfs['end_juxta_after_TM01'] = np.where(utils.isNaN(dfs['start_juxta_after_TM01']) == True, np.nan,
                                                                       dfs['len_query_alignment_sequence'])
                            elif len(list_of_TMDs) > 1:
                                #problem('dfs["TM02_start_in_SW_alignment"] cannot exist yet, because the script iterates through the TMDs one at a time')
                                dfs['end_juxta_after_TM01'] = dfs["TM01_end_in_SW_alignment"] + (
                                (dfs["TM02_start_in_SW_alignment"] - dfs["TM01_end_in_SW_alignment"]) / 2).apply(lambda x: int(x) if not np.isnan(x) else np.nan)
                                # dfs['seq_juxta_after_TM01_in_query'] = dfs[dfs['start_juxta_after_TM01'].notnull()].apply(utils.slice_juxta_after_TMD_in_query, args = (TMD,), axis=1)
                                # dfs['seq_juxta_after_TM01_in_match'] = dfs[dfs['end_juxta_after_TM01'].notnull()].apply(utils.slice_juxta_after_TMD_in_match, args = (TMD,), axis=1)

                        # the analysis is slow, so don't repeat TM01 if there is only one TM helix in the protein
                        if len(list_of_TMDs) > 1:
                            if not TMD == "TM01" and not TMD == last_TMD_of_acc:
                                dfs = juxta_function_1(dfs, TMD)
                                # dfs['start_juxta_after_%s'%TMD] = np.where(utils.isNaN(dfs['TM%.2d_start_in_SW_alignment'%(int(TMD[2:])+1)])==True,np.nan,dfs['%s_end_in_SW_alignment'%TMD])
                                # dfs['end_juxta_before_%s'%TMD] = np.where(dfs["%s_start_in_SW_alignment"%TMD]!=0,dfs["%s_start_in_SW_alignment"%TMD],np.nan)
                                # dfs['end_juxta_after_%s'%TMD] = dfs["%s_end_in_SW_alignment"%TMD]+((dfs["TM%.2d_start_in_SW_alignment"%(int(TMD[2:])+1)]-dfs["%s_end_in_SW_alignment"%TMD])/2).apply(lambda x :int(x) if not np.isnan(x) else np.nan)
                                # dfs['start_juxta_before_%s'%TMD] = np.where(dfs["end_juxta_after_TM%.2d"%(int(TMD[2:])-1)] == dfs['end_juxta_before_%s'%TMD] ,dfs["end_juxta_after_TM%.2d"%(int(TMD[2:])-1)],dfs["end_juxta_after_TM%.2d"%(int(TMD[2:])-1)])
                                # dfs['seq_juxta_after_%s_in_query'%TMD] = dfs[dfs['start_juxta_after_%s'%TMD].notnull()].apply(utils.slice_juxta_after_TMD_in_query, args = (TMD,), axis=1)
                                # dfs['seq_juxta_after_%s_in_match'%TMD] = dfs[dfs['end_juxta_after_%s'%TMD].notnull()].apply(utils.slice_juxta_after_TMD_in_match, args = (TMD,), axis=1)

                            if TMD == last_TMD_of_acc:
                                dfs['start_juxta_before_%s'%TMD] = dfs['end_juxta_after_TM%.2d' % (int(TMD[2:]) - 1)]
                                dfs['end_juxta_before_%s'%TMD] = dfs['%s_start_in_SW_alignment'%TMD]
                                dfs['start_juxta_after_%s'%TMD] = np.where(
                                    dfs['%s_end_in_SW_alignment'%TMD] == dfs['len_query_alignment_sequence'], np.nan,
                                    dfs['%s_end_in_SW_alignment'%TMD])
                                dfs['end_juxta_after_%s'%TMD] = np.where(utils.isNaN(dfs['start_juxta_after_%s'%TMD]) == True, np.nan,
                                                                           dfs['len_query_alignment_sequence'])
                                # dfs['seq_juxta_after_%s_in_query'%TMD] = dfs[dfs['start_juxta_after_%s'%TMD].notnull()].apply(utils.slice_juxta_after_TMD_in_query, args = (TMD,), axis=1)
                                # dfs['seq_juxta_after_%s_in_query'%TMD] = dfs.query_alignment_sequence[int(dfs['start_juxta_after_TM10']):int(dfs['end_juxta_after_TM10'])]
                                # dfs['seq_juxta_after_%s_in_match'%TMD] =

                        last_TMD_of_acc = list_of_TMDs[-1]
                        dfs['seq_juxta_before_%s_in_query'%TMD] = dfs[dfs['start_juxta_before_%s'%TMD].notnull()].apply(
                            utils.slice_juxta_before_TMD_in_query, args=(TMD,), axis=1)
                        dfs['seq_juxta_before_%s_in_match'%TMD] = dfs[dfs['start_juxta_before_%s'%TMD].notnull()].apply(
                            utils.slice_juxta_before_TMD_in_match, args=(TMD,), axis=1)
                        if not TMD == last_TMD_of_acc:
                            dfs['seq_juxta_after_%s_in_query'%TMD] = dfs[dfs['end_juxta_after_%s'%TMD].notnull()].apply(
                                utils.slice_juxta_after_TMD_in_query, args=(TMD,), axis=1)
                            dfs['seq_juxta_after_%s_in_match'%TMD] = dfs[dfs['end_juxta_after_%s'%TMD].notnull()].apply(
                                utils.slice_juxta_after_TMD_in_match, args=(TMD,), axis=1)
                        else:
                            dfs['seq_juxta_after_%s_in_query'%TMD] = np.nan
                            dfs['seq_juxta_after_%s_in_match'%TMD] = np.nan
                            for hit in dfs.index:
                                if not utils.isNaN(dfs['start_juxta_after_%s'%TMD])[hit]:
                                    # altered to .loc rather than ['seq_juxta_after_%s_in_match'%TMD][hit] after SettingWithCopyWarning
                                    dfs.loc[hit, 'seq_juxta_after_%s_in_match'%TMD] = dfs.match_alignment_sequence[hit][
                                                                                        int(dfs.loc[hit, "start_juxta_after_%s"%TMD]):int(
                                                                                            dfs.loc[hit, "end_juxta_after_%s"%TMD])]
                                    dfs.loc[hit, 'seq_juxta_after_%s_in_query'%TMD] = dfs.query_alignment_sequence[hit][
                                                                                        int(dfs.loc[hit, "start_juxta_after_%s"%TMD]):int(
                                                                                            dfs.loc[hit, "end_juxta_after_%s"%TMD])]

                if n_TMDs_w_homol == 0:
                    # there was no data obtained from the csv file, which probably means the original XML file was not properly downloaded
                    # there is no need to continue the script for this protein
                    logging.warning('%s does not contain any valid homologues at all. '
                                    'CSV and/or XML file is damaged. Re-downloading simap homologue XML may be necessary. '
                                    'AAIMON ratios will not be calculated for this protein.' % protein_name)

                if n_TMDs_w_homol > 0:
                    # create a boolean column that describel whether the sequence is above the minimum gapped identity
                    cr_minimum_identity_of_full_protein = set_["cr_minimum_identity_of_full_protein"]
                    dfs['cr_gapped_ident_above_cutoff'] = dfs['FASTA_gapped_identity'] > cr_minimum_identity_of_full_protein

                    '''
                    3) values associated each TMD, such as average AAIMON ratio
                    '''
                    # calculate the nonTMD percentage identity and gaps
                    korbinian.cons_ratio.calc_nonTMD_perc_ident_and_gaps(acc, dfs, set_, df, list_of_TMDs, logging)

                    # calculate AAISMON etc for each TMD
                    for TMD in list_of_TMDs:
                        df, dfs = korbinian.cons_ratio.calc_AAIMON(acc, TMD, dfs, set_, df, logging)

                    homol_cr_ratios_zip = df.loc[acc, 'homol_cr_ratios_zip']
                    AAIMON_hist_path_prefix = df.loc[acc, 'AAIMON_hist_path_prefix']
                    # save histograms for each TMD of that protein, with relative conservation
                    korbinian.cons_ratio.save_hist_AAIMON_ratio_single_protein(dfs, set_, list_of_TMDs, homol_cr_ratios_zip, AAIMON_hist_path_prefix)

                # remove columns to make output csv smaller
                if set_['drop_columns_to_reduce_csv_filesize']:
                    list_cols_to_drop = ['match_alignment_sequence', 'query_alignment_sequence', 'alignment_markup',
                                         'nonTMD_seq_query', 'nonTMD_markup']
                    for col in list_cols_to_drop:
                        if col in dfs.columns:
                            dfs.drop(col, axis=1, inplace=True)
                #dfs.to_csv(df.loc[acc, 'SIMAP_csv_analysed_path'], sep=",", quoting=csv.QUOTE_NONNUMERIC)
                # save dfs with homologues for a single protein, as a single zipped csv
                #utils.save_df_to_csv_zip(dfs, df.loc[acc, 'homol_df_orig_zip'], open_method="w")
                print("important, work out a save method here" / 4)

                df.loc[acc, 'num_hits_with_SW_align_node'] = dfs['hit_contains_SW_node'].value_counts()[True]
                logging.info('num_hits_with_SW_align_node: %s' % df.loc[acc, 'num_hits_with_SW_align_node'])

        df.loc[acc, "create_fasta_or_cons_ratio_single_protein"] = True
        # save to csv after each protein is analysed, incrementally adding the extra data
        df.to_csv(pathdict["list_summary_csv"], sep=",", quoting=csv.QUOTE_NONNUMERIC)

        logging.info("~~~~~run_calculate_AAIMON_ratios is finished ~~~~~")