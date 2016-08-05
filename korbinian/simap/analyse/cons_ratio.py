import ast
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import csv
import os
import tarfile
import korbinian
import korbinian.mtutils as utils
import zipfile

def calculate_AAIMON_ratios(pathdict, set_, logging):
    logging.info('~~~~~~~~~~~~starting calculate_AAIMON_ratios~~~~~~~~~~~~')
    overwrite_prev_calculated_AAIMON_ratios = set_["overwrite_prev_calculated_AAIMON_ratios"]
    df = pd.read_csv(pathdict["list_summary_csv"])

    #filter to remove sequences where no TMDs are found,
    df = df.loc[df['list_of_TMDs'].notnull()]
    # #filter to remove sequences where no TMDs are found (if string)
    # df.loc[df['list_of_TMDs'] != 'nan']
    # #filter to remove sequences where no TMDs are found (if string)
    # df = df.loc[df['list_of_TMDs'] != 'nan']
    #determine if the dataframe already contains some previously analysed data
    dataframe_contains_prev_calc_AAIMON_ratios = True if 'TM01_AAIMON_ratio_mean' in df.columns else False
    logging.info('dataframe_contains_prev_calc_AAIMON_ratios = %s' % dataframe_contains_prev_calc_AAIMON_ratios)
    #iterate over the dataframe for proteins with an existing list_of_TMDs. Note that acc = uniprot accession here.
    for acc in df.loc[df['list_of_TMDs'] != 'nan'].index:
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
                    logging.info('calculate_AAIMON_ratios skipped, prev_calc_AAIMON_ratio_for_this_protein_exists = True for %s' % df.loc[acc, 'uniprot_acc'])
        else:
            #if the settings says to overwrite the data, ignore its existence
            prev_calc_AAIMON_ratio_for_this_protein_exists = False
        if overwrite_prev_calculated_AAIMON_ratios == True or prev_calc_AAIMON_ratio_for_this_protein_exists == False:
            protein_name = df.loc[acc, 'protein_name']
            logging.info('%s' % protein_name)
            # SIMAP_csv_from_XML_tarfile = df.loc[acc, 'SIMAP_csv_from_XML_tarfile']
            # if os.path.exists(SIMAP_csv_from_XML_tarfile):
            #     try:
            #         with tarfile.open(df.loc[acc, 'SIMAP_csv_from_XML_tarfile'], mode='r:gz') as tar:
            #             if df.loc[acc, 'SIMAP_csv_from_XML'] in [tarinfo.name for tarinfo in tar]:
            #                 SIMAP_csv_from_XML_exists = True
            #             else:
            #                 SIMAP_csv_from_XML_exists = False
            #     except (EOFError, tarfile.ReadError, OSError):
            #         #file may be corrupted, if script stopped unexpectedly before compression was finished
            #         logging.info('%s seems to be corrupted. File will be deleted.'
            #                      'Will need to be parsed again next time program is run' %
            #                       df.loc[acc, 'SIMAP_csv_from_XML_tarfile'])
            #         SIMAP_csv_from_XML_exists = False
            #         os.remove(df.loc[acc, 'SIMAP_csv_from_XML_tarfile'])
            # else:
            #     SIMAP_csv_from_XML_exists = False

            SIMAP_csv_from_XML_exists = False
            if os.path.exists(df.loc[acc, 'homol_csv_zip']):
                try:
                    with zipfile.ZipFile(df.loc[acc, 'homol_csv_zip'], "r", zipfile.ZIP_DEFLATED) as openzip:
                        if openzip.namelist()[0][:-4] == ".csv":
                            SIMAP_csv_from_XML_exists = True
                except:
                    #file may be corrupted, if script stopped unexpectedly before compression was finished
                    logging.info('%s seems to be corrupted. File will be deleted.'
                                 'Will need to be parsed again next time program is run' % df.loc[acc, 'homol_csv_zip'])
            if SIMAP_csv_from_XML_exists:
                # run the analysis function, which slices homologues, calculates identity, and saves output files
                analyse_homologues_single_protein(df.loc[acc, 'homol_csv_zip'], acc, protein_name, set_, df, pathdict, logging)

    logging.info('calculate_AAIMON_ratios is finished.')

def juxta_function_1(dfs, TMD):
    dfs['start_juxta_after_%s'%TMD] = np.where(utils.isNaN(dfs['TM%.2d_start_in_SW_alignment'%(int(TMD[2:])+1)])==True,np.nan,dfs['%s_end_in_SW_alignment'%TMD])
    dfs['end_juxta_before_%s'%TMD] = np.where(dfs["%s_start_in_SW_alignment"%TMD]!=0,dfs["%s_start_in_SW_alignment"%TMD],np.nan)
    dfs['end_juxta_after_%s'%TMD] = dfs["%s_end_in_SW_alignment"%TMD]+((dfs["TM%.2d_start_in_SW_alignment"%(int(TMD[2:])+1)]-dfs["%s_end_in_SW_alignment"%TMD])/2).apply(lambda x :int(x) if not np.isnan(x) else np.nan)
    dfs['start_juxta_before_%s'%TMD] = np.where(dfs["end_juxta_after_TM%.2d"%(int(TMD[2:])-1)] == dfs['end_juxta_before_%s'%TMD] ,dfs["end_juxta_after_TM%.2d"%(int(TMD[2:])-1)],dfs["end_juxta_after_TM%.2d"%(int(TMD[2:])-1)])
    return dfs

def analyse_homologues_single_protein(homol_csv_zip, acc, protein_name, set_, df, pathdict, logging):
    # with tarfile.open(SIMAP_csv_from_XML_tarfile, 'r:gz') as tar_in:
    #SIMAP_csv_from_XML_extracted = tar_in.extractfile(df.loc[acc, 'SIMAP_csv_from_XML'])
    dfs = utils.open_df_from_csv_zip(homol_csv_zip)
    # reopen the csv file containing all the homologue data for that particular protein as a pandas dataframe (labelled Data Frame SIMAP, or dfs)
    #dfs = pd.read_csv(SIMAP_csv_from_XML_extracted, sep=",", index_col=0, quoting=csv.QUOTE_NONNUMERIC)
    # create a column with the length of the Smith Waterman alignment (same length for query, markup and match)
    if True in dfs['hit_contains_SW_node']:
        at_least_one_hit_contains_SW_node = True
    else:
        at_least_one_hit_contains_SW_node = False
    if at_least_one_hit_contains_SW_node:
        try:
            dfs['len_query_alignment_sequence'] = dfs['query_alignment_sequence'].str.len()
        except KeyError:
            # dataframe does not contain query_alignment_sequence,
            # which means that the XML file is probably damaged somehow
            logging.warning('SIMAP_csv_from_XML seems to be damaged for %s' % protein_name)
            dfs['query_alignment_sequence'] = ''
            dfs['len_query_alignment_sequence'] = 0
    else:
        dfs['query_alignment_sequence'] = ''
        dfs['len_query_alignment_sequence'] = 0

    ########################################################################################
    #                                                                                      #
    #             Find disallowed words (e.g. synthetic, patent)                           #
    #                                                                                      #
    ########################################################################################
    # add the list of words to the globals, to be accessed by utils.find_disallowed_words
    cr_words_not_allowed_in_description = ast.literal_eval(set_["cr_words_not_allowed_in_description"])
    # collect disallowed words in hit protein description (patent, synthetic, etc)
    dfs['cr_list_disallowed_words_in_descr'] = dfs['uniprot_description'].dropna().apply(utils.find_disallowed_words, args=(cr_words_not_allowed_in_description,))
    # create a boolean column to select hits that do not contain these words in the description
    dfs['cr_disallowed_words_not_in_descr'] = dfs['cr_list_disallowed_words_in_descr'] == '[]'

    # check if there are non-IUPAC amino acids in the sequence (frequently large gaps from NG sequencing data)
    dfs['X_in_match_seq'] = 'X' in dfs['match_alignment_sequence']

    list_of_TMDs = ast.literal_eval(df.loc[acc, 'list_of_TMDs'])
    number_of_TMDs_containing_some_homologue_data = 0

    for TMD in list_of_TMDs:
        ########################################################################################
        #                                                                                      #
        #                                Slice out TMD regions                                 #
        #                                                                                      #
        ########################################################################################
        '''slice the TMD regions out of the alignment markup and match sequences
        the method first finds the indices of the TMD region in the query, and uses these indices to slice
        the filters are not applied yet, so the sequences can be viewed in the csv file for analysis
        '''
        # create regex string for each TMD
        query_TMD_sequence = df.loc[acc, '%s_seq'%TMD]
        # create TMD regex search string (e.g. L.*L.*M.*L.*L.*M.*L.*L.*M.*L.*L.*M.*L.*L.*M.*)
        TMD_regex_ss = utils.create_regex_string(query_TMD_sequence)
        # select to only include data where the XML contained a SW node, and then apply function for a regex search
        query_seqs_ser = dfs.query('hit_contains_SW_node == True')['query_alignment_sequence']

        dfs['%s_start_end_tuple_in_SW_alignment'%TMD] = query_seqs_ser.apply(utils.get_start_and_end_of_TMD_in_query,args=(TMD_regex_ss,))
        '''the output of the regex search is a tuple with three components.
        1) a match (e.g. True),
        2) the start of the match (e.g. 124),
        3) the end of the match (e.g. 144)
        '''
        # for convenience, separate the components into separate columns
        columns_from_regex_output = ['%s_in_SW_alignment'%TMD, '%s_start_in_SW_alignment'%TMD,'%s_end_in_SW_alignment'%TMD]
        # n = the index (1,2,3) in the tuple
        # col = column item in the list (e.g. 'TM09_in_SW_alignment')
        for n, col in enumerate(columns_from_regex_output):
            # first filter to analyse only columns that contain a SW node
            df_match = dfs.query('hit_contains_SW_node == True')
            # add a new column which is named TM01_start, etc, and insert the appropriate integer (start or stop) or bool from the tuple
            dfs[col] = df_match['%s_start_end_tuple_in_SW_alignment'%TMD].apply(lambda x: x[n])

        ########################################################################################
        #                                                                                      #
        #           Define juxtamembrane regions associated with each TMD                      #
        #                                                                                      #
        ########################################################################################
        # convert the tuple of (True, 32, 53) into separate dataframes.
        # http://stackoverflow.com/questions/29550414/how-to-split-column-of-tuples-in-pandas-dataframe
        df_match = pd.DataFrame(dfs['%s_start_end_tuple_in_SW_alignment' % TMD].values.tolist())
        dfs["%s_start_in_SW_alignment" % TMD] = df_match[1]
        dfs["%s_end_in_SW_alignment" % TMD] = df_match[2]

        last_TMD_of_acc = list_of_TMDs[-1]

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
                problem('dfs["TM02_start_in_SW_alignment"] cannot exist yet, because the script iterates through the TMDs one at a time')
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


        ########################################################################################
        #                                                                                      #
        #                Apply function to slice homologues and count gaps                     #
        #                                                                                      #
        ########################################################################################

        # in some cases, there is no data to obtain as the hit_contains_SW_node = False for too many sequences, giving no start_in_SW_alignment
        number_of_rows_containing_data = dfs[dfs['%s_start_in_SW_alignment'%TMD].notnull()].shape[0]
        if number_of_rows_containing_data == 0:
            logging.info('%s does not have any valid homologues for %s. '
                         'Re-downloading simap homologue XML may be necessary.' % (protein_name, TMD))
        if number_of_rows_containing_data != 0:
            number_of_TMDs_containing_some_homologue_data += 1
            # apply the slicing function to the homologues
            dfs = slice_homologues_and_count_gaps(acc, TMD, df, dfs, set_)

    if number_of_TMDs_containing_some_homologue_data == 0:
        # there was no data obtained from the csv file, which probably means the original XML file was not properly downloaded
        # there is no need to continue the script for this protein
        logging.warning('%s does not contain any valid homologues at all. '
                        'CSV and/or XML file is damaged. Re-downloading simap homologue XML may be necessary. '
                        'AAIMON ratios will not be calculated for this protein.' % protein_name)

    if number_of_TMDs_containing_some_homologue_data > 0:
        # create a boolean column that describel whether the sequence is above the minimum gapped identity
        cr_minimum_identity_of_full_protein = set_["cr_minimum_identity_of_full_protein"]
        dfs['cr_gapped_ident_above_cutoff'] = dfs['FASTA_gapped_identity'] > cr_minimum_identity_of_full_protein

        ########################################################################################
        #                                                                                      #
        #                                 nonTMD calculations                                  #
        #                                                                                      #
        ########################################################################################
        # filter to only analyse sequences with tho following:
        # 1) full protein identity above cutoff
        # 2)containing smith waterman alignment in XML file
        # 3) hit description lacks disallowedd words (patent, synthetic, etc)
        dfs_filt = dfs.query('cr_gapped_ident_above_cutoff == True and '
                             'hit_contains_SW_node == True and '
                             'cr_disallowed_words_not_in_descr == True')

        if dfs_filt.shape[0] > 0:

            # check if all tmds are in SW alignment
            # create a list of columns to reindex the DataFrame
            list_columns_TMD_in_SW_alignment = []
            for TMD in list_of_TMDs:
                # TMD found by regex
                list_columns_TMD_in_SW_alignment.append('%s_in_SW_alignment'%TMD)
                # TMD matching useful sequence
                list_columns_TMD_in_SW_alignment.append('%s_in_SW_align_match'%TMD)

            # create a slice of the filtered dataframe that only contains the relevant columns (N.B. copy=False, this will provide a view, not a copy)
            df2 = dfs_filt.reindex(index=dfs.index, columns=list_columns_TMD_in_SW_alignment, copy=False)
            # create a new column in the original dataframe that shows that ALL TMDs have been found in the SW alignment
            dfs['all_tmds_in_SW_alignment'] = df2.dropna().all(axis=1)
            # filter to contain only hits with all tmds
            # create a copy of the original dfs dataframe containing only hits where all tmds are found in the match
            dfs_nonTMD = dfs.loc[dfs['all_tmds_in_SW_alignment'].notnull()].query('all_tmds_in_SW_alignment == True')
            # filter to contain only hits that do not include disallowed words
            dfs_nonTMD = dfs_nonTMD.query('cr_disallowed_words_not_in_descr == True')
            # filter to contain only hits where the index for the TMD is present
            first_TMD_start_index = '%s_start_in_SW_alignment' % list_of_TMDs[0]

            dfs_nonTMD = dfs_nonTMD.loc[dfs_nonTMD[first_TMD_start_index].notnull()]
            dfs_nonTMD['nonTMD_index_tuple_first'] = dfs_nonTMD[first_TMD_start_index].apply(lambda x: (0, int(x)))

            # create start and stop indices for all sections between tmds
            # the start of the last nonTMD section will be the end of the last TMD
            dfs_nonTMD['nonTMD_index_tuple_last0'] = dfs_nonTMD[
                '%s_end_in_SW_alignment' % list_of_TMDs[-1]].dropna().astype('int32')
            # the end of the last nonTMD section will be the end of the full alignment sequence
            dfs_nonTMD['nonTMD_index_tuple_last1'] = dfs_nonTMD['len_query_alignment_sequence'].dropna().astype('int32')
            # join to make a tuple
            # dfs_nonTMD['nonTMD_index_tuple_last'] = dfs_nonTMD[['nonTMD_index_tuple_last0', 'nonTMD_index_tuple_last1']].apply(tuple, axis=1)

            # create the index tuple
            dfs_nonTMD['nonTMD_index_tuple_last'] = dfs_nonTMD.apply(utils.create_indextuple_nonTMD_last, axis=1)

            ########################################################################################
            #                                                                                      #
            #                    create the indices for the nonTMD region                          #
            #                                                                                      #
            ########################################################################################

            # for each TMD EXCEPT the last, which ends at the sequence end, create the indices for the nonTMD region (after the TMD)
            for TM_Nr in range(len(list_of_TMDs) - 1):
                # the TMD is the equivalent item in the list
                TMD = list_of_TMDs[TM_Nr]
                # the next TMD, which contains the end index, is the next item in the list
                next_TMD = list_of_TMDs[TM_Nr + 1]
                # select only the columns in the dataframe that are of interest, and change the data type to integer
                index_columns = ['%s_end_in_SW_alignment'%TMD, '%s_start_in_SW_alignment' % next_TMD]
                dfs_nonTMD[index_columns] = dfs_nonTMD[index_columns].astype('int64')
                # create a tuple containing the indices for the nonTMD sequence regions in between each TMD (middle indices)
                dfs_nonTMD['nonTMD_index_%s'%TMD] = tuple(zip(dfs_nonTMD['%s_end_in_SW_alignment'%TMD],
                                                                dfs_nonTMD['%s_start_in_SW_alignment' % next_TMD]))

            # now join all the indices together to make one tuple of tuples for the non-TMD region
            # dfs_nonTMD = dfs.query('all_tmds_in_SW_alignment == True')
            dfs_nonTMD['nested_tuple_indices_all_nonTMD_regions'] = dfs_nonTMD[['nonTMD_index_tuple_last']].apply(tuple,axis=1)

            # create a view of the dataframe that contains only the desired columns
            # first create a list of the desired columns
            # start with the first tuple, from 0 to the start of the first TMD
            list_of_nonTMD_index_columns = ['nonTMD_index_tuple_first']
            # create a nonTMD region for each of the TMDs (except the last one)
            list_from_TMs = ['nonTMD_index_%s'%TMD2 for TMD2 in list_of_TMDs[:-1]]
            # join lists
            list_of_nonTMD_index_columns = list_of_nonTMD_index_columns + list_from_TMs
            list_of_nonTMD_index_columns += ['nonTMD_index_tuple_last']

            # create the new view by reindexing the dataframe with the list of desired columns
            dfs_tuples = dfs_nonTMD.reindex(index=dfs_nonTMD.index, columns=list_of_nonTMD_index_columns, copy=False)
            # now for convenience, these tuples can be combined together to form one column, with a tuple of tuples
            # first convert all values in each row to a list, excluding the index column
            list_tuple_indices_all_nonTMD_regions = list(dfs_tuples.itertuples(index=False))
            # convert to a series, and reindex with the original index from the dataframe
            tuples_series = pd.Series(list_tuple_indices_all_nonTMD_regions, index=dfs_nonTMD.index)
            # for some reason, the tuples are a pandas object "Pandas(nonTMD_index_tuple_first=(0, 592), nonTMD_index_tuple_last=(615, 618))"
            # convert to simple tuples
            tuples_series = tuples_series.apply(lambda x: tuple(x))

            dfs_nonTMD['nested_tuple_indices_all_nonTMD_regions'] = tuples_series
            # change to a string, in case this solves the weird effect with only the last tuple shown
            dfs_nonTMD['nested_tuple_indices_all_nonTMD_regions'] = dfs_nonTMD['nested_tuple_indices_all_nonTMD_regions'].astype(str)
            # you can test that the original index is maintained as follows:
            # add the series as a new column in the original dataframe. Missing data (when not all TMDs found) will be filled using np.nan
            dfs['nested_tuple_indices_all_nonTMD_regions'] = dfs_nonTMD['nested_tuple_indices_all_nonTMD_regions']

            # filter to remove incomplete sequences
            dfs_nonTMD = dfs_nonTMD.query('all_tmds_in_SW_alignment == True')
            # define the string for slicing as a numpy array
            # use the numpy vectorize function, which effectively applies the function in a for loop (not optimized for speed)
            # dfs_nonTMD['nonTMD_seq_query'] = np.vectorize(utils.slice_with_nested_tuple)(np.array(dfs_nonTMD['query_alignment_sequence']),np.array(dfs_nonTMD['nested_tuple_indices_all_nonTMD_regions']))
            # dfs_nonTMD['nonTMD_markup'] = np.vectorize(utils.slice_with_nested_tuple)(np.array(dfs_nonTMD['alignment_markup']),np.array(dfs_nonTMD['nested_tuple_indices_all_nonTMD_regions']))
            # dfs_nonTMD['nonTMD_seq_match'] = np.vectorize(utils.slice_with_nested_tuple)(np.array(dfs_nonTMD['match_alignment_sequence']),np.array(dfs_nonTMD['nested_tuple_indices_all_nonTMD_regions']))

            ########################################################################################
            #                                                                                      #
            #                  slice out the nonTMD region for each homologue                      #
            #                                                                                      #
            ########################################################################################
            # due to problems with np.vectorize and the pandas methods, slice the sequences one at a time with a simple 'for loop'
            # for each hit, perform the slice
            for hit in dfs_nonTMD.index:
                dfs_nonTMD.loc[hit, 'nonTMD_seq_query'] = utils.slice_with_nested_tuple(
                    dfs_nonTMD.loc[hit, 'query_alignment_sequence'],
                    dfs_nonTMD.loc[hit, 'nested_tuple_indices_all_nonTMD_regions'])
                dfs_nonTMD.loc[hit, 'nonTMD_markup'] = utils.slice_with_nested_tuple(
                    dfs_nonTMD.loc[hit, 'alignment_markup'],
                    dfs_nonTMD.loc[hit, 'nested_tuple_indices_all_nonTMD_regions'])
                dfs_nonTMD.loc[hit, 'nonTMD_seq_match'] = utils.slice_with_nested_tuple(
                    dfs_nonTMD.loc[hit, 'match_alignment_sequence'],
                    dfs_nonTMD.loc[hit, 'nested_tuple_indices_all_nonTMD_regions'])
            # transfer to original dataframe (index should still match the original, partial seqs will be filled with np.nan)
            dfs['nonTMD_seq_query'] = dfs_nonTMD['nonTMD_seq_query']
            dfs['nonTMD_markup'] = dfs_nonTMD['nonTMD_markup']
            dfs['nonTMD_seq_match'] = dfs_nonTMD['nonTMD_seq_match']

            ########################################################################################
            #                                                                                      #
            #                calculate nonTMD len, percent identity, gaps, etc                     #
            #                                                                                      #
            ########################################################################################
            # calculate identical residues in the nonTMD region (simply count the pipes '|' in the markup sequence)
            dfs['nonTMD_num_ident_res'] = dfs['nonTMD_markup'].dropna().apply(lambda x: x.count('|'))
            # calculate similar residues in the nonTMD region (simply count the colons ':' in the markup sequence)
            dfs['nonTMD_num_sim_res'] = dfs['nonTMD_markup'].dropna().apply(lambda x: x.count(':'))
            # add the identical and similar residues together to get the total number of similar + identical residues
            dfs['nonTMD_num_sim_plus_ident_res'] = dfs['nonTMD_num_ident_res'] + dfs['nonTMD_num_sim_res']

            # count the gaps in the nonTMD sequence of the query
            dfs['nonTMD_q_num_gaps'] = dfs['nonTMD_seq_query'].dropna().apply(lambda x: x.count('-'))
            # count the gaps in the nonTMD sequence of the match
            dfs['nonTMD_m_num_gaps'] = dfs['nonTMD_seq_match'].dropna().apply(lambda x: x.count('-'))
            # calculate the length of the nonTMD sequences, which may include gaps
            dfs['len_nonTMD_seq_query'] = dfs['nonTMD_seq_query'].dropna().str.len()
            dfs['len_nonTMD_seq_match'] = dfs['nonTMD_seq_match'].dropna().str.len()
            # calculate the number aligned sequences, excluding gaps (length of query, or length of match, whichever is shorter)
            dfs['len_nonTMD_align'] = dfs[['len_nonTMD_seq_query', 'len_nonTMD_seq_match']].dropna(how='all').min(axis=1)

            # calculate the length of the nonTMD sequence excluding gaps
            dfs['len_nonTMD_q_excl_gaps'] = dfs['len_nonTMD_seq_query'] - dfs['nonTMD_q_num_gaps']
            dfs['len_nonTMD_m_excl_gaps'] = dfs['len_nonTMD_seq_match'] - dfs['nonTMD_m_num_gaps']
            # calculate the lenth of the alignment by finding which seq excl gaps is smaller
            dfs['len_nonTMD_align'] = dfs[['len_nonTMD_q_excl_gaps', 'len_nonTMD_m_excl_gaps']].min(axis=1)

            # calculate the percentage identity of the nonTMD region (number of identical residues divided by the length excluding gaps)
            # used for the Amino Acid Identity : Membranous over Nonmembranous (AAIMON ratio)
            # note that the length = length of the aligned residues excluding gaps
            dfs['nonTMD_perc_ident'] = dfs['nonTMD_num_ident_res'] / dfs['len_nonTMD_align']
            dfs['nonTMD_perc_sim'] = dfs['nonTMD_num_sim_res'] / dfs['len_nonTMD_align']
            dfs['nonTMD_perc_sim_plus_ident'] = dfs['nonTMD_num_sim_plus_ident_res'] / dfs['len_nonTMD_align']
            # calculate the average number of gaps per residue in the nonTMD alignment
            # filter to analyse only sequences that are valid (length > 0)
            dfs_filt_gaps = dfs.loc[dfs['len_nonTMD_q_excl_gaps'] != 0]
            # calculate number of gaps in query AND match
            dfs_filt_gaps['nonTMD_qm_num_gaps'] = dfs_filt_gaps['nonTMD_q_num_gaps'] + dfs_filt_gaps['nonTMD_m_num_gaps']
            # add to simap dataframe
            dfs['nonTMD_qm_num_gaps'] = dfs_filt_gaps['nonTMD_qm_num_gaps']
            # gaps per query residue for both query and match = ((gaps in query + gaps in match)/2))/length of query excluding gaps
            dfs['nonTMD_qm_gaps_per_q_residue'] = dfs_filt_gaps['nonTMD_qm_num_gaps'] / 2 / dfs_filt_gaps['len_nonTMD_q_excl_gaps']

            # SSR ratio calculations take a long time and show no difference to AAIMON
            SSR_ratio_calculations = False
            if SSR_ratio_calculations:
                conduct_ssr_ratio_calculations()

            # create optional filter string, related to X in the full sequence
            cr_X_filt_str = " and X_in_match_seq == False" if set_["cr_X_allowed_in_full_seq"] == False else ""

            # re-filter the original dataframe to create another copy with the desired sequences
            cons_ratio_query_str = 'cr_gapped_ident_above_cutoff == True and ' \
                                   'hit_contains_SW_node == True and ' \
                                   'cr_disallowed_words_not_in_descr == True' \
                                   '{}'.format(cr_X_filt_str)

            dfs_filt = dfs.query(cons_ratio_query_str)
            ########################################################################################
            #                                                                                      #
            #     calculate average nonTMD len, perc ident, etc and add to list of proteins        #
            #                                                                                      #
            ########################################################################################

            '''Calculate average values, add to original dataframe.
               1) values associated with the FASTA output of SIMAP
            '''
            # fasta identity
            df.loc[acc, 'FASTA_ident_mean'] = float('%0.2f' % dfs['FASTA_identity'].mean())
            # number of identical residues in FASTA alignment
            dfs['FASTA_num_ident_res'] = dfs_filt['FASTA_identity'] / 100 * dfs_filt['FASTA_overlap']
            df.loc[acc, 'FASTA_num_ident_res'] = float('%0.2f' % dfs_filt['FASTA_identity'].mean())

            '''2) values associated with the nonTMD region
            '''
            # add the average values regarding the nonTMD region to the original file/dataframe with each protein
            df.loc[acc, 'len_nonTMD_seq_match_mean'] = float('%0.2f' % dfs_filt['len_nonTMD_seq_match'].dropna().mean())
            df.loc[acc, 'nonTMD_perc_ident_mean'] = float('%0.3f' % dfs_filt['nonTMD_perc_ident'].dropna().mean())
            df.loc[acc, 'nonTMD_perc_sim_mean'] = float('%0.3f' % dfs_filt['nonTMD_perc_sim'].dropna().mean())
            df.loc[acc, 'nonTMD_perc_sim_plus_ident_mean'] = float('%0.3f' % dfs_filt['nonTMD_perc_sim_plus_ident'].dropna().mean())
            df.loc[acc, 'len_nonTMD_align_mean'] = float('%0.2f' % dfs_filt['len_nonTMD_align'].dropna().mean())
            df.loc[acc, 'nonTMD_qm_gaps_per_q_residue_mean'] = float('%0.2f' % dfs_filt['nonTMD_qm_gaps_per_q_residue'].dropna().mean())
            logging.info('nonTMD_qm_gaps_per_q_residue : %0.5f' % df.loc[acc, 'nonTMD_qm_gaps_per_q_residue_mean'])

            '''3) values associated each TMD, such as average AAIMON ratio
            '''
            # calculate AAISMON etc for each TMD
            for TMD in list_of_TMDs:
                df, dfs = calc_AAIMON(acc, TMD, dfs, dfs_filt, set_, df, logging)

            save_hist_and_fastA(acc, dfs, dfs_filt, set_, df, list_of_TMDs, logging)

    df.loc[acc, "analyse_homologues_single_protein"] = True
    # save to csv after each protein is analysed, incrementally adding the extra data
    df.to_csv(pathdict["list_summary_csv"], sep=",", quoting=csv.QUOTE_NONNUMERIC)

def slice_homologues_and_count_gaps(acc, TMD, df, dfs, set_):
    len_query_TMD = len(df.loc[acc, '%s_seq'%TMD])
    # use small throwaway functions to slice each TMD from the query, markup and match sequences
    # notnull() removes missing data
    # explanation: dfs['new_column_with_selected_seq'] = dfs[dfs['only_rows_containing_data]].apply(utils.slice_function_that_specifies_columns_with_start_and_stop)
    dfs['%s_SW_query_seq'%TMD] = dfs[dfs['%s_start_in_SW_alignment'%TMD].notnull()].apply(utils.slice_SW_query_TMD_seq,args=(TMD,), axis=1)
    dfs['%s_SW_markup_seq'%TMD] = dfs[dfs['%s_start_in_SW_alignment'%TMD].notnull()].apply(utils.slice_SW_markup_TMD,args=(TMD,), axis=1)
    dfs['%s_SW_match_seq'%TMD] = dfs[dfs['%s_start_in_SW_alignment'%TMD].notnull()].apply(utils.slice_SW_match_TMD_seq,args=(TMD,), axis=1)
    # count the number of gaps in the query and match sequences
    dfs['%s_SW_query_num_gaps'%TMD] = dfs['%s_SW_query_seq'%TMD].dropna().apply(lambda x: x.count('-'))
    dfs['%s_SW_match_num_gaps'%TMD] = dfs['%s_SW_match_seq'%TMD].dropna().apply(lambda x: x.count('-'))
    # calculate the length of the match TMD seq excluding gaps
    dfs['%s_SW_m_seq_len'%TMD] = dfs['%s_SW_match_seq'%TMD].str.len()
    # for the alignment length, take the smallest value from the length of query or match
    # this will exclude gaps from the length in the following calculations, preventing false "low conservation" where the query TMD is much longer than the match TMD)
    # note that for most calculations this is somewhat redundant, because the max number of acceptable gaps in sequence is probable ~2
    dfs['%s_SW_align_len'%TMD] = dfs['%s_SW_m_seq_len'%TMD].apply(lambda x: x if x < len_query_TMD else len_query_TMD)
    # create a boolean column that allows filtering by the accepted number of gaps, according to the settings file
    dfs['%s_cr_SW_query_acceptable_n_gaps'%TMD] = dfs['%s_SW_query_num_gaps'%TMD] <= set_["cr_max_n_gaps_in_query_TMD"]
    dfs['%s_cr_SW_match_acceptable_n_gaps'%TMD] = dfs['%s_SW_match_num_gaps'%TMD] <= set_["cr_max_n_gaps_in_match_TMD"]
    # count identical residues between query and match TMDs by counting the number of pipes in the markup string
    dfs['%s_SW_num_ident_res'%TMD] = dfs['%s_SW_markup_seq'%TMD].dropna().apply(lambda x: x.count('|'))
    dfs['%s_SW_num_sim_res'%TMD] = dfs['%s_SW_markup_seq'%TMD].dropna().apply(lambda x: x.count(':'))
    # check that the TMD seq in match is not just 100% gaps!
    dfs['%s_in_SW_align_match'%TMD] = dfs['%s_SW_num_ident_res'%TMD].dropna() != 0
    dfs['%s_in_SW_align_match'%TMD].fillna(value=False)

    # the percentage identity of that TMD is defined as the number of identical residues (pipes in markup) divided by the length of the the aligned residues (excluding gaps, based on the length of the shortest TMD, either match or query)
    # note that the nonTMD percentage identity is calculated the same way
    dfs['%s_perc_ident'%TMD] = dfs['%s_SW_num_ident_res'%TMD] / dfs['%s_SW_align_len'%TMD]
    # calculate percentage similar residues
    dfs['%s_perc_sim'%TMD] = dfs['%s_SW_num_sim_res'%TMD] / dfs['%s_SW_align_len'%TMD]
    # add together to obtain the percentage similar + identical residues
    dfs['%s_perc_sim_plus_ident'%TMD] = dfs['%s_perc_ident'%TMD] + dfs['%s_perc_sim'%TMD]
    # add to main dataframe
    dfs['%s_perc_ident'%TMD] = dfs['%s_perc_ident'%TMD]
    # calculate the average number of gaps per residue in the TMD alignment
    # (number of gaps)/(length of sequence excluding gaps)
    dfs['%s_SW_q_gaps_per_q_residue'%TMD] = dfs['%s_SW_query_num_gaps'%TMD].dropna() / len_query_TMD
    return dfs

def calc_AAIMON(acc, TMD, dfs, dfs_filt, set_, df, logging):
    # following the general filters, filter to only analyse sequences with TMD identity above cutoff, and a nonTMD_perc_ident above zero ,to avoid a divide by zero error
    min_identity_of_TMD_initial_filter = set_['cr_min_identity_of_TMD_initial_filter']
    dfs_filt_AAIMON = dfs_filt.loc[dfs['%s_perc_ident'%TMD] >= min_identity_of_TMD_initial_filter]
    # avoid a divide by zero error in the unlikely case that there are no_identical_residues_in_alignment
    dfs_filt_AAIMON = dfs_filt_AAIMON.loc[dfs_filt_AAIMON['nonTMD_perc_ident'] != 0]
    # calculate the Amino Acid Identity : Membranous Over Nonmembranous
    dfs_filt_AAIMON['%s_AAIMON_ratio'%TMD] = dfs_filt_AAIMON['%s_perc_ident'%TMD] / dfs_filt_AAIMON['nonTMD_perc_ident']
    # calculate the Amino Acid Similarity : Membranous Over Nonmembranous (AASMON) (includes similarity + identity based on the matrix used in the SW alignment of SIMAP)
    dfs_filt_AAIMON['%s_AASMON_ratio'%TMD] = dfs_filt_AAIMON['%s_perc_sim_plus_ident'%TMD] / dfs_filt_AAIMON[ 'nonTMD_perc_sim_plus_ident']
    df.loc[acc, '%s_SW_q_gaps_per_q_residue_mean'%TMD] = dfs_filt_AAIMON['%s_SW_q_gaps_per_q_residue'%TMD].dropna().mean()
    logging.info('%s_SW_q_gaps_per_q_residue Average: %0.3e' %(TMD, df.loc[acc, '%s_SW_q_gaps_per_q_residue_mean'%TMD]))

    # add to original dataframe with the list of uniprot sequences
    df.loc[acc, '%s_perc_ident_mean'%TMD] = dfs_filt_AAIMON['%s_perc_ident'%TMD].mean()
    df.loc[acc, '%s_perc_sim_mean'%TMD] = dfs_filt_AAIMON['%s_perc_sim'%TMD].mean()
    df.loc[acc, '%s_perc_sim_plus_ident_mean'%TMD] = dfs_filt_AAIMON['%s_perc_sim_plus_ident'%TMD].mean()
    df.loc[acc, '%s_AAIMON_ratio_mean'%TMD] = float(dfs_filt_AAIMON['%s_AAIMON_ratio'%TMD].mean())
    df.loc[acc, '%s_AAIMON_ratio_std'%TMD] = dfs_filt_AAIMON['%s_AAIMON_ratio'%TMD].std()
    df.loc[acc, '%s_AASMON_ratio_mean'%TMD] = dfs_filt_AAIMON['%s_AASMON_ratio'%TMD].mean()
    df.loc[acc, '%s_AASMON_ratio_std'%TMD] = dfs_filt_AAIMON['%s_AASMON_ratio'%TMD].std()
    logging.info('AAIMON MEAN %s: %0.2f' % (TMD, df.loc[acc, '%s_AAIMON_ratio_mean'%TMD]))
    logging.info('AAISON MEAN %s: %0.2f' % (TMD, df.loc[acc, '%s_AASMON_ratio_mean'%TMD]))

    # add to the dataframe with the SIMAP data for that particular protein
    dfs['%s_AAIMON_ratio'%TMD] = dfs_filt_AAIMON['%s_AAIMON_ratio'%TMD]
    dfs['%s_AASMON_ratio'%TMD] = dfs_filt_AAIMON['%s_AASMON_ratio'%TMD]

    dfs['%s_ratio_length_of_TMD_to_rest_of_alignment'%TMD] = dfs['%s_SW_query_seq'%TMD].str.len() / dfs['FASTA_overlap']
    dfs['%s_ratio_length_of_query_TMD_to_rest_of_match_protein'%TMD] = dfs['%s_SW_query_seq'%TMD].str.len() / dfs[ 'len_full_match_seq']

    df.loc[acc, '%s_ratio_length_of_TMD_to_rest_of_alignment_mean'%TMD] = float('%0.2f' % dfs['%s_ratio_length_of_TMD_to_rest_of_alignment'%TMD].dropna().mean())
    df.loc[acc, '%s_ratio_length_of_query_TMD_to_rest_of_match_protein_mean'%TMD] = float('%0.2f' % dfs['%s_ratio_length_of_query_TMD_to_rest_of_match_protein'%TMD].dropna().mean())

    return df, dfs

def save_hist_and_fastA(acc, dfs, dfs_filt, set_, df, list_of_TMDs, logging):
    # use linspace to get a fixid number of points between tha min and the max for the histogram
    # set up evenly distributed bins between the chosen min and max
    # if possible, 1.0 should be in the centre of a bin, to catch cases where a lot of homologues have a ratio that approximates 1
    linspace_binlist = np.linspace(set_["1p_smallest_bin"],
                                   set_["1p_largest_bin"],
                                   set_["1p_number_of_bins"])
    # add 30 as the last bin, to make sure 100% of the data is added to the histogram, including major outliers
    binlist = np.append(linspace_binlist,
                        set_["1p_final_highest_bin"])

    # se default font size for text in the plot
    fontsize = 4
    # use a dictionary to organise the saving of multiple plots in multiple figures, with a certain number of plots per figure
    n_plots_per_fig = 4
    nrows_in_each_fig = 2
    ncols_in_each_fig = 2
    dict_organising_subplots = utils.create_dict_organising_subplots(
        n_plots_per_fig=n_plots_per_fig,
        n_rows=nrows_in_each_fig,
        n_cols=ncols_in_each_fig
    )
    with tarfile.open(df.loc[acc, 'output_tarfile_path'], mode='w:gz') as tar_out:
        # calculate ratio of Amino acid Identity Membranous Over Nonmembranous  (AAIMON ratio)
        for TMD in list_of_TMDs:
            len_query_TMD = len(df.loc[acc, '%s_seq'%TMD])
            # following the general filters, filter to only analyse sequences with TMD identity above cutoff, and a nonTMD_perc_ident above zero ,to avoid a divide by zero error
            min_identity_of_TMD_initial_filter = set_['cr_min_identity_of_TMD_initial_filter']
            #print(dfs_filt['%s_AAIMON_ratio'%TMD])
            #RESET SO THAT IT TAKES IT FROM DFS AGAIN
            #dfs_filt_AAIMON = dfs_filt.loc[dfs['%s_perc_ident'%TMD] >= min_identity_of_TMD_initial_filter]
            # find the TMD number (starting from 1)
            TMD_Nr = list_of_TMDs.index(TMD) + 1
            # use the dictionary to obtain the figure number, plot number in figure, plot indices, etc
            newfig, savefig, fig_nr, plot_nr_in_fig, row_nr, col_nr = dict_organising_subplots[TMD_Nr]
            # if the TMD is the last one, the figure should be saved
            if TMD_Nr == len(list_of_TMDs):
                savefig = True
            # if a new figure should be created (either because the orig is full, or the last TMD is analysed)
            if newfig:
                # create a new figure
                fig, axarr = plt.subplots(nrows=nrows_in_each_fig,
                                          ncols=ncols_in_each_fig)  # sharex=True
            # create numpy array of membranous over nonmembranous conservation ratios (identity)
            hist_data_I = np.array(dfs['%s_AAIMON_ratio'%TMD].dropna())
            # use numpy to create a histogram
            freq_counts, bin_array = np.histogram(hist_data_I, bins=binlist)
            # assuming all of the bins are exactly the same size, make the width of the column equal to 70% of each bin
            col_width = float('%0.3f' % (0.7 * (bin_array[1] - bin_array[0])))
            # when align='center', the central point of the bar in the x-axis is simply the middle of the bins ((bin_0-bin_1)/2, etc)
            centre_of_bar_in_x_axis = (bin_array[:-2] + bin_array[1:-1]) / 2
            # add the final bin, which is physically located just after the last regular bin but represents all higher values
            centre_of_bar_in_x_axis = np.append(centre_of_bar_in_x_axis,
                                                centre_of_bar_in_x_axis[-1] +
                                                centre_of_bar_in_x_axis[0])
            barcontainer_I = axarr[row_nr, col_nr].bar(left=centre_of_bar_in_x_axis,
                                                       height=freq_counts, align='center',
                                                       width=col_width, color="#0489B1",
                                                       alpha=0.5)  # edgecolor='black',
            # create numpy array of membranous over nonmembranous conservation ratios (identity + similarity)
            hist_data_S = np.array(dfs['%s_AASMON_ratio'%TMD].dropna())
            # use numpy to create a histogram
            freq_counts, bin_array = np.histogram(hist_data_S, bins=binlist)
            # create a line graph rather than a bar graph for the AAISON (ident + similarity)
            linecontainer_S = axarr[row_nr, col_nr].plot(centre_of_bar_in_x_axis, freq_counts,
                                                         color="#0101DF", alpha=0.5)
            # other colours that are compatible with colourblind readers: #8A084B Dark red, #B45F04 deep orange, reddish purple #4B088A
            # http://html-color-codes.info/
            # label the x-axis for each plot, based on the TMD
            axarr[row_nr, col_nr].set_xlabel('%s conservation ratio (membranous over nonmembranous)'%TMD,
                                             fontsize=fontsize)
            if savefig:
                # take x-axis min from settings
                xlim_min = set_["1p_smallest_bin"]
                # take x-axis max from settings
                xlim_max = set_["1p_largest_bin"]
                # apply the following formatting changes to all plots in the figure
                for ax in axarr.flat:
                    # set x-axis min
                    ax.set_xlim(xlim_min, xlim_max)
                    # set x-axis ticks
                    # use the slide selection to select every second item in the list as an xtick(axis label)
                    ax.set_xticks([float('%0.1f' % c) for c in centre_of_bar_in_x_axis[::3]])
                    ax.set_ylabel('freq', rotation='vertical', fontsize=fontsize)
                    # change axis font size
                    ax.tick_params(labelsize=fontsize)
                    # create legend?#http://stackoverflow.com/questions/9834452/how-do-i-make-a-single-legend-for-many-subplots-with-matplotlib
                    ax.legend(['AASMON (identity + similarity)', 'AAIMON (identity)'],
                              loc='upper right', fontsize=fontsize)
                    # add background grid
                    ax.grid(True, color='0.75', alpha=0.5)
                # automatically tighten the layout of plots in the figure
                fig.tight_layout()
                # save files
                fig.savefig(df.loc[acc, 'AAIMON_hist_BASENAMEPATH'] + '_%01d.png' % fig_nr,
                            format='png', dpi=200)
                fig.savefig(df.loc[acc, 'AAIMON_hist_BASENAMEPATH'] + '_%01d.pdf' % fig_nr,
                            format='pdf')
                # close figure
                plt.close('all')
                # add to tarfile
                tar_out.add(df.loc[acc, 'AAIMON_hist_BASENAMEPATH'] + '_%01d.png' % fig_nr,
                            arcname=df.loc[acc, 'AAIMON_hist_BASENAME'] + '_%01d.png' % fig_nr)
                tar_out.add(df.loc[acc, 'AAIMON_hist_BASENAMEPATH'] + '_%01d.pdf' % fig_nr,
                            arcname=df.loc[acc, 'AAIMON_hist_BASENAME'] + '_%01d.pdf' % fig_nr)
                # delete files
                # os.remove(df.loc[acc,'AAIMON_hist_BASENAMEPATH'] + '_%01d.png' % fig_nr, format='png', dpi=200)
                os.remove(df.loc[acc, 'AAIMON_hist_BASENAMEPATH'] + '_%01d.pdf' % fig_nr)

            dfs = korbinian.simap.filter_and_save_fastA(df, dfs, acc, TMD, set_, tar_out, logging)

        # remove columns to make output csv smaller
        if set_['drop_columns_to_reduce_csv_filesize']:
            list_cols_to_drop = ['match_alignment_sequence', 'query_alignment_sequence', 'alignment_markup',
                                 'nonTMD_seq_query', 'nonTMD_markup']
            for col in list_cols_to_drop:
                if col in dfs.columns:
                    dfs.drop(col, axis=1, inplace=True)
        dfs.to_csv(df.loc[acc, 'SIMAP_csv_analysed_path'], sep=",", quoting=csv.QUOTE_NONNUMERIC)
        tar_out.add(df.loc[acc, 'SIMAP_csv_analysed_path'], arcname=df.loc[acc, 'SIMAP_csv_analysed'])
        # delete original uncompressed file
        os.remove(df.loc[acc, 'SIMAP_csv_analysed_path'])
        df.loc[acc, 'num_hits_with_SW_align_node'] = dfs['hit_contains_SW_node'].value_counts()[True]
        logging.info('num_hits_with_SW_align_node: %s' % df.loc[acc, 'num_hits_with_SW_align_node'])

def conduct_ssr_ratio_calculations():
    '''  _________________________________________SSR ratio calculations____________________________________________________
    '''

    #                    def calc_score_ss_qTMD(dfs):
    #                        score = sum(utils.score_pairwise(seq1=dfs['%s_SW_query_seq'%TMD], seq2=dfs['%s_SW_query_seq'%TMD],
    #                                                     matrix=aa_sub_matrix,
    #                                                     gap_open_penalty=gap_open_penalty,
    #                                                     gap_extension_penalty=gap_extension_penalty))
    #                        return(score)
    #
    #                    def calc_score_ss_mTMD(dfs):
    #                        score = sum(utils.score_pairwise(seq1=dfs['%s_SW_match_seq'%TMD], seq2=dfs['%s_SW_match_seq'%TMD],
    #                                                     matrix=aa_sub_matrix,
    #                                                     gap_open_penalty=gap_open_penalty,
    #                                                     gap_extension_penalty=gap_extension_penalty))
    #                        return(score)
    #
    #                    def calc_score_qTMD_mTMD(dfs):
    #                        score = sum(utils.score_pairwise(seq1=dfs['%s_SW_query_seq'%TMD], seq2=dfs['%s_SW_match_seq'%TMD],
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
    #                    for gap_open_penalty in range(set_["gap_open_penalty_max"], set_["gap_open_penalty_increment"]):
    #                        #print(gap_open_penalty)
    #                        #for simplicity, give the gap open and gap extend the same value
    #                        gap_extension_penalty = gap_open_penalty
    #                        for matrix_name in list_of_aa_sub_matrices:
    #                            #so long as the matrix is imported into python, eval will convert the matrix name to an object
    #                            aa_sub_matrix = ast.literal_eval(matrix_name)
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
    #                                dfs_nonTMD = dfs_nonTMD.loc[dfs['%s_SW_query_seq'%TMD].notnull()]
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