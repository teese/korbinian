import korbinian
import korbinian.mtutils as utils
import pandas as pd

def calc_nonTMD_perc_ident_and_gaps(acc, dfs, set_, df, list_of_TMDs, logging):

    ########################################################################################
    #                                                                                      #
    #                           nonTMD calculations    [AAIMON]                            #
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
            list_columns_TMD_in_SW_alignment.append('%s_in_SW_alignment' % TMD)
            # TMD matching useful sequence
            list_columns_TMD_in_SW_alignment.append('%s_in_SW_align_match' % TMD)

        # create a slice of the filtered dataframe that only contains the relevant columns (N.B. copy=False, this will provide a view, not a copy)
        df2 = dfs_filt.reindex(index=dfs.index, columns=list_columns_TMD_in_SW_alignment, copy=False)
        # create a new column in the original dataframe that shows that ALL TMDs have been found in the SW alignment
        dfs['all_tmds_in_SW_alignment'] = df2.dropna().all(axis=1)
        # filter to contain only hits with all tmds

        ########################################################################################
        #                                                                                      #
        #                 start processing nonTMD region [AAIMON]                              #
        #                                                                                      #
        ########################################################################################
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
        #            create the indices for the nonTMD region [AAIMON]                         #
        #                                                                                      #
        ########################################################################################

        # for each TMD EXCEPT the last, which ends at the sequence end, create the indices for the nonTMD region (after the TMD)
        for TM_Nr in range(len(list_of_TMDs) - 1):
            # the TMD is the equivalent item in the list
            TMD = list_of_TMDs[TM_Nr]
            # the next TMD, which contains the end index, is the next item in the list
            next_TMD = list_of_TMDs[TM_Nr + 1]
            # select only the columns in the dataframe that are of interest, and change the data type to integer
            index_columns = ['%s_end_in_SW_alignment' % TMD, '%s_start_in_SW_alignment' % next_TMD]
            dfs_nonTMD[index_columns] = dfs_nonTMD[index_columns].astype('int64')
            # create a tuple containing the indices for the nonTMD sequence regions in between each TMD (middle indices)
            dfs_nonTMD['nonTMD_index_%s' % TMD] = tuple(zip(dfs_nonTMD['%s_end_in_SW_alignment' % TMD],
                                                            dfs_nonTMD['%s_start_in_SW_alignment' % next_TMD]))

        # now join all the indices together to make one tuple of tuples for the non-TMD region
        # dfs_nonTMD = dfs.query('all_tmds_in_SW_alignment == True')
        dfs_nonTMD['nested_tuple_indices_all_nonTMD_regions'] = dfs_nonTMD[['nonTMD_index_tuple_last']].apply(tuple, axis=1)

        # create a view of the dataframe that contains only the desired columns
        # first create a list of the desired columns
        # start with the first tuple, from 0 to the start of the first TMD
        list_of_nonTMD_index_columns = ['nonTMD_index_tuple_first']
        # create a nonTMD region for each of the TMDs (except the last one)
        list_from_TMs = ['nonTMD_index_%s' % TMD2 for TMD2 in list_of_TMDs[:-1]]
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
        dfs_nonTMD['nested_tuple_indices_all_nonTMD_regions'] = dfs_nonTMD[
            'nested_tuple_indices_all_nonTMD_regions'].astype(str)
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
        #              slice out the nonTMD region for each homologue [AAIMON]                 #
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
        #         calculate nonTMD len, percent identity, gaps, etc [AAIMON]                   #
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
            korbinian.cons_ratio.ssr.conduct_ssr_ratio_calculations()

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
        #   calculate average nonTMD len, perc ident, etc and add to list of proteins [AAIMON] #
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

        return df, dfs