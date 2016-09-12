import korbinian
import korbinian.mtutils as utils
import pandas as pd

def slice_nonTMD_seqs(dfs, df_nonTMD_sliced, list_of_TMDs):
    ########################################################################################
    #                                                                                      #
    #                           nonTMD calculations    [AAIMON]                            #
    #                                                                                      #
    ########################################################################################
    if dfs.shape[0] > 0:

        # check if all tmds are in SW alignment
        # create a list of columns to reindex the DataFrame
        list_columns_TMD_in_SW_alignment = []
        for TMD in list_of_TMDs:
            # TMD found by regex
            list_columns_TMD_in_SW_alignment.append('%s_in_SW_alignment' % TMD)
            # TMD matching useful sequence
            #list_columns_TMD_in_SW_alignment.append('%s_in_SW_align_match' % TMD)

        # create a slice of the filtered dataframe that only contains the relevant columns (N.B. copy=False, this will provide a view, not a copy)
        df2 = df_nonTMD_sliced.reindex(index=df_nonTMD_sliced.index, columns=list_columns_TMD_in_SW_alignment, copy=False)
        # create a new column in the original dataframe that shows that ALL TMDs have been found in the SW alignment
        df_nonTMD_sliced['all_tmds_in_SW_alignment'] = df2.dropna().all(axis=1)
        # filter to contain only hits with all tmds

        ########################################################################################
        #                                                                                      #
        #                 start processing nonTMD region [AAIMON]                              #
        #                                                                                      #
        ########################################################################################
        # create a copy of the original df_cr dataframe containing only hits where all tmds are found in the match
        #df_cr_nonTMD = df_cr.loc[df_cr['all_tmds_in_SW_alignment'].notnull()].query('all_tmds_in_SW_alignment == True')

        # drop any homologues where not all TMDs were fould in the match
        df_nonTMD_sliced.query('all_tmds_in_SW_alignment == True', inplace=True)
        # filter to contain only hits where the index for the TMD is present
        first_TMD_start_index = '%s_start_in_SW_alignment' % list_of_TMDs[0]

        # FILTER REMOVED! HOPEFULLY STILL WORKS :)
        #df_cr_nonTMD = df_cr_nonTMD.loc[df_cr_nonTMD[first_TMD_start_index].notnull()]
        df_nonTMD_sliced['nonTMD_index_tuple_first'] = df_nonTMD_sliced[first_TMD_start_index].apply(lambda x: (0, int(x)))

        # create start and stop indices for all sections between tmds
        # the start of the last nonTMD section will be the end of the last TMD
        df_nonTMD_sliced['nonTMD_index_tuple_last0'] = df_nonTMD_sliced['%s_end_in_SW_alignment' % list_of_TMDs[-1]].dropna().astype('int32')

        # the end of the last nonTMD section will be the end of the full alignment sequence
        df_nonTMD_sliced['nonTMD_index_tuple_last1'] = df_nonTMD_sliced['len_query_align_seq'].dropna().astype('int32')
        # join to make a tuple
        # df_cr_nonTMD['nonTMD_index_tuple_last'] = df_cr_nonTMD[['nonTMD_index_tuple_last0', 'nonTMD_index_tuple_last1']].apply(tuple, axis=1)

        # create the index tuple
        df_nonTMD_sliced['nonTMD_index_tuple_last'] = df_nonTMD_sliced.apply(utils.create_indextuple_nonTMD_last, axis=1)


        utils.aaa(df_nonTMD_sliced)


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
            df_nonTMD_sliced[index_columns] = df_nonTMD_sliced[index_columns].astype('int64')
            # create a tuple containing the indices for the nonTMD sequence regions in between each TMD (middle indices)
            df_nonTMD_sliced['nonTMD_index_%s' % TMD] = tuple(zip(df_nonTMD_sliced['%s_end_in_SW_alignment' % TMD],df_nonTMD_sliced['%s_start_in_SW_alignment' % next_TMD]))

        # now join all the indices together to make one tuple of tuples for the non-TMD region
        # df_cr_nonTMD = df_cr.query('all_tmds_in_SW_alignment == True')
            df_nonTMD_sliced['nested_tuple_indices_all_nonTMD_regions'] = df_nonTMD_sliced[['nonTMD_index_tuple_last']].apply(tuple, axis=1)

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
        df_cr_tuples = df_nonTMD_sliced.reindex(index=df_nonTMD_sliced.index, columns=list_of_nonTMD_index_columns, copy=False)
        # now for convenience, these tuples can be combined together to form one column, with a tuple of tuples
        # first convert all values in each row to a list, excluding the index column
        list_tuple_indices_all_nonTMD_regions = list(df_cr_tuples.itertuples(index=False))
        # convert to a series, and reindex with the original index from the dataframe
        tuples_series = pd.Series(list_tuple_indices_all_nonTMD_regions, index=df_nonTMD_sliced.index)
        # for some reason, the tuples are a pandas object "Pandas(nonTMD_index_tuple_first=(0, 592), nonTMD_index_tuple_last=(615, 618))"
        # convert to simple tuples
        df_nonTMD_sliced['nested_tuple_indices_all_nonTMD_regions'] = tuples_series.apply(lambda x: tuple(x))
        # change to a string, in case this solves the weird effect with only the last tuple shown
        df_nonTMD_sliced['nested_tuple_indices_all_nonTMD_regions'] = df_nonTMD_sliced['nested_tuple_indices_all_nonTMD_regions'].astype(str)
        # you can test that the original index is maintained as follows:
        
        # filter dfs to only contain the rows of interest as contained by df_TMD_indice (homologues that contain all TMDs)
        dfs = dfs.loc[df_nonTMD_sliced.index,:]
        
        # add the series as a new column in the original dataframe. Missing data (when not all TMDs found) will be filled using np.nan
        dfs['nested_tuple_indices_all_nonTMD_regions'] = df_nonTMD_sliced['nested_tuple_indices_all_nonTMD_regions']

        # # filter to remove incomplete sequences
        # df_cr_nonTMD = df_cr_nonTMD.query('all_tmds_in_SW_alignment == True')
        # define the string for slicing as a numpy array
        # use the numpy vectorize function, which effectively applies the function in a for loop (not optimized for speed)
        # df_cr_nonTMD['nonTMD_seq_query'] = np.vectorize(utils.slice_with_nested_tuple)(np.array(df_cr_nonTMD['query_align_seq']),np.array(df_cr_nonTMD['nested_tuple_indices_all_nonTMD_regions']))
        # df_cr_nonTMD['nonTMD_markup'] = np.vectorize(utils.slice_with_nested_tuple)(np.array(df_cr_nonTMD['align_markup_seq']),np.array(df_cr_nonTMD['nested_tuple_indices_all_nonTMD_regions']))
        # df_cr_nonTMD['nonTMD_seq_match'] = np.vectorize(utils.slice_with_nested_tuple)(np.array(df_cr_nonTMD['match_align_seq']),np.array(df_cr_nonTMD['nested_tuple_indices_all_nonTMD_regions']))

        ########################################################################################
        #                                                                                      #
        #              slice out the nonTMD region for each homologue [AAIMON]                 #
        #                                                                                      #
        ########################################################################################
        # due to problems with np.vectorize and the pandas methods, slice the sequences one at a time with a simple 'for loop'
        # for each hit, perform the slice
        for hit in dfs.index:
            df_nonTMD_sliced.loc[hit, 'nonTMD_seq_query'] = utils.slice_with_nested_tuple(dfs.loc[hit, 'query_align_seq'],dfs.loc[hit, 'nested_tuple_indices_all_nonTMD_regions'])
            df_nonTMD_sliced.loc[hit, 'nonTMD_markup'] = utils.slice_with_nested_tuple(dfs.loc[hit, 'align_markup_seq'],dfs.loc[hit, 'nested_tuple_indices_all_nonTMD_regions'])
            df_nonTMD_sliced.loc[hit, 'nonTMD_seq_match'] = utils.slice_with_nested_tuple(dfs.loc[hit, 'match_align_seq'],dfs.loc[hit, 'nested_tuple_indices_all_nonTMD_regions'])
        # # transfer to original dataframe (index should still match the original, partial seqs will be filled with np.nan)
        # df_cr['nonTMD_seq_query'] = df_cr_nonTMD['nonTMD_seq_query']
        # df_cr['nonTMD_markup'] = df_cr_nonTMD['nonTMD_markup']
        # df_cr['nonTMD_seq_match'] = df_cr_nonTMD['nonTMD_seq_match']
        
        return df_nonTMD_sliced


def calc_nonTMD_perc_ident_and_gaps(acc, df_nonTMD, mean_ser, logging):


        ########################################################################################
        #                                                                                      #
        #         calculate nonTMD len, percent identity, gaps, etc [AAIMON]                   #
        #                                                                                      #
        ########################################################################################
        # calculate identical residues in the nonTMD region (simply count the pipes '|' in the markup sequence)
        df_nonTMD['nonTMD_num_ident_res'] = df_nonTMD['nonTMD_markup'].str.count('|')
        # calculate similar residues in the nonTMD region (simply count the colons ':' in the markup sequence)
        df_nonTMD['nonTMD_num_sim_res'] = df_nonTMD['nonTMD_markup'].str.count(':')
        # add the identical and similar residues together to get the total number of similar + identical residues
        df_nonTMD['nonTMD_num_sim_plus_ident_res'] = df_nonTMD['nonTMD_num_ident_res'] + df_nonTMD['nonTMD_num_sim_res']

        # count the gaps in the nonTMD sequence of the query
        df_nonTMD['nonTMD_q_num_gaps'] = df_nonTMD['nonTMD_seq_query'].str.count('-')
        # count the gaps in the nonTMD sequence of the match
        df_nonTMD['nonTMD_m_num_gaps'] = df_nonTMD['nonTMD_seq_match'].str.count('-')
        # calculate the length of the nonTMD sequences, which may include gaps
        df_nonTMD['len_nonTMD_seq_query'] = df_nonTMD['nonTMD_seq_query'].str.len()
        df_nonTMD['len_nonTMD_seq_match'] = df_nonTMD['nonTMD_seq_match'].str.len()
        # calculate the number aligned sequences, excluding gaps (length of query, or length of match, whichever is shorter)
        df_nonTMD['len_nonTMD_align'] = df_nonTMD[['len_nonTMD_seq_query', 'len_nonTMD_seq_match']].dropna(how='all').min(axis=1)

        # calculate the length of the nonTMD sequence excluding gaps
        df_nonTMD['len_nonTMD_q_excl_gaps'] = df_nonTMD['len_nonTMD_seq_query'] - df_nonTMD['nonTMD_q_num_gaps']
        df_nonTMD['len_nonTMD_m_excl_gaps'] = df_nonTMD['len_nonTMD_seq_match'] - df_nonTMD['nonTMD_m_num_gaps']
        # calculate the lenth of the alignment by finding which seq excl gaps is smaller
        df_nonTMD['len_nonTMD_align'] = df_nonTMD[['len_nonTMD_q_excl_gaps', 'len_nonTMD_m_excl_gaps']].min(axis=1)

        # calculate the percentage identity of the nonTMD region (number of identical residues divided by the length excluding gaps)
        # used for the Amino Acid Identity : Membranous over Nonmembranous (AAIMON ratio)
        # note that the length = length of the aligned residues excluding gaps
        df_nonTMD['nonTMD_perc_ident'] = df_nonTMD['nonTMD_num_ident_res'] / df_nonTMD['len_nonTMD_align']
        df_nonTMD['nonTMD_perc_sim'] = df_nonTMD['nonTMD_num_sim_res'] / df_nonTMD['len_nonTMD_align']
        df_nonTMD['nonTMD_perc_sim_plus_ident'] = df_nonTMD['nonTMD_num_sim_plus_ident_res'] / df_nonTMD['len_nonTMD_align']
        # calculate the average number of gaps per residue in the nonTMD alignment
        # filter to analyse only sequences that are valid (length > 0)
        df_nonTMD_filt_gaps = df_nonTMD.loc[df_nonTMD['len_nonTMD_q_excl_gaps'] != 0]
        # calculate number of gaps in query AND match
        df_nonTMD_filt_gaps['nonTMD_qm_num_gaps'] = df_nonTMD_filt_gaps['nonTMD_q_num_gaps'] + df_nonTMD_filt_gaps['nonTMD_m_num_gaps']
        # add to simap dataframe
        df_nonTMD['nonTMD_qm_num_gaps'] = df_nonTMD_filt_gaps['nonTMD_qm_num_gaps']
        # gaps per query residue for both query and match = ((gaps in query + gaps in match)/2))/length of query excluding gaps
        df_nonTMD['nonTMD_qm_gaps_per_q_residue'] = df_nonTMD_filt_gaps['nonTMD_qm_num_gaps'] / 2 / df_nonTMD_filt_gaps['len_nonTMD_q_excl_gaps']

        # SSR ratio calculations take a long time and show no difference to AAIMON
        SSR_ratio_calculations = False
        if SSR_ratio_calculations:
            korbinian.cons_ratio.ssr.conduct_ssr_ratio_calculations()

        ########################################################################################
        #                                                                                      #
        #   calculate average nonTMD len, perc ident, etc and add to list of proteins [AAIMON] #
        #                                                                                      #
        ########################################################################################

        ''' calculate values associated with the nonTMD region and add to main dataframe
        '''
        # add the average values regarding the nonTMD region to the original file/dataframe with each protein
        mean_ser['len_nonTMD_seq_match_mean'] = float('%0.2f' % df_nonTMD['len_nonTMD_seq_match'].dropna().mean())
        mean_ser['nonTMD_perc_ident_mean'] = float('%0.3f' % df_nonTMD['nonTMD_perc_ident'].dropna().mean())
        mean_ser['nonTMD_perc_sim_mean'] = float('%0.3f' % df_nonTMD['nonTMD_perc_sim'].dropna().mean())
        mean_ser['nonTMD_perc_sim_plus_ident_mean'] = float('%0.3f' % df_nonTMD['nonTMD_perc_sim_plus_ident'].dropna().mean())
        mean_ser['len_nonTMD_align_mean'] = float('%0.2f' % df_nonTMD['len_nonTMD_align'].dropna().mean())
        mean_ser['nonTMD_qm_gaps_per_q_residue_mean'] = float('%0.2f' % df_nonTMD['nonTMD_qm_gaps_per_q_residue'].dropna().mean())
        #logging.info('nonTMD_qm_gaps_per_q_residue : %0.5f' % mean_ser['nonTMD_qm_gaps_per_q_residue_mean'])

        return mean_ser, df_nonTMD