import korbinian.mtutils as utils
import pandas as pd

def slice_TMD_homol_and_count_gaps(acc, TMD, df, dfs, set_, logging, n_TMDs_w_homol):
    ########################################################################################
    #                                                                                      #
    #                Slice out TMD regions [fasta and AAIMON]                              #
    #                                                                                      #
    ########################################################################################
    '''slice the TMD regions out of the alignment markup and match sequences
    the method first finds the indices of the TMD region in the query, and uses these indices to slice
    the filters are not applied yet, so the sequences can be viewed in the csv file for analysis
    '''
    # create regex string for each TMD
    query_TMD_sequence = df.loc[acc, '%s_seq' % TMD]
    # create TMD regex search string (e.g. L.*L.*M.*L.*L.*M.*L.*L.*M.*L.*L.*M.*L.*L.*M.*)
    TMD_regex_ss = utils.create_regex_string(query_TMD_sequence)
    # select to only include data where the XML contained a SW node, and then apply function for a regex search
    dfs_with_SW_node = dfs.query('hit_contains_SW_node == True')
    # create series of query_alignment_seqs
    query_seqs_ser = dfs_with_SW_node['query_align_seq']

    # obtain the bool, start, end of TMD seqs in the match sequences. Add to the new TMD-specific dataframe.
    dfs['%s_start_end_list_in_SW_alignment' % TMD] = query_seqs_ser.apply(utils.get_start_and_end_of_TMD_in_query,
                                                                             args=(TMD_regex_ss,))
    '''the output of the regex search is a list with three components.
    1) a match (e.g. True),
    2) the start of the match (e.g. 124),
    3) the end of the match (e.g. 144)
    '''
    # for convenience, separate the components into separate columns
    columns_from_regex_output = ['%s_in_SW_alignment' % TMD, '%s_start_in_SW_alignment' % TMD, '%s_end_in_SW_alignment' % TMD]
    # n = the index (1,2,3) in the tuple
    # col = column item in the list (e.g. 'TM09_in_SW_alignment')
    for n, col in enumerate(columns_from_regex_output):
        # add a new column which is named TM01_start, etc, and insert the appropriate integer (start or stop) or bool from the tuple
        dfs[col] = dfs['%s_start_end_list_in_SW_alignment' % TMD].dropna().apply(lambda x: x[n])
    ########################################################################################
    #                                                                                      #
    #          Apply function to slice homologues and count gaps [fasta and AAIMON]        #
    #                                                                                      #
    ########################################################################################

    # in some cases, there is no data to obtain as the hit_contains_SW_node = False for too many sequences, giving no start_in_SW_alignment
    number_of_rows_containing_data = dfs[dfs['%s_start_end_list_in_SW_alignment' % TMD].notnull()].shape[0]
    # drop the original listlike from the regex search
    #NOT NECESSARY. DFS WILL NOT BE RETURNED!
    #df_TMD.drop('%s_start_end_list_in_SW_alignment' % TMD, inplace=True, axis=1)
    if number_of_rows_containing_data == 0:
        logging.info('%s does not have any valid homologues for %s. '
                     'Re-downloading simap homologue XML may be necessary.' % (df.loc[acc,"protein_name"], TMD))
    if number_of_rows_containing_data != 0:
        n_TMDs_w_homol += 1
        len_query_TMD = len(df.loc[acc, '%s_seq' % TMD])
        # apply the slicing function to the homologues
        # df_TMD = korbinian.cons_ratio.slice_TMD_homol_and_count_gaps(TMD, len_query_TMD, df_TMD,
        #                                                               set_["cr_max_n_gaps_in_query_TMD"],
        #                                                               set_["cr_max_n_gaps_in_match_TMD"])
        # use small throwaway functions to slice each TMD from the query, markup and match sequences
        # notnull() removes missing data
        # explanation: df_TMD['new_column_with_selected_seq'] = df_TMD[df_TMD['only_rows_containing_data]].apply(utils.slice_function_that_specifies_columns_with_start_and_stop)
        dfs_sel = dfs[dfs['%s_start_in_SW_alignment'%TMD].notnull()]
        # df_TMD['%s_SW_query_seq'%TMD] = dfs_sel.apply(utils.slice_SW_query_TMD_seq,args=(TMD,), axis=1)
        # df_TMD['%s_SW_markup_seq'%TMD] = dfs_sel.apply(utils.slice_SW_markup_TMD,args=(TMD,), axis=1)
        # df_TMD['%s_SW_match_seq'%TMD] = dfs_sel.apply(utils.slice_SW_match_TMD_seq,args=(TMD,), axis=1)

        #create a new dataframe to hold the sliced homol sequences for that TMD, and number of gaps, etc
        df_TMD = dfs_sel[["organism", "description"]].copy()

        # transfer the columns with indices across to the df_TMD
        cols = ['%s_in_SW_alignment' % TMD, '%s_start_in_SW_alignment' % TMD, '%s_end_in_SW_alignment' % TMD]
        for col in cols:
            df_TMD[col] = dfs_sel[col]

        # slice TMDs out of dfs_sel, and save them in the new df_TMD
        df_TMD['%s_SW_query_seq' % TMD] = dfs_sel.apply(utils.slice_SW_query_TMD_seq, args=(TMD,), axis=1)
        df_TMD['%s_SW_markup_seq' % TMD] = dfs_sel.apply(utils.slice_SW_markup_TMD, args=(TMD,), axis=1)
        df_TMD['%s_SW_match_seq' % TMD] = dfs_sel.apply(utils.slice_SW_match_TMD_seq, args=(TMD,), axis=1)

        ########################################################################################
        #                                                                                      #
        #                         FASTA: slice TMD plus surrounding sequence                   #
        #                                                                                      #
        #       NOTE: probably not used in AAIMON, but if it is not here,                      #
        # the full match_alignment_seq needs to be saved for the fasta algorithms somewhere.   #
        #    Putting it here allows these large, memory-chewing columns to be dropped.         #
        #                                                                                      #
        ########################################################################################
        # redefine the number of amino acids before and after the TMD to be inserted into the FastA files
        n_aa_before_tmd = set_["n_aa_before_tmd"]
        n_aa_after_tmd = set_["n_aa_after_tmd"]

        # define the start of theTMD + surrounding sequence
        dfs['%s_start_in_SW_alignment_plus_surr' % TMD] = dfs['%s_start_in_SW_alignment' % TMD] - n_aa_before_tmd
        # replace negative values with zero
        dfs.loc[dfs['%s_start_in_SW_alignment_plus_surr' % TMD] < 0, '%s_start_in_SW_alignment_plus_surr' % TMD] = 0
        # define the end of the TMD + surrounding sequence. In python slicing, this end can be longer than the sequence.
        dfs['%s_end_in_SW_alignment_plus_surr' % TMD] = dfs['%s_end_in_SW_alignment' % TMD] + n_aa_after_tmd

        # select sequences that seem to have a start
        dfs = dfs.loc[dfs['%s_start_in_SW_alignment_plus_surr' % TMD].notnull()]
        # slice out the match seq + the surrounding sequence
        df_TMD['%s_SW_match_seq_plus_surr' % TMD] = dfs.apply(utils.slice_SW_match_TMD_seq_plus_surr, args=(TMD,),axis=1)
        #and the same for the TMD + surrounding sequence, useful to examine the TMD interface
        # NOT DEEMED NECESSARY. WHY WOULD YOU NEED TO SLICE QUERY OR MARKUP + SURROUNDING?
        # df_TMD['%s_SW_query_seq_plus_surr'%TMD] = dfs.apply(utils.slice_SW_query_TMD_seq_plus_surr, args=(TMD,), axis=1)
        # df_TMD['%s_SW_markup_seq_plus_surr'%TMD] = dfs.apply(utils.slice_SW_markup_TMD_plus_surr, args=(TMD,), axis=1)

        # count the number of gaps in the query and match sequences
        df_TMD['%s_SW_query_num_gaps'%TMD] = df_TMD['%s_SW_query_seq'%TMD].str.count("-")
        df_TMD['%s_SW_match_num_gaps'%TMD] = df_TMD['%s_SW_match_seq'%TMD].str.count("-")
        # calculate the length of the match TMD seq excluding gaps
        df_TMD['%s_SW_m_seq_len'%TMD] = df_TMD['%s_SW_match_seq'%TMD].str.len()
        # for the alignment length, take the smallest value from the length of query or match
        # this will exclude gaps from the length in the following calculations, preventing false "low conservation" where the query TMD is much longer than the match TMD)
        # note that for most calculations this is somewhat redundant, because the max number of acceptable gaps in sequence is probable ~2
        # use the mask function (faster than a lambda function) to replace any lengths larger than the query, with the query
        df_TMD['%s_SW_align_len' % TMD] = df_TMD['%s_SW_m_seq_len' % TMD].mask(df_TMD['%s_SW_m_seq_len' % TMD] > len_query_TMD, len_query_TMD)
        # create a boolean column that allows filtering by the accepted number of gaps, according to the settings file
        df_TMD['%s_cr_SW_query_acceptable_n_gaps'%TMD] = df_TMD['%s_SW_query_num_gaps'%TMD] <= set_["cr_max_n_gaps_in_query_TMD"]
        df_TMD['%s_cr_SW_match_acceptable_n_gaps'%TMD] = df_TMD['%s_SW_match_num_gaps'%TMD] <= set_["cr_max_n_gaps_in_match_TMD"]
        # count identical residues between query and match TMDs by counting the number of pipes in the markup string
        df_TMD['%s_SW_num_ident_res'%TMD] = df_TMD['%s_SW_markup_seq'%TMD].str.count('|')
        df_TMD['%s_SW_num_sim_res'%TMD] = df_TMD['%s_SW_markup_seq'%TMD].str.count(':')
        # check that the TMD seq in match is not just 100% gaps!
        df_TMD['%s_in_SW_align_match'%TMD] = df_TMD['%s_SW_num_ident_res'%TMD].dropna() != 0
        df_TMD['%s_in_SW_align_match'%TMD].fillna(value=False)

        # the percentage identity of that TMD is defined as the number of identical residues (pipes in markup) divided by the length of the the aligned residues (excluding gaps, based on the length of the shortest TMD, either match or query)
        # note that the nonTMD percentage identity is calculated the same way
        df_TMD['%s_perc_ident'%TMD] = df_TMD['%s_SW_num_ident_res'%TMD] / df_TMD['%s_SW_align_len'%TMD]
        # calculate percentage similar residues
        df_TMD['%s_perc_sim'%TMD] = df_TMD['%s_SW_num_sim_res'%TMD] / df_TMD['%s_SW_align_len'%TMD]
        # add together to obtain the percentage similar + identical residues
        df_TMD['%s_perc_sim_plus_ident'%TMD] = df_TMD['%s_perc_ident'%TMD] + df_TMD['%s_perc_sim'%TMD]
        # add to main dataframe
        df_TMD['%s_perc_ident'%TMD] = df_TMD['%s_perc_ident'%TMD]
        # calculate the average number of gaps per residue in the TMD alignment
        # (number of gaps)/(length of sequence excluding gaps)
        df_TMD['%s_SW_q_gaps_per_q_residue'%TMD] = df_TMD['%s_SW_query_num_gaps'%TMD].dropna() / len_query_TMD
    return df_TMD

def calc_AAIMON(TMD, df_cr, mean_ser, logging):

    # calculate the Amino Acid Identity : Membranous Over Nonmembranous
    df_cr['%s_AAIMON_ratio'%TMD] = df_cr['%s_perc_ident'%TMD] / df_cr['nonTMD_perc_ident']
    # calculate the Amino Acid Similarity : Membranous Over Nonmembranous (AASMON) (includes similarity + identity based on the matrix used in the SW alignment of SIMAP)
    df_cr['%s_AASMON_ratio'%TMD] = df_cr['%s_perc_sim_plus_ident'%TMD] / df_cr['nonTMD_perc_sim_plus_ident']
    mean_ser['%s_SW_q_gaps_per_q_residue_mean'%TMD] = df_cr['%s_SW_q_gaps_per_q_residue'%TMD].dropna().mean()
    #logging.info('%s_SW_q_gaps_per_q_residue Average: %0.3e' %(TMD, mean_ser['%s_SW_q_gaps_per_q_residue_mean'%TMD]))

    # add to original dataframe with the list of uniprot sequences
    mean_ser['%s_perc_ident_mean'%TMD] = df_cr['%s_perc_ident'%TMD].mean()
    mean_ser['%s_perc_sim_mean'%TMD] = df_cr['%s_perc_sim'%TMD].mean()
    mean_ser['%s_perc_sim_plus_ident_mean'%TMD] = df_cr['%s_perc_sim_plus_ident'%TMD].mean()
    mean_ser['%s_AAIMON_ratio_mean'%TMD] = float(df_cr['%s_AAIMON_ratio'%TMD].mean())
    mean_ser['%s_AAIMON_ratio_std'%TMD] = df_cr['%s_AAIMON_ratio'%TMD].std()
    mean_ser['%s_AASMON_ratio_mean'%TMD] = df_cr['%s_AASMON_ratio'%TMD].mean()
    mean_ser['%s_AASMON_ratio_std'%TMD] = df_cr['%s_AASMON_ratio'%TMD].std()

    # add to the dataframe with the SIMAP data for that particular protein
    df_cr['%s_AAIMON_ratio'%TMD] = df_cr['%s_AAIMON_ratio'%TMD]
    df_cr['%s_AASMON_ratio'%TMD] = df_cr['%s_AASMON_ratio'%TMD]

    df_cr['%s_ratio_length_of_TMD_to_rest_of_alignment'%TMD] = df_cr['%s_SW_query_seq'%TMD].str.len() / df_cr['FASTA_overlap']
    df_cr['%s_ratio_length_of_query_TMD_to_rest_of_match_protein'%TMD] = df_cr['%s_SW_query_seq'%TMD].str.len() / df_cr[ 'len_full_match_seq']

    mean_ser['%s_ratio_length_of_TMD_to_rest_of_alignment_mean'%TMD] = float('%0.2f' % df_cr['%s_ratio_length_of_TMD_to_rest_of_alignment'%TMD].dropna().mean())
    mean_ser['%s_ratio_length_of_query_TMD_to_rest_of_match_protein_mean'%TMD] = float('%0.2f' % df_cr['%s_ratio_length_of_query_TMD_to_rest_of_match_protein'%TMD].dropna().mean())
    return mean_ser, df_cr