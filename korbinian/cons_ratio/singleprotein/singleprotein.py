import korbinian.mtutils as utils

def slice_homologues_and_count_gaps(acc, TMD, df, dfs, set_):
    len_query_TMD = len(df.loc[acc, '%s_seq'%TMD])
    # use small throwaway functions to slice each TMD from the query, markup and match sequences
    # notnull() removes missing data
    # explanation: dfs['new_column_with_selected_seq'] = dfs[dfs['only_rows_containing_data]].apply(utils.slice_function_that_specifies_columns_with_start_and_stop)
    dfs_sel = dfs[dfs['%s_start_in_SW_alignment'%TMD].notnull()]
    # dfs['%s_SW_query_seq'%TMD] = dfs_sel.apply(utils.slice_SW_query_TMD_seq,args=(TMD,), axis=1)
    # dfs['%s_SW_markup_seq'%TMD] = dfs_sel.apply(utils.slice_SW_markup_TMD,args=(TMD,), axis=1)
    # dfs['%s_SW_match_seq'%TMD] = dfs_sel.apply(utils.slice_SW_match_TMD_seq,args=(TMD,), axis=1)
    dfs['%s_SW_query_seq' % TMD] = dfs_sel.apply(utils.slice_SW_query_TMD_seq, args=(TMD,), axis=1)
    dfs['%s_SW_markup_seq' % TMD] = dfs_sel.apply(utils.slice_SW_markup_TMD, args=(TMD,), axis=1)
    dfs['%s_SW_match_seq' % TMD] = dfs_sel.apply(utils.slice_SW_match_TMD_seq, args=(TMD,), axis=1)

    # count the number of gaps in the query and match sequences
    dfs['%s_SW_query_num_gaps'%TMD] = dfs['%s_SW_query_seq'%TMD].str.count("-")
    dfs['%s_SW_match_num_gaps'%TMD] = dfs['%s_SW_match_seq'%TMD].str.count("-")
    # calculate the length of the match TMD seq excluding gaps
    dfs['%s_SW_m_seq_len'%TMD] = dfs['%s_SW_match_seq'%TMD].str.len()
    # for the alignment length, take the smallest value from the length of query or match
    # this will exclude gaps from the length in the following calculations, preventing false "low conservation" where the query TMD is much longer than the match TMD)
    # note that for most calculations this is somewhat redundant, because the max number of acceptable gaps in sequence is probable ~2
    # use the mask function (faster than a lambda function) to replace any lengths larger than the query, with the query
    dfs['%s_SW_align_len' % TMD] = dfs['%s_SW_m_seq_len' % TMD].mask(dfs['%s_SW_m_seq_len' % TMD] > len_query_TMD, len_query_TMD)
    # create a boolean column that allows filtering by the accepted number of gaps, according to the settings file
    dfs['%s_cr_SW_query_acceptable_n_gaps'%TMD] = dfs['%s_SW_query_num_gaps'%TMD] <= set_["cr_max_n_gaps_in_query_TMD"]
    dfs['%s_cr_SW_match_acceptable_n_gaps'%TMD] = dfs['%s_SW_match_num_gaps'%TMD] <= set_["cr_max_n_gaps_in_match_TMD"]
    # count identical residues between query and match TMDs by counting the number of pipes in the markup string
    dfs['%s_SW_num_ident_res'%TMD] = dfs['%s_SW_markup_seq'%TMD].str.count('|')
    dfs['%s_SW_num_sim_res'%TMD] = dfs['%s_SW_markup_seq'%TMD].str.count(':')
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

def calc_AAIMON(acc, TMD, dfs, set_, df, logging):
    # following the general filters, filter to only analyse sequences with TMD identity above cutoff, and a nonTMD_perc_ident above zero ,to avoid a divide by zero error
    min_identity_of_TMD_initial_filter = set_['cr_min_identity_of_TMD_initial_filter']
    dfs_filt_AAIMON = dfs.loc[dfs['%s_perc_ident'%TMD] >= min_identity_of_TMD_initial_filter]
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