import ast
import korbinian.utils as utils

def conduct_ssr_ratio_calculations(dfs, list_of_TMDs, s, list_of_aa_sub_matrices):

    def calc_score_ss_qTMD(dfs):
       score = sum(utils.score_pairwise(seq1=dfs['%s_SW_query_seq'%TMD], seq2=dfs['%s_SW_query_seq'%TMD],
                                    matrix=aa_sub_matrix,
                                    gap_open_penalty=gap_open_penalty,
                                    gap_extension_penalty=gap_extension_penalty))
       return(score)

    def calc_score_ss_mTMD(dfs):
       score = sum(utils.score_pairwise(seq1=dfs['%s_SW_match_seq'%TMD], seq2=dfs['%s_SW_match_seq'%TMD],
                                    matrix=aa_sub_matrix,
                                    gap_open_penalty=gap_open_penalty,
                                    gap_extension_penalty=gap_extension_penalty))
       return(score)

    def calc_score_qTMD_mTMD(dfs):
       score = sum(utils.score_pairwise(seq1=dfs['%s_SW_query_seq'%TMD], seq2=dfs['%s_SW_match_seq'%TMD],
                                    matrix=aa_sub_matrix,
                                    gap_open_penalty=gap_open_penalty,
                                    gap_extension_penalty=gap_extension_penalty))
       return(score)

    def calc_ss_q_nonTMD(dfs):
       score = sum(utils.score_pairwise(seq1=dfs['nonTMD_seq_query'], seq2=dfs['nonTMD_seq_query'],
                                    matrix=aa_sub_matrix,
                                    gap_open_penalty=gap_open_penalty,
                                    gap_extension_penalty=gap_extension_penalty))
       return(score)

    def calc_ss_m_nonTMD(dfs):
       score = sum(utils.score_pairwise(seq1=dfs['nonTMD_seq_match'], seq2=dfs['nonTMD_seq_match'],
                                    matrix=aa_sub_matrix,
                                    gap_open_penalty=gap_open_penalty,
                                    gap_extension_penalty=gap_extension_penalty))
       return(score)

    def calc_q_m_nonTMD(dfs):
       score = sum(utils.score_pairwise(seq1=dfs['nonTMD_seq_query'], seq2=dfs['nonTMD_seq_match'],
                                    matrix=aa_sub_matrix,
                                    gap_open_penalty=gap_open_penalty,
                                    gap_extension_penalty=gap_extension_penalty))
       return(score)

    for gap_open_penalty in range(s["gap_open_penalty_max"], s["gap_open_penalty_increment"]):
       #sys.stdout.write(gap_open_penalty)
       #for simplicity, give the gap open and gap extend the same value
       gap_extension_penalty = gap_open_penalty
       for matrix_name in list_of_aa_sub_matrices:
           #so long as the matrix is imported into python, eval will convert the matrix name to an object
           aa_sub_matrix = ast.literal_eval(matrix_name)
           #update the matrix (unsure what this deos! Taken directly from Stackoverflow)
           aa_sub_matrix.update(((b, a), val) for (a, b), val in list(aa_sub_matrix.items()))
           column_basename = 'sim_ratio_%s_gapo%i' % (matrix_name.replace("'", "")[0:-7], gap_open_penalty)
           #sys.stdout.write(column_name)
           #scores of query and match. Jan used match-query and query-match because in some cases in his algorithm it gave different scores. In this case, the score is independent of the sequence order, and 2*score_qTMD_mTMD can be used instead of score_qTMD_mTMD + score_mTMD_qTMD
           #score_qTMD_mTMD = sum(utils.score_pairwise(seq1=SW_query_TMD_seq, seq2=SW_match_TMD_seq,
           #                         matrix=aa_sub_matrix,
           #                         gap_open_penalty=gap_open_penalty,
           #                         gap_extension_penalty=gap_extension_penalty))
           #sys.stdout.write(score_qTMD_mTMD)
           dfs_nonTMD = dfs.query('"X" not in match_align_seq')
           dfs_nonTMD = dfs_nonTMD.loc[dfs_nonTMD['match_align_seq'].notnull()]
           dfs_nonTMD = dfs_nonTMD.loc[dfs_nonTMD['nonTMD_seq_query'].notnull()]
           dfs_nonTMD = dfs_nonTMD.loc[dfs_nonTMD['match_align_seq'].apply(lambda x : 'X' not in x)]

           #score/self-score ratio of nonTMD query
           dfs_nonTMD[column_basename + '_ss_q_nonTMD'] = dfs_nonTMD.apply(calc_ss_q_nonTMD, axis = 1)
           #score/self-score ratio of match
           dfs_nonTMD[column_basename + '_ss_m_nonTMD'] = dfs_nonTMD.apply(calc_ss_m_nonTMD, axis = 1)
           #scores of query and match. Jan used match-query and query-match because in some cases in his algorithm it gave different scores. In this case, the score is independent of the sequence order, and 2*score_qTMD_mTMD can be used instead of score_qTMD_mTMD + score_mTMD_qTMD
           dfs_nonTMD[column_basename + '_q_m_nonTMD'] = dfs_nonTMD.apply(calc_q_m_nonTMD, axis = 1)
           #calculate the score/selfscore ratio
           dfs_nonTMD[column_basename + '_ssr_nonTMD'] = dfs_nonTMD[column_basename + '_q_m_nonTMD'] * 2 / (dfs_nonTMD[column_basename + '_ss_q_nonTMD'] + dfs_nonTMD[column_basename + '_ss_m_nonTMD'])
           #add to main dataframe
           dfs[column_basename + '_ssr_nonTMD'] = dfs_nonTMD[column_basename + '_ssr_nonTMD']

           for TMD in list_of_TMDs:
               column_name = TMD + '_' + column_basename
               dfs_nonTMD = dfs_nonTMD.loc[dfs['%s_SW_query_seq'%TMD].notnull()]
               #dfs_nonTMD = dfs_nonTMD.loc[dfs['X_in_match_seq'] == False]
               #score/self-score ratio of query
               dfs[column_name + '_ss_qTMD'] = dfs_nonTMD.apply(calc_score_ss_qTMD, axis = 1)
               #score/self-score ratio of match
               dfs[column_name + '_ss_mTMD'] = dfs_nonTMD.apply(calc_score_ss_mTMD, axis = 1)
               #scores of query and match. Jan used match-query and query-match because in some cases in his algorithm it gave different scores. In this case, the score is independent of the sequence order, and 2*score_qTMD_mTMD can be used instead of score_qTMD_mTMD + score_mTMD_qTMD
               dfs[column_name + '_qTMD_mTMD'] = dfs_nonTMD.apply(calc_score_qTMD_mTMD, axis = 1)
               #score/self-score ratio
               dfs[column_name + '_ssrTMD'] = dfs[column_name + '_qTMD_mTMD'] * 2 / (dfs[column_name + '_ss_qTMD'] + dfs[column_name + '_ss_mTMD'])
               #calculate the ssrTMD/ssr_nonTMD
               dfs_filt3 = dfs.loc[dfs[column_name + '_ssrTMD'].notnull()]
               dfs_filt3 = dfs.loc[dfs[column_name + '_ssrTMD'] > 0]
               dfs_filt3 = dfs.loc[dfs[column_basename + '_ssr_nonTMD'] > 0]
               dfs_filt3[column_name + '_ssrTMD_over_nonTMD'] = dfs[column_name + '_ssrTMD'] / dfs[column_basename + '_ssr_nonTMD']
               #add to main dataframe
               dfs[column_name + '_ssrTMD_over_nonTMD'] = dfs_filt3[column_name + '_ssrTMD_over_nonTMD']
    '''  _________________________________________END SSR ratio calculations____________________________________________________
    '''