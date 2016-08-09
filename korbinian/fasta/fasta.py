import ast
import os
import korbinian.mtutils as utils
import zipfile

def filter_and_save_fasta(df, dfs, acc, TMD, set_, logging, zipout_fasta):
    list_fasta_cols = ['fa_ident_above_cutoff', 'FASTA_gapped_identity', "fa_min_identity_of_full_protein", "X_in_match_seq",
                       "match_alignment_sequence", "cr_disallowed_words_not_in_descr", "hit_contains_SW_node"]
    list_fasta_cols_TMD = ['%s_SW_match_seq'%TMD, '%s_start_in_SW_alignment_plus_surr'%TMD, '%s_end_in_SW_alignment_plus_surr'%TMD, '%s_SW_query_num_gaps'%TMD]

    # redefine the number of amino acids before and after the TMD to be inserted into the FastA files
    fa_aa_before_tmd = set_["fa_aa_before_tmd"]
    fa_aa_after_tmd = set_["fa_aa_after_tmd"]

    # define the start of theTMD + surrounding sequence
    dfs['%s_start_in_SW_alignment_plus_surr'%TMD] = dfs['%s_start_in_SW_alignment'%TMD] - fa_aa_before_tmd
    # replace negative values with zero
    dfs.loc[dfs['%s_start_in_SW_alignment_plus_surr'%TMD] < 0, '%s_start_in_SW_alignment_plus_surr'%TMD] = 0
    # define the end of the TMD + surrounding sequence
    dfs['%s_end_in_SW_alignment_plus_surr'%TMD] = dfs['%s_end_in_SW_alignment'%TMD] + fa_aa_after_tmd
    # replace values longer than the actual sequence with the length of the sequence
    # dfs.loc[dfs['%s_end_in_SW_alignment_plus_surr'%TMD] > dfs['len_query_alignment_sequence'],
    #         '%s_end_in_SW_alignment_plus_surr'%TMD] = dfs['len_query_alignment_sequence']

    # and the same for the TMD + surrounding sequence, useful to examine the TMD interface
    #dfs['%s_SW_query_seq_plus_surr'%TMD] = dfs_fa.apply(utils.slice_SW_query_TMD_seq_plus_surr, args=(TMD,), axis=1)
    #dfs['%s_SW_markup_seq_plus_surr'%TMD] = dfs_fa.apply(utils.slice_SW_markup_TMD_plus_surr, args=(TMD,), axis=1)

    dfs['fa_ident_above_cutoff'] = dfs['FASTA_gapped_identity'] > set_["fa_min_identity_of_full_protein"]
    dfs['fa_ident_below_cutoff'] = dfs['FASTA_gapped_identity'] < set_["fa_max_identity_of_full_protein"]

    fa_X_filt_full_str = " and X_in_match_seq == False" if set_["fa_X_allowed_in_full_seq"] == False else ""

    # if "X" is allowed in the full sequence, check if X is in the selected sequence
    if set_["fa_X_allowed_in_full_seq"] == True:
        dfs['X_in_%s'%TMD] = 'X' in dfs['match_alignment_sequence']
        fa_X_allowed_in_sel_seq = set_["fa_X_allowed_in_sel_seq"]
        fa_X_filt_sel_str = " and X_in_%s == False"%TMD if fa_X_allowed_in_sel_seq == False else ""
    else:
        fa_X_filt_sel_str = ""

    # create a boolean column that allows filtering by the accepted number of gaps, according to the settings file
    dfs['%s_fa_SW_query_acceptable_n_gaps'%TMD] = dfs['%s_SW_query_num_gaps'%TMD] <= set_["fa_max_n_gaps_in_query_TMD"]
    dfs['%s_fa_SW_match_acceptable_n_gaps'%TMD] = dfs['%s_SW_match_num_gaps'%TMD] <= set_["fa_max_n_gaps_in_match_TMD"]
    # measure the hydrophobicity of each TMD
    # %timeit 46.6 ms per loop for 325 homologues
    dfs['%s_SW_match_seq_hydro' % TMD] = dfs['%s_SW_match_seq'%TMD].dropna().apply(lambda x: utils.calc_hydrophob(x))

    # check if the list of disallowed proteins is the same.
    fa_search_disallowed_words = False
    if set_["fa_words_not_allowed_in_description"] == set_["cr_words_not_allowed_in_description"]:
        # if the word searching has already been done, simply copy the results across
        if "cr_disallowed_words_not_in_descr" in dfs.columns:
            dfs["fa_disallowed_words_not_in_descr"] = dfs["cr_disallowed_words_not_in_descr"]
        else:
            fa_search_disallowed_words = True
    else:
        fa_search_disallowed_words = True

    # conduct the text searching for disallowed words
    if fa_search_disallowed_words == True:
        fa_words_not_allowed_in_description = ast.literal_eval(set_["fa_words_not_allowed_in_description"])
        # collect disallowed words in hit protein description (patent, synthetic, etc)
        dfs['fa_list_disallowed_words_in_descr'] = dfs['uniprot_description'].dropna().apply(utils.find_disallowed_words, args=(fa_words_not_allowed_in_description,))
        # create a boolean column to select hits that do not contain these words in the description
        dfs['fa_disallowed_words_not_in_descr'] = dfs['fa_list_disallowed_words_in_descr'] == '[]'

    # select sequences that seem to have a start
    dfs_sel2 = dfs.loc[dfs['%s_start_in_SW_alignment_plus_surr'%TMD].notnull()]
    dfs['%s_SW_match_seq_plus_surr'%TMD] = dfs_sel2.apply(utils.slice_SW_match_TMD_seq_plus_surr, args=(TMD,), axis=1)

    '''re-filter the original dataframe to create another copy with the desired sequences
    note that some values were added after filtering in the last round,
    but all were added to the dataframe dfs, not the copy dfs_filt

    fa_query_filt_str = 'fa_ident_above_cutoff == True and ' \ # aa identity above cutoff
                        'fa_ident_below_cutoff == True and '\ # aa identity below cutoff
                         'hit_contains_SW_node == True and '\ # homologue file is not missing data
                         'fa_disallowed_words_not_in_descr == True'\ # not a patent
                         '{TMD}_fa_SW_query_acceptable_n_gaps == True' \ # not too many gaps in query
                         '{TMD}_fa_SW_match_acceptable_n_gaps == True' \ # not too many gaps in match
                         '{TMD}_SW_m_seq_len > 1' \ # smith waterman match sequence length longer than 1
                         '{Xfull}' \ # acceptable number of X in full protein
                         '{Xsel} # acceptable number of X in selection region of protein
    '''

    fa_query_filt_str = 'fa_ident_above_cutoff == True and ' \
                        'fa_ident_below_cutoff == True and '\
                         'hit_contains_SW_node == True and '\
                         'fa_disallowed_words_not_in_descr == True and '\
                         '{TMD}_fa_SW_query_acceptable_n_gaps == True and ' \
                         '{TMD}_fa_SW_match_acceptable_n_gaps == True and ' \
                         '{TMD}_SW_m_seq_len > 1 and ' \
                         '{TMD}_start_in_SW_alignment >= 1' \
                         '{Xfull}' \
                         '{Xsel}'.format(TMD=TMD, Xfull=fa_X_filt_full_str, Xsel=fa_X_filt_sel_str)

    dfs_fa = dfs.query(fa_query_filt_str)

    # refilter to obtain dfs_fa with the sliced sequences
    # timeit result: filtering of 5000 homologues took ~8 ms, deemed worth it as you need to filter before slicing anyway

    dfs_fa = dfs.query(fa_query_filt_str)

    # # and the same for the TMD + surrounding sequence, useful to examine the TMD interface
    # dfs['%s_SW_query_seq_plus_surr'%TMD] = dfs_fa[dfs_fa['%s_start_in_SW_alignment_plus_surr'%TMD].notnull()].apply(
    #     utils.slice_SW_query_TMD_seq_plus_surr, args=(TMD,), axis=1)
    # dfs['%s_SW_markup_seq_plus_surr'%TMD] = dfs_fa[dfs_fa['%s_start_in_SW_alignment_plus_surr'%TMD].notnull()].apply(
    #     utils.slice_SW_markup_TMD_plus_surr, args=(TMD,), axis=1)
    # dfs['%s_SW_match_seq_plus_surr'%TMD] = dfs_fa[dfs_fa['%s_start_in_SW_alignment_plus_surr'%TMD].notnull()].apply(
    #     utils.slice_SW_match_TMD_seq_plus_surr, args=(TMD,), axis=1)

    # start with the same dataframe copy that has filtered for gapped identity, etc (the AAIMAN ratio is not necessary for the FastA saving)
    # remove hits lacking sequence and also remove hits with too many gaps in TMD from either query or match
    #dfs_filt_FastA = dfs_filt.loc[dfs['%s_SW_match_seq'%TMD].notnull()].loc[dfs['%s_cr_SW_query_acceptable_n_gaps'%TMD]].loc[dfs['%s_cr_SW_match_acceptable_n_gaps'%TMD]]
    # setup the file names again. Note that the file should already exist, and the query sequence included.

    if set_["save_fasta_plus_surr"] == True:
        # setup extension string for fastA plus surrounding sequence (interfacial region)
        fasta_savelist = ["", "_plus_surr"]
    else:
        fasta_savelist = [""]

    for s in fasta_savelist:
        #fasta_file = df.loc[acc, 'fasta_file%s_BASENAME'%s] + '%s.fas' % TMD
        fasta_file_path = df.loc[acc, 'fasta_file%s_BASENAMEPATH'%s] + '%s.fas' % TMD
        print(fasta_file_path)
        with open(fasta_file_path, 'w') as f:
            # add the query sequence, if desired
            if set_["add_query_seq"]:
                # add original query seq to fasta file. Note that the first SIMAP hit is excluded below.
                f.write('>00_%s_%s_uniprot_query%s\n%s\n' % (df.loc[acc, 'protein_name'], TMD, s, df.loc[acc, '%s_seq%s'%(TMD,s)]))
            if set_["remove_redundant_seqs"]:
                # select the non-redundant sequences by using the pandas duplicated function
                nr_dfs_fa = dfs_fa.loc[~dfs_fa['%s_SW_match_seq%s'%(TMD,s)].duplicated()]
            else:
                nr_dfs_fa = dfs_fa

            # if gaps are not allowed in the selected sequence, filter to remove before saving
            if set_["fa_X_allowed_in_sel_seq"] == False:
                nr_dfs_fa = nr_dfs_fa.dropna().loc[~nr_dfs_fa['%s_SW_match_seq%s'%(TMD,s)].str.contains("X")]

            if s == "_plus_surr":
                if set_["fa_X_allowed_in_full_seq"] == True:
                    if set_["fa_X_allowed_in_sel_seq"] == False:
                        # note this is actually faster than tha str.contains("X") function
                        X_not_in_seq_ser = nr_dfs_fa['%s_SW_match_seq%s'%(TMD,s)].apply(lambda x : "X" not in x)
                        n_before = nr_dfs_fa.shape[0]
                        nr_dfs_fa = nr_dfs_fa.loc[X_not_in_seq_ser]
                        print("{} seqs removed due to X in seq plus surr".format(n_before - nr_dfs_fa.shape[0]))

            """REMOVED, FIRST HIT IS ALREADY EXCLUDED"""
            # # add the first non-redundant sequence from the homologues, but only if it is not the same as the query
            # if df.loc[acc, '%s_seq%s'%(TMD,s)] != nr_dfs_fa.loc[0, '%s_SW_match_seq%s'%(TMD,s)]:
            #     f.write('>%04d_%s_%s\n%s\n' % (1, str(nr_dfs_fa.loc[0, 'organism'])[:30],
            #                                    str(nr_dfs_fa.loc[0, 'uniprot_description'])[:30],
            #                                    nr_dfs_fa.loc[0, '%s_SW_match_seq%s'%(TMD,s)]))
            # for each hit after the first one, add the sequence to the fastA file
            for row in nr_dfs_fa.index: # formerly excluding the first hit
                # add the original query seq
                f.write('>%04d_%s_%s\n%s\n' % (row, str(dfs.loc[row, 'organism'])[:30],
                                               str(dfs.loc[row, 'uniprot_description'])[:30],
                                               dfs.loc[row, '%s_SW_match_seq%s'%(TMD,s)]))
                # logging.info('saved ' + fasta_file_path)

        # transfer temporary file into zip
        zipout_fasta.write(fasta_file_path, arcname=os.path.basename(fasta_file_path))
        # delete temporary file
        os.remove(fasta_file_path)
        n_fa_saved = int(nr_dfs_fa.shape[0])
        logging.info("TMD{} saved to fasta, {} sequences.".format(s, n_fa_saved))

            # if s == "":
            #     fasta_file_path = r"D:\Schweris\Projects\Xiao\20160728 SIMAP vs HHBLITS fasta" + "\\" + '%s_10gap.fas' % df.loc[acc, 'protein_name']
            #     with open(fasta_file_path, 'w') as f:
            #         # add the query sequence, if desired
            #         if set_["add_query_seq"]:
            #             # add original query seq to fasta file. Note that the first SIMAP hit is excluded below.
            #             f.write('>00_%s_%s_uniprot_query%s\n%s\n' % (df.loc[acc, 'protein_name'], TMD, s, df.loc[acc, '%s_seq%s' % (TMD, s)]))
            #         if set_["remove_redundant_seqs"]:
            #             # select the non-redundant sequences by using the pandas duplicated function
            #             nr_dfs_fa = dfs_fa.loc[~dfs_fa['%s_SW_match_seq%s' % (TMD, s)].duplicated()]
            #         else:
            #             nr_dfs_fa = dfs_fa
            #
            #         for row in nr_dfs_fa.index:  # formerly excluding the first hit
            #             # add the original query seq
            #             f.write('>%04d_%s_%s\n%s\n' % (row, str(dfs.loc[row, 'organism'])[:30],
            #                                            str(dfs.loc[row, 'uniprot_description'])[:30],
            #                                            dfs.loc[row, '%s_SW_match_seq%s' % (TMD, s)]))


    return dfs