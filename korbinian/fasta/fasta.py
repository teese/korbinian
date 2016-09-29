import ast
import csv
import os
import pandas as pd
import korbinian.mtutils as utils
import zipfile

def filter_and_save_fasta(p):

    # df = pd.read_csv(pathdict["list_summary_csv"], sep=",", quoting=csv.QUOTE_NONNUMERIC, index_col=0)
    # #iterate over the dataframe for proteins with an existing list_of_TMDs. acc = uniprot accession.
    # for acc in df.loc[df['list_of_TMDs'].notnull()].loc[df['list_of_TMDs'] != 'nan'].index:

    pathdict, set_, logging = p["pathdict"], p["set_"], p["logging"]
    protein_name = p["protein_name"]
    acc = p["acc"]
    print(acc, end=", ", flush=True)
    # if the homol_df_orig_zip file does not exist, skip that protein
    if not os.path.exists(p['homol_df_orig_zip']):
        logging.info("{} Protein skipped, file not found.".format(p['homol_df_orig_zip']))
        return
    # open the dataframe containing the "match_align_seq" etc for each hit in the homologues
    dfh = utils.open_df_from_pickle_zip(p['homol_df_orig_zip'], filename=os.path.basename(p['homol_df_orig_pickle']), delete_corrupt=True)
    # skip the protein if the file was corrupt
    if dfh.empty:
        logging.info("{} Protein skipped, file deleted as it is possibly corrupt.".format(p['homol_df_orig_zip']))
        return
    #list_of_TMDs = p['list_of_TMDs'].strip("[']").split(", ")
    list_of_TMDs = ast.literal_eval(p['list_of_TMDs'])

    if os.path.isfile(p["fa_fasta_zip"]):
        os.remove(p["fa_fasta_zip"])
    zipout_fasta = zipfile.ZipFile(p["fa_fasta_zip"], mode="a", compression=zipfile.ZIP_DEFLATED)

    ########################################################################################
    #                                                                                      #
    #           Filter based on homol hit properties (non-TMD-specific)                    #
    #                                                                                      #
    ########################################################################################

    fa_X_filt_full_str = " and X_in_match_seq == False" if set_["fa_X_allowed_in_full_seq"] == False else ""

    fa_homol_query_str = 'FASTA_gapped_identity > {min_ident} and ' \
                        'FASTA_gapped_identity < {max_ident} and ' \
                        'hit_contains_SW_node == True and ' \
                        'disallowed_words_not_in_descr == True' \
                        '{Xfull}'.format(Xfull=fa_X_filt_full_str, min_ident=set_["fa_min_identity_of_full_protein"], max_ident=set_["fa_max_identity_of_full_protein"])

    # filter based on the query string
    dfh.query(fa_homol_query_str, inplace=True)

    for TMD in list_of_TMDs:
        # open the dataframe containing the sequences, gap counts, etc for that TMD only
        df_fa = utils.open_df_from_pickle_zip(p['fa_cr_sliced_TMDs_zip'], filename="{}_{}_sliced_df.pickle".format(protein_name, TMD), delete_corrupt=True)
        # filter based on dfh above, for general homologue settings (e.g. % identity of full protein)
        df_fa = df_fa.loc[dfh.index,:]

        ########################################################################################
        #                                                                                      #
        #           Filter and create FastA files with TMDs from homologues                    #
        #                                                                                      #
        ########################################################################################

        list_fasta_cols = ['fa_ident_above_cutoff', 'FASTA_gapped_identity', "fa_min_identity_of_full_protein", "X_in_match_seq",
                           "match_align_seq", "disallowed_words_not_in_descr", "hit_contains_SW_node"]
        list_fasta_cols_TMD = ['%s_SW_match_seq'%TMD, '%s_start_in_SW_alignment_plus_surr'%TMD, '%s_end_in_SW_alignment_plus_surr'%TMD, '%s_SW_query_num_gaps'%TMD]

        # if "X" is allowed in the full sequence, check if X is in the selected sequence
        if set_["fa_X_allowed_in_full_seq"] == True:
            df_fa['X_in_%s'%TMD] = 'X' in df_fa['%s_SW_match_seq'%TMD]
            fa_X_allowed_in_sel_seq = set_["fa_X_allowed_in_sel_seq"]
            fa_X_filt_sel_str = " and X_in_%s == False"%TMD if fa_X_allowed_in_sel_seq == False else ""
        else:
            fa_X_filt_sel_str = ""

        # create a boolean column that allows filtering by the accepted number of gaps, according to the settings file
        df_fa['%s_fa_SW_query_acceptable_n_gaps'%TMD] = df_fa['%s_SW_query_num_gaps'%TMD] <= set_["fa_max_n_gaps_in_query_TMD"]
        df_fa['%s_fa_SW_match_acceptable_n_gaps'%TMD] = df_fa['%s_SW_match_num_gaps'%TMD] <= set_["fa_max_n_gaps_in_match_TMD"]
        # measure the hydrophobicity of each TMD
        # %timeit 46.6 ms per loop for 325 homologues
        df_fa['%s_SW_match_seq_hydro' % TMD] = df_fa['%s_SW_match_seq'%TMD].dropna().apply(lambda x: utils.calc_hydrophob(x))

        '''re-filter the original dataframe to create another copy with the desired sequences
        note that some values were added after filtering in the last round,
        but all were added to the dataframe df_fa, not the copy df_fa_filt

        fa_query_filt_str = 'fa_ident_above_cutoff == True and ' \ # aa identity above cutoff
                            'fa_ident_below_cutoff == True and '\ # aa identity below cutoff
                             'hit_contains_SW_node == True and '\ # homologue file is not missing data
                             'disallowed_words_not_in_descr == True'\ # not a patent
                             '{TMD}_fa_SW_query_acceptable_n_gaps == True' \ # not too many gaps in query
                             '{TMD}_fa_SW_match_acceptable_n_gaps == True' \ # not too many gaps in match
                             '{TMD}_SW_m_seq_len > 1' \ # smith waterman match sequence length longer than 1
                             '{Xfull}' \ # acceptable number of X in full protein
                             '{Xsel} # acceptable number of X in selection region of protein
        '''

        fa_query_filt_str = '{TMD}_fa_SW_query_acceptable_n_gaps == True and ' \
                             '{TMD}_fa_SW_match_acceptable_n_gaps == True and ' \
                             '{TMD}_SW_m_seq_len > 1' \
                             '{Xsel}'.format(TMD=TMD, Xsel=fa_X_filt_sel_str)

        df_fa.query(fa_query_filt_str, inplace=True)

        # # and the same for the TMD + surrounding sequence, useful to examine the TMD interface
        # df_fa['%s_SW_query_seq_plus_surr'%TMD] = df_fa_fa[df_fa_fa['%s_start_in_SW_alignment_plus_surr'%TMD].notnull()].apply(
        #     utils.slice_SW_query_TMD_seq_plus_surr, args=(TMD,), axis=1)
        # df_fa['%s_SW_markup_seq_plus_surr'%TMD] = df_fa_fa[df_fa_fa['%s_start_in_SW_alignment_plus_surr'%TMD].notnull()].apply(
        #     utils.slice_SW_markup_TMD_plus_surr, args=(TMD,), axis=1)
        # df_fa['%s_SW_match_seq_plus_surr'%TMD] = df_fa_fa[df_fa_fa['%s_start_in_SW_alignment_plus_surr'%TMD].notnull()].apply(
        #     utils.slice_SW_match_TMD_seq_plus_surr, args=(TMD,), axis=1)

        # start with the same dataframe copy that has filtered for gapped identity, etc (the AAIMAN ratio is not necessary for the FastA saving)
        # remove hits lacking sequence and also remove hits with too many gaps in TMD from either query or match
        #df_fa_filt_FastA = df_fa_filt.loc[df_fa['%s_SW_match_seq'%TMD].notnull()].loc[df_fa['%s_cr_SW_query_acceptable_n_gaps'%TMD]].loc[df_fa['%s_cr_SW_match_acceptable_n_gaps'%TMD]]
        # setup the file names again. Note that the file should already exist, and the query sequence included.

        if set_["save_fasta_plus_surr"] == True:
            # setup extension string for fastA plus surrounding sequence (interfacial region)
            fasta_savelist = ["", "_plus_surr"]
        else:
            fasta_savelist = [""]

        for s in fasta_savelist:
            #fasta_file = p['fasta_file%s_BASENAME'%s] + '%s.fas' % TMD
            fasta_file_path = p['fasta_file%s_BASENAMEPATH'%s] + '%s.fas' % TMD
            with open(fasta_file_path, 'w') as f:
                # add the query sequence, if desired
                if set_["add_query_seq"]:
                    # add original query seq to fasta file. Note that the first SIMAP hit is excluded below.
                    #if '%s_seq%s'%(TMD,s) in df.columns:
                    if '%s_seq%s' % (TMD, s) in p:
                        f.write('>00_%s_%s_uniprot_query%s\n%s\n' % (p['protein_name'], TMD, s, p['%s_seq%s'%(TMD,s)]))
                if set_["remove_redundant_seqs"]:
                    # select the non-redundant sequences by using the pandas duplicated function
                    nr_dfs_fa = df_fa.loc[~df_fa['%s_SW_match_seq%s'%(TMD,s)].duplicated()]
                else:
                    nr_dfs_fa = df_fa

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
                            n_X_removed = n_before - nr_dfs_fa.shape[0]
                            if n_X_removed > 1:
                                logging.info("{} {} seqs removed due to X in seq plus surr".format(acc, n_X_removed))

                """REMOVED, FIRST HIT IS ALREADY EXCLUDED"""
                # # add the first non-redundant sequence from the homologues, but only if it is not the same as the query
                # if p['%s_seq%s'%(TMD,s)] != nr_dfs_fa.loc[0, '%s_SW_match_seq%s'%(TMD,s)]:
                #     f.write('>%04d_%s_%s\n%s\n' % (1, str(nr_dfs_fa.loc[0, 'organism'])[:30],
                #                                    str(nr_dfs_fa.loc[0, 'description'])[:30],
                #                                    nr_dfs_fa.loc[0, '%s_SW_match_seq%s'%(TMD,s)]))
                # for each hit after the first one, add the sequence to the fastA file
                for row in nr_dfs_fa.index: # formerly excluding the first hit
                    # add the original query seq
                    f.write('>%04d_%s_%s\n%s\n' % (row, str(nr_dfs_fa.loc[row, 'organism'])[:30],
                                                   str(nr_dfs_fa.loc[row, 'description'])[:30],
                                                   nr_dfs_fa.loc[row, '%s_SW_match_seq%s'%(TMD,s)]))
                    # logging.info('saved ' + fasta_file_path)

            # transfer temporary file into zip
            zipout_fasta.write(fasta_file_path, arcname=os.path.basename(fasta_file_path))
            # delete temporary file
            os.remove(fasta_file_path)
            n_fa_saved = int(nr_dfs_fa.shape[0])
            logging.info("{} {}{} saved to fasta, {} sequences.".format(acc, TMD, s, n_fa_saved))

                # if s == "":
                #     fasta_file_path = r"D:\Schweris\Projects\Xiao\20160728 SIMAP vs HHBLITS fasta" + "\\" + '%s_10gap.fas' % p['protein_name']
                #     with open(fasta_file_path, 'w') as f:
                #         # add the query sequence, if desired
                #         if set_["add_query_seq"]:
                #             # add original query seq to fasta file. Note that the first SIMAP hit is excluded below.
                #             f.write('>00_%s_%s_uniprot_query%s\n%s\n' % (p['protein_name'], TMD, s, p['%s_seq%s' % (TMD, s)]))
                #         if set_["remove_redundant_seqs"]:
                #             # select the non-redundant sequences by using the pandas duplicated function
                #             nr_dfs_fa = df_fa.loc[~df_fa['%s_SW_match_seq%s' % (TMD, s)].duplicated()]
                #         else:
                #             nr_dfs_fa = df_fa
                #
                #         for row in nr_dfs_fa.index:  # formerly excluding the first hit
                #             # add the original query seq
                #             f.write('>%04d_%s_%s\n%s\n' % (row, str(dfs.loc[row, 'organism'])[:30],
                #                                            str(dfs.loc[row, 'description'])[:30],
                #                                            dfs.loc[row, '%s_SW_match_seq%s' % (TMD, s)]))

    # close the zipfile
    zipout_fasta.close()
    return acc, True, "0"
