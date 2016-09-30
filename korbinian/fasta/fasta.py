import ast
import csv
import os
import pandas as pd
import korbinian.utils as utils
import zipfile

def filter_and_save_fasta(p):
    """Filters homologues obtained for each TMD, and saves as an unaligned FastA file.

    First filters based on homol hit properties (non-TMD-specific, e.g. percentage identity of full protein)
    For each TMD:
        Filters based on properties of each TMD (for example number of gaps in TMD sequence)


    Parameters
    ----------
    p : dict
        Protein Dictionary. Contains all input settings, sequences and filepaths related to a single protein.
        Protein-specific data is extracted from one row of the the list summary, e.g. List05_summary.csv, which is read as df.
        p also contains the GENERAL korbinian settings and filepaths for that list (pathdict, set_, logging)


        Components
        ----------
        pathdict : dict
            Dictionary of the key paths and files associated with that List number.
        set_ : dict
            Settings dictionary extracted from excel settings file.
        logging : logging.Logger
            Logger for printing to console and/or logfile.
            If multiprocessing == True, logging.info etc will only print to console.
        p : protein-specific dictionary components
            acc, list_of_TMDs, description, TM01_seq, etc

    Returns
    -------
    In all cases, a tuple (str, bool, str) is returned.

    if sucessful:
        return acc, True, "0"
    if not successful:
        return acc, False, "specific warning or reason why protein failed"

    Saved Files and Figures
    -----------------------
    PROTEIN_NAME_fa_fasta.zip (in homol folder, subfolder first two letters of accession)
        PROTEIN_NAME_homol_seq_TM01.fas
        PROTEIN_NAME_homol_seq_plus_surr_TM01.fas
    """
    pathdict, set_, logging = p["pathdict"], p["set_"], p["logging"]
    protein_name = p["protein_name"]
    acc = p["acc"]
    print(acc, end=", ", flush=True)
    # if the homol_df_orig_zip file does not exist, skip that protein
    if not os.path.exists(p['fa_cr_sliced_TMDs_zip']):
        message = "{} skipped, fa_cr_sliced_TMDs_zip not found.".format(acc)
        logging.info(message)
        return acc, False, message
    # open the dataframe containing the "match_align_seq" etc for each hit in the homologues
    dfh = utils.open_df_from_pickle_zip(p['homol_df_orig_zip'], filename=os.path.basename(p['homol_df_orig_pickle']), delete_corrupt=True)
    # skip the protein if the file was corrupt
    if dfh.empty:
        message = "{} Protein skipped, file deleted as it is possibly corrupt.".format(p['homol_df_orig_zip'])
        logging.info(message)
        return acc, False, message
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
    # if none of the homologues match this filter, skip the protein
    if dfh.empty:
        message = "{} skipped. No homologues left after filtering based on FASTA_gapped_identity, etc.".format(acc)
        logging.warning(message)
        return acc, False, message

    for TMD in list_of_TMDs:
        # open the dataframe containing the sequences, gap counts, etc for that TMD only
        df_fa = utils.open_df_from_pickle_zip(p['fa_cr_sliced_TMDs_zip'], filename="{}_{}_sliced_df.pickle".format(protein_name, TMD), delete_corrupt=True)
        if df_fa.empty:
            message = "{} skipped. df_fa of {} is an empty dataframe, due to lack of homologues or improper zipfile.".format(acc, TMD)
            logging.warning(message)
            # skip protein
            return acc, False, message
        # check if dataframes share any indices
        # isdisjoint will return True if two sets have a null intersection.
        no_common_hits = set(dfh.index).isdisjoint(set(df_fa.index))
        if no_common_hits == True:
            # there are no common elements in the homologues filtered for general properties (e.g. %ID of full prot)
            # and the TMD-specific properties (E.g. n_gaps in TMD). SKIP THIS TMD.
            message = "{} {} skipped. df_fa and dfh share no common homologues.".format(acc, TMD)
            logging.warning(message)
            # skip this TMD (other TMDs might be okay??)
            continue
        # filter based on dfh above, for general homologue settings (e.g. % identity of full protein)
        df_fa = df_fa.loc[dfh.index,:]

        ########################################################################################
        #                                                                                      #
        #           Filter and create FastA files with TMDs from homologues                    #
        #                                                                                      #
        ########################################################################################

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

    # close the zipfile
    zipout_fasta.close()
    return acc, True, "0"
