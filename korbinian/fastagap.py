import ast
import os
import korbinian
import korbinian.utils as utils
import sys
import zipfile
from multiprocessing import Pool

def run_fastagap_save(pathdict, s, logging):
    logging.info('~~~~~~~~~~~~           starting run_fastagap_save             ~~~~~~~~~~~~')
    # if multiprocessing is used, log only to the console
    p_dict_logging = logging if s["use_multiprocessing"] != True else utils.Log_Only_To_Console()
    # create list of protein dictionaries to process
    list_p = korbinian.utils.convert_summary_csv_to_input_list(s, pathdict, p_dict_logging)
    # number of processes is the number the settings, or the number of proteins, whichever is smallest
    n_processes = s["multiprocessing_cores"] if s["multiprocessing_cores"] < len(list_p) else len(list_p)

    if s["use_multiprocessing"]:
        with Pool(processes=n_processes) as pool:
            fastagap_list = pool.map(korbinian.fastagap.fastagap_save, list_p)
            # log the list of protein results (e.g. acc, "simap", True) to the actual logfile, not just the console
            for item in fastagap_list:
                logging.info(item)
    else:
        for p in list_p:
            korbinian.fastagap.fastagap_save(p)
    logging.info('~~~~~~~~~~~~           finished run_fastagap_save             ~~~~~~~~~~~~')

def fastagap_save(p):
    """ Saves TMD_plus_surr for homologues that contain gaps, for the fastagap analysis.

    First filters based on homol hit properties (non-TMD-specific, e.g. percentage identity of full protein)
    For each TMD:
        Filters based on properties of each TMD_plus_surr (for example number of gaps in TMD sequence)

    Parameters
    ----------
    p : dict
        Protein Dictionary. Contains all input settings, sequences and filepaths related to a single protein.
        Protein-specific data is extracted from one row of the the list summary, e.g. List05_summary.csv, which is read as df.
        p also contains the GENERAL korbinian settings and filepaths for that list (pathdict, s, logging)


        *Components of p*

        pathdict : dict
            Dictionary of the key paths and files associated with that List number.
        s : dict
            Settings dictionary extracted from excel settings file.
        logging : logging.Logger
            Logger for printing to console and/or logfile.
            If multiprocessing == True, logging.info etc will only print to console.
        p : protein-specific dictionary components
            acc, list_of_TMDs, description, TM01_seq, etc

    Returns
    -------
    In all cases, a tuple (str, bool, str) is returned.

    if successful:
        return acc, True, "0"
    if not successful:
        return acc, False, "specific warning or reason why protein failed"

    Saved Files and Figures
    -----------------------
    PROTEIN_NAME_fa_fastagap.zip (in homol folder, subfolder first two letters of accession)
        PROTEIN_NAME_homol_seq_plus_surr_TM01.fas
    """
    pathdict, s, logging = p["pathdict"], p["s"], p["logging"]
    protein_name = p["protein_name"]
    acc = p["acc"]
    sys.stdout.write("%s, "%acc)
    sys.stdout.flush()
    # if the fa_cr_sliced_TMDs_zip file does not exist, skip that protein
    if not os.path.exists(p['fa_cr_sliced_TMDs_zip']):
        message = "{} skipped, fa_cr_sliced_TMDs_zip not found.".format(acc)
        logging.info(message)
        return acc, False, message
    # if the homol_df_orig_zip file does not exist, skip that protein
    if not os.path.exists(p['homol_df_orig_zip']):
        message = "{} skipped, homol_df_orig_zip not found.".format(acc)
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

    # remove any existing outputfiles
    if os.path.isfile(p["fa_fastagap_zip"]):
        os.remove(p["fa_fastagap_zip"])

    zipout_fasta = zipfile.ZipFile(p["fa_fastagap_zip"], mode="a", compression=zipfile.ZIP_DEFLATED)

    ########################################################################################
    #                                                                                      #
    #           Filter based on homol hit properties (non-TMD-specific)                    #
    #                                                                                      #
    ########################################################################################

    # filtering similar to fasta.py, except that X in not allowed in any match sequence
    fa_homol_query_str ='FASTA_gapped_identity > {min_ident} & ' \
                        'FASTA_gapped_identity < {max_ident} & ' \
                        'hit_contains_SW_node == True & ' \
                        'disallowed_words_not_in_descr == True &' \
                        'X_in_match_seq == False'.format(min_ident=s["fa_min_identity_of_full_protein"], max_ident=s["fa_max_identity_of_full_protein"])

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

        # create a boolean column that allows filtering by the accepted number of gaps, according to the settings file
        #df_fa['%s_fa_SW_query_acceptable_n_gaps'%TMD] = df_fa['%s_SW_query_num_gaps'%TMD] <= s["fa_max_n_gaps_in_query_TMD"]
       # df_fa['%s_fa_SW_match_acceptable_n_gaps'%TMD] = df_fa['%s_SW_match_num_gaps'%TMD] <= s["fa_max_n_gaps_in_match_TMD"]
        # measure the hydrophobicity of each TMD
        # %timeit 46.6 ms per loop for 325 homologues
        df_fa['%s_SW_match_seq_hydro' % TMD] = df_fa['%s_SW_match_seq'%TMD].dropna().apply(lambda x: utils.calc_hydrophob(x))

        above_hessa_cutoff = df_fa['%s_SW_match_seq_hydro' % TMD] > s["gap_max_hydrophilicity_Hessa"]
        vc = above_hessa_cutoff.value_counts()
        if True in vc.index:
            n = vc[True]
            logging.info("{acc} {TMD} {n} homologues have hydrophilicity above_hessa_cutoff of {c}".format(acc=acc, TMD=TMD, n=n, c=s["gap_max_hydrophilicity_Hessa"]))

        min_number_of_gaps = s["fa_min_n_gaps_in_match_TMD_plus_surr"]
        max_number_of_gaps = s["fa_max_n_gaps_in_match_TMD_plus_surr"]

        if '{}_SW_query_seq_plus_surr'.format(TMD) not in df_fa.columns:
            message = "'{}_SW_query_seq_plus_surr'.format(TMD) not in df_fa. Suggest re-running slice script"
            logging.info(message)
            return acc, False, message

        df_fa['{}_SW_match_plus_surr_num_gaps'.format(TMD)] = df_fa['{}_SW_match_seq_plus_surr'.format(TMD)].str.count("-")
        df_fa['{}_SW_query_plus_surr_num_gaps'.format(TMD)] = df_fa['{}_SW_query_seq_plus_surr'.format(TMD)].str.count("-")
        # create string for the pandas.query syntax NOTE: gaps are not mentioned. Assume that gaps are NOT allowed in full protein.
        fa_query_filt_str = '{min_} <= {TMD}_SW_match_plus_surr_num_gaps <= {max_} &' \
                            '{TMD}_SW_match_seq_hydro <= {hydro_limit}'.format(TMD=TMD,min_=min_number_of_gaps,max_=max_number_of_gaps,
                                                                               hydro_limit = s["gap_max_hydrophilicity_Hessa"])

        # filter based on TMD-specific features
        df_fa.query(fa_query_filt_str, inplace=True)

        fastagap_file_path = p['fa_fastagap_base'] + '%s.fas' % TMD
        with open(fastagap_file_path, 'w') as f:
            # add the query sequence, if desired
            if s["add_query_seq"]:
                # add original query seq to fasta file. Note that the first SIMAP hit is excluded below.
                #if '%s_seq%s'%(TMD,s) in df.columns:
                if '%s_seq_plus_surr' % (TMD) in p:
                    f.write('>00_%s_%s_uniprot_query_plus_surr\n%s\n' % (p['protein_name'], TMD, p['%s_seq_plus_surr'%(TMD)]))
            if s["remove_redundant_seqs"]:
                # select the non-redundant sequences by using the pandas duplicated function
                nr_dfs_fa = df_fa.loc[~df_fa['%s_SW_match_seq_plus_surr'%(TMD)].duplicated()]
            else:
                nr_dfs_fa = df_fa

            # if gaps are not allowed in the selected sequence, filter to remove before saving
            if s["fa_X_allowed_in_sel_seq"] == False:
                nr_dfs_fa = nr_dfs_fa.dropna().loc[~nr_dfs_fa['%s_SW_match_seq_plus_surr'%(TMD)].str.contains("X")]

            if s["fa_X_allowed_in_full_seq"] == True:
                if s["fa_X_allowed_in_sel_seq"] == False:
                    # note this is actually faster than tha str.contains("X") function
                    X_not_in_seq_ser = nr_dfs_fa['%s_SW_match_seq_plus_surr'%(TMD)].apply(lambda x : "X" not in x)
                    n_before = nr_dfs_fa.shape[0]
                    nr_dfs_fa = nr_dfs_fa.loc[X_not_in_seq_ser]
                    n_X_removed = n_before - nr_dfs_fa.shape[0]
                    if n_X_removed > 1:
                        logging.info("{} {} seqs removed due to X in seq plus surr".format(acc, n_X_removed))

            """REMOVED, FIRST HIT IS ALREADY EXCLUDED?"""
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
                                               nr_dfs_fa.loc[row, '%s_SW_match_seq_plus_surr'%(TMD)]))

        # transfer temporary file into zip
        zipout_fasta.write(fastagap_file_path, arcname=os.path.basename(fastagap_file_path))
        # delete temporary file
        os.remove(fastagap_file_path)
        n_fa_saved = int(nr_dfs_fa.shape[0])
        logging.info("{} {}_plus_surr saved to fasta, {} sequences.".format(acc, TMD, n_fa_saved))

    # close the zipfile
    zipout_fasta.close()
    return acc, True, "0"
