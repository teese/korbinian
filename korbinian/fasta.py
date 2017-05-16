import ast
import os
import korbinian
import korbinian.utils as utils
import sys
import zipfile
from multiprocessing import Pool
# import debugging tools
from korbinian.utils import pr, pc, pn, aaa

def run_create_fasta(pathdict, s, logging):
    logging.info('~~~~~~~~~~~~         starting filter_and_save_fasta           ~~~~~~~~~~~~')
    # if multiprocessing is used, log only to the console
    p_dict_logging = logging if s["use_multiprocessing"] != True else utils.Log_Only_To_Console()
    # create list of protein dictionaries to process
    list_p = korbinian.utils.convert_summary_csv_to_input_list(s, pathdict, p_dict_logging)
    # number of processes is the number the settings, or the number of proteins, whichever is smallest
    n_processes = s["multiprocessing_cores"] if s["multiprocessing_cores"] < len(list_p) else len(list_p)

    if s["use_multiprocessing"]:
        with Pool(processes=n_processes) as pool:
            fasta_list = pool.map(korbinian.fasta.filter_and_save_fasta, list_p)
            # log the list of protein results (e.g. acc, "simap", True) to the actual logfile, not just the console
            logging.info("fasta_list : {}".format(fasta_list))
    else:
        for p in list_p:
            korbinian.fasta.filter_and_save_fasta(p)
    logging.info('~~~~~~~~~~~~       filter_and_save_fasta is finished          ~~~~~~~~~~~~')

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
    PROTEIN_NAME_fa_fasta.zip (in homol folder, subfolder first two letters of accession)
        PROTEIN_NAME_homol_seq_TM01.fas
        PROTEIN_NAME_homol_seq_plus_surr_TM01.fas
    """
    pathdict, s, logging = p["pathdict"], p["s"], p["logging"]
    protein_name = p["protein_name"]
    acc = p["acc"]
    sys.stdout.write("{}. ".format(acc))
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

    if os.path.isfile(p["fa_fasta_zip"]):
        os.remove(p["fa_fasta_zip"])
    zipout_fasta = zipfile.ZipFile(p["fa_fasta_zip"], mode="a", compression=zipfile.ZIP_DEFLATED)

    ########################################################################################
    #                                                                                      #
    #           Filter based on homol hit properties (non-TMD-specific)                    #
    #                                                                                      #
    ########################################################################################

    fa_X_filt_full_str = " & X_in_match_seq == False" if s["fa_X_allowed_in_full_seq"] == False else ""

    fa_homol_query_str = 'FASTA_gapped_identity > {min_ident} & ' \
                        'FASTA_gapped_identity < {max_ident} & ' \
                        'hit_contains_SW_node == True & ' \
                        'disallowed_words_not_in_descr == True' \
                        '{Xfull}'.format(Xfull=fa_X_filt_full_str, min_ident=s["fa_min_identity_of_full_protein"], max_ident=s["fa_max_identity_of_full_protein"])

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
        if s["fa_X_allowed_in_full_seq"] == True:
            df_fa['X_in_%s'%TMD] = df_fa['%s_SW_match_seq'%TMD].str.contains("X")
            fa_X_allowed_in_sel_seq = s["fa_X_allowed_in_sel_seq"]
            fa_X_filt_sel_str = " & X_in_%s == False"%TMD if fa_X_allowed_in_sel_seq == False else ""
        else:
            fa_X_filt_sel_str = ""

        # create a boolean column that allows filtering by the accepted number of gaps, according to the settings file
        #df_fa['%s_fa_SW_query_acceptable_n_gaps'%TMD] = df_fa['%s_SW_query_num_gaps'%TMD] <= s["fa_max_n_gaps_in_query_TMD"]
       # df_fa['%s_fa_SW_match_acceptable_n_gaps'%TMD] = df_fa['%s_SW_match_num_gaps'%TMD] <= s["fa_max_n_gaps_in_match_TMD"]
        # measure the hydrophobicity of each TMD
        # %timeit 46.6 ms per loop for 325 homologues
        df_fa['%s_SW_match_lipo' % TMD] = df_fa['%s_SW_match_seq'%TMD].dropna().apply(lambda x: utils.calc_lipophilicity(x))

        # create string for the pandas.query syntax
        fa_query_filt_str =  '{TMD}_SW_query_num_gaps <= {fa_max_n_gaps_in_query_TMD} & ' \
                             '{TMD}_SW_match_num_gaps <= {fa_max_n_gaps_in_match_TMD} &' \
                             '{TMD}_SW_match_lipo <= {hydro_limit}' \
                             '{Xsel}'.format(TMD=TMD, Xsel=fa_X_filt_sel_str, fa_max_n_gaps_in_query_TMD=s["fa_max_n_gaps_in_query_TMD"],
                                             fa_max_n_gaps_in_match_TMD=s["fa_max_n_gaps_in_match_TMD"], hydro_limit=s["fa_max_hydrophilicity_Hessa"])

        # filter based on TMD-specific features
        df_fa.query(fa_query_filt_str, inplace=True)

        # setup the file names again. Note that the file should already exist, and the query sequence included.
        if s["save_fasta_plus_surr"] == True:
            fasta_savelist_suffixes = ["", "_plus_surr"]
        elif s["save_fasta_plus_surr"] == False:
            fasta_savelist_suffixes = [""]

        # for either the "_plus_surr" or "" suffix:
        for suff in fasta_savelist_suffixes:
            #fasta_file = p['fasta_file%s_BASENAME'%s] + '%s.fas' % TMD
            fasta_file_path = p['fasta_file%s_BASENAMEPATH'%suff] + '%s.fas' % TMD
            with open(fasta_file_path, 'w') as f:
                # add the query sequence, if desired
                if s["add_query_seq"]:
                    # add original query seq to fasta file. Note that the first SIMAP hit is excluded below.
                    #if '%s_seq%s'%(TMD,s) in df.columns:
                    if '%s_seq%s' % (TMD, suff) in p:
                        f.write('>00_%s_%s_uniprot_query%s\n%s\n' % (p['protein_name'], TMD, suff, p['%s_seq%s'%(TMD,suff)]))
                if s["remove_redundant_seqs"]:
                    # select the non-redundant sequences by using the pandas duplicated function
                    nr_dfs_fa = df_fa.loc[~df_fa['%s_SW_match_seq%s'%(TMD,suff)].duplicated()]
                else:
                    nr_dfs_fa = df_fa

                # if gaps are not allowed in the selected sequence, filter to remove before saving
                if s["fa_X_allowed_in_sel_seq"] == False:
                    nr_dfs_fa = nr_dfs_fa.dropna().loc[~nr_dfs_fa['%s_SW_match_seq%s'%(TMD,suff)].str.contains("X")]

                if suff == "_plus_surr":
                    if s["fa_X_allowed_in_full_seq"] == True:
                        if s["fa_X_allowed_in_sel_seq"] == False:
                            # note this is actually faster than tha str.contains("X") function
                            X_not_in_seq_ser = nr_dfs_fa['%s_SW_match_seq%s'%(TMD,suff)].apply(lambda x : "X" not in x)
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
                                                   nr_dfs_fa.loc[row, '%s_SW_match_seq%s'%(TMD,suff)]))

            # transfer temporary file into zip
            zipout_fasta.write(fasta_file_path, arcname=os.path.basename(fasta_file_path))
            # delete temporary file
            os.remove(fasta_file_path)
            n_fa_saved = int(nr_dfs_fa.shape[0])
            logging.info("{} {}{} saved to fasta, {} sequences.".format(acc, TMD, suff, n_fa_saved))

    # close the zipfile
    zipout_fasta.close()
    return acc, True, "0"
