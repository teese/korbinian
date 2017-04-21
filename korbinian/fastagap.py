import ast
import csv
import korbinian
import korbinian.utils as utils
import zipfile
from multiprocessing import Pool
import pandas as pd
import pickle
import os
import sys
from Bio import AlignIO
import numpy as np

def save_fastagap(pathdict, s, logging):
    """Runs fastagap_save for each protein, using multiprocessing Pool.

    Parameters
    ----------
    pathdict : dict
        Dictionary of the key paths and files associated with that List number.
    s : dict
        Settings dictionary extracted from excel settings file.
    logging : logging.Logger
        Logger for printing to console and/or logfile.
        If multiprocessing == True, logging.info etc will only print to console.

    """

    logging.info('~~~~~~~~~~~~           starting save_fastagap             ~~~~~~~~~~~~')
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
    logging.info('~~~~~~~~~~~~           finished save_fastagap             ~~~~~~~~~~~~')

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
    if os.path.isfile(p["fastagap_zip"]):
        os.remove(p["fastagap_zip"])

    zipout_fastagap = zipfile.ZipFile(p["fastagap_zip"], mode="a", compression=zipfile.ZIP_DEFLATED)

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

    n_valid_homol_dict = {}

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
        df_fa['%s_SW_match_lipo' % TMD] = df_fa['%s_SW_match_seq'%TMD].dropna().apply(lambda x: utils.calc_lipophilicity(x))

        above_hessa_cutoff = df_fa['%s_SW_match_lipo' % TMD] > s["max_lipo_homol"]
        vc = above_hessa_cutoff.value_counts()
        if True in vc.index:
            n = vc[True]
            logging.info("{acc} {TMD} {n}/{m} homologues have hydrophilicity above_hessa_cutoff "
                         "of {c}".format(acc=acc, TMD=TMD, n=n, c=s["max_lipo_homol"],m=df_fa.shape[0]))

        min_number_of_gaps = s["fa_min_n_gaps_in_match_TMD_plus_surr"]
        max_number_of_gaps = s["fa_max_n_gaps_in_match_TMD_plus_surr"]

        if '{}_SW_query_seq_plus_surr'.format(TMD) not in df_fa.columns:
            message = "'{}_SW_query_seq_plus_surr'.format(TMD) not in df_fa. Suggest re-running slice script"
            logging.info(message)
            return acc, False, message

        # count the number of gaps "-" in the sequence
        df_fa['{}_SW_match_plus_surr_num_gaps'.format(TMD)] = df_fa['{}_SW_match_seq_plus_surr'.format(TMD)].str.count("-")
        df_fa['{}_SW_query_plus_surr_num_gaps'.format(TMD)] = df_fa['{}_SW_query_seq_plus_surr'.format(TMD)].str.count("-")
        # calculate the length of the seq plus_surr for the match
        # because of alignment truncations, the surrounding sequence is not necessarily available
        df_fa["{}_SW_match_seq_plus_surr_len".format(TMD)] = df_fa["{}_SW_match_seq_plus_surr".format(TMD)].str.len()
        # calculate the query length plus surrounding sequence (end-start) + number of aa before TMD + number of aa after TMD
        query_len_plus_surr = p["{}_end".format(TMD)] - p["{}_start".format(TMD)] + s["n_aa_before_tmd"] + s["n_aa_after_tmd"]
        # the seq_plus_surr is truncated if the match length (including gaps) is shorter than the query sequence length
        df_fa["{}_plus_surr_is_truncated".format(TMD)] = df_fa["{}_SW_match_seq_plus_surr_len".format(TMD)] < query_len_plus_surr
        # # calculate the actual length of the seq_plus_surr for the homologues (len - n_gaps)
        # df_fa["{}_SW_match_seq_plus_surr_len_excl_gaps".format(TMD)] = df_fa["{}_SW_match_seq_plus_surr_len".format(TMD)] - df_fa['{}_SW_match_plus_surr_num_gaps'.format(TMD)]

        # create string for the pandas.query syntax NOTE: gaps are not mentioned. Assume that gaps are NOT allowed in full protein.
        fa_query_filt_str = '{min_} <= {TMD}_SW_match_plus_surr_num_gaps <= {max_} &' \
                            '{TMD}_SW_match_lipo <= {hydro_limit} & ' \
                            '{TMD}_plus_surr_is_truncated == False'.format(TMD=TMD,min_=min_number_of_gaps,max_=max_number_of_gaps,
                                                                               hydro_limit = s["max_lipo_homol"])

        # filter based on TMD-specific features
        df_fa.query(fa_query_filt_str, inplace=True)

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

        if nr_dfs_fa.empty:
            # there are no valid sequence sequences
            n_valid_homol_dict[TMD] = 0
            continue
        else:
            n_valid_homol_dict[TMD] = nr_dfs_fa.shape[0]
            # REMOVE ANY GAPS. This means that the gap analysis depends on the alignment from ClustalO, not the original SW alignment
            nr_dfs_fa[r'%s_SW_match_seq_plus_surr' % (TMD)] = nr_dfs_fa['%s_SW_match_seq_plus_surr' % (TMD)].str.replace("-", "")

        if n_valid_homol_dict[TMD] != 0:

            """REMOVED, FIRST HIT IS ALREADY EXCLUDED?"""
            # # add the first non-redundant sequence from the homologues, but only if it is not the same as the query
            # if p['%s_seq%s'%(TMD,s)] != nr_dfs_fa.loc[0, '%s_SW_match_seq%s'%(TMD,s)]:
            #     f.write('>%04d_%s_%s\n%s\n' % (1, str(nr_dfs_fa.loc[0, 'organism'])[:30],
            #                                    str(nr_dfs_fa.loc[0, 'description'])[:30],
            #                                    nr_dfs_fa.loc[0, '%s_SW_match_seq%s'%(TMD,s)]))
            # for each hit after the first one, add the sequence to the fastA file

            fastagap_filtered_fasta_path = os.path.normpath(p['fastagap_base'] + '%s_filtered.fas' % TMD)

            if os.path.isfile(fastagap_filtered_fasta_path):
                os.remove(fastagap_filtered_fasta_path)

            with open(fastagap_filtered_fasta_path, 'w') as f:

                # add the query sequence, if desired
                if s["add_query_seq"]:
                    # add original query seq to fasta file. Note that the first SIMAP hit is excluded below.
                    #if '%s_seq%s'%(TMD,s) in df.columns:
                    if '%s_seq_plus_surr'%TMD in p and '%s_seq_plus_surr'%TMD != np.nan:
                        f.write('>00_%s_%s_uniprot_query_plus_surr\n%s\n' % (p['protein_name'], TMD, p['%s_seq_plus_surr'%(TMD)]))
                    else:
                        message = "{} {}_seq_plus_surr is is missing (np.nan or not in database)".format(acc, TMD)
                        logging.warning(message)
                        continue

                for row in nr_dfs_fa.index: # formerly excluding the first hit
                    # add the original query seq
                    f.write('>%04d_%s_%s\n%s\n' % (row, str(nr_dfs_fa.loc[row, 'organism'])[:30],
                                                   str(nr_dfs_fa.loc[row, 'description'])[:30],
                                                   nr_dfs_fa.loc[row, '%s_SW_match_seq_plus_surr'%(TMD)]))


            # transfer temporary file into zip
            zipout_fastagap.write(fastagap_filtered_fasta_path, arcname=os.path.basename(fastagap_filtered_fasta_path))
            n_fa_saved = nr_dfs_fa.shape[0]
            logging.info("{} {}_plus_surr saved to fasta, {} sequences.".format(acc, TMD, n_fa_saved))

            fastagap_aligned_fasta_path = p['fastagap_base'] + '%s_aligned.fas' % TMD
            clustal_path = os.path.normpath(s["clustal_path"])
            """
            cp : path to clustalo executable
            -i filt_fasta : input unaligned file
            --output-order={input-order} : reorders the output sequences so the query is at the top again
            --force forces the overwriting of any existing output files
            """
            string_to_run_omega_for_alignment = "{cp}  -i  {filt_fasta} -o {out} --output-order={inputorder} --force".format(cp=clustal_path, filt_fasta=fastagap_filtered_fasta_path, out=fastagap_aligned_fasta_path, inputorder = "{input-order}")

            command = utils.Command(string_to_run_omega_for_alignment)
            # set the timeout for 10 seconds
            command.run(timeout=1)
            # there is sometimes a problem with "file used by process", probably because the command is not fully finished.
            # Wait 2 seconds, allowing time to save and finish with the alignment file.
            utils.sleep_x_seconds(2, print_stuff=False)

            # delete temporary files
            os.remove(fastagap_filtered_fasta_path)
            # zip and delete aligned fasta file
            if os.path.isfile(fastagap_aligned_fasta_path):
                zipout_fastagap.write(fastagap_aligned_fasta_path, arcname=os.path.basename(fastagap_aligned_fasta_path))
                os.remove(fastagap_aligned_fasta_path)

    # save the dictionary containing the number of valid homologues for each TMD as a csv file, and transfer to zip
    fastagap_n_valid_homol_csv = p['fastagap_base'] + 'n_valid_homol.csv'
    pd.Series(n_valid_homol_dict).to_csv(fastagap_n_valid_homol_csv)
    zipout_fastagap.write(fastagap_n_valid_homol_csv, arcname=os.path.basename(fastagap_n_valid_homol_csv))
    os.remove(fastagap_n_valid_homol_csv)

    # close the zipfile
    zipout_fastagap.close()
    return acc, True, "0"

def run_calc_fastagap_densities(pathdict, s, logging):
    logging.info("~~~~~~~~~~~~         starting calc_fastagap_densities           ~~~~~~~~~~~~")
    df = pd.read_csv(pathdict["list_csv"], sep=",", quoting=csv.QUOTE_NONNUMERIC, index_col=0)
    # get list of accessions that could not be downloaded, and can immediately be excluded
    not_in_homol_db = utils.get_acc_list_from_txt(pathdict["acc_not_in_homol_db_txt"])
    acc_kept = set(df.index) - set(not_in_homol_db)
    # filter to remove proteins not in the homologue database
    df = df.loc[acc_kept, :]
    # remove any proteins from list that do not have a list of TMDs
    df = df.loc[df.list_of_TMDs.notnull()]
    # create a list to contain all gap positions for all proteins
    pos_with_gaps_for_all_TMDs_all_proteins = []

    # iterate over the dataframe for proteins with an existing list_of_TMDs. acc = uniprot accession.
    for acc in df.index:
        protein_name = df.loc[acc, 'protein_name']
        sys.stdout.write("{} ".format(protein_name))
        sys.stdout.flush()
        # define gap output file path
        if not os.path.exists(df.loc[acc, "fastagap_zip"]):
            logging.info("{} skipped. {} does not exist".format(acc, df.loc[acc, "fastagap_zip"]))
            continue
        # open the csv file that contains the number of valid homologue sequences for each TMD
        fastagap_n_valid_homol_csv = os.path.basename(df.loc[acc, 'fastagap_base'] + 'n_valid_homol.csv')

        try:
            # open zipfile for reading, and csv as a handle
            with zipfile.ZipFile(df.loc[acc, "fastagap_zip"], "r", zipfile.ZIP_DEFLATED) as openzip:
                with openzip.open(fastagap_n_valid_homol_csv) as csv_file_handle:
                    # read as pandas dataframe for "number of valid" homologues
                    df_nv = pd.read_csv(csv_file_handle, header=None)

                    df_nv.columns = ["TMD", "n_valid_homol"]
                    df_nv.set_index("TMD", inplace=True)
                    # filter to remove any TMDs with less than 1 homologue
                    df_nv.query("n_valid_homol >= 1", inplace=True)
                    if df_nv.empty:
                        # skip this protein, none of the TMDs have enough homologues
                        continue

                    with zipfile.ZipFile(df.loc[acc, "fastagap_pos_arrays_zip"], "w", zipfile.ZIP_DEFLATED) as zipout:

                        for TMD in df_nv.index:
                            TMD_query_seq = df.loc[acc,'%s_seq'%TMD]
                            fastagap_aligned_fasta_path = df.loc[acc,'fastagap_base'] + '%s_aligned.fas' % TMD
                            # extract file from zip
                            filename = os.path.basename(fastagap_aligned_fasta_path)
                            path = os.path.dirname(fastagap_aligned_fasta_path)
                            try:
                                openzip.extract(filename, path)
                            except PermissionError:
                                # sometimes windows and python like to hang on to a file, and don't let go
                                # rename to an alternative to avoid permission errors
                                fastagap_aligned_fasta_path = fastagap_aligned_fasta_path[:-4] + "_2.fas"
                                filename = filename[:-4] + "_2.fas"
                                openzip.extract(filename, path)
                            # read as a biopython alignment object
                            alignment = AlignIO.read(fastagap_aligned_fasta_path, "fasta")
                            # convert alignment to a pandas dataframe
                            dfal = pd.DataFrame(np.array(alignment))
                            # create TMD regex search string (e.g. L.*L.*M.*L.*L.*M.*L.*L.*M.*L.*L.*M.*L.*L.*M.*)
                            TMD_regex_ss = utils.create_regex_string(TMD_query_seq)
                            # the query sequence in the alignment is the first sequence
                            query_seq_in_alignment = str(alignment[0].seq)
                            # the regex result should be a tuple of bool, start, end, e.g. (True, 30, 50)
                            regex_result = utils.get_start_and_end_of_TMD_in_query(query_seq_in_alignment, TMD_regex_ss)
                            if isinstance(regex_result, list):
                                contains_TMD_bool, start, end = regex_result
                            # extract the query TMD, which may include gaps
                            query_TMD_in_alignment = query_seq_in_alignment[start:end]
                            # the column labels before the TMD start should be -10, -9, -8, -7, -6, -5, -4, -3, -2, -1, 0.0, for example
                            new_cols_before_TM = list(range(-start, 0))
                            # relabel gap positions, so they lie in between the original TMD index
                            # e.g. L0 I1 I2 Q3 L4 F5 A6 T7 C8 F9 L10 A11 S12 -12.5 L13 M14 F15 F16 W17 E18 P19 I20, where the gap between 12 and 13 is given the position 12.5
                            if "-" in query_TMD_in_alignment:
                                # number of gaps (used to adjust the index as you go along)
                                n_gaps = 0
                                # final column labels for the TM residues
                                cols_TM = []
                                for n, aa in enumerate(query_TMD_in_alignment):
                                    # any gaps result in the original position -1 e.g. S12 -13 L14 will be converted to S12 -12.5 L13
                                    n -= n_gaps
                                    if aa == "-":
                                        # gap positions themselves lie in between the original indices of the surrounding aa e.g. S12 -13 will be converted to S12 -12.5
                                        n = n - 0.5 - n_gaps
                                        n_gaps += 1
                                    cols_TM.append(n)
                                    #sys.stdout.write("{}{} ".format(aa,n))
                            else:
                                # if there are no gaps in the query, the indices are the same
                                cols_TM = list(range(len(TMD_query_seq)))
                            # convert to numpy array
                            cols_TM = np.array(cols_TM)
                            # normalise between 0 and 1, for the TMD region
                            cols_TM_norm_to_1 = cols_TM / cols_TM.max()

                            # start the cols after the TM with 2, so they don't interfere with the TMD region, which goes from 0-1
                            # e.g.  2.0   3.0   4.0   5.0   6.0   7.0   8.0   9.0   10.0  11.0  12.0
                            cols_after_TM = list(range(2, len(query_seq_in_alignment) - end + 2))

                            # the new columns (gap positions) are made by combining the before, TM, and after_TM
                            full_new_cols = new_cols_before_TM + list(cols_TM_norm_to_1) + cols_after_TM
                            # replace original columns
                            if dfal.shape[1] == len(full_new_cols):
                                dfal.columns = full_new_cols
                            else:
                                message = "{} {} skipped. full_new_cols with gap indices is not the same length as the dfa1.columns.".format(acc, TMD)
                                logging.info(message)
                                continue
                            # save the dataframe, with original sequences (before converting non-gaps to np.nan)
                            csv_out = fastagap_aligned_fasta_path[:-4] + "pos_array.csv"
                            dfal.to_csv(csv_out)
                            zipout.write(csv_out, arcname=os.path.basename(csv_out))

                            # convert all gaps to True, and non-gaps to np.nan
                            dfal = dfal.replace("-", np.nan).isnull().replace(False, np.nan)
                            # drop all columns that contain only nan (i.e. only non-gaps)
                            dfal.dropna(how="all", axis=1, inplace=True)
                            # the columns are the list of gap positions for this TMD of this protein
                            pos_with_gaps_for_this_TMD = list(dfal.columns)
                            logging.info("{} {} gap pos : {}".format(acc, TMD, pos_with_gaps_for_this_TMD))
                            # add to the large list for all proteins
                            pos_with_gaps_for_all_TMDs_all_proteins += pos_with_gaps_for_this_TMD

                            # remove the non-zipped temporary files
                            os.remove(fastagap_aligned_fasta_path)
                            os.remove(csv_out)

        except (KeyError, pd.io.common.EmptyDataError, zipfile.BadZipFile):
            logging.info("{} {} corrupt zip file, will be removed".format(acc, df.loc[acc, "fastagap_zip"]))
            os.remove(df.loc[acc, "fastagap_zip"])
            continue

    # save in summaries folder as pickle
    with open(pathdict["gap_fastagap_all_pos_pickle"], "wb") as f:
        pickle.dump(pos_with_gaps_for_all_TMDs_all_proteins, f, protocol=pickle.HIGHEST_PROTOCOL)
    logging.info("pos_with_gaps_for_all_TMDs_all_proteins saved as {}".format(pathdict["gap_fastagap_all_pos_pickle"]))