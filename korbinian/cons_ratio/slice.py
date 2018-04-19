import ast
import os
import pickle
import zipfile
import korbinian
import numpy as np
import pandas as pd
from korbinian import utils as utils
from multiprocessing import Pool
import sys
# import debugging tools
from korbinian.utils import pr, pc, pn, aaa

def run_slice_TMDs_from_homologues(pathdict, s, logging):
    """For a list of proteins, slice TMD sequences from homologues and count gaps.

    Parameters
    ----------
    pathdict : dict
        Dictionary of the key paths and files associated with that List number.
    s : dict
        Settings dictionary extracted from excel settings file.
    logging : logging.Logger
        Logger for printing to console and logfile.

    Saved Files and Figures
    -----------------------
    see slice_1_TMD_from_homol
    """

    logging.info('~~~~~~~~~~~~                 starting run_slice_TMDs_from_homologues                ~~~~~~~~~~~~')
    # if multiprocessing is used, log only to the console
    p_dict_logging = logging if s["use_multiprocessing"] != True else utils.Log_Only_To_Console()
    # get list of accessions that could not be downloaded, and can immediately be excluded
    not_in_homol_db = utils.get_acc_list_from_txt(pathdict["acc_not_in_homol_db_txt"])
    # create list of protein dictionaries to process
    list_p = korbinian.utils.convert_summary_csv_to_input_list(s, pathdict, p_dict_logging, list_excluded_acc=not_in_homol_db)

    # number of processes is the number the settings, or the number of proteins, whichever is smallest
    n_processes = s["multiprocessing_cores"] if s["multiprocessing_cores"] < len(list_p) else len(list_p)

    if s["use_multiprocessing"]:
        with Pool(processes=n_processes) as pool:
            slice_list = pool.map(korbinian.cons_ratio.slice.slice_TMD_1_prot_from_homol, list_p)
            # log the list of protein results (e.g. acc, "simap", True) to the actual logfile, not just the console
        logging.info("\nslice_list : {}".format(slice_list))
    else:
        for p in list_p:
            korbinian.cons_ratio.slice.slice_TMD_1_prot_from_homol(p)
    logging.info('\n~~~~~~~~~~~~                 finished run_slice_TMDs_from_homologues                ~~~~~~~~~~~~')

def slice_TMD_1_prot_from_homol(p):
    """ Slices TMDs from homologues, for a single protein in the list.

     - checks that the homol_df_orig_zip file exists with the full homologue sequences
     - if slice_juxtamembrane_regions is chosen, conducts JM slicing (currently not stable)
     - removes any old, existing fa_cr_sliced_TMDs_zip files
     - creates df_nonTMD_sliced Dataframe to hold the sliced nonTMD regions, based on the indices from all the regex searches
        - for each TMD:
             - identifies and slices out the TMD region from the query, markup and match from each SW alignment
             - the dataframe with sliced sequences for each TMD is added to the fa_cr_sliced_TMDs_zip as PROTEIN_NAME_TM01_sliced_df.pickle, PROTEIN_NAME_TM02_sliced_df.pickle, etc
             - adds the indices for each TMD to the df_nonTMD_sliced Dataframe
        - when each TMD is finished :
            df_nonTMD_sliced uses the indices for each TMD to create the indices for the full nonTMD region (korbinian.cons_ratio.slice.slice_nonTMD_seqs)
            df_nonTMD_sliced contains nonTMD region, as one large slice of query, markup or match from alignment, pieces joined end-to-end
            df_nonTMD_sliced is saved in fa_cr_sliced_TMDs_zip as PROTEIN_NAME_nonTMD_sliced_df.pickle

    Parameters
    ----------
    p : dict
        Protein Dictionary. Contains all input settings, sequences and filepaths related to a single protein.
        Protein-specific data is extracted from one row of the the list summary, e.g. List05_summary.csv, which is read as df.
        p also contains the GENERAL korbinian settings and filepaths for that list (pathdict, s, logging)

        Components
        ----------
        pathdict : dict
            Dictionary of the key paths and files associated with that List number.
        s : dict
            Settings dictionary extracted from excel settings file.
        logging : logging.Logger
            Logger for printing to console and/or logfile.
            If multiprocessing == True, logging.info etc will only print to console.
        p : protein-specific dictionary components
            acc, list_of_TMDs, description, TM01_seq, etc

    Dataframes
    ----------
    dfs
        Dataframe for Sequences
        index = hit_num
        columns = md5, FASTA_gapped_identity, hit_contains_SW_node, organism, X_in_match_seq, disallowed_words_not_in_descr, etc
    df_TMD
        Dataframe for 1 TMD, from 1 protein
        index = hit_num
        columns = 'organism', 'description', 'TM01_in_SW_alignment', 'TM01_start_in_SW_alignment', 'TM01_end_in_SW_alignment', 'TM01_SW_query_seq', 'TM01_SW_markup_seq', 'TM01_SW_match_seq', etc
    df_nonTMD_sliced
        Dataframe for all nonTMD region, from 1 protein
        index = hit_num
        columns = 'nested_tuple_indices_all_nonTMD_regions', 'nonTMD_markup', 'nonTMD_seq_match', 'len_query_align_seq', 'TM01_in_SW_alignment', 'TM01_start_in_SW_alignment', 'TM01_end_in_SW_alignment', 'TM02_in_SW_alignment',  etc

    Saved Files and Figures
    -----------------------
    fa_cr_sliced_TMDs_zip
        df_nonTMD_temp_pickle, e.g. A4ARX1_nonTMD_sliced_df.pickle
        TM_temp_pickle, E.g. A4ARX1_TM01_sliced_df.pickle

    Returns
    -------
    In all cases, a tuple (str, bool, str) is returned.

    if successful:
        return acc, True, "0"
    if not successful:
        return acc, False, "specific warning or reason why protein failed"
    """
    pathdict, s, logging = p["pathdict"], p["s"], p["logging"]
    acc = p["acc"]
    sys.stdout.write("{} ".format(acc))
    sys.stdout.flush()
    protein_name = p['protein_name']
    if not os.path.exists(p['homol_df_orig_zip']):
        warning = "{} Protein skipped. File does not exist".format(p['homol_df_orig_zip'])
        logging.info(warning)
        return acc, False, warning

    if utils.file_is_old(p['homol_df_orig_zip'], s["oldest_acceptable_file_date"]):
        os.remove(p['homol_df_orig_zip']),
        message = "{} skipped, file is old and has been deleted".format(acc)
        logging.info(message)
        return acc, False, message

    dfs = utils.open_df_from_pickle_zip(p['homol_df_orig_zip'], delete_corrupt=True)
    if dfs.empty:
        warning = "{} Protein skipped, file deleted as it is possibly corrupt.".format(p['homol_df_orig_zip'])
        logging.info(warning)
        return acc, False, warning

    list_of_TMDs = ast.literal_eval(p['list_of_TMDs'])

    # create a boolean "p_is_multipass" to show whether protein is multipass (>1TMD) or singlepass (1TMD)
    if "TM02" in list_of_TMDs:
        p_is_multipass = True
    else:
        p_is_multipass = False

    # create counter for number of TMDs with some homologue data
    n_TMDs_w_homol = 0
    fa_cr_sliced_TMDs_zip = p['fa_cr_sliced_TMDs_zip']
    if os.path.isfile(fa_cr_sliced_TMDs_zip):
        if s["overwrite_sliced_homologues"] == True:
            # delete any existing sliced zipfile
            os.remove(fa_cr_sliced_TMDs_zip)
        else:
            warning = "{} skipped, output from slice_TMD_1_prot_from_homol already exists".format(acc)
            logging.info(warning)
            # skip this protein
            return acc, False, warning

    utils.make_sure_path_exists(fa_cr_sliced_TMDs_zip, isfile=True)

    # open new zipfile (NOTE, it must be closed later!!)
    with zipfile.ZipFile(fa_cr_sliced_TMDs_zip, mode="a", compression=zipfile.ZIP_DEFLATED) as homol_sliced_zip:

        # get directory for zip (and other temp files to be transferred)
        homol_dir = os.path.dirname(fa_cr_sliced_TMDs_zip)
        # create a specific dataframe to hold the nonTMD region, including indices (True, start, end) of all the TMD segments
        if "len_query_align_seq" not in dfs.columns:
            warning = "{} len_query_align_seq not in columns, protein skipped for slice_TMD_1_prot_from_homol".format(acc)
            logging.warning(warning)
            #skip protein
            return acc, False, warning

        # add the FASTA_gapped_identity and length of the alignment sequence from dfs, to act as the "end" of all the nonTMD regions
        df_nonTMD_sliced = dfs[['len_query_align_seq', 'SW_query_coverage']].copy()

        # start with an empty dataframe, that will be replaced if there is any data to analyse
        df_TMD = pd.DataFrame()
        for TMD in list_of_TMDs:
            query_TMD_sequence = p['%s_seq' % TMD]
            if type(query_TMD_sequence) == float:
                warning = "{} {} query_TMD_sequence is np.nan! skipping protein.".format(acc, TMD)
                logging.warning(warning)
                # skip protein
                return acc, False, warning
            ## SHOULD NOT BE NECESSARY. OMPdb DATABASE NOW FIXED TO AVOID NAN VALUES IN TM_SEQ
            # if isinstance(query_TMD_sequence, float):
            #     warning = "{} {} query_TMD_sequence is a float ({}), probably np.nan.".format(acc, TMD, query_TMD_sequence)
            #     logging.warning(warning)
            #     return acc, False, warning
            df_TMD = korbinian.cons_ratio.slice.slice_1_TMD_from_homol(acc, TMD, query_TMD_sequence, dfs, s, logging)
            if df_TMD.empty:
                warning = "{} {} df_TMD.empty, probably number_of_rows_containing_data == 0".format(acc, TMD, query_TMD_sequence)
                logging.warning(warning)
                # skip TMD, as number_of_rows_containing_data == 0
                # here I really should skip the protein too. It's tempting to use goto :). "from goto import goto" (http://entrian.com/goto/)
                return acc, False, warning
            n_TMDs_w_homol += 1
            # transfer the columns with indices across to the df_nonTMD_sliced
            cols = ['%s_in_SW_alignment' % TMD, '%s_start_in_SW_alignment' % TMD, '%s_end_in_SW_alignment' % TMD]
            for col in cols:
                df_nonTMD_sliced[col] = df_TMD[col]

            TM_temp_pickle = os.path.join(homol_dir, "{}_{}_sliced_df.pickle".format(protein_name, TMD))
            with open(TM_temp_pickle, "wb") as f:
                pickle.dump(df_TMD, f, protocol=pickle.HIGHEST_PROTOCOL)
            homol_sliced_zip.write(TM_temp_pickle, arcname=os.path.basename(TM_temp_pickle))
            os.remove(TM_temp_pickle)
            sys.stdout.write(".")
            sys.stdout.flush()

        if df_TMD.empty:
            # skip protein, as number_of_rows_containing_data == 0 for at least one TMD (or at least the last TMD)
            warning = "{} skipped, number_of_rows_containing_data == 0 for at least one TMD".format(acc)
            logging.info(warning)
            return acc, False, warning

        df_nonTMD_sliced = korbinian.cons_ratio.slice.slice_nonTMD_seqs(dfs, df_nonTMD_sliced, list_of_TMDs)
        if df_nonTMD_sliced.empty:
            warning = "{} df_nonTMD_sliced is empty, probably this means no homologues contain all TMDs".format(acc)
            logging.warning(warning)
            #skip protein
            return acc, False, warning

        if s["slice_juxtamembrane_regions"] == True:
            for TMD in list_of_TMDs:
                ########################################################################################
                #                                                                                      #
                #        Define juxtamembrane regions associated with each TMD  [AAIMON]               #
                #                                                                                      #
                ########################################################################################
                # convert the tuple of (True, 32, 53) into separate dataframes.
                # http://stackoverflow.com/questions/29550414/how-to-split-column-of-tuples-in-pandas-dataframe

                if p_is_multipass:
                    next_TMD = "TM{:02d}".format(int(TMD[2:]) + 1)
                    prev_TMD = "TM{:02d}".format(int(TMD[2:]) - 1)
                    #df_next_TMD = df_TMD = korbinian.cons_ratio.slice.slice_1_TMD_from_homol(acc, next_TMD, query_TMD_sequence, dfs, s, logging)
                    #if TMD != "TM01":
                    #    df_prev_TMD = df_TMD = korbinian.cons_ratio.slice.slice_1_TMD_from_homol(acc, prev_TMD, query_TMD_sequence, dfs, s, logging)

                last_TMD_of_acc = list_of_TMDs[-1]

                if TMD == "TM01":
                    # np.where syntax: np.where(boolean_query, value_if_query_true, value_if_query_false)
                    # @RJ, If TM01_start_in_SW_alignment is not an integer above 0, replaces with np.nan?
                    df_nonTMD_sliced['start_juxta_before_TM01'] = np.where(df_nonTMD_sliced['TM01_start_in_SW_alignment'] > 0, 0, np.nan)
                    # if the TM01_start_in_SW_alignment is 0, there is no JM region N-terminal to the TMD, therefore replace end_juxta_before_TM01 with np.nan, otherwise use TM01_start_in_SW_alignment
                    df_nonTMD_sliced['end_juxta_before_TM01'] = np.where(df_nonTMD_sliced['TM01_start_in_SW_alignment'] == 0, np.nan, df_nonTMD_sliced['TM01_start_in_SW_alignment'])
                    # set the start of the juxta as the end of the TMD
                    df_nonTMD_sliced['start_juxta_after_TM01'] = df_nonTMD_sliced['TM01_end_in_SW_alignment']
                    # if there is only one TMD (search for TM02 rather than measuring length of list, in case of signal peptides)
                    if p_is_multipass:
                        # open up the dataframes of the next and previous TMD
                        # define the end_juxta_after_TM01 as the TM01 end + half of the TM01_to_TM02 JM region
                        # NOTE, due to np.nan this is a float. will be converted to integers later
                        df_nonTMD_sliced['end_juxta_after_TM01'] = df_nonTMD_sliced["TM01_end_in_SW_alignment"] + ((df_nonTMD_sliced["TM02_start_in_SW_alignment"] - df_nonTMD_sliced["TM01_end_in_SW_alignment"]) / 2)

                        # RJ original
                        ## problem('df_nonTMD_sliced["TM02_start_in_SW_alignment"] cannot exist yet, because the script iterates through the TMDs one at a time')
                        # df_nonTMD_sliced['end_juxta_after_TM01'] = df_nonTMD_sliced["TM01_end_in_SW_alignment"] + ((df_nonTMD_sliced["TM02_start_in_SW_alignment"] - df_nonTMD_sliced["TM01_end_in_SW_alignment"]) / 2).apply(lambda x: int(x) if not np.isnan(x) else np.nan)

                        # RJ commented out
                        # df_nonTMD_sliced['seq_juxta_after_TM01_in_query'] = df_nonTMD_sliced[df_nonTMD_sliced['start_juxta_after_TM01'].notnull()].apply(utils.slice_juxta_after_TMD_in_query, args = (TMD,), axis=1)
                        # df_nonTMD_sliced['seq_juxta_after_TM01_in_match'] = df_nonTMD_sliced[df_nonTMD_sliced['end_juxta_after_TM01'].notnull()].apply(utils.slice_juxta_after_TMD_in_match, args = (TMD,), axis=1)

                    else:
                        # if there is only one TMD, TM01 == last_TMD_of_acc
                        # @RJ replace with df_nonTMD_sliced['end_juxta_after_TM01'] = df_nonTMD_sliced['len_query_align_seq'] and use dropna to avoid nans later?
                        df_nonTMD_sliced['end_juxta_after_TM01'] = np.where(utils.isNaN(df_nonTMD_sliced['start_juxta_after_TM01']) == True, np.nan, df_nonTMD_sliced['len_query_align_seq'])

                # the analysis is slow, so don't repeat TM01 if there is only one TM helix in the protein
                if p_is_multipass:
                    if not TMD == "TM01" and not TMD == last_TMD_of_acc:
                        df_nonTMD_sliced = juxta_function_1(df_nonTMD_sliced, TMD)
                        # df_nonTMD_sliced['start_juxta_after_%s'%TMD] = np.where(utils.isNaN(df_nonTMD_sliced['TM%.2d_start_in_SW_alignment'%(int(TMD[2:])+1)])==True,np.nan,df_nonTMD_sliced['%s_end_in_SW_alignment'%TMD])
                        # df_nonTMD_sliced['end_juxta_before_%s'%TMD] = np.where(df_nonTMD_sliced["%s_start_in_SW_alignment"%TMD]!=0,df_nonTMD_sliced["%s_start_in_SW_alignment"%TMD],np.nan)
                        # df_nonTMD_sliced['end_juxta_after_%s'%TMD] = df_nonTMD_sliced["%s_end_in_SW_alignment"%TMD]+((df_nonTMD_sliced["TM%.2d_start_in_SW_alignment"%(int(TMD[2:])+1)]-df_nonTMD_sliced["%s_end_in_SW_alignment"%TMD])/2).apply(lambda x :int(x) if not np.isnan(x) else np.nan)
                        # df_nonTMD_sliced['start_juxta_before_%s'%TMD] = np.where(df_nonTMD_sliced["end_juxta_after_TM%.2d"%(int(TMD[2:])-1)] == df_nonTMD_sliced['end_juxta_before_%s'%TMD] ,df_nonTMD_sliced["end_juxta_after_TM%.2d"%(int(TMD[2:])-1)],df_nonTMD_sliced["end_juxta_after_TM%.2d"%(int(TMD[2:])-1)])
                        # df_nonTMD_sliced['seq_juxta_after_%s_in_query'%TMD] = df_nonTMD_sliced[df_nonTMD_sliced['start_juxta_after_%s'%TMD].notnull()].apply(utils.slice_juxta_after_TMD_in_query, args = (TMD,), axis=1)
                        # df_nonTMD_sliced['seq_juxta_after_%s_in_match'%TMD] = df_nonTMD_sliced[df_nonTMD_sliced['end_juxta_after_%s'%TMD].notnull()].apply(utils.slice_juxta_after_TMD_in_match, args = (TMD,), axis=1)

                    if TMD == last_TMD_of_acc:
                        df_nonTMD_sliced['start_juxta_before_%s' % TMD] = df_nonTMD_sliced['end_juxta_after_%s' % prev_TMD]
                        df_nonTMD_sliced['end_juxta_before_%s' % TMD] = df_nonTMD_sliced['%s_start_in_SW_alignment' % TMD]
                        df_nonTMD_sliced['start_juxta_after_%s' % TMD] = np.where(
                            df_nonTMD_sliced['%s_end_in_SW_alignment' % TMD] == df_nonTMD_sliced['len_query_align_seq'], np.nan,
                            df_nonTMD_sliced['%s_end_in_SW_alignment' % TMD])
                        df_nonTMD_sliced['end_juxta_after_%s' % TMD] = np.where(utils.isNaN(df_nonTMD_sliced['start_juxta_after_%s' % TMD]) == True, np.nan,
                                                                   df_nonTMD_sliced['len_query_align_seq'])
                        # df_nonTMD_sliced['seq_juxta_after_%s_in_query'%TMD] = df_nonTMD_sliced[df_nonTMD_sliced['start_juxta_after_%s'%TMD].notnull()].apply(utils.slice_juxta_after_TMD_in_query, args = (TMD,), axis=1)
                        # df_nonTMD_sliced['seq_juxta_after_%s_in_query'%TMD] = df_nonTMD_sliced.query_align_seq[int(df_nonTMD_sliced['start_juxta_after_TM10']):int(df_nonTMD_sliced['end_juxta_after_TM10'])]
                        # df_nonTMD_sliced['seq_juxta_after_%s_in_match'%TMD] =
                else:
                    # the end_juxta_after_TM01 is already defined, nothing else needs to be done for the single-pass proteins
                    pass

                last_TMD_of_acc = list_of_TMDs[-1]
                index_juxta = df_nonTMD_sliced['start_juxta_before_%s' % TMD].notnull().index
                q = np.array(dfs.loc[index_juxta, "query_align_seq"])
                st = np.array(df_nonTMD_sliced.loc[index_juxta, 'start_juxta_before_%s' % TMD])
                st = st.astype(int)
                en = np.array(df_nonTMD_sliced.loc[index_juxta, 'end_juxta_before_%s' % TMD])
                en = en.astype(int)
                q_sliced = [q[i][st[i]:en[i]] for i in range(len(q))]
                df_nonTMD_sliced['seq_juxta_before_%s_in_query' % TMD] = pd.Series(q_sliced, index=index_juxta)
                m = np.array(dfs.loc[index_juxta, "match_align_seq"])
                m_sliced = [m[i][st[i]:en[i]] for i in range(len(m))]
                df_nonTMD_sliced['seq_juxta_before_%s_in_match' % TMD] = pd.Series(m_sliced, index=index_juxta)

                #df_nonTMD_sliced['seq_juxta_before_%s_in_query' % TMD] = df_nonTMD_sliced[df_nonTMD_sliced['start_juxta_before_%s' % TMD].notnull()].apply(utils.slice_juxta_before_TMD_in_query, args=(TMD,), axis=1)
                #df_nonTMD_sliced['seq_juxta_before_%s_in_match' % TMD] = df_nonTMD_sliced[df_nonTMD_sliced['start_juxta_before_%s' % TMD].notnull()].apply(utils.slice_juxta_before_TMD_in_match, args=(TMD,), axis=1)
                if not TMD == last_TMD_of_acc:
                    index_juxta = df_nonTMD_sliced['end_juxta_after_%s' % TMD].notnull().index
                    st = np.array(df_nonTMD_sliced.loc[index_juxta, 'start_juxta_after_%s' % TMD])
                    st = st.astype(int)
                    en = np.array(df_nonTMD_sliced.loc[index_juxta, 'end_juxta_after_%s' % TMD])
                    en = en.astype(int)
                    q_sliced = [q[i][st[i]:en[i]] for i in range(len(q))]
                    df_nonTMD_sliced['seq_juxta_after_%s_in_query' % TMD] = pd.Series(q_sliced, index=index_juxta)
                    m_sliced = [m[i][st[i]:en[i]] for i in range(len(m))]
                    df_nonTMD_sliced['seq_juxta_after_%s_in_match' % TMD] = pd.Series(m_sliced, index=index_juxta)
                    #df_nonTMD_sliced['seq_juxta_after_%s_in_query' % TMD] = df_nonTMD_sliced[df_nonTMD_sliced['end_juxta_after_%s' % TMD].notnull()].apply(utils.slice_juxta_after_TMD_in_query, args=(TMD,), axis=1)
                    #df_nonTMD_sliced['seq_juxta_after_%s_in_match' % TMD] = df_nonTMD_sliced[df_nonTMD_sliced['end_juxta_after_%s' % TMD].notnull()].apply(utils.slice_juxta_after_TMD_in_match, args=(TMD,), axis=1)
                else:
                    index_juxta = df_nonTMD_sliced['start_juxta_after_%s' % TMD].notnull().index
                    st = np.array(df_nonTMD_sliced.loc[index_juxta, 'start_juxta_after_%s' % TMD])
                    st = st.astype(int)
                    en = np.array(df_nonTMD_sliced.loc[index_juxta, 'end_juxta_after_%s' % TMD])
                    en = en.astype(int)
                    q_sliced = [q[i][st[i]:en[i]] for i in range(len(q))]
                    df_nonTMD_sliced['seq_juxta_after_%s_in_query' % TMD] = pd.Series(q_sliced, index=index_juxta)
                    m_sliced = [m[i][st[i]:en[i]] for i in range(len(m))]
                    df_nonTMD_sliced['seq_juxta_after_%s_in_match' % TMD] = pd.Series(m_sliced, index=index_juxta)
                    # df_nonTMD_sliced['seq_juxta_after_%s_in_query' % TMD] = np.nan
                    # df_nonTMD_sliced['seq_juxta_after_%s_in_match' % TMD] = np.nan
                    # for hit in df_nonTMD_sliced.index:
                    #     if not utils.isNaN(df_nonTMD_sliced['start_juxta_after_%s' % TMD])[hit]:
                    #         # altered to .loc rather than ['seq_juxta_after_%s_in_match'%TMD][hit] after SettingWithCopyWarning
                    #         df_nonTMD_sliced.loc[hit, 'seq_juxta_after_%s_in_match' % TMD] = df_nonTMD_sliced.match_align_seq[hit][int(df_nonTMD_sliced.loc[hit, "start_juxta_after_%s" % TMD]):int(df_nonTMD_sliced.loc[hit, "end_juxta_after_%s" % TMD])]
                    #         df_nonTMD_sliced.loc[hit, 'seq_juxta_after_%s_in_query' % TMD] = df_nonTMD_sliced.query_align_seq[hit][int(df_nonTMD_sliced.loc[hit, "start_juxta_after_%s" % TMD]):int(df_nonTMD_sliced.loc[hit, "end_juxta_after_%s" % TMD])]

        df_nonTMD_temp_pickle = os.path.join(homol_dir, "{}_nonTMD_sliced_df.pickle".format(protein_name))
        with open(df_nonTMD_temp_pickle, "wb") as f:
            pickle.dump(df_nonTMD_sliced, f, protocol=pickle.HIGHEST_PROTOCOL)
        homol_sliced_zip.write(df_nonTMD_temp_pickle, arcname=os.path.basename(df_nonTMD_temp_pickle))
        os.remove(df_nonTMD_temp_pickle)
        return acc, True, "0"

def slice_1_TMD_from_homol(acc, TMD, query_TMD_sequence, dfs, s, logging):
    """For a single protein, slice TMD sequences from homologues and count gaps.

    Slices out the TMD region for each homologue.

    The TMD region is found through a regex match to the query sequence.
        regex search string (e.g. L.*L.*M.*L.*L.*M.*L.*L.*M.*L.*L.*M.*L.*L.*M.*)
        If the full TMD is not in the alignment, it will NOT BE FOUND!

    Note the output of the regex search is a list with three components.
        1) a match (e.g. True),
        2) the start of the match (e.g. 124),
        3) the end of the match (e.g. 144)

    This function also slices out the TMD plus surrounding sequence (usually 10 residues each side, depending on the settings file)

    Parameters
    ----------
    acc : str
        Accession number for protein
    TMD : str
        String denoting transmembrane domain number (e.g. "TM01")
    df : pd.DataFrame
        Dataframe containing the list of proteins for analysis, usually containing data from UniProt, with each row representing a different protein.
    dfs : pd.Dataframe
        Dataframe for Sequences (dfs).
        This is the dataframe containing the full homologue sequences from the BLAST-like data analysis.
        It is stored in in homol_df_orig_zip, along with the csv containing the pretty version of homologues.
    s : dict
        Settings dictionary extracted from excel settings file.
    logging : logging.Logger
        Logger for printing to console and logfile.

    Returns
    -------
    df_TMD : pd.DataFrame
        Dataframe containing the sliced TMD sequence, and various calculated features related to that TMD, for example, the number of gaps.
        Will be saved in PROTEIN_NAME_fa_cr_sliced_TMDs.zip as PROTEIN_NAME_TM01_sliced_df.pickle (for TM01, for example)
        There will be a separate df_TMD, and a separate PROTEIN_NAME_TM01_sliced_df.pickle for each TMD or region.
    """
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
    #query_TMD_sequence = df.loc[acc, '%s_seq' % TMD]
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
    if number_of_rows_containing_data != 0:
        #len_query_TMD = len(df.loc[acc, '%s_seq' % TMD])
        # apply the slicing function to the homologues
        # df_TMD = korbinian.cons_ratio.slice_1_TMD_from_homol(TMD, len_query_TMD, df_TMD,
        #                                                               s["cr_max_n_gaps_in_query_TMD"],
        #                                                               s["cr_max_n_gaps_in_match_TMD"])
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
        n_aa_before_tmd = s["n_aa_before_tmd"]
        n_aa_after_tmd = s["n_aa_after_tmd"]

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
        # slice the query sequence. Could be useful for fasta_gap analysis.
        df_TMD['%s_SW_query_seq_plus_surr'%TMD] = dfs.apply(utils.slice_SW_query_TMD_seq_plus_surr, args=(TMD,), axis=1)

        #and the same for the TMD + surrounding sequence, useful to examine the TMD interface
        # NOT DEEMED NECESSARY. WHY WOULD YOU NEED TO SLICE QUERY OR MARKUP + SURROUNDING?
        # df_TMD['%s_SW_markup_seq_plus_surr'%TMD] = dfs.apply(utils.slice_SW_markup_TMD_plus_surr, args=(TMD,), axis=1)
        ########################################################################################
        #                                                                                      #
        #        count number of gaps in the query and match TMD aligned sequence              #
        #                                                                                      #
        ########################################################################################
        """This is used as a filter in filter_and_save_fasta, therefore is conducted earlier in the slicing function. """
        # count the number of gaps in the query and match sequences
        df_TMD['%s_SW_query_num_gaps' % TMD] = df_TMD['%s_SW_query_seq' % TMD].str.count("-")
        df_TMD['%s_SW_match_num_gaps' % TMD] = df_TMD['%s_SW_match_seq' % TMD].str.count("-")

        ########################################################################################
        #                                                                                      #
        #     calculate the average number of gaps per residue in the TMD alignment            #
        #             (number of gaps)/(length of sequence excluding gaps)                     #
        #                                                                                      #
        ########################################################################################
        df_TMD['%s_SW_q_gaps_per_q_residue' % TMD] = df_TMD['%s_SW_query_num_gaps' % TMD].dropna() / len(query_TMD_sequence)

        # calculate hydrophobicity
        df_TMD['%s_SW_match_lipo' % TMD] = df_TMD['%s_SW_match_seq'%TMD].dropna().apply(lambda x: utils.calc_lipophilicity(x))

    else:
        logging.info('{} does not have any valid homologues for {}. Re-downloading simap homologue XML may be necessary.'.format(acc, TMD))
        df_TMD = pd.DataFrame()

    return df_TMD

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
        if df_nonTMD_sliced.empty:
            # there are no homologues with all TMDs in the Smith Waterman alignment
            # return an empty dataframe
            return pd.DataFrame()

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

def juxta_function_1(dfs, TMD):
    """
    Parameters
    ----------
    dfs : pd.Dataframe
        Dataframe for Sequences (dfs). This is the dataframe containing the full homologue sequences from the BLAST-like data analysis.
    TMD : str
        String denoting transmembrane domain number (e.g. "TM01")

    Returns
    -------

    """
    next_TM = "TM{:02d}".format(int(TMD[2:]) + 1)
    prev_TM = "TM{:02d}".format(int(TMD[2:])-1)

    dfs['start_juxta_after_%s'%TMD] = np.where(utils.isNaN(dfs['{}_start_in_SW_alignment'.format(next_TM)]) == True, np.nan, dfs['%s_end_in_SW_alignment'%TMD])

    dfs['end_juxta_before_%s'%TMD] = np.where(dfs["%s_start_in_SW_alignment"%TMD] != 0, dfs["%s_start_in_SW_alignment"%TMD], np.nan)

    dfs['end_juxta_after_%s'%TMD] = dfs["%s_end_in_SW_alignment"%TMD] + ((dfs["{}_start_in_SW_alignment".format(next_TM)]-dfs["%s_end_in_SW_alignment"%TMD])/2).apply(lambda x :int(x) if not np.isnan(x) else np.nan)

    dfs['start_juxta_before_%s'%TMD] = np.where(dfs["end_juxta_after_{}".format(prev_TM)] == dfs['end_juxta_before_%s'%TMD], dfs["end_juxta_after_{}".format(prev_TM)],dfs["end_juxta_after_{}".format(prev_TM)])
    return dfs