import ast
import os
import pickle
import zipfile

import korbinian
import numpy as np
import pandas as pd
from korbinian import utils as utils

def slice_TMD_homol_and_count_gaps(acc, TMD, query_TMD_sequence, dfs, set_, logging):
    """Slice TMD sequences from homologues and count gaps.

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
    set_ : dict
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
        len_query_TMD = len(query_TMD_sequence)
        # for the alignment length, take the smallest value from the length of query or match
        # this will exclude gaps from the length in the following calculations, preventing false "low conservation" where the query TMD is much longer than the match TMD)
        # note that for most calculations this is somewhat redundant, because the max number of acceptable gaps in sequence is probable ~2
        # use the mask function (faster than a lambda function) to replace any lengths larger than the query, with the query
        df_TMD['%s_SW_align_len' % TMD] = df_TMD['%s_SW_m_seq_len' % TMD].mask(df_TMD['%s_SW_m_seq_len' % TMD] > len_query_TMD, len_query_TMD)
        # create a boolean column that allows filtering by the accepted number of gaps, according to the settings file
        df_TMD['%s_cr_SW_query_acceptable_n_gaps'%TMD] = df_TMD['%s_SW_query_num_gaps'%TMD] <= set_["cr_max_n_gaps_in_query_TMD"]
        df_TMD['%s_cr_SW_match_acceptable_n_gaps'%TMD] = df_TMD['%s_SW_match_num_gaps'%TMD] <= set_["cr_max_n_gaps_in_match_TMD"]
        # count identical residues between query and match TMDs by counting the number of pipes in the markup string
        # NOTE THAT df_TMD['%s_SW_markup_seq'%TMD].str.count('|') DOES NOT WORK, as "|" has a function in regex and needs to be escaped
        df_TMD['%s_SW_num_ident_res'%TMD] = df_TMD['%s_SW_markup_seq'%TMD].str.count('\|')
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
        # calculate the average number of gaps per residue in the TMD alignment
        # (number of gaps)/(length of sequence excluding gaps)
        df_TMD['%s_SW_q_gaps_per_q_residue'%TMD] = df_TMD['%s_SW_query_num_gaps'%TMD].dropna() / len_query_TMD
    else:
        logging.info('%s does not have any valid homologues for %s. '
                     'Re-downloading simap homologue XML may be necessary.' % (acc, TMD))
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

def slice_TMDs_from_homologues(p):
    """ Slices TMDs from homologues

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
    """
    pathdict, set_, logging = p["pathdict"], p["set_"], p["logging"]
    acc = p["acc"]
    # df = pd.read_csv(pathdict["list_summary_csv"], sep=",", quoting=csv.QUOTE_NONNUMERIC, index_col=0, low_memory=False)
    # #iterate over the dataframe for proteins with an existing list_of_TMDs. Note that acc = uniprot accession here.
    # for acc in df.loc[df['list_of_TMDs'].notnull()].loc[df['list_of_TMDs'] != 'nan'].index:
    protein_name = p['protein_name']
    print(acc, end=", ", flush=True)
    if not os.path.exists(p['homol_df_orig_zip']):
        warning = "{} Protein skipped. File does not exist".format(p['homol_df_orig_zip'])
        logging.info(warning)
        return acc, False, warning

    dfs = utils.open_df_from_pickle_zip(p['homol_df_orig_zip'], delete_corrupt=True)
    if dfs.empty:
        warning = "{} Protein skipped, file deleted as it is possibly corrupt.".format(p['homol_df_orig_zip'])
        logging.info(warning)
        return acc, False, warning

    if set_["slice_juxtamembrane_regions"]:
        ########################################################################################
        #                                                                                      #
        #        Define juxtamembrane regions associated with each TMD  [AAIMON]               #
        #                                                                                      #
        ########################################################################################
        # convert the tuple of (True, 32, 53) into separate dataframes.
        # http://stackoverflow.com/questions/29550414/how-to-split-column-of-tuples-in-pandas-dataframe

        # dfs_with_match = dfs['%s_start_end_list_in_SW_alignment' % TMD].dropna()
        # df_match_start_end = pd.DataFrame(dfs_with_match.values.tolist())
        # df_match_start_end.index = dfs_with_match.index
        # dfs["%s_start_in_SW_alignment" % TMD] = df_match_start_end[1]
        # dfs["%s_end_in_SW_alignment" % TMD] = df_match_start_end[2]
        list_of_TMDs = ["temp", "list"]
        TMD = "temp, needs checking"
        last_TMD_of_acc = list_of_TMDs[-1]

        if set_["slice_juxtamembrane_regions"] == True:
            if TMD == "TM01":
                dfs['start_juxta_before_TM01'] = np.where(dfs['TM01_start_in_SW_alignment'] > 0, 0, np.nan)
                dfs['end_juxta_before_TM01'] = np.where(dfs['TM01_start_in_SW_alignment'] == 0, np.nan,
                                                        dfs['TM01_start_in_SW_alignment'])
                dfs['start_juxta_after_TM01'] = dfs['TM01_end_in_SW_alignment']
                if len(list_of_TMDs) == 1:
                    # if there is only one TMD, TM01 == last_TMD_of_acc
                    dfs['end_juxta_after_TM01'] = np.where(utils.isNaN(dfs['start_juxta_after_TM01']) == True, np.nan,
                                                           dfs['len_query_align_seq'])
                elif len(list_of_TMDs) > 1:
                    # problem('dfs["TM02_start_in_SW_alignment"] cannot exist yet, because the script iterates through the TMDs one at a time')
                    dfs['end_juxta_after_TM01'] = dfs["TM01_end_in_SW_alignment"] + (
                        (dfs["TM02_start_in_SW_alignment"] - dfs["TM01_end_in_SW_alignment"]) / 2).apply(lambda x: int(x) if not np.isnan(x) else np.nan)
                    # dfs['seq_juxta_after_TM01_in_query'] = dfs[dfs['start_juxta_after_TM01'].notnull()].apply(utils.slice_juxta_after_TMD_in_query, args = (TMD,), axis=1)
                    # dfs['seq_juxta_after_TM01_in_match'] = dfs[dfs['end_juxta_after_TM01'].notnull()].apply(utils.slice_juxta_after_TMD_in_match, args = (TMD,), axis=1)

            # the analysis is slow, so don't repeat TM01 if there is only one TM helix in the protein
            if len(list_of_TMDs) > 1:
                if not TMD == "TM01" and not TMD == last_TMD_of_acc:
                    dfs = juxta_function_1(dfs, TMD)
                    # dfs['start_juxta_after_%s'%TMD] = np.where(utils.isNaN(dfs['TM%.2d_start_in_SW_alignment'%(int(TMD[2:])+1)])==True,np.nan,dfs['%s_end_in_SW_alignment'%TMD])
                    # dfs['end_juxta_before_%s'%TMD] = np.where(dfs["%s_start_in_SW_alignment"%TMD]!=0,dfs["%s_start_in_SW_alignment"%TMD],np.nan)
                    # dfs['end_juxta_after_%s'%TMD] = dfs["%s_end_in_SW_alignment"%TMD]+((dfs["TM%.2d_start_in_SW_alignment"%(int(TMD[2:])+1)]-dfs["%s_end_in_SW_alignment"%TMD])/2).apply(lambda x :int(x) if not np.isnan(x) else np.nan)
                    # dfs['start_juxta_before_%s'%TMD] = np.where(dfs["end_juxta_after_TM%.2d"%(int(TMD[2:])-1)] == dfs['end_juxta_before_%s'%TMD] ,dfs["end_juxta_after_TM%.2d"%(int(TMD[2:])-1)],dfs["end_juxta_after_TM%.2d"%(int(TMD[2:])-1)])
                    # dfs['seq_juxta_after_%s_in_query'%TMD] = dfs[dfs['start_juxta_after_%s'%TMD].notnull()].apply(utils.slice_juxta_after_TMD_in_query, args = (TMD,), axis=1)
                    # dfs['seq_juxta_after_%s_in_match'%TMD] = dfs[dfs['end_juxta_after_%s'%TMD].notnull()].apply(utils.slice_juxta_after_TMD_in_match, args = (TMD,), axis=1)

                if TMD == last_TMD_of_acc:
                    dfs['start_juxta_before_%s' % TMD] = dfs['end_juxta_after_TM%.2d' % (int(TMD[2:]) - 1)]
                    dfs['end_juxta_before_%s' % TMD] = dfs['%s_start_in_SW_alignment' % TMD]
                    dfs['start_juxta_after_%s' % TMD] = np.where(
                        dfs['%s_end_in_SW_alignment' % TMD] == dfs['len_query_align_seq'], np.nan,
                        dfs['%s_end_in_SW_alignment' % TMD])
                    dfs['end_juxta_after_%s' % TMD] = np.where(utils.isNaN(dfs['start_juxta_after_%s' % TMD]) == True, np.nan,
                                                               dfs['len_query_align_seq'])
                    # dfs['seq_juxta_after_%s_in_query'%TMD] = dfs[dfs['start_juxta_after_%s'%TMD].notnull()].apply(utils.slice_juxta_after_TMD_in_query, args = (TMD,), axis=1)
                    # dfs['seq_juxta_after_%s_in_query'%TMD] = dfs.query_align_seq[int(dfs['start_juxta_after_TM10']):int(dfs['end_juxta_after_TM10'])]
                    # dfs['seq_juxta_after_%s_in_match'%TMD] =

            last_TMD_of_acc = list_of_TMDs[-1]
            dfs['seq_juxta_before_%s_in_query' % TMD] = dfs[dfs['start_juxta_before_%s' % TMD].notnull()].apply(
                utils.slice_juxta_before_TMD_in_query, args=(TMD,), axis=1)
            dfs['seq_juxta_before_%s_in_match' % TMD] = dfs[dfs['start_juxta_before_%s' % TMD].notnull()].apply(
                utils.slice_juxta_before_TMD_in_match, args=(TMD,), axis=1)
            if not TMD == last_TMD_of_acc:
                dfs['seq_juxta_after_%s_in_query' % TMD] = dfs[dfs['end_juxta_after_%s' % TMD].notnull()].apply(
                    utils.slice_juxta_after_TMD_in_query, args=(TMD,), axis=1)
                dfs['seq_juxta_after_%s_in_match' % TMD] = dfs[dfs['end_juxta_after_%s' % TMD].notnull()].apply(
                    utils.slice_juxta_after_TMD_in_match, args=(TMD,), axis=1)
            else:
                dfs['seq_juxta_after_%s_in_query' % TMD] = np.nan
                dfs['seq_juxta_after_%s_in_match' % TMD] = np.nan
                for hit in dfs.index:
                    if not utils.isNaN(dfs['start_juxta_after_%s' % TMD])[hit]:
                        # altered to .loc rather than ['seq_juxta_after_%s_in_match'%TMD][hit] after SettingWithCopyWarning
                        dfs.loc[hit, 'seq_juxta_after_%s_in_match' % TMD] = dfs.match_align_seq[hit][
                                                                            int(dfs.loc[hit, "start_juxta_after_%s" % TMD]):int(
                                                                                dfs.loc[hit, "end_juxta_after_%s" % TMD])]
                        dfs.loc[hit, 'seq_juxta_after_%s_in_query' % TMD] = dfs.query_align_seq[hit][
                                                                            int(dfs.loc[hit, "start_juxta_after_%s" % TMD]):int(
                                                                                dfs.loc[hit, "end_juxta_after_%s" % TMD])]

    # # get length of seq. Previously this was a lambda function.
    # if True in dfs['hit_contains_SW_node'].tolist():
    #     dfs['len_query_align_seq'] = dfs['query_align_seq'].str.len()
    # else:
    #     dfs['query_align_seq'] = np.nan
    #     dfs['len_query_align_seq'] = np.nan

    #     try:
    #         dfs['len_query_align_seq'] = dfs['query_align_seq'].str.len()
    #     except KeyError:
    #         pass
            # dataframe does not contain query_align_seq,
            # which means that the XML file is probably damaged somehow
    #         logging.warning('SIMAP_csv_from_XML seems to be damaged for %s' % protein_name)
    #         dfs['query_align_seq'] = ''
    #         dfs['len_query_align_seq'] = 0

    list_of_TMDs = ast.literal_eval(p['list_of_TMDs'])
    # create counter for number of TMDs with some homologue data
    n_TMDs_w_homol = 0
    fa_cr_sliced_TMDs_zip = p['fa_cr_sliced_TMDs_zip']
    if os.path.isfile(fa_cr_sliced_TMDs_zip):
        if set_["overwrite_sliced_homologues"] == True:
            # delete any existing sliced zipfile
            os.remove(fa_cr_sliced_TMDs_zip)
        else:
            warning = "{} skipped, output from slice_TMDs_from_homologues already exists".format(acc)
            logging.info(warning)
            # skip this protein
            return acc, False, warning

    # open new zipfile (NOTE, it must be closed later!!)
    with zipfile.ZipFile(fa_cr_sliced_TMDs_zip, mode="a", compression=zipfile.ZIP_DEFLATED) as homol_sliced_zip:

        # get directory for zip (and other temp files to be transferred)
        homol_dir = os.path.dirname(fa_cr_sliced_TMDs_zip)
        # create a specific dataframe to hold the nonTMD region, including indices (True, start, end) of all the TMD segments
        if "len_query_align_seq" not in dfs.columns:
            warning = "{} len_query_align_seq not in columns, protein skipped for slice_TMDs_from_homologues".format(acc)
            logging.warning(warning)
            #skip protein
            return acc, False, warning

        # add the FASTA_gapped_identity and length of the alignment sequence from dfs, to act as the "end" of all the nonTMD regions
        df_nonTMD_sliced = dfs[['len_query_align_seq']].copy()
        # start with an empty dataframe, that will be replaced if there is any data to analyse
        df_TMD = pd.DataFrame()
        for TMD in list_of_TMDs:
            query_TMD_sequence = p['%s_seq' % TMD]
            ## SHOULD NOT BE NECESSARY. OMPdb DATABASE NOW FIXED TO AVOID NAN VALUES IN TM_SEQ
            # if isinstance(query_TMD_sequence, float):
            #     warning = "{} {} query_TMD_sequence is a float ({}), probably np.nan.".format(acc, TMD, query_TMD_sequence)
            #     logging.warning(warning)
            #     return acc, False, warning
            df_TMD = korbinian.cons_ratio.singleprotein.slice.slice_TMD_homol_and_count_gaps(acc, TMD, query_TMD_sequence, dfs, set_, logging)
            if df_TMD.empty:
                warning = "{} {} df_TMD.empty, probably number_of_rows_containing_data == 0".format(acc, TMD, query_TMD_sequence)
                logging.warning(warning)
                # skip TMD, as number_of_rows_containing_data == 0
                # here I really should skip the protein too. It's tempting to use goto :). "from goto import goto" (http://entrian.com/goto/)
                continue
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
            #korbinian.cons_ratio.slice_nonTMD_seqs(dfs, df_nonTMD_sliced, list_of_TMDs)
        if df_TMD.empty:
            # skip protein, as number_of_rows_containing_data == 0 for at least one TMD (or at least the last TMD)
            warning = "{} skipped, number_of_rows_containing_data == 0 for at least one TMD".format(acc)
            logging.info(warning)
            return acc, False, warning

        df_nonTMD_sliced = korbinian.cons_ratio.singleprotein.slice.slice_nonTMD_seqs(dfs, df_nonTMD_sliced, list_of_TMDs)

        df_nonTMD_temp_pickle = os.path.join(homol_dir, "{}_nonTMD_sliced_df.pickle".format(protein_name))
        with open(df_nonTMD_temp_pickle, "wb") as f:
            pickle.dump(df_nonTMD_sliced, f, protocol=pickle.HIGHEST_PROTOCOL)
        homol_sliced_zip.write(df_nonTMD_temp_pickle, arcname=os.path.basename(df_nonTMD_temp_pickle))
        os.remove(df_nonTMD_temp_pickle)
        return  acc, True, "0"

        # ########################################################################################
        # #                                                                                      #
        # #                Slice out TMD regions [fasta and AAIMON]                              #
        # #                                                                                      #
        # ########################################################################################
        #
        # # create a new dataframe to hold the sliced homol sequences for that TMD, and number of gaps, etc
        # df_TMD = pd.DataFrame()
        #
        # '''slice the TMD regions out of the alignment markup and match sequences
        # the method first finds the indices of the TMD region in the query, and uses these indices to slice
        # the filters are not applied yet, so the sequences can be viewed in the csv file for analysis
        # '''
        # # create regex string for each TMD
        # query_TMD_sequence = p['%s_seq'%TMD]
        # # create TMD regex search string (e.g. L.*L.*M.*L.*L.*M.*L.*L.*M.*L.*L.*M.*L.*L.*M.*)
        # TMD_regex_ss = utils.create_regex_string(query_TMD_sequence)
        # # select to only include data where the XML contained a SW node, and then apply function for a regex search
        # dfs_with_SW_node = dfs.query('hit_contains_SW_node == True')
        # # create series of query_alignment_seqs
        # query_seqs_ser = dfs_with_SW_node['query_align_seq']
        #
        # # obtain the bool, start, end of TMD seqs in the match sequences. Add to the new TMD-specific dataframe.
        # df_TMD['%s_start_end_list_in_SW_alignment'%TMD] = query_seqs_ser.apply(utils.get_start_and_end_of_TMD_in_query,args=(TMD_regex_ss,))
        # '''the output of the regex search is a list with three components.
        # 1) a match (e.g. True),
        # 2) the start of the match (e.g. 124),
        # 3) the end of the match (e.g. 144)
        # '''
        # # for convenience, separate the components into separate columns
        # columns_from_regex_output = ['%s_in_SW_alignment'%TMD, '%s_start_in_SW_alignment'%TMD,'%s_end_in_SW_alignment'%TMD]
        # # n = the index (1,2,3) in the tuple
        # # col = column item in the list (e.g. 'TM09_in_SW_alignment')
        # for n, col in enumerate(columns_from_regex_output):
        #     # add a new column which is named TM01_start, etc, and insert the appropriate integer (start or stop) or bool from the tuple
        #     df_TMD[col] = df_TMD['%s_start_end_list_in_SW_alignment'%TMD].dropna().apply(lambda x: x[n])
        # # drop the original listlike from the regex search
        # df_TMD.drop('%s_start_end_list_in_SW_alignment' % TMD, inplace=True)
        #
        # ########################################################################################
        # #                                                                                      #
        # #          Apply function to slice homologues and count gaps [fasta and AAIMON]        #
        # #                                                                                      #
        # ########################################################################################
        #
        # # in some cases, there is no data to obtain as the hit_contains_SW_node = False for too many sequences, giving no start_in_SW_alignment
        # number_of_rows_containing_data = df_TMD[df_TMD['%s_start_end_list_in_SW_alignment' % TMD].notnull()].shape[0]
        # if number_of_rows_containing_data == 0:
        #     logging.info('%s does not have any valid homologues for %s. '
        #                  'Re-downloading simap homologue XML may be necessary.' % (protein_name, TMD))
        # if number_of_rows_containing_data != 0:
        #     n_TMDs_w_homol += 1
        #     len_query_TMD = len(p['%s_seq' % TMD])
        #     # apply the slicing function to the homologues
        #     df_TMD = korbinian.cons_ratio.slice_TMD_homol_and_count_gaps(TMD, len_query_TMD, df_TMD, set_["cr_max_n_gaps_in_query_TMD"], set_["cr_max_n_gaps_in_match_TMD"])


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
    dfs['start_juxta_after_%s'%TMD] = np.where(utils.isNaN(dfs['TM%.2d_start_in_SW_alignment'%(int(TMD[2:])+1)])==True,np.nan,dfs['%s_end_in_SW_alignment'%TMD])
    dfs['end_juxta_before_%s'%TMD] = np.where(dfs["%s_start_in_SW_alignment"%TMD]!=0,dfs["%s_start_in_SW_alignment"%TMD],np.nan)
    dfs['end_juxta_after_%s'%TMD] = dfs["%s_end_in_SW_alignment"%TMD]+((dfs["TM%.2d_start_in_SW_alignment"%(int(TMD[2:])+1)]-dfs["%s_end_in_SW_alignment"%TMD])/2).apply(lambda x :int(x) if not np.isnan(x) else np.nan)
    dfs['start_juxta_before_%s'%TMD] = np.where(dfs["end_juxta_after_TM%.2d"%(int(TMD[2:])-1)] == dfs['end_juxta_before_%s'%TMD] ,dfs["end_juxta_after_TM%.2d"%(int(TMD[2:])-1)],dfs["end_juxta_after_TM%.2d"%(int(TMD[2:])-1)])
    return dfs