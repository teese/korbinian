import ast
import pandas as pd
import csv
import itertools
import korbinian
import numpy as np
import os
import pickle
import re
import sys
import korbinian.utils as utils
from multiprocessing import Pool

def run_calculate_gap_densities(pathdict, s, logging):
    """Runs calculate_gap_densities using multiprocessing Pool.

    Uses multiprocessing to generate a separate output file for each protein, which can be "gathered" later
    using the gather_gap_densities function.

    Parameters
    ----------
    pathdict : dict
        Dictionary of the key paths and files associated with that List number.
    s : dict
        Settings dictionary extracted from excel settings file.
    logging : logging.Logger
        Logger for printing to console and/or logfile.
        If multiprocessing == True, messages sent to the logger, e.g. logging.info(message), will only print to console.

    """
    logging.info("~~~~~~~~~~~~          starting calculate_gap_densities            ~~~~~~~~~~~~")
    # if multiprocessing is used, log only to the console
    p_dict_logging = logging if s["use_multiprocessing"] != True else utils.Log_Only_To_Console()
    # set current working directory as the data_dir/homol, where temp files will be saved before moving to zip
    os.chdir(os.path.join(s["data_dir"], "homol"))

    not_in_homol_db = utils.get_acc_list_from_txt(pathdict["acc_not_in_homol_db_txt"])

    # create list of protein dictionaries to process
    list_p = korbinian.utils.convert_summary_csv_to_input_list(s, pathdict, p_dict_logging, list_excluded_acc=not_in_homol_db)
    # number of processes is the number the settings, or the number of proteins, whichever is smallest
    n_processes = s["multiprocessing_cores"] if s["multiprocessing_cores"] < len(list_p) else len(list_p)

    if s["use_multiprocessing"]:
        with Pool(processes=n_processes) as pool:
            gap_return_statement_list = pool.map(korbinian.gap.calculate_gap_densities, list_p)
            # log the list of protein results (e.g. acc, "simap", True) to the actual logfile, not just the console
            for message in gap_return_statement_list:
                logging.info(message)
    else:
        for p in list_p:
            korbinian.gap.calculate_gap_densities(p)
    logging.info("~~~~~~~~~~~~         calculate_gap_densities is finished          ~~~~~~~~~~~~")

def calculate_gap_densities(p):
    """For one protein, calculates the gap positions amongst homologues.

    Based on the scripts from Rimma Jenske.
    Some changes and annotations by Mark, including the adaption to korbinian, and conversion to multiprocessing.
    Note that Rimma's code for the slicing of JM regions has been inserted into "slice.py".

    Parameters
    ----------
    p : dict
        Protein Dictionary. Contains all input settings, sequences and filepaths related to a single protein.
        Protein-specific data is extracted from one row of the the list summary, e.g. List05_summary.csv, which is read as df.
        p also contains the GENERAL korbinian settings and filepaths for that list (pathdict, s, logging)

        Components of p :
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
    """
    pathdict, s, logging = p["pathdict"], p["s"], p["logging"]
    acc = p["acc"]

    ######################################################################################################################
    #                                                                                                                    #
    #                          Define some constants. Skip protein if already done.                                      #
    #                                                                                                                    #
    ######################################################################################################################

    # Maximum number of gaps for tmds to be considered
    allowed_gaps_per_tmd = s["gap_allowed_gaps_per_tmd"]
    # 24 for beta barrel proteins, can be altered if only several TMDs to consider
    max_number_of_tmds = s["max_number_of_tmds"]
    # # iterate through each protein that has a list_of_TMDs
    protein_name = p["protein_name"]
    # define output file path
    gapout_csv_path = "{}_gapout.csv".format(p['homol_base'])

    # The next steps (the main analysis) is only executed, if previous analysis can be overwritten or no analysis has yet been done
    if s["overwrite_previous_gap_analysis"] == False:
        if os.path.isfile(gapout_csv_path):
            message = "{} skipped, gaps already analysed."
            logging.info(message)
            return acc, False, message
        
    logging.info(acc)
    list_of_TMDs = ast.literal_eval(p["list_of_TMDs"])

    if not os.path.isfile(p['homol_df_orig_zip']):
        message = "{} skipped, {} not found.".format(acc, p['homol_df_orig_zip'])
        logging.info(message)
        return acc, False, message

    ######################################################################################################################
    #                                                                                                                    #
    #                     Filter by FASTA_gapped_identity (of full protein), etc.                                        #
    #    Open homol_df_orig_pickle, which contains the full sequence of all homologues, and GappedIdentity, etc          #
    #                                                                                                                    #
    ######################################################################################################################
    dfh = utils.open_df_from_pickle_zip(p['homol_df_orig_zip'], filename=os.path.basename(p['homol_df_orig_pickle']), delete_corrupt=True)
    dfh.index = dfh.index.astype(int)

    if dfh.empty:
        message = "{} Protein skipped, file deleted as it is possibly corrupt.".format(p['homol_df_orig_zip'])
        logging.info(message)
        return acc, False, message

    gap_homol_query_str = 'FASTA_gapped_identity > {min_ident} and ' \
                          'FASTA_gapped_identity < {max_ident} and ' \
                          'hit_contains_SW_node == True and ' \
                          'disallowed_words_not_in_descr == True and ' \
                          'X_in_match_seq == False'.format(min_ident=s["gap_min_identity_of_full_protein"], max_ident=s["gap_max_identity_of_full_protein"])

    # filter based on the query string
    dfh.query(gap_homol_query_str, inplace=True)
    if dfh.empty:
        message = "{} skipped, filtering by gap_homol_query_str did not leave any valid homologues.".format(acc)
        logging.info(message)
        return acc, False, message

    # keep the index of dfh after filtering based on protein-wide filters
    dfh_filt_index = dfh.index

    # remove dfh from memory, as it shouldn't be necessary anymore
    del dfh

    # open the dataframe for nonTMD sliced, including all TM indices, JM indices, and JM sequences
    dfnon = utils.open_df_from_pickle_zip(p['fa_cr_sliced_TMDs_zip'], filename="{}_nonTMD_sliced_df.pickle".format(protein_name), delete_corrupt=True)
    # filter based on dfh above, for general homologue settings (e.g. % identity of full protein)
    try:
        dfnon = dfnon.loc[dfh_filt_index, :]
    except KeyError:
        # use a try/except loop to detect the rare case where there are no valid homologues between both dataframes
        # (running a set operation for all hit_num of all proteins is likely to slow the script, with no great benefit)
        message = "{} skipped, dfnon does not contain any valid homologues."
        logging.info(message)
        return acc, False, message

    ######################################################################################################################
    #                                                                                                                    #
    #              For the Beta-Barrel Dataset, check if the boolean toggle n_term_ec is there. This should be           #
    #               False for all BB proteins,  but you can confirm it by making sure the first residue is labelled      #
    #               as "I" for Inside.                                                                                   #
    #                                                                                                                    #
    ######################################################################################################################

    """ Current code in uniprot_parse
     # information about location of first non-tmd (extracellular or periplasmic/cytoplasmic)
        if len(location_of_non_tmds_in_feature_list) > 0:
            output_dict['loc_start'] = record.features[location_of_non_tmds_in_feature_list[0]][3]
            output_dict['n_term_ec'] = "Extracellular" in output_dict["loc_start"]
        else:
            output_dict['loc_start'] = np.nan
            output_dict['n_term_ec'] = np.nan
    """

    if "n_term_ec" not in p:
        if "Topology" in p:
            # if the first residue is labelled as "inside, I", or "membrane, M"
            if p["Topology"][0] in ["I", "M"]:
                p["n_term_ec"] = False
            elif p["Topology"][0] == "O":
                p["n_term_ec"] = True
            else:
                raise ValueError('p["Topology"][0] not recognized')
        else:
            raise ValueError('n_term_ec not available')
    # in some cases, p["n_term_ec"] is np.nan. Confirm that if empty, it is True.
    n_term_ec = False if p["n_term_ec"] == False else True

    # create empty output dict, to contain all of the lists of gap positions
    gapout_dict = {}

    # for each TMD in the proteins, creates new lists which will contain gap_positions, lists are saved in a column and created again for each tmd
    for tmd in list_of_TMDs:
        sys.stdout.write(".")
        sys.stdout.flush()
        # open the dataframe containing the sequences, gap counts, etc for that TMD only
        df_s1 = utils.open_df_from_pickle_zip(p['fa_cr_sliced_TMDs_zip'], filename="{}_{}_sliced_df.pickle".format(protein_name, tmd), delete_corrupt=True)
        # filter based on dfh above, for general homologue settings (e.g. % identity of full protein)
        df_s1 = df_s1.loc[dfh_filt_index, :]

        """
        cannot filter out rows that contain no gaps, as the JM number of gaps is not yet defined!
        """
        #df_s1["%s_n_gaps_q_and_m"%tmd] = df_s1["%s_SW_query_num_gaps"%tmd] + df_s1.loc["%s_SW_match_num_gaps"%tmd]

        max_before_filter = df_s1["{}_SW_query_num_gaps".format(tmd)].max()
        shape_before = df_s1.shape

        gap_TMD_query_str =  '({TMD}_SW_query_num_gaps <= {allowed_gaps_per_tmd}) & ' \
                             '({TMD}_SW_match_num_gaps <= {allowed_gaps_per_tmd})'.format(TMD=tmd,allowed_gaps_per_tmd=allowed_gaps_per_tmd)

        # filter based on the query string
        df_s1.query(gap_TMD_query_str, inplace=True)

        #len_of_query = len(df_s1["%s_SW_query_seq"%tmd][1]) # Length of first query sequence, which does (usually) not contain any gaps
        # Length of query TM sequence
        len_of_query = len(p["%s_seq" % tmd])
        len_of_query_reversed= ((1/len_of_query)+1) # Reversed length, important if TMD needs to be reversed afterwards
        list_of_gaps_in_tmd = []
        list_of_gaps_intracellular = []
        list_of_gaps_extracellular = []

        df_s1["%s_n_gaps_q_and_m" % tmd] = df_s1["%s_SW_query_num_gaps" % tmd] + df_s1["%s_SW_match_num_gaps" % tmd]

        min_n_gaps_in_TMD = s["gap_min_n_gaps_in_TMD"]
        max_n_gaps_in_TMD = s["gap_max_n_gaps_in_TMD"]
        # the ax number gaps is  the max number for query or match, x2.
        # this still has to be filtered further (query could have 4 gaps, depending on allowed_gaps_per_tmd), but it still cuts down on the number of items in the loop
        max_n_gaps_in_TMD_q_plus_m = min_n_gaps_in_TMD + max_n_gaps_in_TMD
        filt_string = "{mi} <= {TMD}_n_gaps_q_and_m <= {ma}".format(mi = min_n_gaps_in_TMD, TMD=tmd, ma = max_n_gaps_in_TMD_q_plus_m)
        df_s1_TMD = df_s1.query(filt_string)
        # if the query didn't return an empty dataframe
        if df_s1_TMD.shape[0] != 0:
            for hit in df_s1_TMD.index:
                # '''
                # Start of the main gap analysis
                # Code searches for "-" in the TMD sequence and returns the index!! (not the position)
                # '''
                # Following if conditions only refer to gaps in the query!
                # Query gaps are counted as "in between positions", for example: 4,5 refers to a gap between position 4 and 5;
                # if two gaps occur one after another: only one position (between two amino acids is considered)

                # Filter to make sure, that there are 1 or 2 gaps in the query sequence and up to the max allowed gaps in the match
                if (df_s1.loc[hit, "%s_SW_query_num_gaps" % tmd] != 0.0) and (
                    df_s1.loc[hit, "%s_SW_query_num_gaps" % tmd] <= 2.0):  # and (df_s1.loc[hit,"%s_SW_match_num_gaps"%tmd] <= int("%s"%allowed_gaps_per_tmd)):

                    # Stores the endpoints in a temp list; endpoints are used, to switch from python indices to numbers
                    list_of_gaps_per_hit_in_query = [m.start() for m in re.finditer("-", df_s1.loc[hit, "%s_SW_query_seq" % tmd]) if m.start()]

                    # if there are two gaps in the query (and 2 allowed), code checks if they are side by side (difference of 1)
                    # and appends the gap position, else appends both gap positions
                    # 0.5 is substracted in order to make them "in between" position;
                    # if two gaps are observed, 1.5 is substracted from the second one, since the residue positions are moved due to the first gap
                    if len(list_of_gaps_per_hit_in_query) == 2 and allowed_gaps_per_tmd == 2:

                        if list_of_gaps_per_hit_in_query[1] - list_of_gaps_per_hit_in_query[0] == 1:
                            list_of_gaps_in_tmd.append(list_of_gaps_per_hit_in_query[0] - 0.5)

                        else:
                            list_of_gaps_in_tmd.append(list_of_gaps_per_hit_in_query[0] - 0.5)
                            list_of_gaps_in_tmd.append(list_of_gaps_per_hit_in_query[1] - 1.5)

                    # if there is only one gap in query or only one gap is allowed, it appends the first (and only) gap from the list_of_gaps_per_hit_in_query to the list_of_TMDs
                    else:
                        if len(list_of_gaps_per_hit_in_query) == 1:
                            list_of_gaps_in_tmd.append(list_of_gaps_per_hit_in_query[0] - 0.5)

                # Following if conditions only refer to gaps in the match!
                # Query gaps are counted as deletions of positions; for example: 4 refers to a gap on position 4;
                # if two gaps occur one after another, both are considered since two actual amino acids from the original query are deleted
                # Since the gap positions are dependend on the query sequence, query-gap positions in the same alignment have to be considered as well

                # Filter to make sure, that there are 1 or 2 gaps in the match sequence and up to the max allowed gaps in the query
                if (df_s1.loc[hit, "%s_SW_query_num_gaps" % tmd] <= 2.0) and (df_s1.loc[hit, "%s_SW_match_num_gaps" % tmd] <= 2.0) \
                        and (df_s1.loc[hit, "%s_SW_match_num_gaps" % tmd] != 0.0):

                    # It's not sure that the list of hits in query was already determined, maybe there were no gaps, anyway here it is important how many
                    list_of_gaps_per_hit_in_query = [m.start() for m in re.finditer("-", df_s1.loc[hit, "%s_SW_query_seq" % tmd]) if m.start()]
                    list_of_gaps_per_hit_in_match = [m.start() for m in re.finditer("-", df_s1.loc[hit, "%s_SW_match_seq" % tmd]) if m.start()]

                    if len(list_of_gaps_per_hit_in_match) > 0:
                        for n in list(reversed(list_of_gaps_per_hit_in_match)):
                            substracted_value = len([m < n for m in list_of_gaps_per_hit_in_query])
                            list_of_gaps_in_tmd.append(abs(n - substracted_value))

        # create a view that shows only the desired sequences
        dfnon_j = dfnon[["seq_juxta_before_%s_in_query" % tmd, "seq_juxta_before_%s_in_match"%tmd, "seq_juxta_after_%s_in_query"%tmd, "seq_juxta_after_%s_in_match"%tmd]]
        dfnon_j = dfnon_j.dropna(how="all")

        for hit in dfnon_j.index:
            #######
            # Start of the Juxta Consideration
            # In the case of n_term being located intracellular:
            # there are 4 groups: 1. juxta_before_odd_TMDs + 2.juxta_after_even_TMDs 3. Juxta_before_even_TMDs + 4. Juxta_after_odd_TMDs
            # 1 + 2 --> Intracellular
            # 3 + 4 --> Extracellular
            # If the n_term is extracellular, that it's the other way round. 1+2 --> Extracellular 3+4 --> Intracellular
            ### The data will already be flipped in order to align extracellular and intracellular parts, extracellular: + , intracellular: -

            # juxta before_odd_TMDs:
            if "SP01" in list_of_TMDs:
                return ValueError ("The gap analysis is not currently designed for proteins with signal peptides.")

            tmd_int = int(tmd[-2:])  # Integer of TMD number
            if utils.isOdd(tmd_int) == True:  # also für 1 , 3 ...

                # list of gap indices
                # makes sure that the search is done in a string
                if type(dfnon_j.loc[hit,"seq_juxta_before_%s_in_query"%tmd]) == str:

                    list_of_gaps_in_query_before_odd = [m.start()+0.5 for m in re.finditer("-",dfnon_j.loc[hit,"seq_juxta_before_%s_in_query"%tmd][::-1])if m.start()+0.5 < 31]

                    # if one gap is found, code checks location and appends it
                    if len (list_of_gaps_in_query_before_odd)==1:
                        if n_term_ec == False:
                            list_of_gaps_intracellular.append(list_of_gaps_in_query_before_odd[0])
                        else:
                            list_of_gaps_extracellular.append(list_of_gaps_in_query_before_odd[0])
                    # if more than one gap is found, code checks if the gapy are one after another in the query!
                    if len (list_of_gaps_in_query_before_odd)>1.0:
                        following_gap = 0
                        rev_value = list_of_gaps_in_query_before_odd[0]
                        for n in list_of_gaps_in_query_before_odd:
                            if n-following_gap == rev_value:
                                if n_term_ec == False:
                                    list_of_gaps_intracellular.append(n-following_gap)

                                    following_gap = following_gap+1
                                else:
                                    list_of_gaps_extracellular.append(n-following_gap)
                                    following_gap = following_gap+1
                            else:
                                if n_term_ec == False:
                                    list_of_gaps_intracellular.append(n-following_gap)
                                    following_gap = following_gap+1
                                    rev_value = n
                                else:
                                    list_of_gaps_extracellular.append(n-following_gap)
                                    following_gap = following_gap+1
                                    rev_value = n


                    if type(dfnon_j.loc[hit,"seq_juxta_before_%s_in_match"%tmd])== str:
                    # Makes a list of gaps of the match of the odd juxta before the TMD
                        list_of_gaps_in_query_before_odd = [m.start()+1 for m in re.finditer("-",dfnon_j.loc[hit,"seq_juxta_before_%s_in_query"%tmd][::-1])if m.start() < 32]

                        list_of_gaps_in_match_before_odd = [m.start()+1 for m in re.finditer("-",dfnon_j.loc[hit,"seq_juxta_before_%s_in_match"%tmd][::-1])if m.start() < 32]
                        # This step is essential, to control, if there is a gap before, in the query region
                        for n in list(reversed(list_of_gaps_in_match_before_odd)):
                            greater_values = sum(i< n for i in list_of_gaps_in_query_before_odd)
                            if n_term_ec== False:
                                list_of_gaps_intracellular.append(n-greater_values)

                            else:
                                list_of_gaps_extracellular.append(n-greater_values)

                    # juxta after odd TMDs:

                    if type(dfnon_j.loc[hit,"seq_juxta_after_%s_in_query"%tmd])== str:

                        list_of_gaps_in_query_after_odd = [m.start()+0.5 for m in re.finditer("-",dfnon_j.loc[hit,"seq_juxta_after_%s_in_query"%tmd])if m.start()+0.5 < 31]

                        # if one gap is found, code checks location and appends it
                        if len (list_of_gaps_in_query_after_odd)==1:
                            if n_term_ec == False:
                                list_of_gaps_extracellular.append(list_of_gaps_in_query_after_odd[0])
                            else:
                                list_of_gaps_intracellular.append(list_of_gaps_in_query_after_odd[0])

                        # if more than one gap is found, code checks if the gaps are one after another in the query!
                        if len (list_of_gaps_in_query_after_odd)>1.0:
                            following_gap = 0
                            rev_value = list_of_gaps_in_query_after_odd[0]
                            for n in list_of_gaps_in_query_after_odd:
                                if n+following_gap == rev_value:
                                    if n_term_ec == False:
                                        list_of_gaps_extracellular.append(n-following_gap)
                                        following_gap = following_gap+1
                                    else:
                                        list_of_gaps_intracellular.append(n-following_gap)
                                        following_gap = following_gap+1
                                else:
                                    if n_term_ec == False:
                                        list_of_gaps_extracellular.append(n-following_gap)
                                        following_gap = following_gap+1
                                        rev_value = n
                                    else:
                                        list_of_gaps_intracellular.append(n-following_gap)
                                        following_gap = following_gap+1
                                        rev_value = n
                        # # juxta after odd TMDs:
                        # if type(dfnon_j.loc[hit,"seq_juxta_after_TM%.2d_in_query"%tmd_int])== str:
                        #     list_of_gaps_in_query_after_odd = [m.start()+0.5 for m in re.finditer("-",dfnon_j.loc[hit,"seq_juxta_after_TM%.2d_in_query"%tmd_int])if m.start()+0.5 < 31]
                        #
                        #     #if one gap is found, code checks location and appends it
                        #     if len (list_of_gaps_in_query_after_odd)==1:
                        #        if n_term_ec == False:
                        #            list_of_gaps_extracellular.append(list_of_gaps_in_query_after_odd[0])
                        #        else:
                        #            list_of_gaps_intracellular.append(list_of_gaps_in_query_after_odd[0])
                        #     #if more than one gap is found, code checks if the gaps are one after another in the query!
                        #     if len (list_of_gaps_in_query_after_odd)>1.0:
                        #        following_gap = 0
                        #        rev_value = list_of_gaps_in_query_after_odd[0]
                        #        for n in list_of_gaps_in_query_after_odd:
                        #            if n+following_gap == rev_value:
                        #                if n_term_ec == False:
                        #     #list_of_gaps_extracellular.append(n-following_gap)
                        #                    following_gap = following_gap+1
                        #                else:
                        #                    list_of_gaps_intracellular.append(n-following_gap)
                        #                    following_gap = following_gap+1
                        #            else:
                        #                if n_term_ec == False:
                        #                    list_of_gaps_extracellular.append(n-following_gap)
                        #                    following_gap = following_gap+1
                        #                    rev_value = n
                        #                else:
                        #                    list_of_gaps_intracellular.append(n-following_gap)
                        #                    following_gap = following_gap+1
                        #                    rev_value = n
                else:  # for 2,4

                # juxta before even TMDs:

                # makes sure that the search is done in a string
                    if type(dfnon_j.loc[hit,"seq_juxta_before_%s_in_query"%tmd])== str:

                        # list of gap indices
                        list_of_gaps_in_query_before_even = [m.start()+0.5 for m in re.finditer("-",dfnon_j.loc[hit,"seq_juxta_before_%s_in_query"%tmd])if m.start()+0.5 < 31]

                        # if one gap is found, code checks location and appends it
                        if len (list_of_gaps_in_query_before_even)==1:
                            if n_term_ec == False:
                                list_of_gaps_extracellular.append(list_of_gaps_in_query_before_even[0])
                            else:
                                list_of_gaps_intracellular.append(list_of_gaps_in_query_before_even[0])

                        # if more than one gap is found, code checks if the gaps are one after another in the query!
                        if len (list_of_gaps_in_query_before_even)>1.0:
                            following_gap = 0
                            rev_value = list_of_gaps_in_query_before_even[0]
                            for n in list_of_gaps_in_query_before_even:
                                if n+following_gap == rev_value:
                                    if n_term_ec == False:
                                        list_of_gaps_extracellular.append(n-following_gap)
                                        following_gap = following_gap+1
                                    else:
                                        list_of_gaps_intracellular.append(n-following_gap)
                                        following_gap = following_gap+1
                                else:
                                    if n_term_ec == False:
                                        list_of_gaps_extracellular.append(n-following_gap)
                                        following_gap = following_gap+1
                                        rev_value = n
                                    else:
                                        list_of_gaps_intracellular.append(n-following_gap)
                                        following_gap = following_gap+1
                                        rev_value = n

                    if type(dfnon_j.loc[hit,"seq_juxta_before_%s_in_match"%tmd])== str:
                    # Makes a list of gaps of the match of the odd juxta before the TMD
                        list_of_gaps_in_query_before_even = [m.start()+0.5 for m in re.finditer("-",dfnon_j.loc[hit,"seq_juxta_before_%s_in_query"%tmd])if m.start()+0.5 < 31]

                        list_of_gaps_in_match_before_even = [m.start()+1 for m in re.finditer("-",dfnon_j.loc[hit,"seq_juxta_before_%s_in_match"%tmd])if m.start() < 31]

                        for n in list(reversed(list_of_gaps_in_match_before_even)):
                            greater_values = sum(i< n for i in list_of_gaps_in_query_before_even)
                            if n_term_ec== False:
                                list_of_gaps_extracellular.append(n-greater_values)
                            else:
                                list_of_gaps_intracellular.append(n-greater_values)

                    # juxta after even TMDs:

                    if type(dfnon_j.loc[hit,"seq_juxta_after_%s_in_query"%tmd])== str:

                        list_of_gaps_in_query_after_even = [m.start()+0.5 for m in re.finditer("-",dfnon_j.loc[hit,"seq_juxta_after_%s_in_query"%tmd][::-1]) if m.start()+0.5 < 31]

                        # if one gap is found, code checks location and appends it
                        if len (list_of_gaps_in_query_after_even)==1:
                            if n_term_ec == False:
                                list_of_gaps_intracellular.append(list_of_gaps_in_query_after_even[0])
                            else:
                                list_of_gaps_extracellular.append(list_of_gaps_in_query_after_even[0])

                        # if more than one gap is found, code checks if the gaps are one after another in the query!
                        if len (list_of_gaps_in_query_after_even)>1.0:
                            following_gap = 0
                            rev_value = list_of_gaps_in_query_after_even[0]
                            for n in list_of_gaps_in_query_after_even:
                                if n-following_gap == rev_value:
                                    if n_term_ec == False:
                                        list_of_gaps_intracellular.append(n-following_gap)
                                        following_gap = following_gap+1
                                    else:
                                        list_of_gaps_extracellular.append(n-following_gap)
                                        following_gap = following_gap+1
                                else:
                                    if n_term_ec == False:
                                        list_of_gaps_intracellular.append(n-following_gap)
                                        following_gap = following_gap+1
                                        rev_value = n
                                    else:
                                        list_of_gaps_extracellular.append(n-following_gap)
                                        following_gap = following_gap+1
                                        rev_value = n

                                   # Makes a list of gaps of the match of the odd juxta before the TMD

                    if type(dfnon_j.loc[hit,"seq_juxta_after_%s_in_match"%tmd])== str and type(dfnon_j.loc[hit,"seq_juxta_after_%s_in_query"%tmd])==str:

                        list_of_gaps_in_query_after_even = [m.start()+0.5 for m in re.finditer("-",dfnon_j.loc[hit,"seq_juxta_after_%s_in_query"%tmd][::-1])if m.start()+0.5 < 31]
                        list_of_gaps_in_match_after_even = [m.start()+1 for m in re.finditer("-",dfnon_j.loc[hit,"seq_juxta_after_%s_in_match"%tmd][::-1])if m.start() < 31]

                        for n in list(reversed(list_of_gaps_in_match_after_even)):
                            greater_values = sum(i< n for i in list_of_gaps_in_query_after_even)
                            if n_term_ec== False:
                                list_of_gaps_intracellular.append(n-greater_values)

                            else:
                                list_of_gaps_extracellular.append(n-greater_values)
        # sets of lists are created, to assure, that each gapposition contributes only once to the possible gap positions
        unique_list_of_gaps_in_tmd = list(set(list_of_gaps_in_tmd))
        unique_list_of_gaps_intracellular = list(set(list_of_gaps_intracellular))
        unique_list_of_gaps_extracellular = list(set(list_of_gaps_extracellular))

        # Saves the calculated lists into cells in the columns
        gapout_dict["%s_occurring_gaps"%tmd]=str(unique_list_of_gaps_in_tmd)
        gapout_dict["%s_amount_possible_gap_positions"%tmd]=len(unique_list_of_gaps_in_tmd)

        gapout_dict['juxta_%s_intracellular_possible_gap_positions'%tmd] = str(unique_list_of_gaps_intracellular)
        gapout_dict['juxta_%s_extracellular_possible_gap_positions'%tmd] = str(unique_list_of_gaps_extracellular)
        gapout_dict['juxta_%s_intracellular_num_gaps'%tmd] = len(unique_list_of_gaps_intracellular)
        gapout_dict['juxta_%s_exracellular_num_gaps'%tmd] = len(unique_list_of_gaps_extracellular)

    pd.Series(gapout_dict).to_csv(gapout_csv_path)
    # sys.stdout.write("{} gapout_dict is not saved.".format(tmd))
    # sys.stdout.flush()

    return acc, True, "0"
    # # At the end, sets analysed to true, this is important to not overwrite
    # gapout_dict["gaps_analysed"] = "True"
    # # save to csv
    # df.to_csv(pathdict["list_csv"], sep=",", quoting=csv.QUOTE_NONNUMERIC)

def gather_gap_densities(pathdict, s, logging):
    """Gathers the gap density data saved in the gapout_csv for each protein, and processes data for figure creation.

    currently NOT compatible with multiprocessing

    1) excludes proteins whose homologues could not be downloaded
    2) opens gapout_csv files for each protein, appends to df_gap
    3) in the original protein list, extracts the juxtamembrane region from each TMD


    Parameters
    ----------
    pathdict : dict
        Dictionary of the key paths and files associated with that List number.
    s : dict
        Settings dictionary extracted from excel settings file.
    logging : logging.Logger
        Logger for printing to console and/or logfile.
        If multiprocessing == True, logging.info etc will only print to console.

    Returns
    -------

    """
    logging.info("~~~~~~~~~~~~         starting gather_gap_densities           ~~~~~~~~~~~~")

    ######################################################################################################################
    #                                                                                                                    #
    #               Open list of proteins. Exclude proteins whose homologues could not be downloaded.                    #
    #                                                                                                                    #
    ######################################################################################################################

    df = pd.read_csv(pathdict["list_csv"], sep=",", quoting=csv.QUOTE_NONNUMERIC, index_col=0)
    # get list of accessions that could not be downloaded, and can immediately be excluded
    not_in_homol_db = utils.get_acc_list_from_txt(pathdict["acc_not_in_homol_db_txt"])
    acc_kept = set(df.index) - set(not_in_homol_db)
    # filter to remove proteins not in the homologue database
    df = df.loc[acc_kept, :]
    # remove any proteins from list that do not have a list of TMDs
    df = df.loc[df.list_of_TMDs.notnull()]

    ######################################################################################################################
    #                                                                                                                    #
    #                      Combine the gap data from all proteins into a single dataframe                                #
    #                                                                                                                    #
    ######################################################################################################################
    # create an empty dataframe for gathering the various output files
    df_gap = pd.DataFrame()
    # iterate over the dataframe for proteins with an existing list_of_TMDs. acc = uniprot accession.
    for acc in df.index:
        protein_name = df.loc[acc, 'protein_name']
        sys.stdout.write(" {}".format(protein_name))
        sys.stdout.flush()
        # define gap output file path
        gapout_csv_path = "{}_gapout.csv".format(df.loc[acc,'homol_base'])
        if not os.path.exists(gapout_csv_path):
            logging.info("{} {} Protein skipped. File does not exist".format(acc, gapout_csv_path))
            continue
        # open csv as pandas dataframe (note, it was originally a series, and contains only one column and an index)
        gapout_df = pd.read_csv(gapout_csv_path, index_col=0)
        gapout_df.columns = ["value"]
        gapout_df.loc["acc", "value"] = acc
        gapout_df.loc["list_of_TMDs", "value"] = df.loc[acc, "list_of_TMDs"]
        df_gap = pd.concat([df_gap,gapout_df], axis=1)

    # transpose dataframe df_gap
    df_gap = df_gap.T
    df_gap.set_index("acc", inplace=True)

    df_gap.to_csv(pathdict["list_gap_summary_csv"], sep=",", quoting=csv.QUOTE_NONNUMERIC)

    logging.info("~~~~~~~~~~~~         finished gather_gap_densities           ~~~~~~~~~~~~")

    # # test if the dataframe has already been created, otherwise re-open from uniprot csv file
    # if os.path.isfile(pathdict["dfout10_uniprot_gaps"]):
    #     df = pd.read_csv(pathdict["dfout10_uniprot_gaps"], sep=",", quoting=csv.QUOTE_NONNUMERIC, index_col=[0])
    #     logging.info('df loaded from %s' % pathdict["dfout10_uniprot_gaps"])
    # else:
    #     raise FileNotFoundError(
    #         'No gap analysis has been done yet. %s is not found. Please run calculate calculate_gap_densities' % pathdict[
    #             "dfout10_uniprot_gaps"])

    #df = pd.read_csv(pathdict["list_csv"], sep = ",", quoting = csv.QUOTE_NONNUMERIC, index_col = 0)
    ######################################################################################################################
    #                                                                                                                    #
    #               Get juxtamembrane regions from all proteins. Code is modified from slice_TMD_1_prot_from_homol,      #
    #                                where all homologues were sliced for a particular protein                           #
    #                                                                                                                    #
    ######################################################################################################################

    df["list_of_TMDs"] = df.list_of_TMDs.apply(ast.literal_eval)
    df["len_list_of_TMDs"] = df["list_of_TMDs"].apply(lambda x : len(x))
    index_longest_list_TMDs = df["len_list_of_TMDs"].idxmax()
    longest_list_TMDs = df.loc[index_longest_list_TMDs, "list_of_TMDs"]

    df["is_multipass"] = df.list_of_TMDs.apply(lambda x: "TM02" in x)
    # np.where syntax: np.where(boolean_query, value_if_query_true, value_if_query_false)
    # @RJ, If TM01_start_in_SW_alignment is not an integer above 0, replaces with np.nan?
    df['start_juxta_before_TM01'] = np.where(df['TM01_start'] > 0, 0, np.nan)
    # if the TM01_start_in_SW_alignment is 0, there is no JM region N-terminal to the TMD, therefore replace end_juxta_before_TM01 with np.nan, otherwise use TM01_start_in_SW_alignment
    df['end_juxta_before_TM01'] = np.where(df['TM01_start'] == 0, np.nan, df['TM01_start'])
    # the start of the juxtamembrane region after TM01 is the end of TM01
    df['start_juxta_after_TM01'] = df['TM01_end']

    # divide into single-pass (sp) and multipass (mp) dataframes
    sp = df.loc[df["is_multipass"] == False]
    mp = df.loc[df["is_multipass"]]

    if not mp.empty:
        # for multipass proteins, the end after TM01 is the middle between TM01 and TM02
        mp_end_juxta_after_TM01 = mp["TM01_end"] + ((mp["TM02_start"] - mp["TM01_end"]) / 2)
    else:
        # if there is no multipass in the dataset, return an empty series, which will be ignored
        mp_end_juxta_after_TM01 = pd.Series()
    if not sp.empty:
        # for single-pass proteins, the end after TM01 is the end of the sequence
        sp_end_juxta_after_TM01 = sp['seqlen']
    else:
        # if there is no singlepass in the dataset, return an empty series, which will be ignored
        sp_end_juxta_after_TM01 = pd.Series()
    # join the mp and sp together, and add to the main dataframe as integers
    df['end_juxta_after_TM01'] = pd.concat([mp_end_juxta_after_TM01, sp_end_juxta_after_TM01]).astype(int)

    logging.info('~~~~~~~~~~~~starting extraction of JM regions ~~~~~~~~~~~~')

    for acc_nr, acc in enumerate(df.index):
        sys.stdout.write("{} ".format(acc))
        if acc_nr != 0 and acc_nr % 50 == 0:
            sys.stdout.write("\n")
        sys.stdout.flush()

        list_of_TMDs = df.loc[acc, "list_of_TMDs"]
        last_TMD_of_acc = list_of_TMDs[-1]
        for TMD in list_of_TMDs:
            if TMD == "TM01":
                # skip, it's already done
                continue

            next_TM = "TM{:02d}".format(int(TMD[2:]) + 1)
            prev_TM = "TM{:02d}".format(int(TMD[2:])-1)

            """2222222222222nnnnnnnnn|NNNNNNNNN333333333333333aaaaaaaaaaaaaaaaaZ
            2222222222222 = TM02
            333333333333333 = TM03
            nnnnnnnnn = juxta after TM02
            NNNNNNNNN = juxta before TM03
            aaaaaaaaa = juxta after TM03, the last TMD
            Z = final aa, length of sequence, "seqlen"
            """
            # start_juxta_before_TM02 = |
            df.loc[acc, 'start_juxta_before_%s'%TMD] = df.loc[acc, "end_juxta_after_{}".format(prev_TM)]
            # end_juxta_before_TM02 = TM02_start
            df.loc[acc, 'end_juxta_before_%s'%TMD] = df.loc[acc, "%s_start"%TMD]
            # start_juxta_after_TM02 = TM02_end
            df.loc[acc, 'start_juxta_after_%s'%TMD] = df.loc[acc, '%s_end'%TMD]
            if not TMD == last_TMD_of_acc:
                # end_juxta_after_TM02 = halfway between TM02_end and TM03_start = TM02_end + ((TM03_start - TM02_end) / 2)
                df.loc[acc, 'end_juxta_after_%s' % TMD] = df.loc[acc, "%s_end" % TMD] + ((df.loc[acc, "{}_start".format(next_TM)] - df.loc[acc, "%s_end" % TMD]) / 2)
            elif TMD == last_TMD_of_acc:
                # if TM02 is the last TMD
                # end_juxta_after_TM02 = seqlen
                df.loc[acc, 'end_juxta_after_%s' % TMD] = df.loc[acc, "seqlen"]
        for TMD in list_of_TMDs:
            df.loc[acc, 'len_juxta_before_{}'.format(TMD)] = df.loc[acc, 'end_juxta_before_%s'%TMD] - df.loc[acc, 'start_juxta_before_%s'%TMD]
            df.loc[acc, 'len_juxta_after_{}'.format(TMD)] = df.loc[acc, 'end_juxta_after_%s'%TMD] - df.loc[acc, 'start_juxta_after_%s'%TMD]

    logging.info('~~~~~~~~~~~~finished extraction of JM regions ~~~~~~~~~~~~')

    # if TMD == "TM01":
    #     # np.where syntax: np.where(boolean_query, value_if_query_true, value_if_query_false)
    #     # @RJ, If TM01_start_in_SW_alignment is not an integer above 0, replaces with np.nan?
    #     df['start_juxta_before_TM01'] = np.where(df['TM01_start'] > 0, 0, np.nan)
    #     # if the TM01_start_in_SW_alignment is 0, there is no JM region N-terminal to the TMD, therefore replace end_juxta_before_TM01 with np.nan, otherwise use TM01_start_in_SW_alignment
    #     df['end_juxta_before_TM01'] = np.where(df['TM01_start'] == 0, np.nan, df['TM01_start'])
    #     if df.loc[acc, "is_multipass"]:
    #         df['end_juxta_after_TM01'] = df["TM01_end"] + ((df["TM02_start"] - df["TM01_end"]) / 2)
    #     else:
    #         df['end_juxta_after_TM01'] = np.where(utils.isNaN(df['start_juxta_after_TM01']) == True, np.nan, df['len_query_align_seq'])
    # if df.loc[acc, "is_multipass"]:
    #     if not TMD == "TM01" and not TMD == last_TMD_of_acc:
    #         df = juxta_function_orig(df, TMD)
    #
    #     if TMD == last_TMD_of_acc:
    #         df['start_juxta_before_%s' % TMD] = df['end_juxta_after_TM%.2d' % (int(TMD[2:]) - 1)]
    #         df['end_juxta_before_%s' % TMD] = df['%s_start_in_SW_alignment' % TMD]
    #         df['start_juxta_after_%s' % TMD] = np.where(
    #             df['%s_end_in_SW_alignment' % TMD] == df['len_query_align_seq'], np.nan,
    #             df['%s_end_in_SW_alignment' % TMD])
    #         df['end_juxta_after_%s' % TMD] = np.where(utils.isNaN(df['start_juxta_after_%s' % TMD]) == True, np.nan, df['len_query_align_seq'])


    #df_gap = pd.read_csv(pathdict["list_gap_summary_csv"], sep = ",", quoting = csv.QUOTE_NONNUMERIC, index_col = 0)

    ######################################################################################################################
    #                                                                                                                    #
    #                      Remove proteins with no gap info. NOTE, this converts any "number of gaps"                    #
    #                     to NaN, and therefore may excludes some real data (i.e. proteins with no gaps)                 #
    #                                                                                                                    #
    ######################################################################################################################

    init_n_proteins = df_gap.shape[0]
    df_gap.replace(["[]", 0.0], np.nan, inplace=True)
    df_gap.dropna(how="all", axis=0, inplace=True)
    n_proteins_with_gap_info = df_gap.shape[0]
    n_proteins_excluded = init_n_proteins - n_proteins_with_gap_info
    logging.info("init_n_proteins, n_proteins_with_gap_info, n_proteins_excluded = {}, {}, {}".format(init_n_proteins, n_proteins_with_gap_info, n_proteins_excluded))

    # find the maximum number of TMDs amongst the proteins
    n_TMDs_max = int(df["number_of_TMDs"].max())
    TMD_range = range(1, n_TMDs_max + 1)
    #TMD_range_plus_1 = range(1, n_TMDs_max + 2)
    TMD_range_2nd = range(1, n_TMDs_max + 1, 2)

    for num_TMD in range(1, n_TMDs_max + 1):
        # use literal_eval to convert the stringlist to a python list
        df_gap["TM%.2d_occurring_gaps" % num_TMD] = df_gap["TM%.2d_occurring_gaps" % num_TMD].dropna().apply(ast.literal_eval)

    ######################################################################################################################
    #                                                                                                                    #
    #          Flip the gap propensities in the TMDs so that the intracellular side is always at the start (left)        #
    #                                                                                                                    #
    ######################################################################################################################
    num_of_bins_in_tmd_region = s["num_of_bins_in_tmd_region"]
    flipped = []
    not_flipped = []
    # times 2, because TMDs in gap and in query are considered! --> double amount
    #n_TMDs_in_all_proteins = df.loc[df.gaps_analysed == True, "number_of_TMDs"].sum() * 2
    df_gap["list_of_TMDs"] = df_gap["list_of_TMDs"].apply(ast.literal_eval)
    # not sure why TMDs were counted, and not added, as here
    #total_n_TMDs = 0

    logging.info('~~~~~~~~~~~~ starting flip ~~~~~~~~~~~~')
    n_TMDs_in_all_proteins = 0
    for acc in df_gap.index:
        # count the total number of TMDs in all proteins in list, with gap data
        n_TMDs = len(df_gap.loc[acc,"list_of_TMDs"])
        n_TMDs_in_all_proteins += n_TMDs

        # get number of TMDs for that protein
        #n_TMDs = int(df.loc[acc, "number_of_TMDs"])
        #total_n_TMDs += n_TMDs
        #for num_TMD in [1, 3, 5, 7 etc]

        for num_TMD in range(1, n_TMDs + 1, 2):
            sys.stdout.write("num_TMD", num_TMD)
            list_occurring_gaps = df_gap.loc[acc, "TM%.2d_occurring_gaps" % num_TMD]
            sys.stdout.write("list_occurring_gaps", list_occurring_gaps)
            for n in list_occurring_gaps:
                sys.stdout.write ("n {}".format(n))
                sys.stdout.write ((df.loc[acc,"TM%.2d_seq"%num_TMD]) )
                TMD_seq = df.loc[acc, "TM%.2d_seq" % num_TMD]
                TMD_len = len(TMD_seq)
                n_over_TMD_len_times_n_bins = n / (TMD_len - 1) * num_of_bins_in_tmd_region

                if df.loc[acc, "n_term_ec"] == False:
                    not_flipped.append(n_over_TMD_len_times_n_bins)
                if df.loc[acc, "n_term_ec"] == True:
                    flipped.append(num_of_bins_in_tmd_region - n_over_TMD_len_times_n_bins)
            # if not utils.isNaN(df_gap.loc[acc, "TM%.2d_occurring_gaps" % num_TMD]):
            #     for n in ast.literal_eval(df_gap.loc[acc, "TM%.2d_occurring_gaps" % num_TMD]):
            #         sys.stdout.write ("n {}".format(n))
            #         sys.stdout.write ((df.loc[acc,"TM%.2d_seq"%num_TMD]) )
            #         TMD_len = len(df.loc[acc, "TM%.2d_seq" % num_TMD])
            #         n_over_TMD_len_times_n_bins = n / (TMD_len - 1) * num_of_bins_in_tmd_region
            #
            #         if df.loc[acc, "n_term_ec"] == False:
            #             not_flipped.append(n_over_TMD_len_times_n_bins)
            #         if df.loc[acc, "n_term_ec"] == True:
            #             flipped.append(num_of_bins_in_tmd_region - n_over_TMD_len_times_n_bins)


        for num_TMD in range(2, int(df.loc[acc, "number_of_TMDs"]) + 1, 2):
            if not utils.isNaN(df_gap.loc[acc, "TM%.2d_occurring_gaps" % num_TMD]):
                for n in ast.literal_eval(df_gap.loc[acc, "TM%.2d_occurring_gaps" % num_TMD]):
                    # m = np.prod(df.loc[acc,"TM%.2d_len"%num_TMD]*float(n))
                    if df.loc[acc, "n_term_ec"] == False:
                        flipped.append(num_of_bins_in_tmd_region - ((n / (len(df.loc[acc, "TM%.2d_seq" % num_TMD]) - 1)) * num_of_bins_in_tmd_region))

                    if df.loc[acc, "n_term_ec"] == True:
                        not_flipped.append((n / (len(df.loc[acc, "TM%.2d_seq" % num_TMD]) - 1)) * num_of_bins_in_tmd_region)
    logging.info('~~~~~~~~~~~~ finished flip ~~~~~~~~~~~~')


    ######################################################################################################################
    #                                                                                                                    #
    #           Extract all lists of gap positions. Note that this was orig a single line, but since it is fast,         #
    #                          it is better to split it up and make the code more readable.                              #
    #                                                                                                                    #
    ######################################################################################################################

    nested_list_of_gaps_intracellular = []
    nested_list_of_gaps_extracellular = []
    for TMD_nr in TMD_range:
        # for the list of proteins, create a long list of the intracellular gap positions for that TMD (e.g. TM01)
        intracell_gap_pos_ser = df_gap['juxta_TM{:02d}_intracellular_possible_gap_positions'.format(TMD_nr)].dropna().apply(ast.literal_eval)
        extracell_gap_pos_ser = df_gap['juxta_TM{:02d}_extracellular_possible_gap_positions'.format(TMD_nr)].dropna().apply(ast.literal_eval)
        # append as a list, rather than a series
        nested_list_of_gaps_intracellular.append(list(itertools.chain(*intracell_gap_pos_ser.tolist())))
        nested_list_of_gaps_extracellular.append(list(itertools.chain(*extracell_gap_pos_ser.tolist())))

    """ ORIGINAL LIST COMPREHENSION CODE"""
    # nested_list_of_gaps_intracellular = [ast.literal_eval(m) for n in range (1,25) for m in df['juxta_TM%.2d_intracellular_possible_gap_positions'%n].dropna().tolist()]
    #nested_list_of_gaps_intracellular = [ast.literal_eval(m) for n in TMD_range for m in df['juxta_TM%.2d_intracellular_possible_gap_positions' % n].dropna().tolist()]
    # data for extracellular part --> right
    #nested_list_of_gaps_extracellular = [ast.literal_eval(m) for n in TMD_range for m in df['juxta_TM%.2d_extracellular_possible_gap_positions' % n].dropna().tolist()]

    # join all values in all lists together to make a single list of floats, so they can be used to make a histogram
    hist_data_juxta_intracellular = np.array(list(itertools.chain(*nested_list_of_gaps_intracellular)))
    hist_data_juxta_extracellular = np.array(list(itertools.chain(*nested_list_of_gaps_extracellular)))

    min_value = int(abs(hist_data_juxta_intracellular.min()))
    # y-axis_intracell = [(freq_counts_I.tolist()[::-1][n]/frequency_of_position_intracellular(n)) for n in range (1,int(min_value))]

    #n_TMDs_in_all_proteins = len([TMD for acc in df.index for TMD in ast.literal_eval(df.loc[acc,"list_of_TMDs"])if (df.loc[acc,"gaps_analysed"]==True)and(utils.isNaN(df.loc[acc,"list_of_TMDs"]))])*2
    logging.info("n_TMDs_in_all_proteins {}".format(n_TMDs_in_all_proteins))

    list_of_positionfrequency_extra = []
    list_of_positionfrequency_intra = []

    logging.info("~~~~~~~~~~~~        starting list_of_positionfrequency       ~~~~~~~~~~~~")

    for acc in df.index:
        #if df.loc[acc, "gaps_analysed"] == True:
        logging.info(acc)
        if df.loc[acc, "n_term_ec"] == True:
            # no idea why the TMD range should be larger than the actual number of TMDs!! strange. Changed!
            #for n in TMD_range_plus_1:
            for n in TMD_range:
                list_of_positionfrequency_extra.append(df.loc[acc, 'len_juxta_before_TM%.2d' % n].tolist())
                list_of_positionfrequency_intra.append(df.loc[acc, 'len_juxta_after_TM%.2d' % n].tolist())
            for n in TMD_range_2nd:
                list_of_positionfrequency_extra.append(df.loc[acc, 'len_juxta_after_TM%.2d' % n].tolist())
                list_of_positionfrequency_intra.append(df.loc[acc, 'len_juxta_before_TM%.2d' % n].tolist())

        if df.loc[acc, "n_term_ec"] == False:
            # no idea why the TMD range should be larger than the actual number of TMDs!! strange. Changed!
            #for n in TMD_range_plus_1:
            for n in TMD_range:
                list_of_positionfrequency_extra.append(df.loc[acc, 'len_juxta_after_TM%.2d' % n].tolist())
                list_of_positionfrequency_intra.append(df.loc[acc, 'len_juxta_before_TM%.2d' % n].tolist())
            for n in TMD_range_2nd:
                list_of_positionfrequency_extra.append(df.loc[acc, 'len_juxta_before_TM%.2d' % n].tolist())
                list_of_positionfrequency_intra.append(df.loc[acc, 'len_juxta_after_TM%.2d' % n].tolist())

    logging.info("~~~~~~~~~~~~        finished list_of_positionfrequency       ~~~~~~~~~~~~")

    ######################################################################################################################
    #                                                                                                                    #
    #                Save all generated data to a pickle file, to be opened later by the plotting scripts                #
    #                                                                                                                    #
    ######################################################################################################################
    # create a list of all the outpupts from the gap analysis
    # this can be saved as a pickle, without affecting datatypes
    output_data = [flipped, not_flipped, hist_data_juxta_intracellular, hist_data_juxta_extracellular, min_value, list_of_positionfrequency_extra, list_of_positionfrequency_intra, n_TMDs_in_all_proteins]

    # save the output_data list (includes lists and ints) as a pickle file, to be opened later for plotting
    # e.g. path = "D:\Databases\summaries\03\List03_gap_data.pickle"
    with open(pathdict["gap_data_pickle"], "wb") as pkl_file:
        pickle.dump(output_data, pkl_file, -1)

    logging.info('~~~~~~~~~~~~ finished list_of_positionfrequency_xxxx ~~~~~~~~~~~~')
    logging.info("~~~~~~~~~~~~        gather_gap_densities is finished         ~~~~~~~~~~~~")

