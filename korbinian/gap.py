import ast
import pandas as pd
import numpy as np
import csv
import os
import re
import sys
import korbinian.utils as utils
import zipfile

def calculate_gap_densities(pathdict, s, logging):
    logging.info("~~~~~~~~~~~~          starting calculate_gap_densities            ~~~~~~~~~~~~")

    # If script previously has been run, continues with proteins not being processed yet, or overwrites previous gap analysis
    overwrite_previous_gap_analysis = s["overwrite_previous_gap_analysis"]

    # Maximum number of gaps for tmds to be considered
    allowed_gaps_per_tmd = s["gap_allowed_gaps_per_tmd"]

    # 24 for beta barrel proteins, can be altered if only several TMDs to consider
    max_number_of_tmds = s["max_number_of_tmds"]

    df = pd.read_csv(pathdict["list_summary_csv"], sep = ",", quoting = csv.QUOTE_NONNUMERIC, index_col = 0)

    print("n_term_ec in columns = ", "n_term_ec" in df.columns)

    # iterate through each protein that has a list_of_TMDs
    for acc in df.loc[df['list_of_TMDs'].notnull()].loc[df['list_of_TMDs'] != 'nan'].index:
        protein_name = df.loc[acc, "protein_name"]

        # The next steps (the main analysis) is only executed, if previous analysis can be overwritten or no analysis has yet been done
        if (overwrite_previous_gap_analysis == True) or (df.loc[acc,"gaps_analysed"] != True):
            logging.info(acc)
            list_of_TMDs = ast.literal_eval(df.loc[acc,"list_of_TMDs"])

            dfh = utils.open_df_from_pickle_zip(df.loc[acc,'homol_df_orig_zip'], filename=os.path.basename(df.loc[acc,'homol_df_orig_pickle']), delete_corrupt=True)
            if dfh.empty:
                message = "{} Protein skipped, file deleted as it is possibly corrupt.".format(df.loc[acc,'homol_df_orig_zip'])
                logging.info(message)
                return acc, False, message

            gap_homol_query_str = 'FASTA_gapped_identity > {min_ident} and ' \
                                  'FASTA_gapped_identity < {max_ident} and ' \
                                  'hit_contains_SW_node == True and ' \
                                  'disallowed_words_not_in_descr == True and ' \
                                  'X_in_match_seq == False'.format(min_ident=s["gap_min_identity_of_full_protein"], max_ident=s["gap_max_identity_of_full_protein"])

            # filter based on the query string
            dfh.query(gap_homol_query_str, inplace=True)

            # keep the index of dfh after filtering based on protein-wide filters
            dfh_filt_index = dfh.index

            # remove dfh from memory, as it shouldn't be necessary anymore
            del dfh

            # open the dataframe for nonTMD sliced, including all TM indices, JM indices, and JM sequences
            dfnon = utils.open_df_from_pickle_zip(df.loc[acc, 'fa_cr_sliced_TMDs_zip'], filename="{}_nonTMD_sliced_df.pickle".format(protein_name), delete_corrupt=True)
            # filter based on dfh above, for general homologue settings (e.g. % identity of full protein)
            dfnon = dfnon.loc[dfh_filt_index, :]

            if df.loc[acc,"n_term_ec"] == False:
                reverse_tmd = False
            else:
                reverse_tmd = True
            print (reverse_tmd)

            # for each TMD in the proteins, creates new lists which will contain gappositions, lists are saved in a column and created again for each tmd
            for tmd in list_of_TMDs:
                sys.stdout.write(".")
                sys.stdout.flush()

                # open the dataframe containing the sequences, gap counts, etc for that TMD only
                df_s1 = utils.open_df_from_pickle_zip(df.loc[acc, 'fa_cr_sliced_TMDs_zip'], filename="{}_{}_sliced_df.pickle".format(protein_name, tmd), delete_corrupt=True)
                # filter based on dfh above, for general homologue settings (e.g. % identity of full protein)
                df_s1 = df_s1.loc[dfh_filt_index, :]

                # gap_TMD_query_str =  '{min_gaps} < {TMD}_SW_query_num_gaps < {max_gaps} and ' \
                #                      '{min_gaps} < {TMD}_SW_match_num_gaps < {max_gaps}'.format(TMD=tmd, min_gaps=s["gap_min_n_gaps_in_TMD"], max_gaps=s["gap_max_n_gaps_in_TMD"])

                gap_TMD_query_str =  '({TMD}_SW_query_num_gaps < {allowed_gaps_per_tmd}) & ' \
                                     '({TMD}_SW_match_num_gaps < {allowed_gaps_per_tmd})'.format(TMD=tmd,allowed_gaps_per_tmd=allowed_gaps_per_tmd)

                # filter based on the query string
                df_s1.query(gap_TMD_query_str, inplace=True)

                #len_of_query = len(df_s1["%s_SW_query_seq"%tmd][1]) # Length of first query sequence, which does (usually) not contain any gaps
                # Length of query TM sequence
                len_of_query = len(df.loc[acc, "%s_seq" % tmd])
                len_of_query_reversed= ((1/len_of_query)+1) # Reversed length, important if TMD needs to be reversed afterwards
                list_of_gaps_in_tmd = []
                list_of_gaps_intracellular = []
                list_of_gaps_extracellular = []

                for hit in df_s1.index:
                #'''
                #Start of the main gap analysis
                #Code searches for "-" in the TMD sequence and returns the index!! (not the position)
                #'''
                    # Following if conditions only refer to gaps in the query!
                    # Query gaps are counted as "in between positions", for example: 4,5 refers to a gap between position 4 and 5;
                    # if two gaps occur one after another: only one position (between two amino acids is considered)

                    # Filter to make sure, that there are 1 or 2 gaps in the query sequence and up to the max allowed gaps in the match
                    if (df_s1.loc[hit,"%s_SW_query_num_gaps"%tmd] != 0.0) and (df_s1.loc[hit,"%s_SW_query_num_gaps"%tmd] <= 2.0):# and (df_s1.loc[hit,"%s_SW_match_num_gaps"%tmd] <= int("%s"%allowed_gaps_per_tmd)):

                        # Stores the endpoints in a temp list; endpoints are used, to switch from python indices to numbers
                        list_of_gaps_per_hit_in_query = [m.start() for m in re.finditer("-",df_s1.loc[hit,"%s_SW_query_seq"%tmd]) if m.start()]

                        # if there are two gaps in the query (and 2 allowed), code checks if they are side by side (difference of 1)
                        # and appends the gap position, else appends both gap positions
                        # 0.5 is substracted in order to make them "in between" position;
                        #if two gaps are observed, 1.5 is substracted from the second one, since the residue positions are moved due to the first gap
                        if len(list_of_gaps_per_hit_in_query) == 2 and allowed_gaps_per_tmd==2:

                            if list_of_gaps_per_hit_in_query[1]- list_of_gaps_per_hit_in_query[0] ==1:
                                list_of_gaps_in_tmd.append(list_of_gaps_per_hit_in_query[0]-0.5)

                            else:
                                list_of_gaps_in_tmd.append(list_of_gaps_per_hit_in_query[0]-0.5)
                                list_of_gaps_in_tmd.append(list_of_gaps_per_hit_in_query[1]-1.5)

                        # if there is only one gap in query or only one gap is allowed, it appends the first (and only) gap from the list_of_gaps_per_hit_in_query to the list_of_TMDs
                        else:
                            if len(list_of_gaps_per_hit_in_query) == 1:
                                list_of_gaps_in_tmd.append (list_of_gaps_per_hit_in_query[0]-0.5)


                    # Following if conditions only refer to gaps in the match!
                    # Query gaps are counted as deletions of positions; for example: 4 refers to a gap on position 4;
                    # if two gaps occur one after another, both are considered since two actual amino acids from the original query are deleted
                    # Since the gap positions are dependend on the query sequence, query-gap positions in the same alignment have to be considered as well

                    # Filter to make sure, that there are 1 or 2 gaps in the match sequence and up to the max allowed gaps in the query
                    if (df_s1.loc[hit,"%s_SW_query_num_gaps"%tmd] <=2.0) and (df_s1.loc[hit,"%s_SW_match_num_gaps"%tmd] <= 2.0)\
                        and (df_s1.loc[hit,"%s_SW_match_num_gaps"%tmd] != 0.0):

                        # It's not sure that the list of hits in query was already determined, maybe there were no gaps, anyway here it is important how many
                        list_of_gaps_per_hit_in_query = [m.start() for m in re.finditer("-",df_s1.loc[hit,"%s_SW_query_seq"%tmd]) if m.start()]
                        list_of_gaps_per_hit_in_match = [m.start() for m in re.finditer("-",df_s1.loc[hit,"%s_SW_match_seq"%tmd])if m.start()]
                        #print (list_of_gaps_per_hit_in_query)
                        #print (list_of_gaps_per_hit_in_match)

                        if len(list_of_gaps_per_hit_in_match)>0:
                            for n in list(reversed(list_of_gaps_per_hit_in_match)):
                                substracted_value = len([m<n for m in list_of_gaps_per_hit_in_query])
                                list_of_gaps_in_tmd.append(abs(n-substracted_value))

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
                        if type(dfnon.loc[hit,"seq_juxta_before_%s_in_query"%tmd]) == str:


                            list_of_gaps_in_query_before_odd = [m.start()+0.5 for m in re.finditer("-",dfnon.loc[hit,"seq_juxta_before_%s_in_query"%tmd][::-1])if m.start()+0.5 < 31]

                            # if one gap is found, code checks location and appends it
                            if len (list_of_gaps_in_query_before_odd)==1:
                                if reverse_tmd == False:
                                    list_of_gaps_intracellular.append(list_of_gaps_in_query_before_odd[0])
                                else:
                                    list_of_gaps_extracellular.append(list_of_gaps_in_query_before_odd[0])
                            # if more than one gap is found, code checks if the gapy are one after another in the query!
                            if len (list_of_gaps_in_query_before_odd)>1.0:
                                following_gap = 0
                                rev_value = list_of_gaps_in_query_before_odd[0]
                                for n in list_of_gaps_in_query_before_odd:
                                    if n-following_gap == rev_value:
                                        if reverse_tmd == False:
                                            list_of_gaps_intracellular.append(n-following_gap)

                                            following_gap = following_gap+1
                                        else:
                                            list_of_gaps_extracellular.append(n-following_gap)
                                            following_gap = following_gap+1
                                    else:
                                        if reverse_tmd == False:
                                            list_of_gaps_intracellular.append(n-following_gap)
                                            following_gap = following_gap+1
                                            rev_value = n
                                        else:
                                            list_of_gaps_extracellular.append(n-following_gap)
                                            following_gap = following_gap+1
                                            rev_value = n


                            if type(dfnon.loc[hit,"seq_juxta_before_%s_in_match"%tmd])== str:
                            # Makes a list of gaps of the match of the odd juxta before the TMD
                                list_of_gaps_in_query_before_odd = [m.start()+1 for m in re.finditer("-",dfnon.loc[hit,"seq_juxta_before_%s_in_query"%tmd][::-1])if m.start() < 32]

                                list_of_gaps_in_match_before_odd = [m.start()+1 for m in re.finditer("-",dfnon.loc[hit,"seq_juxta_before_%s_in_match"%tmd][::-1])if m.start() < 32]
                                # This step is essential, to control, if there is a gap before, in the query region
                                for n in list(reversed(list_of_gaps_in_match_before_odd)):
                                    greater_values = sum(i< n for i in list_of_gaps_in_query_before_odd)
                                    if reverse_tmd== False:
                                        list_of_gaps_intracellular.append(n-greater_values)

                                    else:
                                        list_of_gaps_extracellular.append(n-greater_values)

                            # juxta after odd TMDs:

                            if type(dfnon.loc[hit,"seq_juxta_after_%s_in_query"%tmd])== str:

                                list_of_gaps_in_query_after_odd = [m.start()+0.5 for m in re.finditer("-",dfnon.loc[hit,"seq_juxta_after_%s_in_query"%tmd])if m.start()+0.5 < 31]

                                # if one gap is found, code checks location and appends it
                                if len (list_of_gaps_in_query_after_odd)==1:
                                    if reverse_tmd == False:
                                        list_of_gaps_extracellular.append(list_of_gaps_in_query_after_odd[0])
                                    else:
                                        list_of_gaps_intracellular.append(list_of_gaps_in_query_after_odd[0])

                                # if more than one gap is found, code checks if the gaps are one after another in the query!
                                if len (list_of_gaps_in_query_after_odd)>1.0:
                                    following_gap = 0
                                    rev_value = list_of_gaps_in_query_after_odd[0]
                                    for n in list_of_gaps_in_query_after_odd:
                                        if n+following_gap == rev_value:
                                            if reverse_tmd == False:
                                                list_of_gaps_extracellular.append(n-following_gap)
                                                following_gap = following_gap+1
                                            else:
                                                list_of_gaps_intracellular.append(n-following_gap)
                                                following_gap = following_gap+1
                                        else:
                                            if reverse_tmd == False:
                                                list_of_gaps_extracellular.append(n-following_gap)
                                                following_gap = following_gap+1
                                                rev_value = n
                                            else:
                                                list_of_gaps_intracellular.append(n-following_gap)
                                                following_gap = following_gap+1
                                                rev_value = n
                                # # juxta after odd TMDs:
                                # if type(dfnon.loc[hit,"seq_juxta_after_TM%.2d_in_query"%tmd_int])== str:
                                #     list_of_gaps_in_query_after_odd = [m.start()+0.5 for m in re.finditer("-",dfnon.loc[hit,"seq_juxta_after_TM%.2d_in_query"%tmd_int])if m.start()+0.5 < 31]
                                #
                                #     #if one gap is found, code checks location and appends it
                                #     if len (list_of_gaps_in_query_after_odd)==1:
                                #        if reverse_tmd == False:
                                #            list_of_gaps_extracellular.append(list_of_gaps_in_query_after_odd[0])
                                #        else:
                                #            list_of_gaps_intracellular.append(list_of_gaps_in_query_after_odd[0])
                                #     #if more than one gap is found, code checks if the gaps are one after another in the query!
                                #     if len (list_of_gaps_in_query_after_odd)>1.0:
                                #        following_gap = 0
                                #        rev_value = list_of_gaps_in_query_after_odd[0]
                                #        for n in list_of_gaps_in_query_after_odd:
                                #            if n+following_gap == rev_value:
                                #                if reverse_tmd == False:
                                #     #list_of_gaps_extracellular.append(n-following_gap)
                                #                    following_gap = following_gap+1
                                #                else:
                                #                    list_of_gaps_intracellular.append(n-following_gap)
                                #                    following_gap = following_gap+1
                                #            else:
                                #                if reverse_tmd == False:
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
                            if type(dfnon.loc[hit,"seq_juxta_before_%s_in_query"%tmd])== str:

                                # list of gap indices
                                list_of_gaps_in_query_before_even = [m.start()+0.5 for m in re.finditer("-",dfnon.loc[hit,"seq_juxta_before_%s_in_query"%tmd])if m.start()+0.5 < 31]

                                # if one gap is found, code checks location and appends it
                                if len (list_of_gaps_in_query_before_even)==1:
                                    if reverse_tmd == False:
                                        list_of_gaps_extracellular.append(list_of_gaps_in_query_before_even[0])
                                    else:
                                        list_of_gaps_intracellular.append(list_of_gaps_in_query_before_even[0])

                                # if more than one gap is found, code checks if the gapy are one after another in the query!
                                if len (list_of_gaps_in_query_before_even)>1.0:
                                    following_gap = 0
                                    rev_value = list_of_gaps_in_query_before_even[0]
                                    for n in list_of_gaps_in_query_before_even:
                                        if n+following_gap == rev_value:
                                            if reverse_tmd == False:
                                                list_of_gaps_extracellular.append(n-following_gap)
                                                following_gap = following_gap+1
                                            else:
                                                list_of_gaps_intracellular.append(n-following_gap)
                                                following_gap = following_gap+1
                                        else:
                                            if reverse_tmd == False:
                                                list_of_gaps_extracellular.append(n-following_gap)
                                                following_gap = following_gap+1
                                                rev_value = n
                                            else:
                                                list_of_gaps_intracellular.append(n-following_gap)
                                                following_gap = following_gap+1
                                                rev_value = n

                            if type(dfnon.loc[hit,"seq_juxta_before_%s_in_match"%tmd])== str:
                            # Makes a list of gaps of the match of the odd juxta before the TMD
                                list_of_gaps_in_query_before_even = [m.start()+0.5 for m in re.finditer("-",dfnon.loc[hit,"seq_juxta_before_%s_in_query"%tmd])if m.start()+0.5 < 31]

                                list_of_gaps_in_match_before_even = [m.start()+1 for m in re.finditer("-",dfnon.loc[hit,"seq_juxta_before_%s_in_match"%tmd])if m.start() < 31]

                                for n in list(reversed(list_of_gaps_in_match_before_even)):
                                    greater_values = sum(i< n for i in list_of_gaps_in_query_before_even)
                                    if reverse_tmd== False:
                                        list_of_gaps_extracellular.append(n-greater_values)
                                    else:
                                        list_of_gaps_intracellular.append(n-greater_values)



                            # juxta after even TMDs:

                            if type(dfnon.loc[hit,"seq_juxta_after_%s_in_query"%tmd])== str:

                                list_of_gaps_in_query_after_even = [m.start()+0.5 for m in re.finditer("-",dfnon.loc[hit,"seq_juxta_after_%s_in_query"%tmd][::-1]) if m.start()+0.5 < 31]

                                # if one gap is found, code checks location and appends it
                                if len (list_of_gaps_in_query_after_even)==1:
                                    if reverse_tmd == False:
                                        list_of_gaps_intracellular.append(list_of_gaps_in_query_after_even[0])
                                    else:
                                        list_of_gaps_extracellular.append(list_of_gaps_in_query_after_even[0])

                                # if more than one gap is found, code checks if the gaps are one after another in the query!
                                if len (list_of_gaps_in_query_after_even)>1.0:
                                    following_gap = 0
                                    rev_value = list_of_gaps_in_query_after_even[0]
                                    for n in list_of_gaps_in_query_after_even:
                                        if n-following_gap == rev_value:
                                            if reverse_tmd == False:
                                                list_of_gaps_intracellular.append(n-following_gap)
                                                following_gap = following_gap+1
                                            else:
                                                list_of_gaps_extracellular.append(n-following_gap)
                                                following_gap = following_gap+1
                                        else:
                                            if reverse_tmd == False:
                                                list_of_gaps_intracellular.append(n-following_gap)
                                                following_gap = following_gap+1
                                                rev_value = n
                                            else:
                                                list_of_gaps_extracellular.append(n-following_gap)
                                                following_gap = following_gap+1
                                                rev_value = n

                                           # Makes a list of gaps of the match of the odd juxta before the TMD

                            if type(dfnon.loc[hit,"seq_juxta_after_%s_in_match"%tmd])== str and type(dfnon.loc[hit,"seq_juxta_after_%s_in_query"%tmd])==str:

                                list_of_gaps_in_query_after_even = [m.start()+0.5 for m in re.finditer("-",dfnon.loc[hit,"seq_juxta_after_%s_in_query"%tmd][::-1])if m.start()+0.5 < 31]
                                list_of_gaps_in_match_after_even = [m.start()+1 for m in re.finditer("-",dfnon.loc[hit,"seq_juxta_after_%s_in_match"%tmd][::-1])if m.start() < 31]

                                for n in list(reversed(list_of_gaps_in_match_after_even)):
                                    greater_values = sum(i< n for i in list_of_gaps_in_query_after_even)
                                    if reverse_tmd== False:
                                        list_of_gaps_intracellular.append(n-greater_values)

                                    else:
                                        list_of_gaps_extracellular.append(n-greater_values)
                # sets of lists are created, to assure, that each gapposition contributes only once to the possible gap positions
                unique_list_of_gaps_in_tmd = list(set(list_of_gaps_in_tmd))
                unique_list_of_gaps_intracellular = list(set(list_of_gaps_intracellular))
                unique_list_of_gaps_extracellular = list(set(list_of_gaps_extracellular))

                # Saves the calculated lists into cells in the columns
                df.loc[acc,"%s_occuring_gaps"%tmd]=str(unique_list_of_gaps_in_tmd)
                df.loc[acc,"%s_amount_possible_gappositions"%tmd]=len(unique_list_of_gaps_in_tmd)

                df.loc[acc,'juxta_%s_intracellular_possible_gappositions'%tmd] = str(unique_list_of_gaps_intracellular)
                df.loc[acc,'juxta_%s_extracellular_possible_gappositions'%tmd] = str(unique_list_of_gaps_extracellular)
                df.loc[acc,'juxta_%s_intracellular_num_gaps'%tmd] = len(unique_list_of_gaps_intracellular)
                df.loc[acc,'juxta_%s_exracellular_num_gaps'%tmd] = len(unique_list_of_gaps_extracellular)
            # At the end, sets analysed to true, this is important to not overwrite
            df.loc[acc,"gaps_analysed"] = "True"
            logging.info("--Analysed")
            logging.info("~~~~~~~~~~~~         calculate_gap_densities is finished          ~~~~~~~~~~~~")
            # save to csv
            df.to_csv(pathdict["list_summary_csv"], sep=",", quoting=csv.QUOTE_NONNUMERIC)

        else:
            logging.info("Gap analysis for %s already done" %acc)