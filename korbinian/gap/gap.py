import ast
import pandas as pd
import numpy as np
import csv
import os
import re
import tarfile
import korbinian.mtutils as utils

def calculate_gap_densities(pathdict, settingsdict, logging):
    logging.info('~~~~~~~~~~~~starting calculate_gap_densities~~~~~~~~~~~~')

    # If script previously has been run, continues with proteins not being processed yet, or overwrites previous gap analysis
    overwrite_previous_gap_analysis = settingsdict["variables"]["analyse.calculate_gap_densities.overwrite_previous_gap_analysis"]

    # Maximum number of gaps for tmds to be considered
    allowed_gaps_per_tmd = settingsdict["variables"]["analyse.calculate_gap_densities.allowed_gaps_per_tmd"]

    # 24 for beta barrel proteins, can be altered if only several TMDs to consider
    max_number_of_tmds = settingsdict["variables"]["analyse.calculate_gap_densities.max_number_of_tmds"]

    if os.path.isfile(pathdict["dfout10_uniprot_gaps"]):
        df = pd.read_csv(pathdict["dfout10_uniprot_gaps"], sep=",", quoting=csv.QUOTE_NONNUMERIC, index_col=[0])
        logging.info('df loaded from %s' % pathdict["dfout10_uniprot_gaps"])
    elif os.path.isfile(pathdict["dfout08_simap_AAIMON"]):
        df = pd.read_csv(pathdict["dfout08_simap_AAIMON"], sep=",", quoting=csv.QUOTE_NONNUMERIC, index_col=[0])
        logging.info('no gap file found, df loaded from %s' % pathdict["dfout08_simap_AAIMON"])
        # If no previous analysis had been done, uniprot csv is opened and additional columns are created
        df["gaps_analysed"]= np.nan     # Column, containing true or false
        for n in range (1,max_number_of_tmds+1):
            df["TM%.2d_occuring_gaps"%n]=np.nan     # List of unique gappositions, which occur in the tmd
            df["TM%.2d_amount_possible_gappositions"%n]=np.nan
            df['total_amount_of_TM%.2d'%n] = np.nan     # How often TMD is considered
            df['juxta_TM%.2d_intracellular_possible_gappositions'%n] = np.nan   # List of unique gappositions, which occur in the intracellular loop-part/end
            df['juxta_TM%.2d_extracellular_possible_gappositions'%n] = np.nan   # List of unique gappositions, which occur in the intracellular loop-part/end
            df['juxta_TM%.2d_intracellular_num_gaps'%n] = np.nan
            df['juxta_TM%.2d_exracellular_num_gaps'%n] = np.nan
            df['len_juxta_TM%.2d_intracellular'%n] = np.nan
            df['len_juxta_TM%.2d_extracellular'%n] = np.nan
    else:
        raise IOError("df is not in memory, and neither %s nor %s are found" % (pathdict["dfout10_uniprot_gaps"],pathdict["dfout08_simap_AAIMON"]))

    for acc in df.index:

        protein_name = df.loc[acc,'A2_protein_name']
        # The next steps (the main analysis) is only executed, if previous analysis can be overwritten or no analysis has yet been done
        if (overwrite_previous_gap_analysis == True) or (df.loc[acc,"gaps_analysed"] != True):
            logging.info("%s"%acc)
            list_of_TMDs = ast.literal_eval(df.loc[acc,"list_of_TMDs"])

            # Checks if outputfiles (tar) exist
            if os.path.exists(df.loc[acc,'output_tarfile_path']):

            # opens the analysed csv for each protein and loads it into a dataframe
                with tarfile.open(df.loc[acc,'output_tarfile_path'], mode= 'r:gz')as tar:

                # checks if the analysed file exists, otherwise prints that it does not exist
                    if '%s_analysed.csv' % protein_name in tar.getnames():

                    # loads file into analysed csv
                        analysed_csv = tar.extractfile('%s_analysed.csv' % protein_name)
                        analysed_df = pd.read_csv(analysed_csv,low_memory=False,index_col=[0])

                        # checks if first amino acid is located inside (or periplasmatic) or outside, returns a boolean, true or false
                        # if first residue is located inside, every even tmd (tmd2,tmd4,tmd6...) is reversed, otherwise every odd tmd is reversed
                        # output is a boolean for each tmd, depending on the number and on the first amino acid

                        if df.loc[acc,"n_term_ec"] == False:
                            reverse_tmd = False
                        else:
                            reverse_tmd = True
                        print (reverse_tmd)


                        # for each TMD in the proteins, creates new lists which will contain gappositions, lists are saved in a column and created again for each tmd
                        for tmd in list_of_TMDs:
                            print(tmd)
                            tmd_int = int(tmd[-2:]) # Integer of TMD number
                            len_of_query = len(analysed_df["%s_SW_query_seq"%tmd][1]) # Length of first query sequence, which does (usually) not contain any gaps
                            len_of_query_reversed= ((1/len_of_query)+1) # Reversed length, important if TMD needs to be reversed afterwards
                            list_of_gaps_in_tmd = []
                            list_of_gaps_intracellular = []
                            list_of_gaps_extracellular = []


                            for hit in analysed_df.index:

                            #'''
                            #Start of the main gap analysis
                            #Code searches for "-" in the TMD sequence and returns the index!! (not the position)
                            #'''

                                # Following if conditions only refer to gaps in the query!
                                # Query gaps are counted as "in between positions", for example: 4,5 refers to a gap between position 4 and 5;
                                # if two gaps occur one after another: only one position (between two amino acids is considered)

                                # Filter to make sure, that there are 1 or 2 gaps in the query sequence and up to the max allowed gaps in the match
                                if (analysed_df.loc[hit,"%s_SW_query_num_gaps"%tmd] != 0.0) and (analysed_df.loc[hit,"%s_SW_query_num_gaps"%tmd] <= 2.0)\
                                    and (analysed_df.loc[hit,"%s_SW_match_num_gaps"%tmd] <= int("%s"%allowed_gaps_per_tmd)):

                                    # Stores the endpoints in a temp list; endpoints are used, to switch from python indices to numbers
                                    list_of_gaps_per_hit_in_query = [m.start() for m in re.finditer("-",analysed_df.loc[hit,"%s_SW_query_seq"%tmd]) if m.start()]
                                    print (list_of_gaps_per_hit_in_query)

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
                                if (analysed_df.loc[hit,"%s_SW_query_num_gaps"%tmd] <=2.0) and (analysed_df.loc[hit,"%s_SW_match_num_gaps"%tmd] <= 2.0)\
                                    and (analysed_df.loc[hit,"%s_SW_match_num_gaps"%tmd] != 0.0):

                                    # It's not sure that the list of hits in query was already determined, maybe there were no gaps, anyway here it is important how many
                                    list_of_gaps_per_hit_in_query = [m.start() for m in re.finditer("-",analysed_df.loc[hit,"%s_SW_query_seq"%tmd]) if m.start()]
                                    list_of_gaps_per_hit_in_match = [m.start() for m in re.finditer("-",analysed_df.loc[hit,"%s_SW_match_seq"%tmd])if m.start()]
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

                                if utils.isOdd(tmd_int)==True:  # also für 1 , 3 ...


                                    # makes sure that the search is done in a string
                                    if type(analysed_df.loc[hit,"seq_juxta_before_TM%.2d_in_query"%tmd_int])== str:

                                        # list of gap indices
                                        list_of_gaps_in_query_before_odd = [m.start()+0.5 for m in re.finditer("-",analysed_df.loc[hit,"seq_juxta_before_TM%.2d_in_query"%tmd_int][::-1])if m.start()+0.5 < 31]

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


                                        if type(analysed_df.loc[hit,"seq_juxta_before_TM%.2d_in_match"%tmd_int])== str:
                                        # Makes a list of gaps of the match of the odd juxta before the TMD
                                            list_of_gaps_in_query_before_odd = [m.start()+1 for m in re.finditer("-",analysed_df.loc[hit,"seq_juxta_before_TM%.2d_in_query"%tmd_int][::-1])if m.start() < 32]

                                            list_of_gaps_in_match_before_odd = [m.start()+1 for m in re.finditer("-",analysed_df.loc[hit,"seq_juxta_before_TM%.2d_in_match"%tmd_int][::-1])if m.start() < 32]
                                            # This step is essential, to control, if there is a gap before, in the query region
                                            for n in list(reversed(list_of_gaps_in_match_before_odd)):
                                                greater_values = sum(i< n for i in list_of_gaps_in_query_before_odd)
                                                if reverse_tmd== False:
                                                    list_of_gaps_intracellular.append(n-greater_values)

                                                else:
                                                    list_of_gaps_extracellular.append(n-greater_values)

# juxta after odd TMDs:

                                        if type(analysed_df.loc[hit,"seq_juxta_after_TM%.2d_in_query"%tmd_int])== str:

                                            list_of_gaps_in_query_after_odd = [m.start()+0.5 for m in re.finditer("-",analysed_df.loc[hit,"seq_juxta_after_TM%.2d_in_query"%tmd_int])if m.start()+0.5 < 31]

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
     # juxta after odd TMDs:


#                               if type(analysed_df.loc[hit,"seq_juxta_after_TM%.2d_in_query"%tmd_int])== str:

                         #                   list_of_gaps_in_query_after_odd = [m.start()+0.5 for m in re.finditer("-",analysed_df.loc[hit,"seq_juxta_after_TM%.2d_in_query"%tmd_int])if m.start()+0.5 < 31]

                                            # if one gap is found, code checks location and appends it
                        #                    if len (list_of_gaps_in_query_after_odd)==1:
                        #                        if reverse_tmd == False:
                        #                            list_of_gaps_extracellular.append(list_of_gaps_in_query_after_odd[0])
                        #                        else:
                        #                            list_of_gaps_intracellular.append(list_of_gaps_in_query_after_odd[0])

                                            # if more than one gap is found, code checks if the gaps are one after another in the query!
                        #                    if len (list_of_gaps_in_query_after_odd)>1.0:
                         #                       following_gap = 0
                         #                       rev_value = list_of_gaps_in_query_after_odd[0]
                          #                      for n in list_of_gaps_in_query_after_odd:
                          #                          if n+following_gap == rev_value:
                          #                              if reverse_tmd == False:
                          ##                                  list_of_gaps_extracellular.append(n-following_gap)
                          #                                  following_gap = following_gap+1
                          #                              else:
                          #                                  list_of_gaps_intracellular.append(n-following_gap)
                          #                                  following_gap = following_gap+1
                         #                           else:
                         #                               if reverse_tmd == False:
                         #                                   list_of_gaps_extracellular.append(n-following_gap)
                        #                                    following_gap = following_gap+1
                         #                                   rev_value = n
                         #                               else:
                         #                                   list_of_gaps_intracellular.append(n-following_gap)
                         #                                   following_gap = following_gap+1
                         #                                   rev_value = n

                                    else:  # for 2,4

                                    # juxta before even TMDs:

                                    # makes sure that the search is done in a string
                                        if type(analysed_df.loc[hit,"seq_juxta_before_TM%.2d_in_query"%tmd_int])== str:

                                            # list of gap indices
                                            list_of_gaps_in_query_before_even = [m.start()+0.5 for m in re.finditer("-",analysed_df.loc[hit,"seq_juxta_before_TM%.2d_in_query"%tmd_int])if m.start()+0.5 < 31]

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

                                        if type(analysed_df.loc[hit,"seq_juxta_before_TM%.2d_in_match"%tmd_int])== str:
                                        # Makes a list of gaps of the match of the odd juxta before the TMD
                                            list_of_gaps_in_query_before_even = [m.start()+0.5 for m in re.finditer("-",analysed_df.loc[hit,"seq_juxta_before_TM%.2d_in_query"%tmd_int])if m.start()+0.5 < 31]

                                            list_of_gaps_in_match_before_even = [m.start()+1 for m in re.finditer("-",analysed_df.loc[hit,"seq_juxta_before_TM%.2d_in_match"%tmd_int])if m.start() < 31]

                                            for n in list(reversed(list_of_gaps_in_match_before_even)):
                                                greater_values = sum(i< n for i in list_of_gaps_in_query_before_even)
                                                if reverse_tmd== False:
                                                    list_of_gaps_extracellular.append(n-greater_values)
                                                else:
                                                    list_of_gaps_intracellular.append(n-greater_values)



                                        # juxta after even TMDs:

                                        if type(analysed_df.loc[hit,"seq_juxta_after_TM%.2d_in_query"%tmd_int])== str:

                                            list_of_gaps_in_query_after_even = [m.start()+0.5 for m in re.finditer("-",analysed_df.loc[hit,"seq_juxta_after_TM%.2d_in_query"%tmd_int][::-1]) if m.start()+0.5 < 31]

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

                                        if type(analysed_df.loc[hit,"seq_juxta_after_TM%.2d_in_match"%tmd_int])== str and type(analysed_df.loc[hit,"seq_juxta_after_TM%.2d_in_query"%tmd_int])==str:

                                            list_of_gaps_in_query_after_even = [m.start()+0.5 for m in re.finditer("-",analysed_df.loc[hit,"seq_juxta_after_TM%.2d_in_query"%tmd_int][::-1])if m.start()+0.5 < 31]
                                            list_of_gaps_in_match_after_even = [m.start()+1 for m in re.finditer("-",analysed_df.loc[hit,"seq_juxta_after_TM%.2d_in_match"%tmd_int][::-1])if m.start() < 31]

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
                        with open(pathdict["dfout10_uniprot_gaps"], 'w') as csv_out:
                            df.to_csv(csv_out, sep=",", quoting=csv.QUOTE_NONNUMERIC)


                    else:
                        logging.info("Analysed csv for %s does not exist" %acc)
            else:
                logging.info("Output file for %s does not exist" %acc)

        else:
            logging.info("Gap analysis for %s already done" %acc)