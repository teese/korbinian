import korbinian
import numpy as np
import os
import pandas as pd
# import debugging tools
from korbinian.utils import pr, pc, pn, aaa

def calc_AAIMON(TMD, df_cr, len_query_TMD):
    """Calculates the amino acid identity, membranous over nonmembranous (AAIMON) ratio for each homologue, and the average for all homologues of that protein.

    TM01_AAIMON = TM01_perc_ident / nonTMD_perc_ident

    Note that there are several ways of calculating the TM percentage identity, and the nonTMD percentage identity.
    The number of identical residues is easy:
        Number of identical residues = number of pipes in markup

    Percentage identity = Number of identical residues / length

    HOWEVER. The LENGTH can be calculated in different ways.
        - length of query excluding gaps
        - length of query including gaps
        - length of match excluding gaps
        - length of alignment (length of query excluding gaps + number gaps in query + number gaps in match)

    Parameters
    ----------
    TMD : str
        String denoting transmembrane domain number (e.g. "TM01")
    df_cr : pd.DataFrame
        Dataframe with conservation ratios for a particular TMD (or region).
    len_query_TMD : int
        Length of the query TMD sequence.
    mean_ser : pd.Series
        Series containing the mean values for all homologues of a single protein.
        Will be saved as a csv.
        The csv files for each protein will be gathered to create a single dataframe.

    Returns
    -------
    mean_ser : pd.Series
        Returns mean_ser with extra entries
    df_cr : pd.DataFrame
        Returns the dataframe with the extra calculated AAIMON ratios.
    """
    ########################################################################################
    #                                                                                      #
    #              count number of identical (and/or similar) residues in                  #
    #                  the TMD region of interest, for each homologue                      #
    #                                                                                      #
    ########################################################################################
    # count identical residues between query and match TMDs by counting the number of pipes in the markup string
    # NOTE THAT df_cr['%s_SW_markup_seq'%TMD].str.count('|') DOES NOT WORK, as "|" has a function in regex and needs to be escaped
    df_cr['%s_SW_num_ident_res' % TMD] = df_cr['%s_SW_markup_seq' % TMD].str.count('\|')
    df_cr['%s_SW_num_sim_res' % TMD] = df_cr['%s_SW_markup_seq' % TMD].str.count(':')

    ########################################################################################
    #                                                                                      #
    #     calculate the length of the alignment of the TMD region, including gaps          #
    #          (should equal the len_query_TMD + n_gaps_in_q + n_gaps_in_m)                #
    #                                                                                      #
    ########################################################################################
    df_cr['%s_SW_align_len' % TMD] = df_cr['%s_SW_match_seq' % TMD].str.len()

    ##############################################################################################################
    #                                                                                                            #
    #            calculate the length of the alignment of the TMD region, EXCLUDING GAPS                         #
    #    TM01_SW_align_len_excl_gaps = TM01_SW_align_len - TM01_SW_query_num_gaps - TM01_SW_match_num_gaps       #
    #                                                                                                            #
    ##############################################################################################################
    df_cr['%s_SW_align_len_excl_gaps' % TMD] = df_cr['%s_SW_align_len' % TMD] - df_cr['%s_SW_query_num_gaps' % TMD] - df_cr['%s_SW_match_num_gaps' % TMD]
    
    ########################################################################################
    #                                                                                      #
    #              calculate the percentage identity of the TMD region                     #
    #       TM01_perc_ident = TM01_SW_num_ident_res / TM01_SW_align_len_excl_gaps          #
    #                                                                                      #
    ########################################################################################
    # the percentage identity of that TMD is defined as the number of identical residues (pipes in markup) 
    # divided by the length of the the aligned residues (excluding gaps)
    # note that the nonTMD percentage identity is calculated the same way
    df_cr['%s_perc_ident' % TMD] = df_cr['%s_SW_num_ident_res' % TMD] / df_cr['%s_SW_align_len_excl_gaps' % TMD]
    # calculate percentage similar residues
    df_cr['%s_perc_sim' % TMD] = df_cr['%s_SW_num_sim_res' % TMD] / df_cr['%s_SW_align_len_excl_gaps' % TMD]
    # add together to obtain the percentage similar + identical residues
    df_cr['%s_perc_sim_plus_ident' % TMD] = df_cr['%s_perc_ident' % TMD] + df_cr['%s_perc_sim' % TMD]

    ########################################################################################
    #                                                                                      #
    #          calculate Amino Acid Identity : Membranous Over Nonmembranous               #
    #             (number of gaps)/(length of sequence excluding gaps)                     #
    #                                                                                      #
    ########################################################################################
    df_cr['%s_AAIMON'%TMD] = df_cr['%s_perc_ident'%TMD] / df_cr['nonTMD_perc_ident']
    # calculate the Amino Acid Similarity : Membranous Over Nonmembranous (AASMON) (includes similarity + identity based on the matrix used in the SW alignment of SIMAP)
    df_cr['%s_AASMON'%TMD] = df_cr['%s_perc_sim_plus_ident'%TMD] / df_cr['nonTMD_perc_sim_plus_ident']
    # calculate the AAIMON normalised by the random_AA_identity, to exclude identity due to lipophilicity
    df_cr['%s_AAIMON_n' % TMD] = df_cr['%s_AAIMON' % TMD] / df_cr['norm_factor']

    ########################################################################################
    #                                                                                      #
    #  calculate ratio of length of TMD to length of nonTMD excl gaps & full match seq     #
    #                                                                                      #
    ########################################################################################
    df_cr['%s_ratio_len_TMD_to_len_nonTMD'%TMD] = len_query_TMD / df_cr['nonTMD_SW_align_len_excl_gaps']
    df_cr['%s_ratio_len_TMD_to_len_full_match_seq'%TMD] = len_query_TMD / df_cr['len_full_match_seq']

    return df_cr

def calc_nonTMD_perc_ident_and_gaps(df_nonTMD, mean_ser, len_nonTMD_orig_q, rand_nonTM):
    """Calculate the nonTMD percentage identity and gaps.

    Parameters
    ----------
    df_nonTMD : pd.DataFrame
        Dataframe containing the nonTMD sequences, % identity, number of gaps, etc.
        columns : nonTMD_perc_ident etc
        index : homologue number
    mean_ser : pd.Series
        Series containing the mean values for all homologues of a single protein.
        Will be saved as a csv.
        The csv files for each protein will be gathered to create a single dataframe.
    len_nonTMD_orig_q : int
        Length of the nonTMD region in the original UniProt or OMPdb sequence.
    rand_nonTM : float
        Random identity of the nonTM region, based on the AA propensity of nonTM in the dataset.

    Returns
    -------
    df_nonTMD, mean_ser (updated dataframe and dictionary)
    """
    # get the length of the nonTMD region
    df_nonTMD['nonTMD_len'] = df_nonTMD['nonTMD_seq_match'].str.len()

    ########################################################################################
    #                                                                                      #
    #              count number of identical (and/or similar) residues in                  #
    #                      the nonTMD region, for each homologue                           #
    #                                                                                      #
    ########################################################################################
    # calculate identical residues, simply count the pipes '|' in the markup sequence)
    # NOTE THAT df_nonTMD['nonTMD_markup'].str.count('|') DOES NOT WORK, as "|" has a function in regex and needs to be escaped
    df_nonTMD['nonTMD_num_ident_res'] = df_nonTMD['nonTMD_markup'].str.count('\|')
    # calculate similar residues in the nonTMD region (simply count the colons ':' in the markup sequence)
    df_nonTMD['nonTMD_num_sim_res'] = df_nonTMD['nonTMD_markup'].str.count(':')
    # add the identical and similar residues together to get the total number of similar + identical residues
    df_nonTMD['nonTMD_num_sim_plus_ident_res'] = df_nonTMD['nonTMD_num_ident_res'] + df_nonTMD['nonTMD_num_sim_res']

    ########################################################################################
    #                                                                                      #
    #     count number of gaps in the query and match nonTMD aligned sequence              #
    #                                                                                      #
    ########################################################################################
    # count the gaps in the nonTMD sequence of the query
    df_nonTMD['nonTMD_q_num_gaps'] = df_nonTMD['nonTMD_seq_query'].str.count('-')
    # count the gaps in the nonTMD sequence of the match
    df_nonTMD['nonTMD_m_num_gaps'] = df_nonTMD['nonTMD_seq_match'].str.count('-')

    ########################################################################################
    #                                                                                      #
    #    calculate the length of the alignment of the nonTMD alignment, including gaps     #
    #                                                                                      #
    ########################################################################################
    # calculate the length of the nonTMD sequences, which may include gaps
    df_nonTMD['nonTMD_SW_align_len'] = df_nonTMD['nonTMD_seq_query'].str.len()
    mean_ser['nonTMD_SW_align_len_mean'] = float('%0.2f' % df_nonTMD['nonTMD_SW_align_len'].mean()) # dropna(). removed
    # to avoid divide by 0 errors, any homologues with no nonTMD region should be excluded
    # DEPRECATED, seems to be redundant due to the code below (df_nonTMD['nonTMD_SW_align_len_excl_gaps'] = df_nonTMD['nonTMD_SW_align_len_excl_gaps'].replace(0,np.nan))
    #df_nonTMD['nonTMD_SW_align_len'] = df_nonTMD['nonTMD_SW_align_len'].replace(0, np.nan)

    ##############################################################################################################
    #                                                                                                            #
    #            calculate the length of the alignment of the nonTMD region, EXCLUDING GAPS                      #
    #       nonTMD_SW_align_len_excl_gaps = nonTMD_SW_align_len - nonTMD_q_num_gaps - nonTMD_m_num_gaps          #
    #                                                                                                            #
    ##############################################################################################################
    df_nonTMD['nonTMD_SW_align_len_excl_gaps'] = df_nonTMD['nonTMD_SW_align_len'] - df_nonTMD['nonTMD_q_num_gaps'] - df_nonTMD['nonTMD_m_num_gaps']
    # to avoid divide by 0 errors, any homologues with only gaps in the nonTMD region should be excluded
    df_nonTMD['nonTMD_SW_align_len_excl_gaps'] = df_nonTMD['nonTMD_SW_align_len_excl_gaps'].replace(0,np.nan)
    # add mean to output series
    mean_ser['nonTMD_SW_align_len_excl_gaps_mean'] = float('%0.2f' % df_nonTMD['nonTMD_SW_align_len_excl_gaps'].mean()) # dropna(). removed
    # len_nonTMD_orig_q = 100  # to be replaced with a real number!
    df_nonTMD['len_nonTMD_orig_q_minus_nonTMD_SW_align_len_excl_gaps'] = len_nonTMD_orig_q - df_nonTMD['nonTMD_SW_align_len_excl_gaps']
    df_nonTMD['perc_nonTMD_coverage'] = df_nonTMD['nonTMD_SW_align_len_excl_gaps'] / len_nonTMD_orig_q
    mean_ser['len_nonTMD_orig_q_minus_nonTMD_SW_align_len_excl_gaps_mean'] = float('%0.2f' % df_nonTMD['len_nonTMD_orig_q_minus_nonTMD_SW_align_len_excl_gaps'].mean()) # dropna(). removed

    ########################################################################################
    #                                                                                      #
    #            calculate the percentage identity of the nonTMD region                    #
    #      nonTMD_perc_ident = nonTMD_num_ident_res / nonTMD_SW_align_len_excl_gaps        #
    #                                                                                      #
    ########################################################################################
    # calculate the percentage identity of the nonTMD region (number of identical residues divided by the length excluding gaps)
    df_nonTMD['nonTMD_perc_ident_alignable_region'] = df_nonTMD['nonTMD_num_ident_res'] / df_nonTMD['nonTMD_SW_align_len_excl_gaps']

    ########################################################################################
    #                                                                                      #
    #                       adjust_for_trunc_ends_in_local_alignment                       #
    #       nonTMD % identity = pf + r(1 - f)                                              #
    #       where p = percentage identity of alignable nonTM sequence                      #
    #       f = fraction of alignable region (e.g. 0.93)                                   #
    #       r = random identity (e.g. 5.8%, 0.058)                                         #
    #                                                                                      #
    ########################################################################################
    adjust_for_trunc_ends_in_local_alignment = True
    if adjust_for_trunc_ends_in_local_alignment:
        df_nonTMD['nonTMD_perc_ident'] = df_nonTMD['nonTMD_perc_ident_alignable_region'] * df_nonTMD['SW_query_coverage'] + rand_nonTM * (1 - df_nonTMD['SW_query_coverage'])
    else:
        # older, deprecated method
        df_nonTMD['nonTMD_perc_ident'] = df_nonTMD['nonTMD_perc_ident_alignable_region']

    # drop any rows where there is no nonTMD_perc_ident
    df_nonTMD.dropna(subset=["nonTMD_perc_ident"], inplace=True)

    # note that the similarity ratio won't be adjusted for truncated ends! Not currently in use for any calculations or protein comparisons.
    df_nonTMD['nonTMD_perc_sim'] = df_nonTMD['nonTMD_num_sim_res'] / df_nonTMD['nonTMD_SW_align_len_excl_gaps']
    df_nonTMD['nonTMD_perc_sim_plus_ident'] = df_nonTMD['nonTMD_perc_ident_alignable_region'] +  df_nonTMD['nonTMD_perc_sim']

    # add to output dictionary with mean values for all homologues
    mean_ser['nonTMD_perc_ident_mean'] = float('%0.5f' % df_nonTMD['nonTMD_perc_ident'].mean())
    mean_ser['nonTMD_perc_sim_mean'] = float('%0.3f' % df_nonTMD['nonTMD_perc_sim'].mean())
    mean_ser['nonTMD_perc_sim_plus_ident_mean'] = float('%0.3f' % df_nonTMD['nonTMD_perc_sim_plus_ident'].mean())

    ########################################################################################
    #                                                                                      #
    #               calculate the percentage gaps of the nonTMD region                     #
    #                                                                                      #
    ########################################################################################
    df_nonTMD['perc_gaps_nonTMD_SW_align'] = (df_nonTMD['nonTMD_q_num_gaps'] + df_nonTMD['nonTMD_m_num_gaps']) / df_nonTMD['nonTMD_SW_align_len']
    mean_ser['perc_gaps_nonTMD_SW_align_mean'] = float('%0.2f' % df_nonTMD['perc_gaps_nonTMD_SW_align'].mean())

    ########################################################################################
    #                                                                                      #
    # OLD: SSR ratio calculations take a long time and show no real difference to AASMON   #
    #                                                                                      #
    ########################################################################################
    SSR_ratio_calculations = False
    if SSR_ratio_calculations:
        korbinian.cons_ratio.ssr.conduct_ssr_ratio_calculations(df_nonTMD, None, None, None)

    return df_nonTMD, mean_ser