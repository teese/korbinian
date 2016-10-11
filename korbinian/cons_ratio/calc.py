import korbinian
import numpy as np
import pandas as pd

def calc_AAIMON(TMD, df_cr, len_query_TMD, mean_ser):
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
    #        count number of gaps in the query and match TMD aligned sequence              #
    #                                                                                      #
    ########################################################################################
    """This is used as a filter in filter_and_save_fasta, therefore is conducted earlier in the slicing function. """
    ## count the number of gaps in the query and match sequences
    #df_cr['%s_SW_query_num_gaps' % TMD] = df_cr['%s_SW_query_seq' % TMD].str.count("-")
    #df_cr['%s_SW_match_num_gaps' % TMD] = df_cr['%s_SW_match_seq' % TMD].str.count("-")

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
    df_cr['%s_perc_sim' % TMD] = df_cr['%s_SW_num_sim_res' % TMD] / df_cr['%s_SW_align_len' % TMD]
    # add together to obtain the percentage similar + identical residues
    df_cr['%s_perc_sim_plus_ident' % TMD] = df_cr['%s_perc_ident' % TMD] + df_cr['%s_perc_sim' % TMD]
    # add to output dictionary with mean values for all homologues
    mean_ser['%s_perc_ident_mean' % TMD] = df_cr['%s_perc_ident' % TMD].mean()
    mean_ser['%s_perc_sim_mean' % TMD] = df_cr['%s_perc_sim' % TMD].mean()
    mean_ser['%s_perc_sim_plus_ident_mean' % TMD] = df_cr['%s_perc_sim_plus_ident' % TMD].mean()

    ########################################################################################
    #                                                                                      #
    #          calculate Amino Acid Identity : Membranous Over Nonmembranous               #
    #             (number of gaps)/(length of sequence excluding gaps)                     #
    #                                                                                      #
    ########################################################################################
    df_cr['%s_AAIMON_ratio'%TMD] = df_cr['%s_perc_ident'%TMD] / df_cr['nonTMD_perc_ident']
    # calculate the Amino Acid Similarity : Membranous Over Nonmembranous (AASMON) (includes similarity + identity based on the matrix used in the SW alignment of SIMAP)
    df_cr['%s_AASMON_ratio'%TMD] = df_cr['%s_perc_sim_plus_ident'%TMD] / df_cr['nonTMD_perc_sim_plus_ident']
    # add to output dictionary with mean values for all homologues
    mean_ser['%s_AAIMON_ratio_mean' % TMD] = float(df_cr['%s_AAIMON_ratio' % TMD].mean())
    mean_ser['%s_AAIMON_ratio_std' % TMD] = df_cr['%s_AAIMON_ratio' % TMD].std()
    mean_ser['%s_AASMON_ratio_mean' % TMD] = df_cr['%s_AASMON_ratio' % TMD].mean()
    mean_ser['%s_AASMON_ratio_std' % TMD] = df_cr['%s_AASMON_ratio' % TMD].std()

    ########################################################################################
    #                                                                                      #
    #     calculate the average number of gaps per residue in the TMD alignment            #
    #             (number of gaps)/(length of sequence excluding gaps)                     #
    #                                                                                      #
    ########################################################################################
    df_cr['%s_SW_q_gaps_per_q_residue' % TMD] = df_cr['%s_SW_query_num_gaps' % TMD].dropna() / len_query_TMD
    mean_ser['%s_SW_q_gaps_per_q_residue_mean' % TMD] = df_cr['%s_SW_q_gaps_per_q_residue' % TMD].dropna().mean()

    ########################################################################################
    #                                                                                      #
    #  calculate ratio of length of TMD to length of nonTMD excl gaps & full match seq     #
    #                                                                                      #
    ########################################################################################
    df_cr['%s_ratio_len_TMD_to_len_nonTMD'%TMD] = len_query_TMD / df_cr['nonTMD_SW_align_len_excl_gaps']
    df_cr['%s_ratio_len_TMD_to_len_full_match_seq'%TMD] = len_query_TMD / df_cr['len_full_match_seq']
    # add to output dictionary with mean values for all homologues
    mean_ser['%s_ratio_len_TMD_to_len_nonTMD_mean'%TMD] = float('%0.2f' % df_cr['%s_ratio_len_TMD_to_len_nonTMD'%TMD].dropna().mean())
    mean_ser['%s_ratio_len_TMD_to_len_full_match_seq_mean'%TMD] = float('%0.2f' % df_cr['%s_ratio_len_TMD_to_len_full_match_seq'%TMD].dropna().mean())

    ## calculate the length of the match TMD seq excluding gaps
    #df_cr['%s_SW_m_seq_len' % TMD] = df_cr['%s_SW_match_seq' % TMD].str.len()
    #len_query_TMD = len(query_TMD_sequence)
    # for the alignment length, take the smallest value from the length of query or match
    # this will exclude gaps from the length in the following calculations, preventing false "low conservation" where the query TMD is much longer than the match TMD)
    # note that for most calculations this is somewhat redundant, because the max number of acceptable gaps in sequence is probable ~2
    # use the mask function (faster than a lambda function) to replace any lengths larger than the query, with the query
    #df_cr['%s_SW_align_len' % TMD] = df_cr['%s_SW_m_seq_len' % TMD].mask(df_cr['%s_SW_m_seq_len' % TMD] > len_query_TMD, len_query_TMD)
    # # create a boolean column that allows filtering by the accepted number of gaps, according to the settings file
    # df_cr['%s_cr_SW_query_acceptable_n_gaps' % TMD] = df_cr['%s_SW_query_num_gaps' % TMD] <= s["cr_max_n_gaps_in_query_TMD"]
    # df_cr['%s_cr_SW_match_acceptable_n_gaps' % TMD] = df_cr['%s_SW_match_num_gaps' % TMD] <= s["cr_max_n_gaps_in_match_TMD"]
    # check that the TMD seq in match is not just 100% gaps!
    # df_cr['%s_in_SW_align_match' % TMD] = df_cr['%s_SW_num_ident_res' % TMD].dropna() != 0
    # df_cr['%s_in_SW_align_match' % TMD].fillna(value=False)
    #logging.info('%s_SW_q_gaps_per_q_residue Average: %0.3e' %(TMD, mean_ser['%s_SW_q_gaps_per_q_residue_mean'%TMD]))
    # add to original dataframe with the list of uniprot sequences
    return df_cr, mean_ser

def calc_nonTMD_perc_ident_and_gaps(df_nonTMD, mean_ser):
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

    Returns
    -------
    df_nonTMD, mean_ser (updated dataframe and dictionary)
    """
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
    mean_ser['nonTMD_SW_align_len_mean'] = float('%0.2f' % df_nonTMD['nonTMD_SW_align_len'].dropna().mean())
    # to avoid divide by 0 errors, any homologues with no nonTMD region should be excluded
    df_nonTMD['nonTMD_SW_align_len'] = df_nonTMD['nonTMD_SW_align_len'].replace(0, np.nan)

    ##############################################################################################################
    #                                                                                                            #
    #            calculate the length of the alignment of the nonTMD region, EXCLUDING GAPS                      #
    #       nonTMD_SW_align_len_excl_gaps = nonTMD_SW_align_len - nonTMD_q_num_gaps - nonTMD_m_num_gaps          #
    #                                                                                                            #
    ##############################################################################################################
    df_nonTMD['nonTMD_SW_align_len_excl_gaps'] = df_nonTMD['nonTMD_SW_align_len'] - df_nonTMD['nonTMD_q_num_gaps'] - df_nonTMD['nonTMD_m_num_gaps']
    mean_ser['nonTMD_SW_align_len_excl_gaps_mean'] = float('%0.2f' % df_nonTMD['nonTMD_SW_align_len_excl_gaps'].dropna().mean())
    # to avoid divide by 0 errors, any homologues with only gaps in the nonTMD region should be excluded
    df_nonTMD['nonTMD_SW_align_len_excl_gaps'] = df_nonTMD['nonTMD_SW_align_len_excl_gaps'].replace(0,np.nan)

    ########################################################################################
    #                                                                                      #
    #            calculate the percentage identity of the nonTMD region                    #
    #      nonTMD_perc_ident = nonTMD_num_ident_res / nonTMD_SW_align_len_excl_gaps        #
    #                                                                                      #
    ########################################################################################
    # calculate the percentage identity of the nonTMD region (number of identical residues divided by the length excluding gaps)
    df_nonTMD['nonTMD_perc_ident'] = df_nonTMD['nonTMD_num_ident_res'] / df_nonTMD['nonTMD_SW_align_len_excl_gaps']
    df_nonTMD['nonTMD_perc_sim'] = df_nonTMD['nonTMD_num_sim_res'] / df_nonTMD['nonTMD_SW_align_len_excl_gaps']
    df_nonTMD['nonTMD_perc_sim_plus_ident'] = df_nonTMD['nonTMD_perc_ident'] +  df_nonTMD['nonTMD_perc_sim']
    # add to output dictionary with mean values for all homologues
    mean_ser['nonTMD_perc_ident_mean'] = float('%0.3f' % df_nonTMD['nonTMD_perc_ident'].dropna().mean())
    mean_ser['nonTMD_perc_sim_mean'] = float('%0.3f' % df_nonTMD['nonTMD_perc_sim'].dropna().mean())
    mean_ser['nonTMD_perc_sim_plus_ident_mean'] = float('%0.3f' % df_nonTMD['nonTMD_perc_sim_plus_ident'].dropna().mean())

    ########################################################################################
    #                                                                                      #
    #               calculate the percentage gaps of the nonTMD region                     #
    #                                                                                      #
    ########################################################################################
    df_nonTMD['perc_gaps_nonTMD_SW_align'] = (df_nonTMD['nonTMD_q_num_gaps'] - df_nonTMD['nonTMD_m_num_gaps']) / df_nonTMD['nonTMD_SW_align_len']
    mean_ser['perc_gaps_nonTMD_SW_align_mean'] = float('%0.2f' % df_nonTMD['perc_gaps_nonTMD_SW_align'].dropna().mean())

    """ LEGACY CODE TO BE DELETED LATER """
    # calculate the length of the nonTMD sequence excluding gaps
    #df_nonTMD['len_nonTMD_q_excl_gaps'] = df_nonTMD['nonTMD_SW_align_len'] - df_nonTMD['nonTMD_q_num_gaps']
    #df_nonTMD['len_nonTMD_m_excl_gaps'] = df_nonTMD['len_nonTMD_seq_match'] - df_nonTMD['nonTMD_m_num_gaps']
    # calculate the length of the alignment by finding which seq excl gaps is smaller
    # NOTE, THIS IS CURRENTLY TOO SHORT, GIVING NONTMD IDENTITIES ALWAYS ABOVE 1.0. NEEDS TO BE FIXED.
    #df_nonTMD['len_nonTMD_align'] = df_nonTMD[['len_nonTMD_q_excl_gaps', 'len_nonTMD_m_excl_gaps']].min(axis=1)
    # calculate the average number of gaps per residue in the nonTMD alignment
    # # filter to analyse only sequences that are valid (length > 0)
    # df_nonTMD_filt_gaps = df_nonTMD.loc[df_nonTMD['len_nonTMD_q_excl_gaps'] != 0]
    # # calculate number of gaps in query AND match
    # df_nonTMD_filt_gaps['nonTMD_qm_num_gaps'] = df_nonTMD_filt_gaps['nonTMD_q_num_gaps'] + df_nonTMD_filt_gaps['nonTMD_m_num_gaps']
    # # add to simap dataframe
    # df_nonTMD['nonTMD_qm_num_gaps'] = df_nonTMD_filt_gaps['nonTMD_qm_num_gaps']
    # # gaps per query residue for both query and match = ((gaps in query + gaps in match)/2))/length of query excluding gaps
    # df_nonTMD['nonTMD_qm_gaps_per_q_residue'] = df_nonTMD_filt_gaps['nonTMD_qm_num_gaps'] / 2 / df_nonTMD_filt_gaps['len_nonTMD_q_excl_gaps']

    ########################################################################################
    #                                                                                      #
    #    OLD: SSR ratio calculations take a long time and show no difference to AAIMON     #
    #                                                                                      #
    ########################################################################################
    SSR_ratio_calculations = False
    if SSR_ratio_calculations:
        korbinian.cons_ratio.ssr.conduct_ssr_ratio_calculations(df_nonTMD, None, None, None)

    return df_nonTMD, mean_ser