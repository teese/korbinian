import korbinian
import pandas as pd

def calc_AAIMON(TMD, df_cr, mean_ser):
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
    mean_ser : dict
        Dictionary containing the mean values for all homologues of a single protein.
        Will be saved as a csv.
        The csv files for each protein will be gathered to create a single dataframe.
    logging : logging.Logger
        Logger for printing to console and logfile.

    Returns
    -------
    mean_ser : dict
        Returns mean_ser with extra entries
    df_cr : pd.DataFrame
        Returns the dataframe with the extra calculated AAIMON ratios.
    """
    # calculate the Amino Acid Identity : Membranous Over Nonmembranous
    df_cr['%s_AAIMON_ratio'%TMD] = df_cr['%s_perc_ident'%TMD] / df_cr['nonTMD_perc_ident']
    # calculate the Amino Acid Similarity : Membranous Over Nonmembranous (AASMON) (includes similarity + identity based on the matrix used in the SW alignment of SIMAP)
    df_cr['%s_AASMON_ratio'%TMD] = df_cr['%s_perc_sim_plus_ident'%TMD] / df_cr['nonTMD_perc_sim_plus_ident']
    mean_ser['%s_SW_q_gaps_per_q_residue_mean'%TMD] = df_cr['%s_SW_q_gaps_per_q_residue'%TMD].dropna().mean()
    #logging.info('%s_SW_q_gaps_per_q_residue Average: %0.3e' %(TMD, mean_ser['%s_SW_q_gaps_per_q_residue_mean'%TMD]))

    # add to original dataframe with the list of uniprot sequences
    mean_ser['%s_perc_ident_mean'%TMD] = df_cr['%s_perc_ident'%TMD].mean()
    mean_ser['%s_perc_sim_mean'%TMD] = df_cr['%s_perc_sim'%TMD].mean()
    mean_ser['%s_perc_sim_plus_ident_mean'%TMD] = df_cr['%s_perc_sim_plus_ident'%TMD].mean()
    mean_ser['%s_AAIMON_ratio_mean'%TMD] = float(df_cr['%s_AAIMON_ratio'%TMD].mean())
    mean_ser['%s_AAIMON_ratio_std'%TMD] = df_cr['%s_AAIMON_ratio'%TMD].std()
    mean_ser['%s_AASMON_ratio_mean'%TMD] = df_cr['%s_AASMON_ratio'%TMD].mean()
    mean_ser['%s_AASMON_ratio_std'%TMD] = df_cr['%s_AASMON_ratio'%TMD].std()

    df_cr['%s_ratio_length_of_TMD_to_rest_of_alignment'%TMD] = df_cr['%s_SW_query_seq'%TMD].str.len() / df_cr['FASTA_overlap']
    df_cr['%s_ratio_length_of_query_TMD_to_rest_of_match_protein'%TMD] = df_cr['%s_SW_query_seq'%TMD].str.len() / df_cr[ 'len_full_match_seq']

    mean_ser['%s_ratio_length_of_TMD_to_rest_of_alignment_mean'%TMD] = float('%0.2f' % df_cr['%s_ratio_length_of_TMD_to_rest_of_alignment'%TMD].dropna().mean())
    mean_ser['%s_ratio_length_of_query_TMD_to_rest_of_match_protein_mean'%TMD] = float('%0.2f' % df_cr['%s_ratio_length_of_query_TMD_to_rest_of_match_protein'%TMD].dropna().mean())
    return mean_ser, df_cr

def calc_nonTMD_perc_ident_and_gaps(df_nonTMD, mean_ser):
    """Calculate the nonTMD percentage identity and gaps.

    Parameters
    ----------
    df_nonTMD : pd.DataFrame
        Dataframe containing the nonTMD sequences, % identity, number of gaps, etc.
        columns : nonTMD_perc_ident etc
        index : homologue number
    mean_ser : dict
        Dictionary containing the mean values for that protein, from all homologues.
        Will be saved as a separate file for each protein, and gathered later to form a full dataset.

    Returns
    -------
    df_nonTMD, mean_ser (updated dataframe and dictionary)
    """

    ########################################################################################
    #                                                                                      #
    #         calculate nonTMD len, percent identity, gaps, etc [AAIMON]                   #
    #                                                                                      #
    ########################################################################################
    # calculate identical residues in the nonTMD region (simply count the pipes '|' in the markup sequence)
    # NOTE THAT df_nonTMD['nonTMD_markup'].str.count('|') DOES NOT WORK, as "|" has a function in regex and needs to be escaped
    df_nonTMD['nonTMD_num_ident_res'] = df_nonTMD['nonTMD_markup'].str.count('\|')
    # calculate similar residues in the nonTMD region (simply count the colons ':' in the markup sequence)
    df_nonTMD['nonTMD_num_sim_res'] = df_nonTMD['nonTMD_markup'].str.count(':')
    # add the identical and similar residues together to get the total number of similar + identical residues
    df_nonTMD['nonTMD_num_sim_plus_ident_res'] = df_nonTMD['nonTMD_num_ident_res'] + df_nonTMD['nonTMD_num_sim_res']

    # count the gaps in the nonTMD sequence of the query
    df_nonTMD['nonTMD_q_num_gaps'] = df_nonTMD['nonTMD_seq_query'].str.count('-')
    # count the gaps in the nonTMD sequence of the match
    df_nonTMD['nonTMD_m_num_gaps'] = df_nonTMD['nonTMD_seq_match'].str.count('-')
    # calculate the length of the nonTMD sequences, which may include gaps
    df_nonTMD['len_nonTMD_seq_query'] = df_nonTMD['nonTMD_seq_query'].str.len()
    df_nonTMD['len_nonTMD_seq_markup'] = df_nonTMD['nonTMD_markup'].str.len()
    df_nonTMD['len_nonTMD_seq_match'] = df_nonTMD['nonTMD_seq_match'].str.len()

    # calculate the number aligned sequences, excluding gaps (length of query, or length of match, whichever is shorter)
    # NOTE, THIS IS NEEDS TO BE CHECKED.
    #df_nonTMD['len_nonTMD_align'] = df_nonTMD[['len_nonTMD_seq_query', 'len_nonTMD_seq_match']].dropna(how='all').min(axis=1)
    df_nonTMD['len_nonTMD_align'] = df_nonTMD['len_nonTMD_seq_query']


    # calculate the length of the nonTMD sequence excluding gaps
    df_nonTMD['len_nonTMD_q_excl_gaps'] = df_nonTMD['len_nonTMD_seq_query'] - df_nonTMD['nonTMD_q_num_gaps']
    df_nonTMD['len_nonTMD_m_excl_gaps'] = df_nonTMD['len_nonTMD_seq_match'] - df_nonTMD['nonTMD_m_num_gaps']
    # calculate the length of the alignment by finding which seq excl gaps is smaller
    # NOTE, THIS IS CURRENTLY TOO SHORT, GIVING NONTMD IDENTITIES ALWAYS ABOVE 1.0. NEEDS TO BE FIXED.
    #df_nonTMD['len_nonTMD_align'] = df_nonTMD[['len_nonTMD_q_excl_gaps', 'len_nonTMD_m_excl_gaps']].min(axis=1)

    # calculate the percentage identity of the nonTMD region (number of identical residues divided by the length excluding gaps)
    # used for the Amino Acid Identity : Membranous over Nonmembranous (AAIMON ratio)
    # note that the length = length of the aligned residues excluding gaps
    df_nonTMD['nonTMD_perc_ident'] = df_nonTMD['nonTMD_num_ident_res'] / df_nonTMD['len_nonTMD_align']
    df_nonTMD['nonTMD_perc_sim'] = df_nonTMD['nonTMD_num_sim_res'] / df_nonTMD['len_nonTMD_align']
    df_nonTMD['nonTMD_perc_sim_plus_ident'] = df_nonTMD['nonTMD_num_sim_plus_ident_res'] / df_nonTMD['len_nonTMD_align']
    # calculate the average number of gaps per residue in the nonTMD alignment
    # filter to analyse only sequences that are valid (length > 0)
    df_nonTMD_filt_gaps = df_nonTMD.loc[df_nonTMD['len_nonTMD_q_excl_gaps'] != 0]
    # calculate number of gaps in query AND match
    df_nonTMD_filt_gaps['nonTMD_qm_num_gaps'] = df_nonTMD_filt_gaps['nonTMD_q_num_gaps'] + df_nonTMD_filt_gaps['nonTMD_m_num_gaps']
    # add to simap dataframe
    df_nonTMD['nonTMD_qm_num_gaps'] = df_nonTMD_filt_gaps['nonTMD_qm_num_gaps']
    # gaps per query residue for both query and match = ((gaps in query + gaps in match)/2))/length of query excluding gaps
    df_nonTMD['nonTMD_qm_gaps_per_q_residue'] = df_nonTMD_filt_gaps['nonTMD_qm_num_gaps'] / 2 / df_nonTMD_filt_gaps['len_nonTMD_q_excl_gaps']

    # SSR ratio calculations take a long time and show no difference to AAIMON
    SSR_ratio_calculations = False
    if SSR_ratio_calculations:
        korbinian.cons_ratio.ssr.conduct_ssr_ratio_calculations()

    ########################################################################################
    #                                                                                      #
    #   calculate average nonTMD len, perc ident, etc and add to list of proteins [AAIMON] #
    #                                                                                      #
    ########################################################################################

    ''' calculate values associated with the nonTMD region and add to main dataframe
    '''
    # add the average values regarding the nonTMD region to the original file/dataframe with each protein
    mean_ser['len_nonTMD_seq_match_mean'] = float('%0.2f' % df_nonTMD['len_nonTMD_seq_match'].dropna().mean())
    mean_ser['nonTMD_perc_ident_mean'] = float('%0.3f' % df_nonTMD['nonTMD_perc_ident'].dropna().mean())
    mean_ser['nonTMD_perc_sim_mean'] = float('%0.3f' % df_nonTMD['nonTMD_perc_sim'].dropna().mean())
    mean_ser['nonTMD_perc_sim_plus_ident_mean'] = float('%0.3f' % df_nonTMD['nonTMD_perc_sim_plus_ident'].dropna().mean())
    mean_ser['len_nonTMD_align_mean'] = float('%0.2f' % df_nonTMD['len_nonTMD_align'].dropna().mean())
    mean_ser['nonTMD_qm_gaps_per_q_residue_mean'] = float('%0.2f' % df_nonTMD['nonTMD_qm_gaps_per_q_residue'].dropna().mean())
    #logging.info('nonTMD_qm_gaps_per_q_residue : %0.5f' % mean_ser['nonTMD_qm_gaps_per_q_residue_mean'])
    return mean_ser, df_nonTMD