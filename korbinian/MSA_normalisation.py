import csv
import pandas as pd
import numpy as np
import string
import numpy
import random

############################################################################################
#                                                                                          #
#                        calculation of aa propensity                                      #
#                                  based on dataset                                        #
#                                                                                          #
############################################################################################

def cal_aa_propensity_TM_nonTM (df, TM_col='TM01_seq', nonTM_col='nonTMD_seq'):
    """Calculation of amino acid propensity for TM and non-TM region in dataset.

    Parameters
    ----------
    df : pd.DataFrame
        dataframe which contains the TM and non-TM sequences for each protein

    Returns
    -------
    prob_table : pd.DataFrame
        show the aa propensity in TM and non-TM region, respectively
        index is the AA
        columns are the input columns plus aap (e.g. "TM01_seq" + "_aap")
    """

    # create a string to hold all TM segments from all proteins
    massive_string_tm = ""
    for seq in df[TM_col]:
        if type(seq) == str:
            massive_string_tm += seq

    # create a string to hold all non-TM segments from all proteins
    massive_string_ntm = ""
    for seq in df[nonTM_col]:
        if type(seq) == str:
            massive_string_ntm += seq

    # calculate aa propensity in TM region
    prob_table_tm = calc_aa_propensity(massive_string_tm)
    # calculate aa propensity in non-TM region
    prob_table_ntm = calc_aa_propensity(massive_string_ntm)
    # merge the two table into one dataframe
    prob_table = pd.concat([prob_table_tm, prob_table_ntm], axis=1)
    # rename the columns to match the content, with the orig name plus "amino acid propensity"
    prob_table.columns = [TM_col + "_aap", nonTM_col + "_aap"]

    return prob_table


def calc_aa_propensity(seq):
    """calculate aa propensity for each residue in a particular sequence.

    Parameters
    ----------
    seq : string
        TM or non-TM sequence

    Returns
    -------
    freq_ser : pd.Series
        Series containing corresponding aa propensity
    """

    # count absolute number of each residue in the input string
    aa_dict = {}
    for aa in seq:
        aa_dict[aa] = 0
    for aa in seq:
        aa_dict[aa] += 1

    # create a dictionary to hold the propensity of each residue
    prop_dict = {}
    for aa in aa_dict:
        prop_dict[aa] = aa_dict[aa] / len(seq)

    # turn the dictionary into a pd.DataFrame and transpose it so that the residue names are listed in index
    freq_ser = pd.Series(prop_dict)
    freq_ser.index.name = "freq"
    return freq_ser


def cal_random_aa_ident(prob_table, ident=0.9, seq_len=10000, number_seq=500, subset_num=12):
    """Calculation of random amino acid identity based on a particular amino acid propensity.

    Protein regions with a limited aa propensity (e.g. transmembrane regions) have a measurable amino
    acid identity EVEN IN NON-HOMOLOGUES. This is referred to here as the random amino acid identity.
    This formula takes the aa propensity of a sequence or dataset as an input, and calculates the random aa identity.

    Parameters
    ----------
    prob_table : pd.DataFrame
        aa propensity for TM and non-TM region. obtained from the function cal_aa_propensity

    ident: float
        overall identity of randomly created sequence cluster. Default value: 0.9

    seq_len: integer
        length of randomly created sequences. To achieve a more plausible result using randomisation method,
        greater values (> 5000) are recommended. Defalut value: 10,000

    number_seq: integer
        number of aligned sequences. Larger values are recommended. Default value: 500

    subset_num: integer
        currently not in use.

    Returns
    -------
    TM_back_mutation_rate,  nonTM_back_mutation_rate : tuple
        = random identity TM, random identity non-TM
    """

    # extract aa pool and propensity for TM region in form of np.array
    TM_probtable = prob_table["TM"].dropna()
    TM_aa_probabilities = np.array(TM_probtable)
    aa_pool = np.array(TM_probtable.index)

    # calculate number of residues that need to be replaced based on the desired percentage identity. This applies to both
    # TM and non-TM region
    number_mutations = int(seq_len*(1 - ident))

    # generate random sequences, extract the original reference sequence and the sequence cluster
    TM_matrix = generate_ran_seq(seq_len, number_seq, number_mutations, subset_num, aa_pool, TM_aa_probabilities)
    TM_orig_seq = TM_matrix[0]
    TM_mat = TM_matrix[1]

    # make a list of residues in each position (each column in MSA)
    pos = []
    tr_mat = []
    for i in range(seq_len):
        for seq in TM_mat:
            pos.append(seq[i])
        pos = ''.join(pos)
        tr_mat.append(pos)
        pos = []

    # count amino acid frequency for each position in a MSA
    cons_list = []
    for pos in tr_mat:
        cons = calc_aa_propensity(pos)
        cons_list.append(cons)

    # turn the list into a pd.DataFrame
    cons_df = pd.concat(cons_list)
    # len_df = len(cons_df.columns)
    # use the original reference sequence as index to find the conserved residue at each position (column)
    cons_df.index = list(TM_orig_seq)

    # find the conserved residue (which mathches the coresponding one in reference seqeunce) and its propensity at each position
    aa_prop_list = []
    # to avoid redundant counting, use .unique
    for aa in cons_df.index.unique():
        if type(cons_df[aa][aa]) == pd.Series:
            observed_cons_1_row = list(cons_df.loc[aa, aa])
            aa_prop_list += observed_cons_1_row
        else:
            aa_prop_list.append(cons_df.loc[aa, aa])

    # calculation of the average idenity (conservation) at each position
    observed_mean_cons_rate_all_pos = np.array(aa_prop_list).mean()
    # calculate the random identity in TM region in form of back mutation" rate", which represents the fraction of mutation which have
    # resulted in the same aa residue as in the original reference sequence
    TM_back_mutation_rate = (observed_mean_cons_rate_all_pos - ident ) / (1 - ident)


    ### the same procedure for non-Tm region ###

    # extract aa pool and propensity for non-TM region in form of np.array
    nonTM_probtable = prob_table["non-TM"].dropna()
    nonTM_aa_probabilities = np.array(nonTM_probtable)
    aa_pool = np.array(nonTM_probtable.index)

    # generate random sequences, extract the original reference sequence and the sequence cluster
    nonTM_matrix = generate_ran_seq(seq_len, number_seq, number_mutations, subset_num, aa_pool, nonTM_aa_probabilities)
    nonTM_orig_seq = nonTM_matrix[0]
    nonTM_mat = nonTM_matrix[1]

    # make a list of residues in each position (each column in MSA)
    pos = []
    tr_mat = []
    for i in range(seq_len):
        for seq in nonTM_mat:
            pos.append(seq[i])
        pos = ''.join(pos)
        tr_mat.append(pos)
        pos = []

    # count amino acid frequency for each position in a MSA
    cons_list = []
    for pos in tr_mat:
        cons = calc_aa_propensity(pos)
        cons_list.append(cons)

    # turn the list into a pd.DataFrame
    cons_df = pd.concat(cons_list)
    len_df = len(cons_df.columns)
    #cons_df2 = cons_df.copy()
    # use the original reference sequence as index to find the conserved residue at each position (column)
    cons_df.index = list(nonTM_orig_seq)

    # find the conserved residue (which mathches the coresponding one in reference seqeunce) and its propensity at each position
    aa_prop_list = []
    for aa in cons_df.index.unique():
        if type(cons_df[aa][aa]) == pd.Series:
            observed_cons_1_row = list(cons_df.loc[aa, aa])
            aa_prop_list += observed_cons_1_row
        else:
            aa_prop_list.append(cons_df.loc[aa, aa])

    # calculation of the average idenity (conservation) at each position
    observed_cons_rate_all_pos = np.array(aa_prop_list).mean()
    # calculate the random identity in non-TM region in form of back mutation" rate", which represents the fraction of mutation which have
    # resulted in the same aa residue as in the original reference sequence
    nonTM_back_mutation_rate = (observed_cons_rate_all_pos - ident) / (1 - ident )

    return TM_back_mutation_rate,  nonTM_back_mutation_rate


def generate_ran_seq(seq_len, number_seq, number_mutations, subset_num, aa_pool, aa_probabilities):
    """Generation of sequence cluster using randomisation method

    Parameters
    ----------
    seq_len: integer
        length of randomly created sequences. To achieve a more plausible result using randomisation method,
        greater values (> 5000) are recommended. Defalut value: 10,000

    number_seq: integer
        number of aligned sequences. Larger values are recommended. Default value: 500

    number_mutations: integer
        number of residues that will be randomly replaced in a sequence. precalculated in the function cal_random_aa_ident

    subset_num: integer
        currently not in use.

    aa_pool: np.array
        array of amino acid from which residues will be chosen randomly to create a sequence or to replace a residue

    aa_probabilities: np.array
        array of propensities. The order should match that of the amino acid in aa_ppol

    Returns
    -------
    ori_seq: string
        original sequence as reference

    seq_matrix: list
        seqeunce cluster that are created by randomly replacing predetermined number of residues in the reference sequence
    """

    # seq_list = []
    # sublist = ''.join(np.random.choice(aa_pool, p=aa_probabilities) for _ in range(subset_num))
    # subdict = { my_key: prob_table[my_key] for my_key in sublist }
    # pick_list = []
    # for key, prob in subdict.items():
    #    pick_list.extend([key] * int((prob * 100)))

    # generate a reference sequence based on the aa propensity of TM or non-TM region
    ori_seq = "".join(np.random.choice(aa_pool, p=aa_probabilities) for _ in range(int(seq_len)))

    # generate sequence cluster by randomly replacing predetermined number of residues in reference seq
    seq_matrix = []
    # firstly, choose a set of positions whoose aa will be replaced
    for n in range(number_seq):
        inds = [i for i, _ in enumerate(ori_seq) if not ori_seq.isspace()]
        sam = random.sample(inds, number_mutations)
        lst = list(ori_seq)
        # based on aa propensity, replace the residue at each chosen position
        for ind in sam:
            lst[ind] = np.random.choice(aa_pool, p=aa_probabilities)
        new_seq = "".join(lst)

        # append each new sequence to the seq_matrix
        seq_matrix.append(new_seq)

    return ori_seq, seq_matrix



############################################################################################
#                                                                                          #
#                        calculation of normalisation factor                               #
#                                  based on aa identity                                    #
#                                                                                          #
############################################################################################

def cal_MSA_ident_n_factor(ident_tm, rand_tm, rand_ntm):
    """Calculation of the MSA identity normalisation factor

    Parameters
    ----------
    ident_tm: float
        the observed average identity of TM region in your MSA which needs to be normalised

    rand_tm: float
        random identity in TM region, calculated based on your dataset using radomisation method (cal_random_aa_ident)

    rand_ntm: float
        random identity in non-TM region, calculated based on your dataset using radomisation method (cal_random_aa_ident)

    Returns
    -------
    n_factor: float
        normalisation factor which will be applied to your observed TM identity

    TM_ident_n: float
        normalised TM identity for MSA

    Example:
    observed average ident_tm = 0,78, rand_tm = 0.126, rand_ntm = 0.059
    calculated real_cons = 0.748
    calculated ident_ntm = 0.763
    calculated n_factor = 0.78/0.763 = 1.022
    """

    # calculation of real conservation rate based on the random identity in TM region
    real_cons = (ident_tm - rand_tm)/(1 - rand_tm)

    # calculation of the observed conservation in non-TM region
    ident_ntm = (1 - real_cons)*rand_ntm + real_cons

    #calculation of normalisation factor
    n_factor = ident_tm/ident_ntm

    # normalise the TM identity
    TM_ident_n = ident_tm/n_factor
    return 'normalisation factor: %.3f' %n_factor, 'normalised overall TM identity: %.3f' %TM_ident_n




