import csv
import pandas as pd
import numpy as np
import string
import numpy
import random
import sys

############################################################################################
#                                                                                          #
#                        calculation of aa propensity                                      #
#                                  based on dataset                                        #
#                                                                                          #
############################################################################################

def calc_aa_propensity_TM_nonTM(df, TM_col='TM01_seq', nonTM_col='nonTMD_seq'):
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
    massive_string_TM = ""
    for seq in df[TM_col]:
        if type(seq) == str:
            massive_string_TM += seq

    # create a string to hold all non-TM segments from all proteins
    massive_string_nonTM = ""
    for seq in df[nonTM_col]:
        if type(seq) == str:
            massive_string_nonTM += seq

    # calculate aa propensity in TM region
    TM_aa_propensity_ser = calc_aa_propensity(massive_string_TM)
    # calculate aa propensity in non-TM region
    nonTM_aa_propensity_ser = calc_aa_propensity(massive_string_nonTM)
    # merge the two table into one dataframe
    aa_propensity_TM_nonTM_df = pd.concat([TM_aa_propensity_ser, nonTM_aa_propensity_ser], axis=1)
    # rename the columns to match the content, with the orig name plus "amino acid propensity"
    aa_propensity_TM_nonTM_df.columns = [TM_col + "_aap", nonTM_col + "_aap"]

    return aa_propensity_TM_nonTM_df


def calc_aa_propensity(seq):
    """calculate aa propensity for each residue in a particular sequence.

    Parameters
    ----------
    seq : string
        TM or non-TM sequence

    Returns
    -------
    aa_prop_ser : pd.Series
        Series containing corresponding aa propensity
    """

    # count absolute number of each residue in the input string
    number_each_aa_dict = {}

    all_aa = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    # create an dictionary of the numbers {"A" : 57, "C" : 5, ...} etc
    for aa in all_aa:
        number_each_aa_dict[aa] = seq.count(aa)

    # create a dictionary to hold the propensity of each residue
    aa_propensity_dict = {}
    length = len(seq)
    for aa in number_each_aa_dict:
        aa_propensity_dict[aa] = number_each_aa_dict[aa] / length

    # turn the dictionary into a pd.Series
    aa_prop_ser = pd.Series(aa_propensity_dict)
    # normalise so that all the aa propensities add up to 1.0
    # this is important if "X" or "U" is in the sequences
    aa_prop_norm_ser = aa_prop_ser / aa_prop_ser.sum()
    aa_prop_norm_ser.index.name = "freq"
    return aa_prop_norm_ser


def calc_random_aa_ident(aa_prop_ser, seq_len=1000, number_seq=1000, ident=0.7):
    """Calculation of random amino acid identity based on a particular amino acid propensity.

    Protein regions with a limited aa propensity (e.g. transmembrane regions) have a measurable amino
    acid identity EVEN IN NON-HOMOLOGUES. This is referred to here as the random amino acid identity.
    This formula takes the aa propensity of a sequence or dataset as an input, and calculates the random aa identity.

    Parameters
    ----------
    aa_prop_ser : pd.Series
        aa propensity for a particular sequenc or dataset.
        pandas series with all 20 aa as the index, and a proportion (0.08, 0.09 etc) as the data.
        Typically obtained from the function calc_aa_propensity

    ident: float
        desired overall identity of randomly created sequence matrix.

    seq_len: integer
        length of randomly created sequences. To achieve a more plausible result using randomisation method,
        greater values (> 5000) are recommended. Defalut value: 10,000

    number_seq: integer
        number of aligned sequences. Larger values are recommended. Default value: 500

    subset_num: integer
        currently not in use.

    Returns
    -------
    random_aa_identity : float
        random identity due to limited aa propensity
        effectively the average back mutation rate for all amino acids
    """
    # extract aa array and propensity array
    aa_propensities = np.array(aa_prop_ser)
    aa_arr = np.array(aa_prop_ser.index)

    # calculate number of residues that need to be replaced based on the desired percentage identity.
    number_mutations = int(np.round(seq_len*(1 - ident)))

    # generate random sequences, extract the original reference sequence and the sequence cluster
    orig_and_mut_seqs = generate_random_seq(seq_len, number_seq, number_mutations, aa_arr, aa_propensities)
    # extract the original sequence, of which the matrix are variants
    orig_seq = orig_and_mut_seqs[0]
    aa_prop_orig_seq = calc_aa_propensity(orig_seq)
    aa_in_orig_seq_list = list(aa_prop_orig_seq.loc[aa_prop_orig_seq > 0].index)

    # extract the matrix of mutated sequences, slightly different from the orig_seq
    mut_seqs_matrix = orig_and_mut_seqs[1]
    # make a list of residues in each position (each column in MSA)
    list_of_columnwise_strings = []
    for i in range(mut_seqs_matrix.shape[1]):
        """ joins up everything in column
        orig seq : G   I   L   I
        mut1       G   I   L   I
        mut2       G   V   L   I
        mut3       G   I   L   P

        G : GGG
        I : IVI
        L : LLL
        etc.
        """
        # takes one column, joins all aa into a single string
        string_joined_aa_at_that_pos = "".join(mut_seqs_matrix[:, i])
        # adds that string to a list
        list_of_columnwise_strings.append(string_joined_aa_at_that_pos)

    # count amino acid frequency for each position in a MSA
    columnwise_aa_propensities_df = pd.DataFrame()

    # iterate through length of orig_seq
    for n in range(seq_len):
        # in the matrix, the amino acids at that positions can be extracted from the nested list created previously
        string_aa_in_matrix_at_that_pos = list_of_columnwise_strings[n]
        # create a series of amino acid propensities from that nested list (of course, mostly the aa is the original one)
        aa_prop_ser = calc_aa_propensity(string_aa_in_matrix_at_that_pos)
        # add the amino acid propensities as a new column in the dataframe, with the orig_aa number as the column name
        columnwise_aa_propensities_df[n] = aa_prop_ser
    # replace the orig_aa numbers as column names with the orig_aa itself
    columnwise_aa_propensities_df.columns = list(orig_seq)
    """
    columnwise_aa_propensities_df
        index = all 20 amino acids
        columns = all orig aa in sequence
        content = aa propensities for that aa for that position

        orig seq : G        A           L           I
        A          0.0      0.91       0.05         0.02
        C          0.0      0.01       0.04         0.02
        D          0.0      0.0        0.05         0.01
        etc..

    """
    all_aa = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

    # find the conserved residue (which matches the corresponding one in reference sequence) and its propensity at each position
    all_perc_orig_aa_in_matrix_list = []
    # to avoid redundant counting, use .unique
    for aa in all_aa:
        if aa in aa_in_orig_seq_list:
            # if the orig seq contains more than one of the aa of interest
            if type(columnwise_aa_propensities_df.loc[aa,aa]) == pd.Series:
                """
                e.g.

                orig seq : G   I   L   I
                mut1       G   I   L   I
                mut2       G   V   L   I
                mut3       G   I   L   P
                           1  0.75  1  0.75

                for I (aa occurs twice in orig, perc_orig_aa_all_rows_for_that_aa) : [0.75, 0.75]

                for G (aa occurs once, perc_orig_aa_1_row) : 1.0

                """
                perc_orig_aa_all_rows_for_that_aa = list(columnwise_aa_propensities_df.loc[aa, aa])
                # add all items in list to large list of perc_orig_aa
                all_perc_orig_aa_in_matrix_list += perc_orig_aa_all_rows_for_that_aa
            else:
                perc_orig_aa_1_row = columnwise_aa_propensities_df.loc[aa, aa]
                # add that percentage to large list of perc_orig_aa
                all_perc_orig_aa_in_matrix_list.append(perc_orig_aa_1_row)

    # calculation of the average identity (conservation) at each position
    observed_mean_cons_rate_all_pos = np.array(all_perc_orig_aa_in_matrix_list).mean()
    # calculate the random identity in TM region in form of back mutation" rate", which represents the fraction of mutation which have
    # resulted in the same aa residue as in the original reference sequence
    random_aa_identity = (observed_mean_cons_rate_all_pos - ident ) / (1 - ident)
    return random_aa_identity


def generate_random_seq(seq_len, number_seq, number_mutations, list_all_20_aa, probabilities_all_20_aa):
    """Generation of sequence cluster using randomisation method

    Parameters
    ----------
    seq_len: int
        length of randomly created sequences. To achieve a more plausible result using randomisation method,
        greater values (> 5000) are recommended. Defalut value: 10,000

    number_seq: int
        number of aligned sequences. Larger values are recommended. Default value: 500

    number_mutations: int
        number of residues that will be randomly replaced in a sequence. precalculated in the function calc_random_aa_ident

    subset_num: int
        currently not in use.

    list_all_20_aa: np.array
        array of amino acid from which residues will be chosen randomly to create a sequence or to replace a residue

    probabilities_all_20_aa: np.array
        array of propensities. The order should match that of the amino acid in aa_ppol

    Returns
    -------
    ori_seq: str
        original sequence as reference

    seq_matrix: np.array
        sequence cluster that are created by randomly replacing predetermined number of residues in the reference sequence
    """

    # seq_list = []
    # sublist = ''.join(np.random.choice(list_all_20_aa, p=probabilities_all_20_aa) for _ in range(subset_num))
    # subdict = { my_key: prob_table[my_key] for my_key in sublist }
    # pick_list = []
    # for key, prob in subdict.items():
    #    pick_list.extend([key] * int((prob * 100)))

    # generate a reference sequence based on the aa propensity of TM or non-TM region

    orig_seq = "".join(np.random.choice(list_all_20_aa, p=probabilities_all_20_aa) for _ in range(int(seq_len)))

    # generate sequence cluster by randomly replacing predetermined number of residues in reference seq
    seq_matrix = []
    # firstly, choose a set of positions whoose aa will be replaced
    for n in range(number_seq):
        if n % 10 == 0:
            sys.stdout.write(".")
            if n !=0 and n % 300 == 0:
                sys.stdout.write(" please have patience, I'm still calculating \n")
        sys.stdout.flush()

        # create indices (list of positions)
        inds = list(range(seq_len))
        # number of mutations is calculated beforehand. E.g. if ident=0.9, seqlen=100, number_mutations = 10)
        # create a sample of positions to mutate, e.g. [77, 81, 18, 46, 42, 53, 65, 2, 89, 69, ..... and so on
        list_of_aa_positions_to_be_mutated = random.sample(inds, number_mutations)
        orig_seq_as_list = list(orig_seq)
        # based on aa propensity, replace the residue at each chosen position
        for pos in list_of_aa_positions_to_be_mutated:
            orig_seq_as_list[pos] = np.random.choice(list_all_20_aa, p=probabilities_all_20_aa)
        seq_incl_mutations = "".join(orig_seq_as_list)

        # append each new sequence to the seq_matrix
        seq_matrix.append(list(seq_incl_mutations))

    seq_matrix = np.array(seq_matrix)

    return orig_seq, seq_matrix


############################################################################################
#                                                                                          #
#                        calculation of normalisation factor                               #
#                                  based on aa identity                                    #
#                                                                                          #
############################################################################################

def calc_MSA_ident_n_factor(observed_perc_ident_full_seq, rand_perc_ident_TM, rand_perc_ident_nonTM, proportion_seq_TM_residues=0.3):
    """Calculation of the MSA identity normalisation factor



    For this formula, we assume most proteins are multi-pass, and that approximately 30% of the
    residues are TM residues. Therefore a rand_30TM_70nonTM can be calculated, that roughly
    gives the random identity for the full protein.

    rand_30TM_70nonTM = 0.3 * rand_perc_ident_TM + 0.7 * rand_perc_ident_nonTM

    Parameters
    ----------
    observed_perc_ident_full_seq: float
        the observed average identity of TM region in your MSA which needs to be normalised

    rand_perc_ident_TM: float
        random identity in TM region, calculated based on your dataset using radomisation method (calc_random_aa_ident)

    rand_perc_ident_nonTM: float
        random identity in non-TM region, calculated based on your dataset using radomisation method (calc_random_aa_ident)

    proportion_seq_TM_residues : float
        proportion of the sequence length that is the TM region
        To roughly calculate the observed percentage identity of the TM region from the full percentage
        identity, it is necessary to estimate the percentage length of the TM region.
        For the single-pass human dataset this is 0.0681 (6.8% TM region)
        For the multi-pass human dataset this is 0.330 (34% TM region)
        For the non-redundant beta-barrel dataset this is 0.348 (35% TM region)

    Returns
    -------
    n_factor: float
        normalisation factor which will be applied to your observed TM identity

    TM_ident_n: float
        normalised TM identity for MSA

    Example:
    observed_perc_ident_TM = 0,78, rand_perc_ident_TM = 0.126, rand_perc_ident_nonTM = 0.059
    calculated real_perc_identity = 0.748
    calculated observed_perc_ident_nonTM = 0.763
    calculated n_factor = 0.78/0.763 = 1.022
    """
    # calculate proportion of length of full sequence that is nonTM
    proportion_seq_nonTM_residues = 1 - proportion_seq_TM_residues
    # random percentage identity of the full protein, assuming 30% TM region and 70% nonTM region
    rand_perc_ident_full_protein = proportion_seq_TM_residues * rand_perc_ident_TM + proportion_seq_nonTM_residues * rand_perc_ident_nonTM

    # calculation of real conservation rate based on the random identity in TM region
    # solved for R from observed_perc_ident_full_seq = real_perc_identity + (1-real_perc_identity)*rand_perc_ident_full_protein
    # as usual, we assume that the unobserved conservation is a proportion of the observed_changes (i.e. (1-real_perc_identity))
    # and that this proportion is exactly the rand_perc_ident_full_protein * real_changes
    real_perc_identity = (observed_perc_ident_full_seq - rand_perc_ident_full_protein)/(1 - rand_perc_ident_full_protein)

    # from the estimated real_perc_identity of the full protein, calculate the observed percentage identity for the TM region
    observed_perc_ident_TM = (1 - real_perc_identity)*rand_perc_ident_TM + real_perc_identity
    # from the estimated real_perc_identity of the full protein, calculate the observed percentage identity for the nonTM region
    observed_perc_ident_nonTM = (1 - real_perc_identity)*rand_perc_ident_nonTM + real_perc_identity

    #calculation of normalisation factor
    # for randomised sequences, the aa propensity is the ONLY factor giving an effect
    # therefore the ratio of the observed identities gives the normalisation factor
    MSA_TM_nonTM_aa_ident_norm_factor = observed_perc_ident_TM/observed_perc_ident_nonTM

    #sys.stdout.write('\nnormalisation factor: %.3f' %MSA_TM_nonTM_aa_ident_norm_factor)

    return MSA_TM_nonTM_aa_ident_norm_factor

def OLD_calc_MSA_ident_n_factor(observed_perc_ident_TM, rand_perc_ident_TM, rand_perc_ident_nonTM):
    """Calculation of the MSA identity normalisation factor

    To roughly calculate the observed percentage identity of the TM region from the full percentage
    identity, it is necessary to estimate the percentage length of the TM region.
    For the single-pass human dataset this is 0.0681 (6.8% TM region)
    For the multi-pass human dataset this is 0.330 (34% TM region)
    For the non-redundant beta-barrel dataset this is 0.348 (35% TM region)

    For this formula, we assume most proteins are multi-pass, and that approximately 30% of the
    residues are TM residues. Therefore a rand_30TM_70nonTM can be calculated, that roughly
    gives the random identity for the full protein.

    rand_30TM_70nonTM = 0.3 * rand_perc_ident_TM + 0.7 * rand_perc_ident_nonTM

    Parameters
    ----------
    observed_perc_ident_full_seq: float
        the observed average identity of TM region in your MSA which needs to be normalised

    rand_perc_ident_TM: float
        random identity in TM region, calculated based on your dataset using radomisation method (calc_random_aa_ident)

    rand_perc_ident_nonTM: float
        random identity in non-TM region, calculated based on your dataset using radomisation method (calc_random_aa_ident)

    Returns
    -------
    n_factor: float
        normalisation factor which will be applied to your observed TM identity

    TM_ident_n: float
        normalised TM identity for MSA

    Example:
    observed_perc_ident_TM = 0,78, rand_perc_ident_TM = 0.126, rand_perc_ident_nonTM = 0.059
    calculated real_perc_identity = 0.748
    calculated observed_perc_ident_nonTM = 0.763
    calculated n_factor = 0.78/0.763 = 1.022
    """
    # calculation of real conservation rate based on the random identity in TM region
    # solved for R from observed_perc_ident_full_seq = real_perc_identity + (1-real_perc_identity)*rand_perc_ident_full_protein
    # as usual, we assume that the unobserved conservation is a proportion of the observed_changes (i.e. (1-real_perc_identity))
    # and that this proportion is exactly the rand_perc_ident_full_protein * real_changes
    real_perc_identity_TM = (observed_perc_ident_TM - rand_perc_ident_TM)/(1 - rand_perc_ident_TM)

    # # from the estimated real_perc_identity of the full protein, calculate the observed percentage identity for the TM region
    # observed_perc_ident_TM = (1 - real_perc_identity_TM)*rand_perc_ident_TM + real_perc_identity_TM
    # from the estimated real_perc_identity of the full protein, calculate the observed percentage identity for the nonTM region
    observed_perc_ident_nonTM = (1 - real_perc_identity_TM)*rand_perc_ident_nonTM + real_perc_identity_TM

    #calculation of normalisation factor
    # for randomised sequences, the aa propensity is the ONLY factor giving an effect
    # therefore the ratio of the observed identities gives the normalisation factor
    MSA_TM_nonTM_aa_ident_norm_factor = observed_perc_ident_TM/observed_perc_ident_nonTM

    #sys.stdout.write('\nnormalisation factor: %.3f' %MSA_TM_nonTM_aa_ident_norm_factor)

    return MSA_TM_nonTM_aa_ident_norm_factor
