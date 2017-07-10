import argparse
import korbinian
import pandas as pd
import numpy as np
import random
import sys
# import debugging tools
from korbinian.utils import pr, pc, pn, aaa

def calc_aa_propensity_from_csv_col(seq_list_csv_in, aa_prop_csv_out, col_name, sep=","):
    """Calculation of amino acid propensity for TM and non-TM region in dataset.

    Parameters
    ----------
    seq_list_csv_in : csv
        input csv file which contains sequences from region of interest (e.g. TM and nonTM regions), normally as comma separated values

    aa_prop_csv_out: csv
        output csv file which contains the aa propensity for the region of interest

    col_name: str
        specify which column should be used (e.g. TM01_seq or nonTM_seq). This should contain sequences
        from the interested protein region.

    sep: str
        data format. Default: comma separated file
    """

    # open csv
    df = pd.read_csv(seq_list_csv_in, sep=sep)
    # extract column of interest, and drop empty rows
    ser = df[col_name].dropna()

    # create a string to hold segments from all proteins in list
    massive_string_all_prot = ""
    for seq in ser:
        if type(seq) == str:
            massive_string_all_prot += seq

    # calculate aa propensity in region of interest
    aa_propensity_ser = calc_aa_propensity(massive_string_all_prot)
    # save aa propensity series to output csv file
    aa_propensity_ser.to_csv(aa_prop_csv_out, sep="\t")

def calc_aa_propensity_TM_nonTM(df, TM_col='TM01_seq', nonTM_col='nonTMD_seq'):
    """Calculation of amino acid propensity for TM and non-TM region in dataset.

    Parameters
    ----------
    df : pd.DataFrame
        dataframe which contains the TM and non-TM sequences for each protein

    TM_col: str
        column that contains TM sequences

    nonTM_col: str
        column that contains non-TM sequences

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
    aa_prop_norm_ser : pd.Series
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
    # name the index column
    aa_prop_norm_ser.index.name = "freq"
    return aa_prop_norm_ser

def calc_random_aa_ident_from_ser(aa_prop_ser):
    """ Calculates the random AA identity from an input series of AA propensity values.

    Parameters
    ----------
    aa_prop_ser : pd.Series
        Series of AA propensity values
        index : A, C, D, E etc
        values : 0.04, 0.06, .20 etc

    Returns
    -------
    calc_rand_ident : float
        Calculated random identity value, e.g. 0.092.
    """
    all_aa = "ACDEFGHIKLMNPQRSTVWY"
    # check that the index has 20 aa
    if not sorted(list(aa_prop_ser.index)) == sorted(list(all_aa)):
        raise ValueError("Index must be comprised of 20 amino acids, in capitals, single-letter notation.")
    # multiply the probabilities by themselves
    # if 25% of the residues are leucines, and the probability of them staying leucines is 25%
    # that means that 0.25 * 0.25 is the random identity of leucine
    # calculate this for all AA and simply sum!
    # this gives the random identity overall!
    aa_prop_ser = aa_prop_ser**2
    calc_rand_ident = aa_prop_ser.sum()
    return calc_rand_ident

def calc_random_aa_ident(aa_prop_csv_in, rand_seq_ident_csv_out):
    """Mathematical calculation of random amino acid identity based on a particular amino acid propensity.

    Protein regions with a limited aa propensity (e.g. transmembrane regions) have a measurable amino
    acid identity EVEN IN NON-HOMOLOGUES. This is referred to here as the random amino acid identity.
    This formula takes the aa propensity of a sequence or dataset as an input, and calculates the random aa identity.

    Parameters
    ----------
    aa_prop_csv_in: csv
        input csv file containing aa propensity for a particular sequence or dataset.
        Typically obtained from the function calc_aa_propensity_from_csv_col

    rand_seq_ident_csv_out: csv
        Output csv file contaning calculated random aa identity (due to limited aa propensity), and all the input values

    """
    # open csv into a pandas series, normally with all 20 aa as the index, and a proportion (0.08, 0.09 etc) as the data.
    aa_prop_ser = pd.Series.from_csv(aa_prop_csv_in, sep="\t")
    # calculate the random identity mathematically
    random_aa_identity = calc_random_aa_ident_from_ser(aa_prop_ser)
    # create a output series to contain all the output information
    output_ser = pd.Series()
    output_ser["random_sequence_identity_output"] = random_aa_identity
    aa_prop_ser.index = aa_prop_ser.index + "_input"
    output_ser = pd.concat([output_ser, aa_prop_ser])
    # save the setries as csv file
    output_ser.to_csv(rand_seq_ident_csv_out, sep="\t")
    sys.stdout.write("calc_random_aa_ident is finished\n")

def calc_random_aa_ident_via_randomisation(aa_prop_csv_in, rand_seq_ident_csv_out, seq_len=1000, number_seq=1000, ident=0.0, multiprocessing_mode=False):
    """Calculation of random amino acid identity based on a particular amino acid propensity.

    Protein regions with a limited aa propensity (e.g. transmembrane regions) have a measurable amino
    acid identity EVEN IN NON-HOMOLOGUES. This is referred to here as the random amino acid identity.
    This formula takes the aa propensity of a sequence or dataset as an input, and calculates the random aa identity.

    Parameters
    ----------
    aa_prop_csv_in: csv
        input csv file containing aa propensity for a particular sequence or dataset.
        Typically obtained from the function calc_aa_propensity_from_csv_col

    rand_seq_ident_csv_out: csv
        outout csv file contaning calculated random aa identity (due to limited aa propensity), and all the input values

    seq_len: int
        length of randomly created sequences. To achieve a more plausible result using randomisation method,
        greater values (> 5000) are recommended. Defalut value: 1000

    number_seq: int
        number of aligned sequences. Larger values are recommended. Default value: 1000

    ident: float
        desired overall identity of randomly created sequence matrix. This will not affect the random aa identity,
        but smaller values might increase the precision of the calculation. Default value: 0.0 (completely randomised)

    Returns
    -------
    Only returns values when multiprocessing_mode = True (no file is saved)
    random_aa_identity : float
        Random amino acid identity, based on the limited AA propensity in the dataset.
    output_ser = pd.Series
        Output pandas series consisting of the random_aa_identity, and the input aa_prop_ser derived from aa_prop_csv_in.
    """

    # open csv into a pandas series, normally with all 20 aa as the index, and a proportion (0.08, 0.09 etc) as the data.
    aa_prop_ser = pd.Series.from_csv(aa_prop_csv_in, sep="\t")

    # extract aa array and propensity array
    aa_propensities = np.array(aa_prop_ser)
    aa_arr = np.array(aa_prop_ser.index)

    # calculate number of residues that need to be replaced based on the desired percentage identity.
    number_mutations = int(np.round(seq_len*(1 - ident)))

    # generate random sequences, extract the original reference sequence and the sequence cluster
    orig_and_mut_seqs = generate_random_seq(seq_len, number_seq, number_mutations, aa_arr, aa_propensities)
    # sys.stdout.write("\n")
    # sys.stdout.flush()
    # extract the original sequence, of which the matrix are variants
    orig_seq = orig_and_mut_seqs[0]
    # calculate aa propensity and find all used aa in the orig_seq
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
        aa_prop_ser_at_that_pos = calc_aa_propensity(string_aa_in_matrix_at_that_pos)
        # add the amino acid propensities as a new column in the dataframe, with the orig_aa number as the column name
        columnwise_aa_propensities_df[n] = aa_prop_ser_at_that_pos
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

    # create a output series to contain all the output information
    output_ser = pd.Series()
    output_ser["random_sequence_identity_output"] = random_aa_identity
    aa_prop_ser.index = aa_prop_ser.index + "_input"
    output_ser = pd.concat([output_ser, aa_prop_ser])
    if multiprocessing_mode == False:
        # save the setries as csv file
        output_ser.to_csv(rand_seq_ident_csv_out, sep="\t")
    else:
        return random_aa_identity, output_ser


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
        # sys.write something to show that the programming is still running
        if n != 0 and n % 200 == 0:
            sys.stdout.write(".")
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

    # convert the seq_matrix into a np.array to ease further steps (slicing columns)
    seq_matrix = np.array(seq_matrix)

    return orig_seq, seq_matrix

def calc_unobserved_AA_ident(obs_aa_ident_full_protein, rand_ident_region1, rand_ident_region2, fraction_region1_residues, region):
    """Calculate the unobserved amino acid identity from a pairwise or multiple sequence alignment.

    The aim here is to calculate "e" as described in the associated research manuscript.
        e = rx
            x is first calculated using a previously described function
            r is taken from the input variables
            e is returned

    Parameters
    ----------
    obs_aa_ident_full_protein : float
        Observed amino acid identity of the full pairwise alignment
        "observed" in that it always includes some position that have had a substitution to the same AA (e.g. Leu to Leu), which are unobserved.
        Typically excluding gaps, which are heavily biased by splice variants, etc, and are excluded from AAIMON calculations.
        E.g. 0.6, for 60% amino acid identity.
    rand_TM : float
        Random identity of the TM region.
        e.g. 0.124 for 12.4%
    rand_nonTM : float
        Random identity of the nonTM region.
        e.g. 0.059 for 5.9%
    fraction_TM_residues : float
        Fraction of TM residues in protein sequence, e.g. 0.1 (10% TM, 90% nonTM)
        Used to estimate the real underlying AA substitution rate from the observed AA subst. rate.
    region: int
        Region of interest, either 1 or 2.
        Determines which random identity is used for the calculation of e.

    Returns
    -------
    e : float
        Unobserved substitution rate.
        Amino acid identity due solely to the restricted subset of amino acids.

    Usage
    -----
    obs_aa_ident_full_protein = 0.55
    rand_ident_region1 = 0.124
    rand_ident_region2 = 0.058
    fraction_region1_residues = 0.33
    e = calc_unobserved_AA_ident(obs_aa_ident_full_protein, rand_ident_region1, rand_ident_region2, fraction_region1_residues, region)
    """

    # calculate the real, underlying substitution rate in the full protein
    x = korbinian.cons_ratio.norm.calc_real_underlying_subst_rate(obs_aa_ident_full_protein, rand_ident_region1, rand_ident_region2, fraction_region1_residues)
    # get the random identity of the appropriate region, either region1 or region2 (e.g. TM or nonTM)
    if region == 1:
        r = rand_ident_region1
    else:
        r = rand_ident_region2
    # the unobserved substitution rate, e, (a.ka. unobserved AA identity) is simply the random identity
    # multiplied by the real, underlying substitution rate
    e = r * x

    return e

def normalise_AA_ident(obs_aa_ident, obs_aa_ident_full_protein, rand_ident_region1, rand_ident_region2, fraction_region1_residues, region):
    """Normalise the amino acid identity to remove conservation associated with limited amino acid propensities,
    by subtracting the unobserved amino acid identity.

    Parameters
    ----------
    obs_aa_ident
    obs_aa_ident_full_protein
    rand_ident_region1
    rand_ident_region2
    fraction_region1_residues
    region

    Returns
    -------
    a : float
        Normalised amino acid identity.
        Original amino acid identity minus the predicted identity attributable to the
        limited subset of amino acids.

    Usage
    -----
    import numpy as np
    from korbinian.MSA_normalisation import subtract_unobserved_AA_ident
    list_of_aa_ident_TM = np.array([0.48, 0.56, 0.67, 0.35])
    list_of_aa_ident_nonTM = np.array([0.44, 0.53, 0.62, 0.36])
    obs_aa_ident_full_protein = 0.30
    rand_ident_region1 = 0.09
    rand_ident_region2 = 0.058
    fraction_region1_residues = 0.33
    norm_list_of_aa_ident_TM = subtract_unobserved_AA_ident(list_of_aa_ident_TM, obs_aa_ident_full_protein, rand_ident_region1,
                                                            rand_ident_region2, fraction_region1_residues, region=1)
    norm_list_of_aa_ident_nonTM = subtract_unobserved_AA_ident(list_of_aa_ident_nonTM, obs_aa_ident_full_protein, rand_ident_region1,
                                                               rand_ident_region2, fraction_region1_residues, region=2)


    """

    # get the real underlying substitution rate (x) and therefore the unobserved substitution rate, e
    e = calc_unobserved_AA_ident(obs_aa_ident_full_protein, rand_ident_region1, rand_ident_region2, fraction_region1_residues, region)
    # simply subtract the AA identity attributable to the limited subset of AA
    a = obs_aa_ident - e

    return a


############################################################################################
#                                                                                          #
#                        calculation of normalisation factor                               #
#                                  based on aa identity                                    #
#                                                                                          #
############################################################################################

def DEPRECATED_calc_MSA_ident_n_factor(observed_perc_ident_full_seq, rand_perc_ident_TM, rand_perc_ident_nonTM, proportion_seq_TM_residues=0.3):
    """DEPRECATED. Use Calculation of the MSA identity normalisation factor

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
    MSA_aa_ident_norm_factor: float
        normalisation factor which will be applied to your observed TM identity

    Example:
    observed_perc_ident_TM = 0.78, rand_perc_ident_TM = 0.126, rand_perc_ident_nonTM = 0.059
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
    MSA_aa_ident_norm_factor = observed_perc_ident_TM/observed_perc_ident_nonTM

    #sys.stdout.write('\nnormalisation factor: %.3f' %MSA_TM_nonTM_aa_ident_norm_factor)

    return MSA_aa_ident_norm_factor

def create_matrix_artificial_homologues(aa_prop_ser, seq_len, number_seq, number_mutations):
    """ Create a matrix (array) of artificially generated homologue sequences.

    1) creates an original random sequence, based on a given aa propensity
    2) creates an array of "homologues"
    3) in each "homologue", replaces some original aa with a randomly chosen AA, based on the AA propensity

    Parameters
    ----------
    aa_prop_ser : pd.Series
        Amino acid propensity series
            index = A, C, D etc
            values = 0.05, 0.06, 0.14, etc
    seq_len : int
        Length of generated sequences.
    number_seq : int
        Numbef of sequences to be generated
    number_mutations : int
        number of mutations to introduce into the "homologues"

    Returns
    -------
    orig_seq : str
        original randomly generated sequence
    matrix : list
        list of strings of artificial homologues
    """

    # create the original template sequence
    orig_seq = "".join(np.random.choice(aa_prop_ser.index, p=aa_prop_ser) for _ in range(int(seq_len)))

    matrix = []
    for n in range(number_seq):
        # create indices for each AA in orig sequence
        inds = list(range(seq_len))
        # choose a random sample of AA to mutation
        sam = random.sample(inds, number_mutations)
        # convert orig sequence to a list
        seq_list = list(orig_seq)
        # for each index in the random sample, replace the AA with a random AA
        for ind in sam:
            seq_list[ind] = np.random.choice(aa_prop_ser.index, p=aa_prop_ser)
        # # join to make a new sequence
        # new_seq = "".join(seq_list)
        # # append to the matrix of "homologues"
        # matrix.append(new_seq)
        matrix.append(seq_list)
    # convert to a numpy array
    matrix = np.array(matrix)

    return orig_seq, matrix


def count_aa_freq(seq):
    # function to calculate aa propensity for each residue
    # input should be a string
    # output is a pd.DataFrame
    aa_dict = {}
    for aa in seq:
        aa_dict[aa] = 0
    for aa in seq:
        aa_dict[aa] += 1

    prop_dict = {}
    for aa in aa_dict:
        # sys.stdout.write(aa)
        prop_dict[aa] = aa_dict['%s' % aa] / len(seq)

    df = pd.Series(prop_dict)
    df = pd.DataFrame(df, columns=['freq'])
    df = df.transpose()
    return df


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


############################################################################################
#                                                                                          #
#                        Using argparse to enable usage from                               #
#                                  command line                                            #
#                                                                                          #
############################################################################################

# create a parser object to read user inputs from the command line
parser = argparse.ArgumentParser()
# add command-line options
parser.add_argument("-f", "--function",
                    required=True,
                    help=r"Function to be run. Choices are calc_aa_prop, calc_rand_aa_ident, or calc_unobserved_AA_ident")
parser.add_argument("-i", "--input",
                    default=None,
                    help=r'Full path of input file.'
                         r'E.g. "C:\Path\to\your\file.xlsx"')
parser.add_argument("-o",  "--output",
                    default=None,
                    help=r'Full path of output file.'
                         r'E.g. "C:\Path\to\your\file.xlsx"')
parser.add_argument("-c", "--column_name",
                    default=None,
                    help='Column name in input file that should be used for analysis.')
# INPUT for random identity via the randomisation method is not really necessary
# parser.add_argument("-l", "--length",
#                     default=1000,
#                     help='Sequence length for calc_rand_aa_ident.')
# parser.add_argument("-n", "--number_seq",
#                     default=1000,
#                     help='Number of sequences for calc_rand_aa_ident.')
# parser.add_argument("-d", "--ident_in_matrix",
#                     default=0.7,
#                     help='Amino acid identity in mutation matrix for calc_rand_aa_ident.')
parser.add_argument("-x", "--full_length_identity",
                    default=None,
                    help='Average amino acid identity of full sequences in alignment for calc_MSA_n_factor.')
parser.add_argument("-a", "--rand_aa_ident_1",
                    default=None,
                    help='Random aa identity for region 1 (e.g. transmembrane).')
parser.add_argument("-b", "--rand_aa_ident_2",
                    default=None,
                    help='Random aa identity for region 2 (e.g. non-transmembrane).')
parser.add_argument("-af", "--fraction_of_region1_in_full_protein",
                    default=0.3,
                    help="""Average fraction of sequence that is from region 1 (e.g. fract of residues that are
                    transmembrane residues) for use in calc_MSA_n_factor.""")
parser.add_argument("-r", "--region",
                    default=1,
                    help="""Either 1 or 2. Region of position of interest in sequence.""")


# if MSA_normalisation.py is run as the main python script, obtain the options from the command line.
if __name__ == '__main__':
    #sys.stdout.write("\nFor help, run \npython MSA_normalisation.py -h\n")

    # obtain command-line arguments
    args = parser.parse_args()

    if "h" in args:
        parser.print_help()
        sys.stdout.write("")

    if args.function == "calc_aa_prop":
        # write a message if any of the necessary arguments are missing
        if None in [args.input, args.output, args.column_name]:
            sys.stdout.write("Argument missing. Please run as follows:\npython MSA_normalisation.py "
                             "-f calc_aa_prop -i INPUTFILE -o OUTPUTFILE -c COLUMN_NAME")
        else:
            # use the supplied inputs from the command line to run calc_aa_propensity_from_csv_col
            calc_aa_propensity_from_csv_col(seq_list_csv_in=args.input, aa_prop_csv_out=args.output, col_name=args.column_name)
            #sys.stdout.write("\nFinished. Output file : {}".format(args.output))
    elif args.function == "calc_rand_aa_ident":
        # write a message if any of the necessary arguments are missing
        if None in [args.input, args.output]:
            sys.stdout.write("Input or output file path is missing. Please run as follows:\npython MSA_normalisation.py "
                             "-f calc_rand_aa_ident -i INPUTFILE -o OUTPUTFILE -s/-l/-n/-d OPTIONAL_ARGUMENTS")
        else:
            # use the supplied inputs from the command line to run calc_random_aa_ident
            calc_random_aa_ident(aa_prop_csv_in=args.input, rand_seq_ident_csv_out=args.output)
            #sys.stdout.write("\nFinished. Output file : {}".format(args.output))
    elif args.function == "calc_unobserved_AA_ident":
        # write a message if any of the necessary arguments are missing
        if None in [args.full_length_identity, args.rand_aa_ident_1, args.rand_aa_ident_2]:
            sys.stdout.write("Input argument missing. Please run as follows:\npython MSA_normalisation.py "
                             "-f calc_n_factor -x full_length_identity -a rand_aa_ident_1 -b rand_aa_ident_2 "
                             "-o OPTIONAL_OUTPUT_FILE -af OPTIONAL_fraction_of_region1_in_full_protein")
        else:
            # use the supplied inputs from the command line to run calc_unobserved_AA_ident
            full_length_identity = float(args.full_length_identity)
            rand_aa_ident_1 = float(args.rand_aa_ident_1)
            rand_aa_ident_2 = float(args.rand_aa_ident_2)
            fraction_of_region1_in_full_protein = float(args.fraction_of_region1_in_full_protein)
            region = float(args.region)
            unobserved_AA_ident = calc_unobserved_AA_ident(obs_aa_ident_full_protein=full_length_identity,
                                               rand_ident_region1=rand_aa_ident_1,
                                               rand_ident_region2=rand_aa_ident_2,
                                               fraction_region1_residues=fraction_of_region1_in_full_protein,
                                               region=region)


            # if the user has supplied the path to an output file, add the output aa_prop_norm_factor_output and the inputs to a csv and save
            if args.output is not None:
                output_ser = pd.Series()
                output_ser["unobserved_AA_ident"] = unobserved_AA_ident
                output_ser["full_length_identity_input"] = full_length_identity
                output_ser["region"] = region
                output_ser["rand_aa_ident_region1_input"] = rand_aa_ident_1
                output_ser["rand_aa_ident_region2_input"] = rand_aa_ident_2
                output_ser["fraction_of_region1_in_full_protein"] = fraction_of_region1_in_full_protein
                output_ser.to_csv(args.output, sep="\t")
                #sys.stdout.write("\nFinished. Output file : {}".format(args.output))
            else:
                # there is no output file, print the result on the screen
                sys.stdout.write("\nUnobserved amino acid identity =\t{:0.4f}\n".format(unobserved_AA_ident))
                sys.stdout.write("\nTo remove the AA identity associated with the limited subset of amino acids,\n"
                                 "subtract the unobserved amino acid identity from the original observed AA identity.\n")
    else:
        raise ValueError('Function to be run "{}" is not recognised. '
                        'Accepted values are calc_aa_prop, calc_rand_aa_ident, or calc_n_factor'.format(args.f))