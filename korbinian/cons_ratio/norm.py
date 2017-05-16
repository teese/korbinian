import korbinian
import matplotlib.pylab as pylab
import matplotlib.pyplot as plt
import os
import numpy as np
import random
import sys
# import debugging tools
from korbinian.utils import pr, pc, pn, aaa

def calc_real_underlying_subst_rate(obs_aa_ident_full_protein, rand_ident_region1, rand_ident_region2, fraction_region1_residues=0.3):
    """Calculation of the real, underlying substitution rate in a protein with two regions of differing random identity.

    Estimation of the underlying substitution rate is the first step in calculating a normalisation factor,
    to remove the region1/region2 conservation attributable to the restricted subset of amino acids.

    Please see the associated article manuscript (in preparation) for a full explanation of the method.

    Example calculation:
    region1 = transmembrane(TM)
    region2 = non-transmembrane(nonTM)
    i = obs_aa_ident_full_protein = 0.60  #(i.e. 60% amino acid identity, 40% substitutions)
    fraction_TM_residues = 0.1
    # the observed aa subst rate is 1 - obs_aa_ident_full_protein
    b = 1 - i = 1 - 0.60 = 0.4
    # proportion of seq that is Membranous
    m = fraction_TM_residues = 0.1
    # proportion of seq that is Soluble
    s = 1 - m = 1 - 0.1 = 0.90      # i.e. 90% of the residues are nonTM
    # random identity of TM region
    t = rand_ident_region1 = rand_TM = 0.12
    # random identity of nonTM region
    n = rand_ident_region2 = rand_nonTM = 0.06
    # real underlying aa subst rate for full protein
    # solved from b = mx - mtx + sx - snx
    x = b / ((m*-t) + m - n*s + s) = 0.428

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

    Returns
    -------
    x : float
        Amino acid propensity normalisation factor.
        Real, underlying substitution rate for the full protein.

    Usage
    -----
    # calculate the TM/nonTM conservation ratio for a particular query and homologue
    TM_over_nonTM_cons_ratio = percentage_identity_of_TMD / percentage_identity_of_nonTMD

    # calculate the amino acid propensity normalisation factor.
    obs_aa_ident_full_protein, rand_TM, rand_nonTM, fraction_TM_residues = 0.60, 0.136, 0.055, 0.10
    real_underlying_subst_rate = calc_real_underlying_subst_rate(obs_aa_ident_full_protein, rand_TM, rand_nonTM, fraction_TM_residues)

    aa_prop_norm_factor = calc_aa_prop_norm_factor(real_underlying_subst_rate, rand_TM, rand_nonTM)

    # to normalise, divide the AAIMON by the aa_prop_norm_factor
    TM_over_nonTM_cons_ratio_normalised = TM_over_nonTM_cons_ratio / aa_prop_norm_factor
    """
    # the oBserved aa subst rate is 1 - obs_aa_ident_full_protein
    b = 1 - obs_aa_ident_full_protein
    # proportion of seq that is Membranous
    m = fraction_region1_residues
    # proportion of seq that is Soluble
    s = 1 - m
    # random identity of Tm region
    t = rand_ident_region1
    # random identity of NonTM region
    n = rand_ident_region2
    # real underlying aa subst rate for full protein
    # solved from b = mx - mtx + sx - snx
    x = b / ((m * -t) + m - n * s + s)

    return x


def calc_aa_prop_norm_factor(obs_aa_ident_full_protein, rand_ident_region1, rand_ident_region2, fraction_region1_residues=0.3):
    """Calculates the amino acid propensity normalisation factor for region1/region2 conservation ratios
    calculated from pairwise alignments.

    TM/nonTM conservation ratios increase for far homologues, due to the "random" conservation
    attributable to the amino acid propensities (low number of amino acids that are hydrophobic
    enough to enter membranes, resulting in increased amino acid identity amongst homologues).

    The TM/nonTM ratios can be normalised to account for this, based on:
        a) the "random" identity of the TM region (in Oberai 2009, calculated as 0.136. Effectively,
        this suggests that the random alignment of TM seqs will give an amino acid identity of 13.6%)
        b) the "random" identity of the nonTM region (generally 1/20, 0.05. In Oberai 2009, 0.055)
        c) the fraction of TM residues in the full protein
        c) the percentage amino acid identity of the homologue (full protein AA identity)

    This function takes these values as inputs, and returns the normalisation factor, which corresponds to the
    TM/nonTM identity of random alignments, where the only difference is due to aa propensity.

    How can the randTM and rand_nonTM be calculated?
    This correlates to the percentage of "unobserved" substitutions within each region, which
    result in exactly the same amino acid, due to the restricted AA propensity. There are several
    ways that this could be calculated.
        - for AA identity, by taking the sum of the square of the AA propensities for each AA
        - by aligning non-homologous
        - by creating and aligning random sequences, based on a particular amino acid propensity
        - plugging the AA propensity into the scoring formula, for e.g. Entropy
    Note: it is also very likely that this can be calculated based on dS/dN data.

    Example calculation:
    obs_aa_ident_full_protein, rand_TM, rand_nonTM, fraction_TM_residues = 0.60, 0.12, 0.06, 0.1
    # the oBserved aa subst rate is 1 - obs_aa_ident_full_protein
    b = 1 - obs_aa_ident_full_protein = 1 - 0.60 = 0.4
    # proportion of seq that is Membranous
    m = fraction_TM_residues = 0.1
    # proportion of seq that is Soluble
    s = 1 - m = 1 - 0.1 = 0.90
    # random identity of Tm region
    t = rand_TM = 0.12
    # random identity of NonTM region
    n = rand_nonTM = 0.06
    # real underlying aa subst rate for full protein
    # solved from b = mx - mtx + sx - snx
    x = b / ((m*-t) + m - n*s + s) = 0.428

    # since we only want the ratios within TM and nonTM, let m = 1 and s = 1
    unobserved_aa_subst_rate_TM = x * t = 0.428 * 0.12 = 0.051
    unobserved_aa_subst_rate_nonTM = x * n = 0.428 * 0.06 = 0.0256
    # The real aa ident = 1 - real aa subst. rate
    real_aa_ident_full_protein = 1 - x = 0.572
    # observed AA conservation for TM or nonTM
    # = real AA identity, plus a proportion of unobserved AA identity
    obs_aa_cons_TM = real_aa_ident_full_protein + unobserved_aa_subst_rate_TM = 0.572 + 0.051 = 0.623
    obs_aa_cons_nonTM = real_aa_ident_full_protein + unobserved_aa_subst_rate_nonTM = 0.572 + 0.0256 = 0.5976

    # artificial AAIMON, if AA propensity is the only underlying factor
    AAIMON = obs_aa_cons_TM / obs_aa_cons_nonTM = 0.623 / 0.5976 = 1.043

    aa_prop_norm_factor = AAIMON

    This implies that with homologues of 60% amino acid identity, the baseline AAIMON ratio will be 1.043,
    simply due to the reduced amino acid propensity of the hydrophobic residues.

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

    Returns
    -------
    aa_prop_norm_factor : float
        Amino acid propensity normalisation factor.
        Equivalent to the AAIMON ratio of randomised sequences, which have exactly the same
        real underlying substitution rate, and differ in only AA propensity.

    Usage
    -----
    # calculate the Amino Acid Identity: Membranous Over Nonmembranous (AAIMON) ratio
    AAIMON = percentage_identity_of_TMD / percentage_identity_of_nonTMD

    # calculate the amino acid propensity normalisation factor.
    obs_aa_ident_full_protein, rand_TM, rand_nonTM, fraction_TM_residues = 0.60, 0.136, 0.055, 0.10
    aa_prop_norm_factor = calc_aa_prop_norm_factor(obs_aa_ident_full_protein, rand_TM, rand_nonTM, fraction_TM_residues)

    # to normalise, divide the AAIMON by the aa_prop_norm_factor
    AAIMON_normalised = AAIMON_orig / aa_prop_norm_factor
    """
    # calculate the real, underlying substitution rate
    x = calc_real_underlying_subst_rate(obs_aa_ident_full_protein, rand_ident_region1, rand_ident_region2, fraction_region1_residues)

    # define the random identity values
    t = rand_ident_region1
    n = rand_ident_region2

    # since we only want the ratios within TM and nonTM, let m = 1 and s = 1
    m = 1
    s = 1
    # as described in the illustrations in the supplementary data,
    # the unobserved substitutions are the product of the percentage_region, random_ident, and real substitution rate
    unobserved_aa_subst_rate_TM = m * t * x
    unobserved_aa_subst_rate_nonTM = s * n * x

    # The real aa ident = 1 - real aa subst. rate
    real_aa_ident_full_protein = 1 - x
    # observed AA conservation for TM or nonTM
    # Equals the real AA identity, plus a proportion of unobserved AA identity
    obs_aa_cons_TM = real_aa_ident_full_protein + unobserved_aa_subst_rate_TM
    obs_aa_cons_nonTM = real_aa_ident_full_protein + unobserved_aa_subst_rate_nonTM

    # artificial AAIMON ratio, if AA propensity is the only underlying factor
    region1_over_region2_conservation_ratio = obs_aa_cons_TM / obs_aa_cons_nonTM
    # the normalisation factor is defined as the conservation ratio attributable purely to the
    # restricted subset of amino acids
    aa_prop_norm_factor = region1_over_region2_conservation_ratio

    return aa_prop_norm_factor



def old2_calc_AAIMON_aa_prop_norm_factor(obs_aa_ident_full_protein, rand_TM, rand_nonTM, fraction_TM_residues=0.3):
    """OLD METHOD 02.
    Calculates the amino acid propensity normalisation factor for AAIMON ratios calculated from pairwise alignments.

    TM/nonTM conservation ratios increase for far homologues, due to the "random" conservation
    attributable to the amino acid propensities (low number of amino acids that are hydrophobic
    enough to enter membranes, resulting in increased amino acid identity amongst homologues).

    The TM/nonTM ratios can be normalised to account for this, based on:
        a) the "random" identity of the TM region (in Oberai 2009, calculated as 0.136. Effectively,
        this suggests that the random alignment of TM seqs will give an amino acid identity of 13.6%)
        b) the "random" identity of the nonTM region (generally 1/20, 0.05. In Oberai 2009, 0.055)
        c) the fraction of TM residues in the full protein
        c) the percentage amino acid identity of the homologue (full protein AA identity)

    This function takes these values as inputs, and returns the normalisation factor, which corresponds to the
    TM/nonTM identity of random alignments, where the only difference is due to aa propensity.

    How can the randTM and rand_nonTM be calculated?
    This correlates to the percentage of "unobserved" substitutions within each region, which
    result in exactly the same amino acid, due to the restricted AA propensity. There are several
    ways that this could be calculated.
        - for AA identity, by taking the sum of the square of the AA propensities for each AA
        - by aligning non-homologous
        - by creating and aligning random sequences, based on a particular amino acid propensity
        - plugging the AA propensity into the scoring formula, for e.g. Entropy
    Note: it is also very likely that this can be calculated based on dS/dN data.

    Example calculation:
    obs_aa_ident_full_protein, rand_TM, rand_nonTM, fraction_TM_residues = 0.60, 0.12, 0.06, 0.1
    # the oBserved aa subst rate is 1 - obs_aa_ident_full_protein
    b = 1 - obs_aa_ident_full_protein = 1 - 0.60 = 0.4
    # proportion of seq that is Membranous
    m = fraction_TM_residues = 0.1
    # proportion of seq that is Soluble
    s = 1 - m = 1 - 0.1 = 0.90
    # random identity of Tm region
    t = rand_TM = 0.12
    # random identity of NonTM region
    n = rand_nonTM = 0.06
    # real underlying aa subst rate for full protein
    # solved from b = mx - mtx + sx - snx
    x = b / ((m*-t) + m - n*s + s) = 0.428

    # since we only want the ratios within TM and nonTM, let m = 1 and s = 1
    unobserved_aa_subst_rate_TM = x * t = 0.428 * 0.12 = 0.051
    unobserved_aa_subst_rate_nonTM = x * n = 0.428 * 0.06 = 0.0256
    # The real aa ident = 1 - real aa subst. rate
    real_aa_ident_full_protein = 1 - x = 0.572
    # observed AA conservation for TM or nonTM
    # = real AA identity, plus a proportion of unobserved AA identity
    obs_aa_cons_TM = real_aa_ident_full_protein + unobserved_aa_subst_rate_TM = 0.572 + 0.051 = 0.623
    obs_aa_cons_nonTM = real_aa_ident_full_protein + unobserved_aa_subst_rate_nonTM = 0.572 + 0.0256 = 0.5976

    # artificial AAIMON, if AA propensity is the only underlying factor
    AAIMON = obs_aa_cons_TM / obs_aa_cons_nonTM = 0.623 / 0.5976 = 1.043

    aa_prop_norm_factor = AAIMON

    This implies that with homologues of 60% amino acid identity, the baseline AAIMON ratio will be 1.043,
    simply due to the reduced amino acid propensity of the hydrophobic residues.

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

    Returns
    -------
    aa_prop_norm_factor : float
        Amino acid propensity normalisation factor.
        Equivalent to the AAIMON ratio of randomised sequences, which have exactly the same
        real underlying substitution rate, and differ in only AA propensity.

    Usage
    -----
    # calculate the Amino Acid Identity: Membranous Over Nonmembranous (AAIMON) ratio
    AAIMON = percentage_identity_of_TMD / percentage_identity_of_nonTMD

    # calculate the amino acid propensity normalisation factor.
    obs_aa_ident_full_protein, rand_TM, rand_nonTM, fraction_TM_residues = 0.60, 0.136, 0.055, 0.10
    aa_prop_norm_factor = calc_AAIMON_aa_prop_norm_factor(obs_aa_ident_full_protein, rand_TM, rand_nonTM, fraction_TM_residues)

    # to normalise, divide the AAIMON by the aa_prop_norm_factor
    AAIMON_normalised = AAIMON_orig / aa_prop_norm_factor
    """
    # the oBserved aa subst rate is 1 - obs_aa_ident_full_protein
    b = 1 - obs_aa_ident_full_protein
    # proportion of seq that is Membranous
    m = fraction_TM_residues
    # proportion of seq that is Soluble
    s = 1 - m
    # random identity of Tm region
    t = rand_TM
    # random identity of NonTM region
    n = rand_nonTM
    # real underlying aa subst rate for full protein
    # solved from b = mx - mtx + sx - snx
    x = b / ((m * -t) + m - n * s + s)

    # since we only want the ratios within TM and nonTM, let m = 1 and s = 1
    m = 1
    s = 1
    unobserved_aa_subst_rate_TM = m * t * x
    unobserved_aa_subst_rate_nonTM = s * n * x
    # The real aa ident = 1 - real aa subst. rate
    real_aa_ident_full_protein = 1 - x
    # observed AA conservation for TM or nonTM
    # Equals the real AA identity, plus a proportion of unobserved AA identity
    obs_aa_cons_TM = real_aa_ident_full_protein + unobserved_aa_subst_rate_TM
    obs_aa_cons_nonTM = real_aa_ident_full_protein + unobserved_aa_subst_rate_nonTM

    # artificial AAIMON ratio, if AA propensity is the only underlying factor
    AAIMON = obs_aa_cons_TM / obs_aa_cons_nonTM

    aa_prop_norm_factor = AAIMON

    return aa_prop_norm_factor


def OLD_calc_AAIMON_aa_prop_norm_factor(aa_ident, rand_TM, rand_nonTM):
    """OLD method 01, assuming TM and nonTM can have different underlying real AA subst. rates.

    Use NEW method instead!

    Calculates the amino acid propensity normalisation factor for homologues of a particular amino acid identity.

    TM/nonTM conservation ratios increase for far homologues, due to the "random" conservation
    attributable to the amino acid propensities (low number of amino acids that are hydrophobic
    enough to enter membranes, resulting in increased amino acid identity amongst homologues).

    The TM/nonTM ratios can be normalised to account for this, based on:
        a) the "random" identity of the TM region (in Oberai 2009, calculated as 0.136. Effectively,
        this suggests that the random alignment of TM seqs will give an amino acid identity of 13.6%)
        b) the "random" identity of the nonTM region (generally 1/20, 0.05. In Oberai 2009, 0.055)
        c) the percentage amino acid identity of the homologue (full protein AA identity)

    This function takes these three values as inputs, and returns the normalisation factor.

    How can the randTM and rand_nonTM be calculated?
    This correlates to the percentage of "unobserved" substitutions within each region, which
    result in exactly the same amino acid, due to the restricted AA propensity. There are several
    ways that this could be calculated.
        - by aligning non-homologous
        - by creating and aligning random sequences, based on a particular amino acid propensity
        - using the formula for entropy, and the propensity of each amino acid.
    Note: it is also very likely that this can be calculated based on dS/dN data.

    Example calculation:
    aa_ident, rand_TM, rand_nonTM = 0.60, 0.136, 0.055
    observed_changes = 1 - aa_ident = 1 - 0.6 = 0.4 (represented from now on as 40%)
    real_changes_TM = observed_changes / (1 - rand_TM) = 40% / (1-0.136) = 46.3%
    real_changes_nonTM = observed_changes / (1 - rand_nonTM) = 40% / (1-0.055) = 42.32%
    unobserved_changes_TM = real_changes_TM - observed_changes = 46.3% - 40% = 6.3%
    unobserved_changes_nonTM = real_changes_nonTM - observed_changes = 42.32% - 40% = 2.3%
    real_conserv_TM = aa_ident - unobserved_changes_TM = 60% - 6.3% = 53.7%
    real_conserv_nonTM = aa_ident - unobserved_changes_nonTM = 60% - 2.3% = 57.7%
    aa_prop_norm_factor = real_conserv_nonTM / real_conserv_TM =  57.7% / 53.7% = 1.074

    This implies that with homologues of 60% amino acid identity, the baseline AAIMON ratio will be 1.074,
    simply due to the reduced amino acid propensity of the hydrophobic residues.

    Parameters
    ----------
    aa_ident : float
        Amino acid identity (typically excluding gaps, which are heavily biased by splice variants, etc)
        E.g. 0.6, for 60% amino acid identity.
    rand_TM : float
        Random identity of the TM region.
    rand_nonTM : float
        Random identity of the nonTM region.

    Returns
    -------
    aa_prop_norm_factor : float
        Amino acid propensity normalisation factor.

    Usage
    -----
    # calculate the Amino Acid Identity: Membranous Over Nonmembranous (AAIMON) ratio
    AAIMON = percentage_identity_of_TMD / percentage_identity_of_nonTMD

    # calculate the amino acid propensity normalisation factor.
    aa_ident, rand_TM, rand_nonTM = 0.60, 0.136, 0.055
    aa_prop_norm_factor = calc_AAIMON_aa_prop_norm_factor(aa_ident, rand_TM, rand_nonTM)

    # to normalise, divide the AAIMON by the aa_prop_norm_factor
    AAIMON_normalised = AAIMON / aa_prop_norm_factor
    """
    # the proportion of residues where changes have occurred
    observed_changes = 1 - aa_ident
    # the real number of underlying mutations, assuming that a proportion (relative to rand_TM or rand_nonTM) is not
    # visible, as it has mutated to the same amino acid.
    real_changes_TM = observed_changes / (1 - rand_TM)
    real_changes_nonTM = observed_changes / (1 - rand_nonTM)
    # calculation of the unobserved changes (difference between observed and real changes)
    unobserved_changes_TM = real_changes_TM - observed_changes
    unobserved_changes_nonTM = real_changes_nonTM - observed_changes
    # calculation of real conservation (observed conservation - unobserved changes)
    real_conserv_TM = aa_ident - unobserved_changes_TM
    real_conserv_nonTM = aa_ident - unobserved_changes_nonTM
    # calculate normalisation factor
    aa_prop_norm_factor = real_conserv_nonTM / real_conserv_TM

    return aa_prop_norm_factor


def save_graph_for_normalized_AAIMON(acc, AAIMON, norm_AAIMON, aa_ident, zipout, protein_name):
    fig, ax = plt.subplots()
    perc_ident = aa_ident*100
    ax.scatter(perc_ident, norm_AAIMON)
#    ax.set_xticklabels(per_ident, rotation=45)
    ax.set_ylabel('normalised AAIMON', fontsize=14)
    pylab.ylim([0.6, 1.5])

    ax2 = ax.twinx()
    ax2.scatter(perc_ident, AAIMON, color='red', marker='^', alpha=0.5)
    pylab.ylim([0.6, 1.5])
    ax.set_xlabel('% identity of homologue', fontsize=14)
    xlim_min = 40
    xlim_max = 100
    ax.set_xlim(xlim_min, xlim_max)

#    fig.savefig(r'C:\Users\coffee oder tee\Dropbox\Undergraduate\Bachelor thesis\my_documents\plot_folder\%s_normalised_AAIMON.png' % acc, format="png", dpi=200)
    fig.savefig(protein_name + '_AAIMON_normalisation.png', format = 'png', dpi = 200)
    plt.close('all')
    zipout.write(protein_name + '_AAIMON_normalisation.png', arcname=os.path.basename(protein_name) + '_AAIMON_normalisation.png')
    os.remove(protein_name + '_AAIMON_normalisation.png')


def get_perc_ident_random_pairwise_alignment(aa_prop_ser, seq_len, number_mutations):
    """ Obtains the percentage identity between an original random sequence, and a partially randomised "homologue".

    Parameters
    ----------
    aa_prop_ser : pd.Series
        Amino acid propensity of both randomised orig seq, and aa substitution within "homologue"
        index = A, C, D, E, etc
        values = 0.05, 0.01, 0.02, etc
    seq_len : int
        Sequence length of randomised sequences
    number_mutations : int
        Number of mutations to be introduced within the "homologue"

    Returns
    -------
    perc_ident : float
        Percentage identity between random sequence, and partially randomised homologue
    """
    # create the original template sequence
    orig_seq = np.array([np.random.choice(aa_prop_ser.index, p=aa_prop_ser) for _ in range(int(seq_len))])
    # create indices for each AA in orig sequence
    inds = list(range(seq_len))
    # choose a random sample of AA to mutate
    sam = random.sample(inds, number_mutations)
    # convert orig sequence to a list
    homol_seq = list(orig_seq)
    # for each index in the random sample, replace the AA with a random AA
    for ind in sam:
        homol_seq[ind] = np.random.choice(aa_prop_ser.index, p=aa_prop_ser)
    # convert from list to numpy array
    homol_seq = np.array(homol_seq)
    # count how many residues are identical by making True/False array, and counting the "True" values
    aa_is_conserved_bool_array = orig_seq == homol_seq
    n_cons_aa = np.count_nonzero(aa_is_conserved_bool_array)
    # get percentage identity
    perc_ident = n_cons_aa / seq_len
    return perc_ident


def create_list_random_AAIMON_ratios(n_mutations_array, aa_prop_TM, aa_prop_nonTM, seq_len):
    """ Uses randomised pairwise alignments with "homologues" to calculates AAIMON ratios.
    Parameters
    ----------
    n_mutations_array : np.ndarray
        Array containing the number of mutations for each pairwise homologue alignment.
        E.g. [1,2,3,4...700], for up to 700 mutations in a sequence of length 1000 residues
    aa_prop_TM : pd.Series
        Amino acid propensity of TM region
        index = A, C, D, E, etc
        values = 0.05, 0.01, 0.02, etc
    aa_prop_nonTM : pd.Series
        Amino acid propensity of nonTM region
        index = A, C, D, E, etc
        values = 0.05, 0.01, 0.02, etc
    seq_len : int
        Sequence length of randomised sequences

    Returns
    -------
    AAIMON_list : list
        List of AAIMON ratios between random sequences, and partially randomised "homologues" at the desired
        evolutionary distance, based on the array of the number of mutations.
    """
    AAIMON_list = []
    # iterate through each number of mutations (1/1000, 2/1000 ..... 700/1000)
    for n, number_mutations in enumerate(n_mutations_array):
        # get a random pairwise alignment for the TM region
        perc_ident_TM = korbinian.cons_ratio.norm.get_perc_ident_random_pairwise_alignment(aa_prop_TM, seq_len, number_mutations)
        # get a random pairwise alignment for the nonTM region
        perc_ident_nonTM = korbinian.cons_ratio.norm.get_perc_ident_random_pairwise_alignment(aa_prop_nonTM, seq_len, number_mutations)
        # calculate artificial AAIMON
        AAIMON = perc_ident_TM / perc_ident_nonTM
        # append to AAIMON_list
        AAIMON_list.append(AAIMON)
        # calculate the percentage aa substitutions for this position
        perc_aa_subst = number_mutations / seq_len
        # if n % 100 == 0:
        #     sys.stdout.write("\n{}, {:.03f}, {:.03f}, {:.03f}, {:.03f}\n".format(number_mutations, perc_aa_subst, perc_ident_TM, perc_ident_nonTM, AAIMON))
        sys.stdout.flush()

    return AAIMON_list