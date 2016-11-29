import matplotlib.pylab as pylab
import matplotlib.pyplot as plt
import os

def calc_AAIMON_aa_prop_norm_factor(aa_ident, rand_TM, rand_nonTM):
    """Calculates the amino acid propensity normalisation factor for homologues of a particular amino acid identity.

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


