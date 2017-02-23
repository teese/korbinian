import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
import pickle
import korbinian
import sys

def calc_AAIMON_aa_prop_norm_factor_assum_equal_real_aa_sub(observed_aa_ident_full_protein, rand_TM, rand_nonTM, proportion_seq_TM_residues=0.3):
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
    # the oBserved aa subst rate is 1 - observed_aa_ident_full_protein
    b = 1 - observed_aa_ident_full_protein
    # proportion of seq that is Membranous
    m = proportion_seq_TM_residues
    # proportion of seq that is Soluble
    s = 1 - m
    # random identity of Tm region
    t = rand_TM
    # random identity of NonTM region
    n = rand_nonTM
    # real underlying aa subst rate for full protein
    # solved from b = mx - mtx + sx - snx
    x = b / ((m*-t) + m - n*s + s)

    # print("x", x)
    #
    # obs_aa_sub_TM = m*x - m*t*x
    # obs_aa_sub_nonTM = s*x - s*n*x
    # print("obs_aa_sub_TM", obs_aa_sub_TM)
    # print("obs_aa_sub_nonTM", obs_aa_sub_nonTM)
    #
    #
    # obs_cons_TM = 1 - obs_aa_sub_TM
    # obs_cons_nonTM = 1 - obs_aa_sub_nonTM
    #
    # print("obs_cons_TM", obs_cons_TM)
    # print("obs_cons_nonTM", obs_cons_nonTM)
    #
    # AAIMON = obs_cons_TM / obs_cons_nonTM
    # print("AAIMON", AAIMON)
    #
    # aa_prop_norm_factor = AAIMON
    # print("aa_prop_norm_factor", aa_prop_norm_factor)





    # obs_aa_sub_TM =  m - m*t
    # obs_aa_sub_nonTM = s - s*n
    # print("obs_aa_sub_TM", obs_aa_sub_TM)
    # print("obs_aa_sub_nonTM", obs_aa_sub_nonTM)

    #print("------------------------------------------------------")
    # since we only want the ratios within TM and nonTM, let m = 1 and s = 1
    unobserved_aa_subst_rate_TM = x * t
    unobserved_aa_subst_rate_nonTM = x * n
    # real aa ident = 1 - real aa subst. rate
    real_aa_ident_full_protein = 1 - x
    # observed AA conservation for TM or nonTM
    # = real AA identity, plus a proportion of unobserved AA identity
    obs_aa_cons_TM = real_aa_ident_full_protein + unobserved_aa_subst_rate_TM
    obs_aa_cons_nonTM = real_aa_ident_full_protein + unobserved_aa_subst_rate_nonTM

    #print("obs_aa_cons_TM", obs_aa_cons_TM)
    #print("obs_aa_cons_nonTM", obs_aa_cons_nonTM)
    # artificial AAIMON, if AA propensity is the only underlying factor
    AAIMON = obs_aa_cons_TM / obs_aa_cons_nonTM
    #print("AAIMON", AAIMON)

    aa_prop_norm_factor = AAIMON
    #print("aa_prop_norm_factor", aa_prop_norm_factor)


    # # calculate proportion of length of full sequence that is nonTM
    # proportion_seq_nonTM_residues = 1 - proportion_seq_TM_residues
    # # random percentage identity of the full protein, assuming 30% TM region and 70% nonTM region
    # rand_perc_ident_full_protein = proportion_seq_TM_residues * rand_TM + proportion_seq_nonTM_residues * rand_nonTM
    # # the proportion of residues where changes have occurred
    # observed_changes = 1 - aa_ident
    # # # the real number of underlying mutations, assuming that a proportion (relative to rand_TM or rand_nonTM) is not
    # # # visible, as it has mutated to the same amino acid.
    # # real_changes_full_prot = observed_changes / (1 - rand_TM)
    # # real_changes_nonTM = observed_changes / (1 - rand_nonTM)
    #
    # real_changes_full_prot = observed_changes / (1 - rand_perc_ident_full_protein)
    # # calculation of the unobserved changes (difference between observed and real changes)
    # unobserved_changes_TM = real_changes_TM - observed_changes
    # unobserved_changes_nonTM = real_changes_nonTM - observed_changes
    # # calculation of real conservation (observed conservation - unobserved changes)
    # real_conserv_TM = aa_ident - unobserved_changes_TM
    # real_conserv_nonTM = aa_ident - unobserved_changes_nonTM
    # calculate normalisation factor
    #aa_prop_norm_factor = real_conserv_nonTM / real_conserv_TM

    return aa_prop_norm_factor

##########parameters#############
list_number = 1
data_dir = r"D:\Databases"
repeat_randomisation = False

seq_len = 1000
max_num_positions_mutated = 700
n_mutations_array = np.arange(0, max_num_positions_mutated)
real_perc_aa_subst_array = n_mutations_array / seq_len
real_perc_aa_ident_array = 1 - real_perc_aa_subst_array

########################################################################################
#                                                                                      #
#                     Setup paths for input and output files                           #
#                                                                                      #
########################################################################################
# List_rand_TM = r"D:\Databases\summaries\01\List01_rand\List01_rand_TM.csv"
List_rand_TM = os.path.join(data_dir, "summaries\{ln:02d}\List{ln:02d}_rand\List{ln:02d}_rand_TM.csv".format(ln=list_number))
List_rand_nonTM = os.path.join(data_dir, "summaries\{ln:02d}\List{ln:02d}_rand\List{ln:02d}_rand_nonTM.csv".format(ln=list_number))
fig_out = os.path.join(data_dir, "summaries\{ln:02d}\List{ln:02d}_rand\List{ln:02d}_fig_compare_AAIMON_norm_old_new.png".format(ln=list_number))

########################################################################################
#                                                                                      #
#                 Get the AA propensity series for TM and nonTM                        #
#                                                                                      #
########################################################################################

aa_prop_TM = pd.Series.from_csv(List_rand_TM, sep="\t")
# the first line is the random identity. Extract and delete.
rand_TM = aa_prop_TM["random_sequence_identity_output"]
del aa_prop_TM["random_sequence_identity_output"]
aa_prop_TM.index = aa_prop_TM.index.str[:1]

aa_prop_nonTM = pd.Series.from_csv(List_rand_nonTM, sep="\t")
# the first line is the random identity. Extract and delete.
rand_nonTM = aa_prop_nonTM["random_sequence_identity_output"]
del aa_prop_nonTM["random_sequence_identity_output"]
aa_prop_nonTM.index = aa_prop_nonTM.index.str[:1]

print("rand_TM", rand_TM)
print("rand_nonTM", rand_nonTM)


rand_TM = 0.13
rand_nonTM = 0.05
proportion_seq_TM_residues = 0.5

#observed_aa_ident_full_protein = 0.4

observed_aa_ident_full_protein_list = np.linspace(1.0,0.3, 100)
output_norm_list_old = []
output_norm_list_new = []


for observed_aa_ident_full_protein in observed_aa_ident_full_protein_list:
    old = korbinian.cons_ratio.norm.calc_AAIMON_aa_prop_norm_factor(observed_aa_ident_full_protein, rand_TM, rand_nonTM)
    #print("old", old)
    output_norm_list_old.append(old)

    new = calc_AAIMON_aa_prop_norm_factor_assum_equal_real_aa_sub(observed_aa_ident_full_protein, rand_TM, rand_nonTM, proportion_seq_TM_residues)
    output_norm_list_new.append(new)
    #print("new", new)
    #print()


colour_dict = korbinian.utils.create_colour_lists()

obs_aa_sub_full_prot = 1- np.array(observed_aa_ident_full_protein_list)

fig, ax = plt.subplots()
ax.scatter(obs_aa_sub_full_prot, output_norm_list_old, color=colour_dict["TUM_colours"]['TUM5'], label="old")
ax.scatter(obs_aa_sub_full_prot, output_norm_list_new, color=colour_dict["TUM_colours"]['TUM1'], label="new")
ax.legend(loc="upper left")
ax.set_xlabel("observed aa subst. full protein")
ax.set_ylabel("norm factor (AAIMON calculated for random seq)")
fig.savefig(fig_out)