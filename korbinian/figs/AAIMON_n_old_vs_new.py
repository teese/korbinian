import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
import pickle
import korbinian
import sys

# def calc_aa_prop_norm_factor_assum_equal_real_aa_sub(obs_aa_ident_full_protein, rand_TM, rand_nonTM, fraction_TM_residues=0.3):
#     """Calculates the amino acid propensity normalisation factor for AAIMON ratios calculated from pairwise alignments.
#
#     TM/nonTM conservation ratios increase for far homologues, due to the "random" conservation
#     attributable to the amino acid propensities (low number of amino acids that are hydrophobic
#     enough to enter membranes, resulting in increased amino acid identity amongst homologues).
#
#     The TM/nonTM ratios can be normalised to account for this, based on:
#         a) the "random" identity of the TM region (in Oberai 2009, calculated as 0.136. Effectively,
#         this suggests that the random alignment of TM seqs will give an amino acid identity of 13.6%)
#         b) the "random" identity of the nonTM region (generally 1/20, 0.05. In Oberai 2009, 0.055)
#         c) the fraction of TM residues in the full protein
#         c) the percentage amino acid identity of the homologue (full protein AA identity)
#
#     This function takes these values as inputs, and returns the normalisation factor, which corresponds to the
#     TM/nonTM identity of random alignments, where the only difference is due to aa propensity.
#
#     How can the randTM and rand_nonTM be calculated?
#     This correlates to the percentage of "unobserved" substitutions within each region, which
#     result in exactly the same amino acid, due to the restricted AA propensity. There are several
#     ways that this could be calculated.
#         - by aligning non-homologous
#         - by creating and aligning random sequences, based on a particular amino acid propensity
#         - using the formula for entropy, and the propensity of each amino acid.
#     Note: it is also very likely that this can be calculated based on dS/dN data.
#
#     Example calculation:
#     obs_aa_ident_full_protein, rand_TM, rand_nonTM, fraction_TM_residues = 0.60, 0.12, 0.06, 0.1
#     # the oBserved aa subst rate is 1 - obs_aa_ident_full_protein
#     b = 1 - obs_aa_ident_full_protein = 1 - 0.60 = 0.4
#     # proportion of seq that is Membranous
#     m = fraction_TM_residues = 0.1
#     # proportion of seq that is Soluble
#     s = 1 - m = 1 - 0.1 = 0.90
#     # random identity of Tm region
#     t = rand_TM = 0.12
#     # random identity of NonTM region
#     n = rand_nonTM = 0.06
#     # real underlying aa subst rate for full protein
#     # solved from b = mx - mtx + sx - snx
#     x = b / ((m*-t) + m - n*s + s) = 0.428
#
#     # since we only want the ratios within TM and nonTM, let m = 1 and s = 1
#     unobserved_aa_subst_rate_TM = x * t = 0.428 * 0.12 = 0.051
#     unobserved_aa_subst_rate_nonTM = x * n = 0.428 * 0.06 = 0.0256
#     # The real aa ident = 1 - real aa subst. rate
#     real_aa_ident_full_protein = 1 - x = 0.572
#     # observed AA conservation for TM or nonTM
#     # = real AA identity, plus a proportion of unobserved AA identity
#     obs_aa_cons_TM = real_aa_ident_full_protein + unobserved_aa_subst_rate_TM = 0.572 + 0.051 = 0.623
#     obs_aa_cons_nonTM = real_aa_ident_full_protein + unobserved_aa_subst_rate_nonTM = 0.572 + 0.0256 = 0.5976
#
#     # artificial AAIMON, if AA propensity is the only underlying factor
#     AAIMON = obs_aa_cons_TM / obs_aa_cons_nonTM = 0.623 / 0.5976 = 1.043
#
#     aa_prop_norm_factor = AAIMON
#
#     This implies that with homologues of 60% amino acid identity, the baseline AAIMON ratio will be 1.043,
#     simply due to the reduced amino acid propensity of the hydrophobic residues.
#
#     Parameters
#     ----------
#     obs_aa_ident_full_protein : float
#         Observed amino acid identity of the full pairwise alignment
#         "observed" in that it always includes some position that have had a substitution to the same AA (e.g. Leu to Leu), which are unobserved.
#         Typically excluding gaps, which are heavily biased by splice variants, etc, and are excluded from AAIMON calculations.
#         E.g. 0.6, for 60% amino acid identity.
#     rand_TM : float
#         Random identity of the TM region.
#     rand_nonTM : float
#         Random identity of the nonTM region.
#     fraction_TM_residues : float
#         Fraction of TM residues in protein sequence, e.g. 0.1 (10% TM, 90% nonTM)
#         Used to estimate the real underlying AA substitution rate from the observed AA subst. rate.
#
#     Returns
#     -------
#     aa_prop_norm_factor : float
#         Amino acid propensity normalisation factor.
#         Equivalent to the AAIMON ratio of randomised sequences, which have exactly the same
#         real underlying substitution rate, and differ in only AA propensity.
#
#     Usage
#     -----
#     # calculate the Amino Acid Identity: Membranous Over Nonmembranous (AAIMON) ratio
#     AAIMON = percentage_identity_of_TMD / percentage_identity_of_nonTMD
#
#     # calculate the amino acid propensity normalisation factor.
#     obs_aa_ident_full_protein, rand_TM, rand_nonTM, fraction_TM_residues = 0.60, 0.136, 0.055, 0.10
#     aa_prop_norm_factor = calc_AAIMON_aa_prop_norm_factor(obs_aa_ident_full_protein, rand_TM, rand_nonTM, fraction_TM_residues)
#
#     # to normalise, divide the AAIMON by the aa_prop_norm_factor
#     AAIMON_normalised = AAIMON_orig / aa_prop_norm_factor
#     """
#     # the oBserved aa subst rate is 1 - obs_aa_ident_full_protein
#     b = 1 - obs_aa_ident_full_protein
#     # proportion of seq that is Membranous
#     m = fraction_TM_residues
#     # proportion of seq that is Soluble
#     s = 1 - m
#     # random identity of Tm region
#     t = rand_TM
#     # random identity of NonTM region
#     n = rand_nonTM
#     # real underlying aa subst rate for full protein
#     # solved from b = mx - mtx + sx - snx
#     x = b / ((m*-t) + m - n*s + s)
#
#     # since we only want the ratios within TM and nonTM, let m = 1 and s = 1
#     m = 1
#     s = 1
#     unobserved_aa_subst_rate_TM = m * t * x
#     unobserved_aa_subst_rate_nonTM = s * n * x
#     # The real aa ident = 1 - real aa subst. rate
#     real_aa_ident_full_protein = 1 - x
#     # observed AA conservation for TM or nonTM
#     # Equals the real AA identity, plus a proportion of unobserved AA identity
#     obs_aa_cons_TM = real_aa_ident_full_protein + unobserved_aa_subst_rate_TM
#     obs_aa_cons_nonTM = real_aa_ident_full_protein + unobserved_aa_subst_rate_nonTM
#
#     # artificial AAIMON ratio, if AA propensity is the only underlying factor
#     AAIMON = obs_aa_cons_TM / obs_aa_cons_nonTM
#
#     aa_prop_norm_factor = AAIMON
#
#     return aa_prop_norm_factor

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

sys.stdout.write("rand_TM {}".format(rand_TM))
sys.stdout.write("rand_nonTM {}".format(rand_nonTM))


rand_TM = 0.13
rand_nonTM = 0.05
fraction_TM_residues = 0.5

#obs_aa_ident_full_protein = 0.4

obs_aa_ident_full_protein_list = np.linspace(1.0,0.3, 100)
output_norm_list_old = []
output_norm_list_new = []


for obs_aa_ident_full_protein in obs_aa_ident_full_protein_list:
    old = korbinian.cons_ratio.norm.calc_aa_prop_norm_factor(obs_aa_ident_full_protein, rand_TM, rand_nonTM)
    #sys.stdout.write("old", old)
    output_norm_list_old.append(old)

    new = calc_aa_prop_norm_factor_assum_equal_real_aa_sub(obs_aa_ident_full_protein, rand_TM, rand_nonTM, fraction_TM_residues)
    output_norm_list_new.append(new)
    #sys.stdout.write("new", new)
    #sys.stdout.write()


colour_dict = korbinian.utils.create_colour_lists()

obs_aa_sub_full_prot = 1- np.array(obs_aa_ident_full_protein_list)

fig, ax = plt.subplots()
ax.scatter(obs_aa_sub_full_prot, output_norm_list_old, color=colour_dict["TUM_colours"]['TUM5'], label="old")
ax.scatter(obs_aa_sub_full_prot, output_norm_list_new, color=colour_dict["TUM_colours"]['TUM1'], label="new")
ax.legend(loc="upper left")
ax.set_xlabel("observed aa subst. full protein")
ax.set_ylabel("norm factor (AAIMON calculated for random seq)")
fig.savefig(fig_out)

sys.stdout.write(calc_aa_prop_norm_factor_assum_equal_real_aa_sub(0.60, 0.12, 0.06, 0.1))