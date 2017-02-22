import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
import pickle
import korbinian
import sys

##########parameters#############
list_number = 1
data_dir = r"D:\Databases"
repeat_randomisation = False

seq_len = 1000
max_num_positions_mutated = 600
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
pickle_with_artificial_AAIMONs = os.path.join(data_dir, "summaries\{ln:02d}\List{ln:02d}_rand\List{ln:02d}_rand_AAIMONs.pickle".format(ln=list_number))
fig_AAIMON_vs_perc_aa_sub = os.path.join(data_dir, "summaries\{ln:02d}\List{ln:02d}_rand\List{ln:02d}_rand_AAIMON_vs_aa_sub.png".format(ln=list_number))

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

# creating random alignments is slow. Once carried out, save the data and re-use for graphing.
if not os.path.isfile(pickle_with_artificial_AAIMONs) or repeat_randomisation == True:

    ########################################################################################
    #                                                                                      #
    #              Create pairwise alignments of artificial homologue and calc % ident     #
    #                                                                                      #
    ########################################################################################

    nested_list_AAIMONs_3_replicates = []

    # make 3 replicates for each % identity (% real aa subst)
    for i in range(3):

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
            sys.stdout.write(".")
            if n % 60 == 0:
                sys.stdout.write("\n{}, {:.03f}, {:.03f}, {:.03f}, {:.03f}\n".format(number_mutations, perc_aa_subst, perc_ident_TM, perc_ident_nonTM, AAIMON))
            sys.stdout.flush()
        # add each list to a nested list (3 replicates)
        nested_list_AAIMONs_3_replicates.append(AAIMON_list)

    # save to pickle file, so this does not need to be repeated for graphing
    with open(pickle_with_artificial_AAIMONs, "wb") as pkl_file:
        pickle.dump(nested_list_AAIMONs_3_replicates, pkl_file, -1)


########################################################################################
#                                                                                      #
#                             Normalisation and plotting                               #
#                                                                                      #
########################################################################################
# re-open pickle file created earlier
with open(pickle_with_artificial_AAIMONs, "rb") as pkl_file:
    nested_list_AAIMONs_3_replicates = pickle.load(pkl_file)

fig, ax = plt.subplots()
color_nonnorm = "#EE762C"
#color_norm = "#0076B8"
colour_dict = korbinian.utils.create_colour_lists()
color_norm = colour_dict["TUM_colours"]['TUM5']
color_norm_line = "#53A7D5"

proportion_seq_TM_residues = 0.5
# calculate proportion of length of full sequence that is nonTM
proportion_seq_nonTM_residues = 1 - proportion_seq_TM_residues
# random percentage identity of the full protein, assuming 30% TM region and 70% nonTM region
rand_perc_ident_full_protein = proportion_seq_TM_residues * rand_TM + proportion_seq_nonTM_residues * rand_nonTM
print("rand_perc_ident_full_protein", rand_perc_ident_full_protein)

# The unobserved aa_substitition rate is a proportion of the real underlying subst rate

#unobserv_aa_sub_rate = real_perc_aa_subst_array * rand_TM
unobserv_aa_sub_rate = real_perc_aa_subst_array * rand_perc_ident_full_protein

obs_aa_sub_rate = real_perc_aa_subst_array - unobserv_aa_sub_rate
observed_perc_aa_ident_array = 1- obs_aa_sub_rate

print("unobserv_aa_sub_rate", unobserv_aa_sub_rate[::70])
print("obs_aa_sub_rate", obs_aa_sub_rate[::70])
print("observed_perc_aa_ident_array", observed_perc_aa_ident_array[::70])


vfunc = np.vectorize(korbinian.cons_ratio.norm.calc_AAIMON_aa_prop_norm_factor)

norm_factor_array = vfunc(observed_perc_aa_ident_array, rand_TM, rand_nonTM)

#norm_factor_array_obs_test = vfunc(observed_perc_aa_ident_array, rand_TM, rand_nonTM)
print("norm_factor_array", norm_factor_array[::70])


AAIMON_list = []
obs_aa_sub_rate_3rep = []
norm_factor_array_3rep = []
for AAIMON_sublist in nested_list_AAIMONs_3_replicates:
    AAIMON_list += AAIMON_sublist
    obs_aa_sub_rate_3rep = np.append(obs_aa_sub_rate_3rep, obs_aa_sub_rate)
    norm_factor_array_3rep = np.append(norm_factor_array_3rep, norm_factor_array)

norm_AAIMONs = AAIMON_list / norm_factor_array_3rep
ax.scatter(obs_aa_sub_rate_3rep*100, norm_AAIMONs, color=color_norm, s=3, alpha=1)
ax.scatter(obs_aa_sub_rate_3rep*100, AAIMON_list, color=color_nonnorm, s=3, alpha=0.8)


ax.set_ylabel("AAIMON")
ax.set_xlabel("% AA substitutions in full protein")
ax.set_xlim(0, 60)

fig.savefig(fig_AAIMON_vs_perc_aa_sub)

print("mean", norm_AAIMONs[-100:].mean())