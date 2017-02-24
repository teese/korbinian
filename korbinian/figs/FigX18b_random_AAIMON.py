import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
import pickle
import korbinian
import sys
from multiprocessing import Pool



##########parameters#############
list_number = 1
data_dir = r"D:\Databases"
repeat_randomisation = False

seq_len = 2000
max_num_positions_mutated = 1600
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

def create_list_random_AAIMON_ratios_mp (input):
    n_mutations_array, aa_prop_TM, aa_prop_nonTM, seq_len = input
    AAIMON_list = korbinian.cons_ratio.norm.create_list_random_AAIMON_ratios(n_mutations_array, aa_prop_TM, aa_prop_nonTM, seq_len)
    return AAIMON_list

input_mp = (n_mutations_array, aa_prop_TM, aa_prop_nonTM, seq_len)
input_mp_list = [input_mp, input_mp, input_mp,input_mp, input_mp, input_mp,input_mp, input_mp, input_mp, input_mp]

if __name__ == "__main__":

    # creating random alignments is slow. Once carried out, save the data and re-use for graphing.
    if not os.path.isfile(pickle_with_artificial_AAIMONs) or repeat_randomisation == True:

        with Pool(processes=10) as p:
            nested_list_AAIMONs_3_replicates = p.map(create_list_random_AAIMON_ratios_mp, input_mp_list)

        # ########################################################################################
        # #                                                                                      #
        # #              Create pairwise alignments of artificial homologue and calc % ident     #
        # #                                                                                      #
        # ########################################################################################
        #
        # nested_list_AAIMONs_3_replicates = []
        #
        # # make 3 replicates for each % identity (% real aa subst)
        # for i in range(3):
        #     input_mp = n_mutations_array, aa_prop_TM, aa_prop_nonTM, seq_len
        #     AAIMON_list = create_list_random_AAIMON_ratios_mp (input_mp)
        #     # add each list to a nested list (3 replicates)
        #     nested_list_AAIMONs_3_replicates.append(AAIMON_list)

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

    # # OLD NORM FACTOR
    vfunc_old = np.vectorize(korbinian.cons_ratio.norm.OLD_calc_AAIMON_aa_prop_norm_factor)
    norm_factor_array_old = vfunc_old(observed_perc_aa_ident_array, rand_TM, rand_nonTM)

    # NEW NORM FACTOR
    #
    vfunc = np.vectorize(korbinian.cons_ratio.norm.calc_AAIMON_aa_prop_norm_factor)
    norm_factor_array = vfunc(observed_perc_aa_ident_array, rand_TM, rand_nonTM, proportion_seq_TM_residues)

    #norm_factor_array_obs_test = vfunc(observed_perc_aa_ident_array, rand_TM, rand_nonTM)
    # print("norm_factor_array", norm_factor_array[::70])
    #
    #
    AAIMON_list = []
    obs_aa_sub_rate_3rep = []
    norm_factor_array_3rep_new = []
    norm_factor_array_3rep_old = []
    for AAIMON_sublist in nested_list_AAIMONs_3_replicates:
        AAIMON_list += AAIMON_sublist
        obs_aa_sub_rate_3rep = np.append(obs_aa_sub_rate_3rep, obs_aa_sub_rate)
        norm_factor_array_3rep_old = np.append(norm_factor_array_3rep_old, norm_factor_array_old)
        norm_factor_array_3rep_new = np.append(norm_factor_array_3rep_new, norm_factor_array)

    norm_AAIMONs_old = AAIMON_list / norm_factor_array_3rep_old
    norm_AAIMONs_new = AAIMON_list / norm_factor_array_3rep_new
    ax.scatter(obs_aa_sub_rate_3rep*100, norm_AAIMONs_old, color="r", s=3, alpha=1, label="after norm., old method")
    ax.scatter(obs_aa_sub_rate_3rep*100, norm_AAIMONs_new, color=colour_dict["TUM_colours"]['TUM5'], s=3, alpha=1, label="after norm., new method")
    ax.scatter(obs_aa_sub_rate_3rep*100, AAIMON_list, color=color_nonnorm, s=3, alpha=0.8, label="non-normalised")
    ax.scatter(obs_aa_sub_rate_3rep*100, norm_factor_array_3rep_old, s=1, color="r")
    ax.scatter(obs_aa_sub_rate_3rep*100, norm_factor_array_3rep_new, s=1, color=colour_dict["TUM_colours"]['TUM5'])

    z = np.polyfit(obs_aa_sub_rate_3rep*100, AAIMON_list, 13)
    p = np.poly1d(z)
    fitted_values_y = p(obs_aa_sub_rate*100)
    ax.plot(obs_aa_sub_rate*100, fitted_values_y, color="0.5", label="polyfit", linewidth=2)

    z = np.polyfit(obs_aa_sub_rate_3rep*100, norm_AAIMONs_new, 13)
    p = np.poly1d(z)
    fitted_values_y = p(obs_aa_sub_rate*100)
    ax.plot(obs_aa_sub_rate*100, fitted_values_y, color=colour_dict["TUM_colours"]['TUM3'], label="polyfit to new norm. data", linewidth=2)

    z = np.polyfit(obs_aa_sub_rate_3rep*100, norm_AAIMONs_old, 13)
    p = np.poly1d(z)
    fitted_values_y = p(obs_aa_sub_rate*100)
    ax.plot(obs_aa_sub_rate*100, fitted_values_y, color="r", label="polyfit to old norm. data", linewidth=2)



    ax.set_ylabel("AAIMON")
    ax.set_xlabel("% observed AA substitutions in full protein")
    ax.set_xlim(0, 80)
    ax.set_title("Random AAIMON ratios can be normalised to ~1.0")
    ax.legend(loc="upper left")
    fig.savefig(fig_AAIMON_vs_perc_aa_sub)


    old_last = list(norm_AAIMONs_old[-100:]) + list(norm_AAIMONs_old[900:1000]) + list(norm_AAIMONs_old[1900:2000])
    new_last = list(norm_AAIMONs_new[-100:]) + list(norm_AAIMONs_new[900:1000]) + list(norm_AAIMONs_new[1900:2000])


    print("mean", np.array(old_last).mean())
    print("mean", np.array(new_last).mean())
