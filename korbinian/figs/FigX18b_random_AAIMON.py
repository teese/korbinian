import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
import pickle
import korbinian
import sys
from multiprocessing import Pool

##########parameters#############
list_number = 2
data_dir = r"/Volumes/Musik/Databases"
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
List_rand_TM = os.path.normpath(os.path.join(data_dir, "summaries/{ln:02d}/List{ln:02d}_rand/List{ln:02d}_rand_TM.csv".format(ln=list_number)))
List_rand_nonTM = os.path.normpath(os.path.join(data_dir, "summaries/{ln:02d}/List{ln:02d}_rand/List{ln:02d}_rand_nonTM.csv".format(ln=list_number)))
pickle_with_artificial_AAIMONs = os.path.normpath(os.path.join(data_dir, "summaries/{ln:02d}/List{ln:02d}_rand/List{ln:02d}_rand_AAIMONs.pickle".format(ln=list_number)))
fig_AAIMON_vs_perc_aa_sub_png = os.path.normpath(os.path.join(data_dir, "summaries/{ln:02d}/List{ln:02d}_rand/List{ln:02d}_rand_AAIMON_vs_aa_sub.png".format(ln=list_number)))
fig_AAIMON_vs_perc_aa_sub_pdf = os.path.normpath(os.path.join(data_dir, "summaries/{ln:02d}/List{ln:02d}_rand/List{ln:02d}_rand_AAIMON_vs_aa_sub.pdf".format(ln=list_number)))


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

sys.stdout.write("\nrand_TM {}".format(rand_TM))
sys.stdout.write("\nrand_nonTM {}".format(rand_nonTM))


########################################################################################
#                                                                                      #
#        Using multiprocessing, create AAIMONs from randomised sequences               #
#        THIS IS VERY SLOW!!!!                                                         #
#        For plotting, don't forget to change repeat_randomisation to False            #
#                                                                                      #
########################################################################################


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
            nested_list_AAIMONs_10_replicates = p.map(create_list_random_AAIMON_ratios_mp, input_mp_list)

        # save to pickle file, so this does not need to be repeated for graphing
        with open(pickle_with_artificial_AAIMONs, "wb") as pkl_file:
            pickle.dump(nested_list_AAIMONs_10_replicates, pkl_file, -1)

    ########################################################################################
    #                                                                                      #
    #                             Normalisation and plotting                               #
    #                                                                                      #
    ########################################################################################
    # re-open pickle file created earlier
    with open(pickle_with_artificial_AAIMONs, "rb") as pkl_file:
        nested_list_AAIMONs_10_replicates = pickle.load(pkl_file)

    fig, ax = plt.subplots()

    # set plotting parameters
    plt.style.use('seaborn-whitegrid')
    colour_dict = korbinian.utils.create_colour_lists()
    color_nonnorm = colour_dict["TUM_colours"]['TUM1']
    color_norm = colour_dict["TUM_colours"]['TUM2']
    color_norm_line = "#53A7D5"
    color_norm_factor_line = "k"
    alpha = 1
    s = 2
    marker = 'x'
    linewidths = 0.3
    fontsize = 14


    ########################################################################################
    #                                                                                      #
    #         calc. observed_perc_aa_ident_array, used as input for norm factor calc.      #
    #                                                                                      #
    ########################################################################################
    # set the proportion of the sequence that is TM
    # for typical randomisation, this is always 0.5
    proportion_seq_TM_residues = 0.5
    # calculate proportion of length of full sequence that is nonTM
    proportion_seq_nonTM_residues = 1 - proportion_seq_TM_residues
    # random percentage identity of the full protein, assuming 30% TM region and 70% nonTM region
    rand_perc_ident_full_protein = proportion_seq_TM_residues * rand_TM + proportion_seq_nonTM_residues * rand_nonTM
    sys.stdout.write("\nrand_perc_ident_full_protein {}".format(rand_perc_ident_full_protein))

    # The unobserved aa_substitition rate is a proportion of the real underlying subst rate

    #unobserv_aa_sub_rate = real_perc_aa_subst_array * rand_TM
    unobserv_aa_sub_rate = real_perc_aa_subst_array * rand_perc_ident_full_protein

    obs_aa_sub_rate = real_perc_aa_subst_array - unobserv_aa_sub_rate
    observed_perc_aa_ident_array = 1- obs_aa_sub_rate

    sys.stdout.write("\nunobserv_aa_sub_rate {}".format(unobserv_aa_sub_rate[::70]))
    sys.stdout.write("\nobs_aa_sub_rate {}".format(obs_aa_sub_rate[::70]))
    sys.stdout.write("\nobserved_perc_aa_ident_array {}".format(observed_perc_aa_ident_array[::70]))

    ########################################################################################
    #                                                                                      #
    #                   calc. normalisation factors (old and new methods)                  #
    #                                                                                      #
    ########################################################################################

    # OLD NORM FACTOR
    vfunc_old = np.vectorize(korbinian.cons_ratio.norm.OLD_calc_AAIMON_aa_prop_norm_factor)
    norm_factor_array_old = vfunc_old(observed_perc_aa_ident_array, rand_TM, rand_nonTM)

    # NEW NORM FACTOR
    vfunc = np.vectorize(korbinian.cons_ratio.norm.calc_AAIMON_aa_prop_norm_factor)
    norm_factor_array = vfunc(observed_perc_aa_ident_array, rand_TM, rand_nonTM, proportion_seq_TM_residues)


    ########################################################################################
    #                                                                                      #
    #      convert the array nested_list_AAIMONs_10_replicates to a single list            #
    #       orig array:
    #             [1 2 . ],[1 2 . ],[1 2 .  and so on
    #       new single list:
    #               [1 2 . . . . 1 2 . . . . 1 2 .]
    #       this saves plotting the scatter points 10 times, for each sublist
    #                                                                                      #
    ########################################################################################

    AAIMON_list_10rep = []
    obs_aa_sub_rate_10rep = []
    norm_factor_array_10rep_new = []
    norm_factor_array_10rep_old = []
    # iterate through each sublist [1 2 . ],[1 2 . ],[
    for AAIMON_sublist in nested_list_AAIMONs_10_replicates:
        # add datapoints  to end of growing list
        AAIMON_list_10rep += AAIMON_sublist
        obs_aa_sub_rate_10rep = np.append(obs_aa_sub_rate_10rep, obs_aa_sub_rate)
        norm_factor_array_10rep_old = np.append(norm_factor_array_10rep_old, norm_factor_array_old)
        norm_factor_array_10rep_new = np.append(norm_factor_array_10rep_new, norm_factor_array)

    # divide all AAIMONs by the normalisation factor
    norm_AAIMONs_old = AAIMON_list_10rep / norm_factor_array_10rep_old
    norm_AAIMONs_new = AAIMON_list_10rep / norm_factor_array_10rep_new

    # plot the data after normalisation first (not as important as the orig at the beginning)
    ax.scatter(obs_aa_sub_rate_10rep*100, norm_AAIMONs_new, color=color_norm, s=s, alpha=alpha, marker=marker, linewidths=linewidths, label="after normalisation")
    # plot the original data now, so it is clearly seen
    ax.scatter(obs_aa_sub_rate_10rep*100, AAIMON_list_10rep, color=color_nonnorm, s=s, alpha=alpha, marker=marker, linewidths=linewidths, label="before normalisation")

    # plot the AAIMON normalisation factor
    # this is the AAIMON expected for completely random sequences, which differ only in their AA propensity
    # plot only the first replicate (0: 1600)
    norm_curve_x = obs_aa_sub_rate_10rep[:max_num_positions_mutated] * 100
    norm_curve_y = norm_factor_array_10rep_new[:max_num_positions_mutated]
    ax.plot(norm_curve_x, norm_curve_y, color=color_norm_factor_line, linewidth=2, label="calculated normalisation factor")


    ########################################################################################
    #                                                                                      #
    #                        STUFF THAT IS NOT PLOTTED ANYMORE                             #
    #                                                                                      #
    ########################################################################################

    # SCATTER PLOTS SHOWING RESULTS FROM OLD NORMALISATION METHOD
    #ax.scatter(obs_aa_sub_rate_10rep*100, norm_AAIMONs_old, color="r", s=3, alpha=1, label="after norm., old method")
    #ax.scatter(obs_aa_sub_rate_10rep*100, norm_factor_array_10rep_old, s=1, color="r")

    # POLYNOMIAL FIT TO THE ORIGINAL DATA (NOT SHOWN)
    z = np.polyfit(obs_aa_sub_rate_10rep*100, AAIMON_list_10rep, 13)
    p = np.poly1d(z)
    fitted_values_y = p(obs_aa_sub_rate*100)
    #ax.plot(obs_aa_sub_rate*100, fitted_values_y, color="0.5", label="polyfit", linewidth=2)

    # POLYNOMIAL FIT TO THE NEW NORM DATA (NOT SHOWN)
    z = np.polyfit(obs_aa_sub_rate_10rep*100, norm_AAIMONs_new, 13)
    p = np.poly1d(z)
    fitted_values_y = p(obs_aa_sub_rate*100)
    #ax.plot(obs_aa_sub_rate*100, fitted_values_y, color=colour_dict["TUM_colours"]['TUM3'], label="polynomial fit", linewidth=2)

    # POLYNOMIAL FIT TO THE OLD NORM DATA (NOT SHOWN)
    z = np.polyfit(obs_aa_sub_rate_10rep*100, norm_AAIMONs_old, 13)
    p = np.poly1d(z)
    fitted_values_y = p(obs_aa_sub_rate*100)
    #ax.plot(obs_aa_sub_rate*100, fitted_values_y, color="r", label="polyfit to old norm. data", linewidth=2)

    ########################################################################################
    #                                                                                      #
    #                                 AXIS LABELS, ETC                                     #
    #                                                                                      #
    ########################################################################################
    ax.set_ylim(0.9, 1.3)
    ax.set_ylabel("AAIMON", fontsize=fontsize+2)
    ax.set_xlabel("% observed AA substitutions in full protein", fontsize=fontsize+2)
    ax.set_xlim(0, 75)

    # re-order the legend
    handles, labels = ax.get_legend_handles_labels()
    handles = [handles[2], handles[0], handles[1]]
    labels = [labels[2], labels[0], labels[1]]
    legend = ax.legend(handles, labels, loc='upper left', frameon=True, scatterpoints=25, fontsize=fontsize)
    #legend.legendHandles[0]._sizes = [20]
    #legend.legendHandles[2]._sizes = [20]
    ax.tick_params(labelsize=fontsize)
    plt.tight_layout()

    # save fig, e.g. "D:\Databases\summaries\01\List01_rand\List01_rand_AAIMON_vs_aa_sub.png"
    fig.savefig(fig_AAIMON_vs_perc_aa_sub_png, dpi=400)
    fig.savefig(fig_AAIMON_vs_perc_aa_sub_pdf)

    # check that the final normalised values are close to 1.0
    old_last = list(norm_AAIMONs_old[-100:]) + list(norm_AAIMONs_old[900:1000]) + list(norm_AAIMONs_old[1900:2000])
    new_last = list(norm_AAIMONs_new[-100:]) + list(norm_AAIMONs_new[900:1000]) + list(norm_AAIMONs_new[1900:2000])
    sys.stdout.write("\nmean, old norm method {}".format(np.array(old_last).mean()))
    sys.stdout.write("\nmean, new norm method {}".format(np.array(new_last).mean()))

