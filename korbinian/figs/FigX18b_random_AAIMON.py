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

seq_len = 100
max_num_positions_mutated = 60
n_mutations_list = np.arange(0, max_num_positions_mutated)
perc_aa_subst_list = n_mutations_list / seq_len
# List_rand_TM = r"D:\Databases\summaries\01\List01_rand\List01_rand_TM.csv"
List_rand_TM = os.path.join(data_dir, "summaries\{ln:02d}\List{ln:02d}_rand\List{ln:02d}_rand_TM.csv".format(ln=list_number))
List_rand_nonTM = os.path.join(data_dir, "summaries\{ln:02d}\List{ln:02d}_rand\List{ln:02d}_rand_nonTM.csv".format(ln=list_number))
pickle_with_artificial_AAIMONs = os.path.join(data_dir, "summaries\{ln:02d}\List{ln:02d}_rand\List{ln:02d}_rand_AAIMONs.pickle".format(ln=list_number))
fig_AAIMON_vs_perc_aa_sub = os.path.join(data_dir, "summaries\{ln:02d}\List{ln:02d}_rand\List{ln:02d}_rand_AAIMON_vs_aa_sub.png".format(ln=list_number))

if not os.path.isfile(pickle_with_artificial_AAIMONs) or repeat_randomisation == True:
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


    ########################################################################################
    #                                                                                      #
    #              Create pairwise alignment of artificial homologue and calc % ident      #
    #                                                                                      #
    ########################################################################################

    nested_list_AAIMONs_3_replicates = []

    for i in range(3):
        AAIMON_list = []
        for n, number_mutations in enumerate(n_mutations_list):
            perc_ident_TM = korbinian.cons_ratio.norm.get_perc_ident_random_pairwise_alignment(aa_prop_TM, seq_len, number_mutations)
            perc_ident_nonTM = korbinian.cons_ratio.norm.get_perc_ident_random_pairwise_alignment(aa_prop_nonTM, seq_len, number_mutations)
            AAIMON = perc_ident_TM / perc_ident_nonTM
            AAIMON_list.append(AAIMON)
            perc_aa_subst = number_mutations / seq_len
            sys.stdout.write(".")
            if n % 60 == 0:
                sys.stdout.write("\n{}, {:.03f}, {:.03f}, {:.03f}, {:.03f}\n".format(number_mutations, perc_aa_subst, perc_ident_TM, perc_ident_nonTM, AAIMON))
            sys.stdout.flush()
        nested_list_AAIMONs_3_replicates.append(AAIMON_list)


    with open(pickle_with_artificial_AAIMONs, "wb") as pkl_file:
        pickle.dump(nested_list_AAIMONs_3_replicates, pkl_file, -1)

with open(pickle_with_artificial_AAIMONs, "rb") as pkl_file:
    nested_list_AAIMONs_3_replicates = pickle.load(pkl_file)

fig, ax = plt.subplots()
color_nonnorm = "#EE762C"
color_norm = "#0076B8"
color_norm_line = "#53A7D5"

for AAIMON_list in nested_list_AAIMONs_3_replicates:
    ax.scatter(perc_aa_subst_list, AAIMON_list, color=color_nonnorm, s=5)

ax.set_ylabel("AAIMON")
ax.set_xlabel("% AA substitutions in full protein")
ax.set_xlim(0)

fig.savefig(fig_AAIMON_vs_perc_aa_sub)