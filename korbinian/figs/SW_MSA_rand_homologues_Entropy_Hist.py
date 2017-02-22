import numpy as np
import random
import pandas as pd
import korbinian

##########parameters#############

seq_len = 10000
number_seq = 50
number_mutations = 2000
subset_num = 12
ident = 100 * (seq_len - number_mutations) / seq_len

List_rand_TM = r"D:\Databases\summaries\01\List01_rand\List01_rand_TM.csv"
#List_rand_TM = r"D:\Databases\summaries\03\List03_rand\List03_rand_TM.csv"

aa_prop_ser = pd.Series.from_csv(List_rand_TM, sep="\t")
# the first line is the random identity. Extract and delete.
randTM = aa_prop_ser["random_sequence_identity_output"]
del aa_prop_ser["random_sequence_identity_output"]
aa_prop_ser.index = aa_prop_ser.index.str[:1]

########################################################################################
#                                                                                      #
#                 Create a matrix of artificial homologues                             #
#                                                                                      #
########################################################################################

orig_seq, matrix = korbinian.MSA_normalisation.create_matrix_artificial_homologues(aa_prop_ser, seq_len, number_seq, number_mutations)
print(len(matrix), len(matrix[0]))

########################################################################################
#                                                                                      #
#              DO NOT USE. Analysis by AA doesn't work well here                       #
#                                                                                      #
########################################################################################
analyse_by_aa = False
if analyse_by_aa:
    print(orig_seq)
    all_aa = "ACDEFGHIKLMNPQRSTVWY"
    print(len(all_aa))
    perc_ident_ser = pd.Series()

    for aa in all_aa:
        orig_seq_array = np.array(list(orig_seq))
        if aa in orig_seq_array:
            # get indices in orig seq
            list_indices_in_orig_seq_with_this_aa = np.where(orig_seq_array == aa)[0]
            #print(aa)
            #print(list_indices_in_orig_seq_with_this_aa)
            # get those positions in homologues
            rows_in_matrix_with_this_aa = matrix[:,list_indices_in_orig_seq_with_this_aa]

            # # convert to a 1D array of residues
            # rows_in_matrix_with_this_aa_flattened = rows_in_matrix_with_this_aa.flatten()
            # # get value counts
            # values, counts = np.unique(rows_in_matrix_with_this_aa_flattened, return_counts=True)
            # ser_counts = pd.Series(counts, index=values)
            # n_aa_with_orig = ser_counts[aa]
            # perc_ident = n_aa_with_orig / len(rows_in_matrix_with_this_aa_flattened)
            # #print(perc_ident)
            # perc_ident_ser[aa] = perc_ident

            rows_in_matrix_with_this_aa

        else:
            # if this aa is not in the orig seq (e.g. Cysteine), leave value empty
            perc_ident_ser[aa] = np.nan

    print(perc_ident_ser)

########################################################################################
#                                                                                      #
#             Take each column (for pos 1, e.g. Leu) and make a string of all          #
#              aa in that column (e.g. LLLILLLCLLLGLLLLW)                              #
#              By counting the aa in that string, it is possible to get the            #
#               % identity at that position (e.g. 0.89, 89%)                           #
#                                                                                      #
########################################################################################

matrix_df = pd.DataFrame(matrix)

cons_list = []

for pos in matrix_df.columns:
    string_seq_at_that_pos = "".join(matrix_df[pos])
    cons = korbinian.MSA_normalisation.count_aa_freq(string_seq_at_that_pos)
#    print(cons)
    cons_list.append(cons)

########################################################################################
#                                                                                      #
#             Join all these AA percentages into one big dataframe                     #
#            index = each AA in orig_seq (L,W,F,I,L etc)                               #
#            columns = each AA in aa_list (A,C,D,E,F etc)                              #
#   values =  % of the aa in that position (column) in artificial alignment (e.g. 0.89)#
#                                                                                      #
########################################################################################
cons_df = pd.concat(cons_list)
len_df = len(cons_df.columns)
cons_df.index = list(orig_seq)

aa_prop_list = []

########################################################################################
#                                                                                      #
#             Go through each aa separately, and collect the positions where the       #
#           orig AA is L, and the percentage of L in final AA in "homologues"          #
#           since index = orig_aa and columns = ACDE etc, this is as simple as         #
#           df.loc["L", "L"] for all the leucines                                      #
#           if the result is a list, add all to list of % identities, otherwise        #
#           add the single value for that AA at that position                          #
#                                                                                      #
########################################################################################
for aa in cons_df.index.unique():
#    print(cons_df2.loc[aa, aa])
#    print(type(cons_df2.loc[aa, aa]))
    if type(cons_df[aa][aa]) == pd.Series:
        observed_cons_1_row = list(cons_df.loc[aa, aa])
        aa_prop_list += observed_cons_1_row
    else:
        aa_prop_list.append(cons_df.loc[aa, aa])
#print(aa_prop_list)

########################################################################################
#                                                                                      #
#             Calc overall cons rate, and back mutation rate                           #
#                                                                                      #
########################################################################################
observed_cons_rate_all_pos = np.array(aa_prop_list).mean() * 100
back_mutation_rate = 100 * (observed_cons_rate_all_pos - ident)/(100 - ident)

#print(len(aa_prop_list))

print("observed_cons_rate_all_pos", observed_cons_rate_all_pos)
print("back_mutation_rate", back_mutation_rate)

########################################################################################
#                                                                                      #
#             Convert perc_identity to Entropy for each AA [(f(AA) * log2(f(AA))]      #
#                                                                                      #
########################################################################################
for aa in cons_df.columns[:len_df]:
    cons_df['%s_frac' %aa] = cons_df['%s' %aa]* np.log2(cons_df['%s' %aa])
# print(cons_df.columns)

frac_df = cons_df.iloc[:,len_df:]
frac_df['entropy'] = - frac_df.sum(axis=1)

print('mean entropy =',frac_df['entropy'].mean(axis=0))


import matplotlib.pyplot as plt
#%matplotlib inline

# the histogram of the data
n, bins, patches = plt.hist(frac_df['entropy'],12, normed=1, facecolor='lightblue', alpha=0.9)

plt.xlabel('entropy')
plt.ylabel('freq')
plt.title(r'%d aa identity for TM region' %ident)
plt.annotate(r'mean entropy %.3f' %frac_df['entropy'].mean(axis=0),xy =(0.65,0.9),xycoords='axes fraction')
plt.annotate(r'seq len %s'%seq_len, xy =(0.05,0.9),xycoords='axes fraction')
plt.annotate(r'number of seq %s'%number_seq, xy =(0.05,0.85),xycoords='axes fraction')
plt.annotate(r'mean observed identity %.2f'%observed_cons_rate_all_pos, xy =(0.05,0.80),xycoords='axes fraction')
plt.annotate(r'random aa identity %.2f'%back_mutation_rate, xy =(0.05,0.75),xycoords='axes fraction')
#plt.grid(True)

figpath = List_rand_TM[:-4] + "FigX18b_entropy.png"
plt.savefig(figpath)