import numpy as np
import random
import pandas as pd
import korbinian

##########parameters#############

seq_len = 100
number_seq = 50
number_mutations = 20
subset_num = 12
ident = 100 * (seq_len - number_mutations) / seq_len

List01_rand_TM = r"D:\Databases\summaries\01\List01_rand\List01_rand_TM.csv"

aa_prop_ser = pd.Series.from_csv(List01_rand_TM, sep="\t")
# the first line is the random identity. Extract and delete.
randTM = aa_prop_ser["random_sequence_identity_output"]
del aa_prop_ser["random_sequence_identity_output"]
aa_prop_ser.index = aa_prop_ser.index.str[:1]


orig_seq, matrix = korbinian.MSA_normalisation.create_matrix_artificial_homologues(aa_prop_ser, seq_len, number_seq, number_mutations)
print(len(matrix), len(matrix[0]))

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

# aa = np.array(list(orig_seq))
# bb = np.where(aa=="U")
# print(bb)

