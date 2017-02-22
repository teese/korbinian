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

#orig_seq, matrix = korbinian.MSA_normalisation.create_matrix_artificial_homologues(aa_prop_ser, seq_len, number_seq, number_mutations)
#print(len(matrix), len(matrix[0]))

# create the original template sequence
orig_seq = np.array([np.random.choice(aa_prop_ser.index, p=aa_prop_ser) for _ in range(int(seq_len))])

print(orig_seq)

# create indices for each AA in orig sequence
inds = list(range(seq_len))
# choose a random sample of AA to mutation
sam = random.sample(inds, number_mutations)
# convert orig sequence to a list
seq_list = list(orig_seq)
# for each index in the random sample, replace the AA with a random AA
for ind in sam:
    seq_list[ind] = np.random.choice(aa_prop_ser.index, p=aa_prop_ser)
# # join to make a new sequence
# new_seq = "".join(seq_list)
# # append to the matrix of "homologues"

