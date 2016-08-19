import csv
import korbinian.mtutils as utils
import pandas as pd
import numpy as np

omp_nr_fasta = r"D:\Databases\OMPdb\OMPdb.30"
# OMPdb_summary_nr_csv = omp_nr_fasta + "_flatfiles_nr.csv"
# OMPdb_summary_csv_with_TM_seqs = omp_nr_fasta + "_flatfiles_nr_with_seqs.csv"


OMPdb_summary_nr_csv = "D:\Databases\OMPdb\OMPdb_Selected_by_potential_IDs.csv"
OMPdb_summary_csv_with_TM_seqs = OMPdb_summary_nr_csv[:-4] + "_with_seqs.csv"

df_KW = pd.read_csv(OMPdb_summary_nr_csv)

# get sequence length
df_KW["seqlen"] = df_KW["Sequence"].str.len()

# Function which returns list of all M-indices
def getting_membrane_indices(Topo_data):
    m_list = [i for i, topology in enumerate(Topo_data) if topology == "M"]  #find(Topo_data)
    return m_list

# Creating new column M_indices, which contains the indices of Ms
df_KW["M_indices"] = df_KW.Topology.apply(getting_membrane_indices)

# Converting empty entries to NaN
df_KW["M_indices"] = df_KW.M_indices.apply(lambda x: np.nan if x ==[] else x)

num_proteins_BEFORE_dropping_those_without_mem_indices = df_KW.shape[0]

# Extracting entries to a new Dataframe
df_KW = df_KW[df_KW.M_indices.notnull()]

num_proteins_AFTER_dropping_those_without_mem_indices = df_KW.shape[0]

# Filter, cutting of Coverages under 85% & Creating new Index
df_KW = df_KW.loc[df_KW["Coverage(%)"] >= 85]

# df_KW.index = range(1,len(df_KW["Uniprot"])+1)
num_proteins_AFTER_dropping_those_with_coverage_below_85 = df_KW.shape[0]

# Creating new list (nested list)
nested_list_of_membrane_borders = []

# Function which filters out start and end-points
def check_for_border(m_indices):
    for n in range(0, len(m_indices) - 1):
        if m_indices[n] + 1 != m_indices[n + 1]:
            m_borders.append(m_indices[n])
            m_borders.append(m_indices[n + 1])

# Filling nest with lists of start and end-points
for n in df_KW.M_indices:
    m_borders = []
    m_borders.append(n[0])
    check_for_border(n)
    m_borders.append(n[-1])
    nested_list_of_membrane_borders.append(m_borders)

array_membrane_borders = np.array(nested_list_of_membrane_borders)
for subarray in array_membrane_borders:
    #print(subarray[::2] = subarray[::2]*10)
    subarray = np.array(subarray)
    subarray[1::2] = subarray[1::2] + 1
nested_list_of_membrane_borders_python_indexstyle = array_membrane_borders.tolist()

# Creating new column, which contains start and end-points
df_KW["Membrane_Borders"] = nested_list_of_membrane_borders_python_indexstyle

# Creating new column, which contains the Amoung of TMDS
df_KW["number_of_TMDs"] = df_KW.Membrane_Borders.apply(lambda x: len(x) / 2)

# Filter, filters out, if less than 8 or more than 24 TMDs
df_KW["number_of_TMDs"] = df_KW["number_of_TMDs"].apply(lambda x: int(x) if 8 <= x <= 24 else np.nan)

num_proteins_AFTER_dropping_those_without_TMs_between_8_and_24 = df_KW.shape[0]

# Creating new dataframe without nan
df_KW = df_KW[df_KW["number_of_TMDs"].notnull()]

df_KW["TM_indices"] = df_KW["Membrane_Borders"].apply(lambda x : tuple(zip(x[::2], x[1::2])))

# create a list of [TM01, TM02, TM03, etc.
long_list_of_TMDs = []
for i in range(1,50):
    long_list_of_TMDs.append("TM{:02d}".format(i))

# for the .set_value function, set dtype as object
df_KW["list_of_TMDs"] = ""
df_KW["list_of_TMDs"].astype(object)

for row in df_KW.index:
    # get nested tuple of TMDs
    nested_tup_TMs = df_KW.loc[row, "TM_indices"]
    # slice long list of TMD names to get an appropriate list for that protein [TM01, TM02, TM03, etc.
    len_nested_tup_TMs = len(nested_tup_TMs)
    list_of_TMDs = long_list_of_TMDs[:len_nested_tup_TMs]
    # add that list to the dataframe (could also be added as a stringlist, but that's irritating somehow)
    df_KW.set_value(row, "list_of_TMDs", list_of_TMDs)
    # set seq for slicing
    full_seq = df_KW.loc[row, "Sequence"]
    #topology = df_KW.loc[row, "Topology"]
    # iterate through all the TMDs of that protein, slicing out the sequences
    for i in range(len(list_of_TMDs)):
        if i % 10 == 0:
            print(".",end="")
        TMD = list_of_TMDs[i]
        tup = nested_tup_TMs[i]
        df_KW.loc[row, TMD + "_start"] = tup[0]
        df_KW.loc[row, TMD + "_end"] = tup[1]
        df_KW.loc[row, TMD + "_seq"] = utils.slice_with_listlike(full_seq, tup)
        #df_KW.loc[row, TMD + "_top"] = utils.slice_with_listlike(topology, tup)

# nested_list_borders = df_KW["Membrane_Borders"].values.tolist()
# import itertools
# # Function, which converts Start and End-Lists to own Columns and adds None, if no Data
# def border_transpose( x ):
#     tups = itertools.zip_longest( *x )
#     return [ list(t) for t in tups ]
#
# nested_list_borders_T = border_transpose(nested_list_borders)
#
# for n in range(0,24):
#     df_KW["TM%g_start"%(n+1)] = nested_list_borders_T[int(n*2)]
#     df_KW["TM%g_end"%(n+1)]= nested_list_borders_T[int(n*2+1)]
#
# df_KW.iloc[:5,10:]
#
#
# # Adding a column to the Dataframe, with elements containing lists of existing TMDs
#
# list_of_potential_TMDS = ['TM1', 'TM2', 'TM3', 'TM4', 'TM5', 'TM6', 'TM7', 'TM8', 'TM9', 'TM10', 'TM11', 'TM12', 'TM13',
#                           'TM14', 'TM15', 'TM16', 'TM17', 'TM18', 'TM19', 'TM20', 'TM21', 'TM22', 'TM23', 'TM24']
#
# nested_list_of_TMDs = []
#
# # Adds TMD to existing_TMDs if there's a startpoint
# for n in range(1, len(df_KW["Uniprot"]) + 1):
#     existing_TMDs = [potential_TMD for potential_TMD in list_of_potential_TMDS if
#                      np.isnan(df_KW.loc[n, "%s_start" % potential_TMD]) == False]
#     nested_list_of_TMDs.append(existing_TMDs)

# # Transferring nested list to column
# df_KW["List_of_TMDs"] = nested_list_of_TMDs

# df_KW.to_csv("Dataframe with List of TMDs and Membrane Borders.csv")

# # Converting df_KW to strings, in order to save them as excelfile
# df_KW5 = pd.DataFrame(df_KW, dtype=str)
# 
# df_KW5.to_excel("Dataframe with List of TMDs and Membrane Borders.xlsx")

# Sequence of sliced Membranedomains: TM1_seq, TM2_seq...

# nested list, which will contain a list per entry with the sliced sequences
# list_of_sliced_sequences = []
# 
# for m in range (1,len(df_KW["Sequence"])+1):
#     list_of_sliced_sequences_per_entry =[df_KW.loc[m,"Sequence"][int(df_KW.loc[m,"%s_start"%n]):int(df_KW.loc[m,"%s_end" %n])+1] for n in df_KW["List_of_TMDs"][m]]
#     list_of_sliced_sequences.append(list_of_sliced_sequences_per_entry)
# 
# # inverses nested list, new nested list contains a list per TMD-slice
# nested_list_of_slices_per_TMD = border_transpose(list_of_sliced_sequences)
# 
# # Creates new columns "TM1_seq", "TM2_seq" ...
# for n in range(0, 24):
#     df_KW["TM%g_seq" % (n + 1)] = nested_list_of_slices_per_TMD[n]
# 
# df_KW.to_csv("Dataframe with List of TMDs and Membrane Borders and Sliced Sequences.csv")
# 
# # Converting df_KW to strings, in order to save them as excelfile
# df_KW5 = pd.DataFrame(df_KW, dtype=str)
# 
# df_KW5.to_excel("Dataframe with List of TMDs and Membrane Borders and Sliced Sequences.xlsx")

# Creating raw Startpoints (not considering overlapping!), if there is no Value or Value below 0 -> np.nan
# (especially if there is no TMD at all)
# for n in range (1,25):
#     df_KW["start_surrounding_seq_in_query_TM%g" %n] = df_KW["TM%g_start"%n].apply(lambda x: x-10 if x>9 else np.nan)
#     df_KW["end_surrounding_seq_in_query_TM%g" %n] = df_KW["TM%g_end"%n].apply(lambda x: x+10)
#
# # Exchanging NaN with 0 in Sourrounding Start 1, because here, it can be actually 0
# df_KW.start_surrounding_seq_in_query_TM1= df_KW.start_surrounding_seq_in_query_TM1.fillna(value=0)
#
# # Considering the overlapping with TMDs, the sourrounding regions have to be alterd
# for n in range(1, len(df_KW["Uniprot"]) + 1):
#     for m in range(1, 24):
#         if df_KW["end_surrounding_seq_in_query_TM%g" % m][n] >= df_KW["TM%g_start" % (m + 1)][n]:
#             df_KW["end_surrounding_seq_in_query_TM%g" % m][n] = df_KW["TM%g_start" % (m + 1)][n] - 1
#             # else:
#             #   df_KW["end_surrounding_seq_in_query_TM%g"%m][n] = df_KW["end_surrounding_seq_in_query_TM%g"%m][n]
#         if df_KW["start_surrounding_seq_in_query_TM%g" % (m + 1)][n] <= df_KW["TM%g_end" % (m)][n]:
#             df_KW["start_surrounding_seq_in_query_TM%g" % (m + 1)][n] = df_KW["TM%g_end" % (m)][n] + 1
#         # else:
#         #    df_KW["start_surrounding_seq_in_query_TM%g"%(m+1)][n]= df_KW["start_surrounding_seq_in_query_TM%g"%(m+1)][n]
#         if df_KW["end_surrounding_seq_in_query_TM24"][n] > df_KW["len_Sequence"][n]:
#             df_KW["end_surrounding_seq_in_query_TM24"][n] = df_KW["len_Sequence"][n] - 1

max_num_TMDs = df_KW["number_of_TMDs"].max()

# fa_aa_before_tmd = set_["fa_aa_before_tmd"]
# fa_aa_after_tmd = set_["fa_aa_after_tmd"]
fa_aa_before_tmd = 10
fa_aa_after_tmd = 10

# currently the loop is run for each TMD, based on the sequence with the most TMDs
for i in range(1, int(max_num_TMDs) + 1):
    TMD = 'TM%02d' % i
    df_KW = utils.get_indices_TMD_plus_surr_for_summary_file(df_KW, TMD, fa_aa_before_tmd, fa_aa_after_tmd)

df_KW.to_csv(OMPdb_summary_csv_with_TM_seqs, sep=",", quoting=csv.QUOTE_NONNUMERIC)

utils.aaa(df_KW)