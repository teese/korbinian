# import csv
# import korbinian.mtutils as utils
# import pandas as pd
# import numpy as np
#
# omp_nr_fasta = r"D:\Databases\OMPdb\OMPdb.30"
# # OMPdb_summary_nr_csv = omp_nr_fasta + "_flatfiles_nr.csv"
# # OMPdb_summary_csv_with_TM_seqs = omp_nr_fasta + "_flatfiles_nr_with_seqs.csv"
#
# OMPdb_summary_nr_csv = "D:\Databases\OMPdb\OMPdb_Selected_by_potential_IDs.csv"
# OMPdb_summary_csv_with_TM_seqs = OMPdb_summary_nr_csv[:-4] + "_with_seqs.csv"
#
#
#
# def get_omp_TM_indices_and_slice_from_summary_table(OMPdb_summary_nr_csv, OMPdb_summary_csv_with_TM_seqs):
#     df_KW = pd.read_csv(OMPdb_summary_nr_csv)
#
#     # get sequence length
#     df_KW["seqlen"] = df_KW["Sequence"].str.len()
#
#     # Function which returns list of all M-indices
#     def getting_membrane_indices(Topo_data):
#         m_list = [i for i, topology in enumerate(Topo_data) if topology == "M"]  #find(Topo_data)
#         return m_list
#
#     # Creating new column M_indices, which contains the indices of Ms
#     df_KW["M_indices"] = df_KW.Topology.apply(getting_membrane_indices)
#
#     # Converting empty entries to NaN
#     df_KW["M_indices"] = df_KW.M_indices.apply(lambda x: np.nan if x ==[] else x)
#
#     num_proteins_BEFORE_dropping_those_without_mem_indices = df_KW.shape[0]
#
#     # Extracting entries to a new Dataframe
#     df_KW = df_KW[df_KW.M_indices.notnull()]
#
#     num_proteins_AFTER_dropping_those_without_mem_indices = df_KW.shape[0]
#
#     # Filter, cutting of Coverages under 85% & Creating new Index
#     df_KW = df_KW.loc[df_KW["Coverage(%)"] >= 85]
#
#     # df_KW.index = range(1,len(df_KW["Uniprot"])+1)
#     num_proteins_AFTER_dropping_those_with_coverage_below_85 = df_KW.shape[0]
#
#     # Creating new list (nested list)
#     nested_list_of_membrane_borders = []
#
#     # Function which filters out start and end-points
#     def check_for_border(m_indices):
#         for n in range(0, len(m_indices) - 1):
#             if m_indices[n] + 1 != m_indices[n + 1]:
#                 m_borders.append(m_indices[n])
#                 m_borders.append(m_indices[n + 1])
#
#     # Filling nest with lists of start and end-points
#     for n in df_KW.M_indices:
#         m_borders = []
#         m_borders.append(n[0])
#         check_for_border(n)
#         m_borders.append(n[-1])
#         nested_list_of_membrane_borders.append(m_borders)
#
#     array_membrane_borders = np.array(nested_list_of_membrane_borders)
#     for subarray in array_membrane_borders:
#         #print(subarray[::2] = subarray[::2]*10)
#         subarray = np.array(subarray)
#         subarray[1::2] = subarray[1::2] + 1
#     nested_list_of_membrane_borders_python_indexstyle = array_membrane_borders.tolist()
#
#     # Creating new column, which contains start and end-points
#     df_KW["Membrane_Borders"] = nested_list_of_membrane_borders_python_indexstyle
#
#     # Creating new column, which contains the Amoung of TMDS
#     df_KW["number_of_TMDs"] = df_KW.Membrane_Borders.apply(lambda x: len(x) / 2)
#
#     # Filter, filters out, if less than 8 or more than 24 TMDs
#     df_KW["number_of_TMDs"] = df_KW["number_of_TMDs"].apply(lambda x: int(x) if 8 <= x <= 24 else np.nan)
#
#     num_proteins_AFTER_dropping_those_without_TMs_between_8_and_24 = df_KW.shape[0]
#
#     # Creating new dataframe without nan
#     df_KW = df_KW[df_KW["number_of_TMDs"].notnull()]
#
#     df_KW["TM_indices"] = df_KW["Membrane_Borders"].apply(lambda x : tuple(zip(x[::2], x[1::2])))
#
#     # create a list of [TM01, TM02, TM03, etc.
#     long_list_of_TMDs = []
#     for i in range(1,50):
#         long_list_of_TMDs.append("TM{:02d}".format(i))
#
#     # for the .set_value function, set dtype as object
#     df_KW["list_of_TMDs"] = ""
#     df_KW["list_of_TMDs"].astype(object)
#
#     for row in df_KW.index:
#         # get nested tuple of TMDs
#         nested_tup_TMs = df_KW.loc[row, "TM_indices"]
#         # slice long list of TMD names to get an appropriate list for that protein [TM01, TM02, TM03, etc.
#         len_nested_tup_TMs = len(nested_tup_TMs)
#         list_of_TMDs = long_list_of_TMDs[:len_nested_tup_TMs]
#         # add that list to the dataframe (could also be added as a stringlist, but that's irritating somehow)
#         df_KW.set_value(row, "list_of_TMDs", list_of_TMDs)
#         # set seq for slicing
#         full_seq = df_KW.loc[row, "Sequence"]
#         #topology = df_KW.loc[row, "Topology"]
#         # iterate through all the TMDs of that protein, slicing out the sequences
#         for i in range(len(list_of_TMDs)):
#             if i % 10 == 0:
#                 print(".",end="")
#             TMD = list_of_TMDs[i]
#             tup = nested_tup_TMs[i]
#             df_KW.loc[row, TMD + "_start"] = tup[0]
#             df_KW.loc[row, TMD + "_end"] = tup[1]
#             df_KW.loc[row, TMD + "_seq"] = utils.slice_with_listlike(full_seq, tup)
#             #df_KW.loc[row, TMD + "_top"] = utils.slice_with_listlike(topology, tup)
#
#     max_num_TMDs = df_KW["number_of_TMDs"].max()
#
#     # fa_aa_before_tmd = set_["fa_aa_before_tmd"]
#     # fa_aa_after_tmd = set_["fa_aa_after_tmd"]
#     fa_aa_before_tmd = 10
#     fa_aa_after_tmd = 10
#
#     # currently the loop is run for each TMD, based on the sequence with the most TMDs
#     for i in range(1, int(max_num_TMDs) + 1):
#         TMD = 'TM%02d' % i
#         df_KW = utils.get_indices_TMD_plus_surr_for_summary_file(df_KW, TMD, fa_aa_before_tmd, fa_aa_after_tmd)
#
#     df_KW.to_csv(OMPdb_summary_csv_with_TM_seqs, sep=",", quoting=csv.QUOTE_NONNUMERIC)