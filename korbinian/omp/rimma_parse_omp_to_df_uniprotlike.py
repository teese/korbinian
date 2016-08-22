# import korbinian.mtutils as utils
#
# ############################################################
# #
# #        from Ipython notebook, Getting main DataFrame
# #
# ############################################################
#
# import pandas as pd
# import numpy as np
#
# main_data = pd.read_csv(r"C:\Users\Rima\Documents\IPython Notebooks\save_file.csv")
#
# df_KW = pd.DataFrame(main_data, columns= ["Uniprot","Family","Gene_Name","Organism","NCBI_TaxID","Coverage(%)","Sequence","len_Sequence","Topology_Reli","Topology"])
#
# # Function which returns list of all M-indices
# def getting_membrane_indices(Topo_data):
#     m_list = [i for i, topology in enumerate(Topo_data) if topology == "M"]  #find(Topo_data)
#     return m_list
#
# # Creating new column M_indices, which contains the indices of Ms
# df_KW["M_indices"] = df_KW.Topology.apply(getting_membrane_indices)
#
# # Converting empty entries to NaN
# df_KW.M_indices = df_KW.M_indices.apply(lambda x: np.nan if x ==[] else x)
#
# # Extracting entries to a new Dataframe
# df_KW1 = df_KW[df_KW.M_indices.notnull()]
#
#
# # Filter, cutting of Coverages under 85% & Creating new Index
# df_KW2 = df_KW1.loc[df_KW1["Coverage(%)"]>=85]
# df_KW2.index = range(1,len(df_KW2["Uniprot"])+1)
#
# df_KW.iloc[:5,7:]
#
# # Creating new list (nested list)
# nested_list_of_membrane_borders = []
#
# # Function which filters out start and end-points
# def check_for_border(m_indices):
#     for n in range(0, len(m_indices) - 1):
#         if m_indices[n] + 1 != m_indices[n + 1]:
#             m_borders.append(m_indices[n])
#             m_borders.append(m_indices[n + 1])
#
#
# # Filling nest with lists of start and end-points
# for n in df_KW2.M_indices:
#     m_borders = []
#     m_borders.append(n[0])
#     check_for_border(n)
#     m_borders.append(n[-1])
#     nested_list_of_membrane_borders.append(m_borders)
#
# # Creating new column, which contains start and end-points
# df_KW2["Membrane_Borders"] = nested_list_of_membrane_borders
#
# # Creating new column, which contains the Amoung of TMDS
# df_KW2["TMD_Amount"] = df_KW2.Membrane_Borders.apply(lambda x: len(x) / 2)
#
# # Filter, filters out, if less than 8 or more than 24 TMDs
# df_KW2["TMD_Amount"] = df_KW2["TMD_Amount"].apply(lambda x: int(x) if 8 <= x <= 24 else np.nan)
#
# df_KW2.iloc[:5, 7:]
#
#
#
#
#
# # Creating new dataframe without nan
# df_KW3=df_KW2[df_KW2["TMD_Amount"].notnull()]
#
# df_KW3.index = range(1,len(df_KW3["Uniprot"])+1)
#
# b = df_KW3["Membrane_Borders"].values.tolist()
#
# df_KW3.iloc[:5,7:]
#
#
#
#
#
# import itertools
#
# # Function, which converts Start and End-Lists to own Columns and adds None, if no Data
# def border_transpose( x ):
#     tups = itertools.zip_longest( *x )
#     return [ list(t) for t in tups ]
#
# a = border_transpose(b)
#
# df_KW4 = df_KW3.copy()
#
# for n in range(0,24):
#     df_KW4["TM%g_start"%(n+1)] = a[int(n*2)]
#     df_KW4["TM%g_end"%(n+1)]= a[int(n*2+1)]
#
# df_KW4.iloc[:5,10:]
#
#
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
# for n in range(1, len(df_KW3["Uniprot"]) + 1):
#     existing_TMDs = [potential_TMD for potential_TMD in list_of_potential_TMDS if
#                      np.isnan(df_KW4.loc[n, "%s_start" % potential_TMD]) == False]
#     nested_list_of_TMDs.append(existing_TMDs)
#
# # Transferring nested list to column
# df_KW4["List_of_TMDs"] = nested_list_of_TMDs
#
# df_KW4.to_csv("Dataframe with List of TMDs and Membrane Borders.csv")
#
# # Converting df_KW3 to strings, in order to save them as excelfile
# df_KW5 = pd.DataFrame(df_KW3, dtype=str)
#
# df_KW5.to_excel("Dataframe with List of TMDs and Membrane Borders.xlsx")
#
#
#
#
# # Sequence of sliced Membranedomains: TM1_seq, TM2_seq...
#
# # nested list, which will contain a list per entry with the sliced sequences
# list_of_sliced_sequences = []
#
# for m in range (1,len(df_KW4["Sequence"])+1):
#     list_of_sliced_sequences_per_entry =[df_KW4.loc[m,"Sequence"][int(df_KW4.loc[m,"%s_start"%n]):int(df_KW4.loc[m,"%s_end" %n])+1] for n in df_KW4["List_of_TMDs"][m]]
#     list_of_sliced_sequences.append(list_of_sliced_sequences_per_entry)
#
# # inverses nested list, new nested list contains a list per TMD-slice
# nested_list_of_slices_per_TMD = border_transpose(list_of_sliced_sequences)
#
# # Creates new columns "TM1_seq", "TM2_seq" ...
# for n in range(0, 24):
#     df_KW4["TM%g_seq" % (n + 1)] = nested_list_of_slices_per_TMD[n]
#
# df_KW4.to_csv("Dataframe with List of TMDs and Membrane Borders and Sliced Sequences.csv")
#
# # Converting df_KW3 to strings, in order to save them as excelfile
# df_KW5 = pd.DataFrame(df_KW4, dtype=str)
#
# df_KW5.to_excel("Dataframe with List of TMDs and Membrane Borders and Sliced Sequences.xlsx")
#
#
#
# # Creating raw Startpoints (not considering overlapping!), if there is no Value or Value below 0 -> np.nan
# # (especially if there is no TMD at all)
# for n in range (1,25):
#     df_KW4["start_surrounding_seq_in_query_TM%g" %n] = df_KW4["TM%g_start"%n].apply(lambda x: x-10 if x>9 else np.nan)
#     df_KW4["end_surrounding_seq_in_query_TM%g" %n] = df_KW4["TM%g_end"%n].apply(lambda x: x+10)
#
# # Exchanging NaN with 0 in Sourrounding Start 1, because here, it can be actually 0
# df_KW4.start_surrounding_seq_in_query_TM1= df_KW4.start_surrounding_seq_in_query_TM1.fillna(value=0)
#
# # Considering the overlapping with TMDs, the sourrounding regions have to be alterd
# for n in range(1, len(df_KW4["Uniprot"]) + 1):
#     for m in range(1, 24):
#         if df_KW4["end_surrounding_seq_in_query_TM%g" % m][n] >= df_KW4["TM%g_start" % (m + 1)][n]:
#             df_KW4["end_surrounding_seq_in_query_TM%g" % m][n] = df_KW4["TM%g_start" % (m + 1)][n] - 1
#             # else:
#             #   df_KW4["end_surrounding_seq_in_query_TM%g"%m][n] = df_KW4["end_surrounding_seq_in_query_TM%g"%m][n]
#         if df_KW4["start_surrounding_seq_in_query_TM%g" % (m + 1)][n] <= df_KW4["TM%g_end" % (m)][n]:
#             df_KW4["start_surrounding_seq_in_query_TM%g" % (m + 1)][n] = df_KW4["TM%g_end" % (m)][n] + 1
#         # else:
#         #    df_KW4["start_surrounding_seq_in_query_TM%g"%(m+1)][n]= df_KW4["start_surrounding_seq_in_query_TM%g"%(m+1)][n]
#         if df_KW4["end_surrounding_seq_in_query_TM24"][n] > df_KW4["len_Sequence"][n]:
#             df_KW4["end_surrounding_seq_in_query_TM24"][n] = df_KW4["len_Sequence"][n] - 1
#
# df_KW4.to_csv("Testi.csv")
