import pandas as pd
import numpy as np
import pickle
import matplotlib.pyplot as plt
import sys

#excel_settings_file = r"D:\Dropbox\korbinian\dropbox_korbinian_settings\korbinian_run_settings_schweris_mp.xlsx"

# extract (FROM SETTINGS) the number of amino acids before and after the TMD
n_aa_before_tmd = 30
n_aa_after_tmd = 30
n_aa_excluded_from_edge = 2
# define number of bins in TM region (~10 for betabarrel, ~21 for single-pass)
n_bins_in_TM_region = 8

# Betabarrel, 10 res each side, hessa30, gap2 in TMD, gap1-6 in TMD_plus_surr
gap_all_pos_path = r"D:\Databases\summaries\03\20161104_BB_fastagap_60-80\List03_gap_all_pos.pickle"
# singlepass, 30 res each side, hessa30, gap2 in TMD, gap1-6 in TMD_plus_surr
gap_all_pos_path = r"D:\Databases\summaries\01\excl_SiPe\20161107_40-100_Hessa30_gap2\List01_gap_all_pos.pickle"
# multipass, 30 res each side, hessa30, gap2 in TMD, gap1-6 in TMD_plus_surr
gap_all_pos_path = r"D:\Databases\summaries\02\excl_SiPe\20161106_MP_excl_SiPe_40-100_Hessa30_gap2\List02_gap_all_pos.pickle"
# Betabarrel, 30 res each side, hessa30, gap2 in TMD, gap1-6 in TMD_plus_surr
gap_all_pos_path = r"D:\Databases\summaries\03\20161107_40-100_Hessa30_gap2\List03_gap_all_pos.pickle"

#gap_all_pos_path = r"D:\Databases\summaries\01\20161104_fastagap_60-80\List01_gap_all_pos.pickle"
#gap_all_pos_path = r"D:\Databases\summaries\01\20161104_fastagap_40-100\List01_gap_all_pos.pickle"

########################################################################################
#                                                                                      #
#                       Load a massive list of gap positions                           #
#                        pos -10 = 10 res N-term of TMD                                #
#                        pos 0 = first res of TMD                                      #
#                        pos 1 = lost res of TMD                                       #
#                        pos 10 = 10th res C-term of TMD                               #
#                                                                                      #
########################################################################################

# open pickle, giving a python list
with open(gap_all_pos_path, "rb") as pkl:
    pos_with_gaps_for_all_TMDs_all_proteins = pickle.load(pkl)

pos_arr = np.array(pos_with_gaps_for_all_TMDs_all_proteins)
pos_arr = pos_arr[pos_arr > -200]
sys.stdout.write("number of gap positions in list", len(pos_with_gaps_for_all_TMDs_all_proteins))
sys.stdout.write("array shape after excluding infinite values", pos_arr.shape)

pd.options.display.float_format = '{:,.2f}'.format
# define the number of aa surrounding, to be used in bin creation
n_aa_surr_left = 37
n_aa_surr_right = n_aa_surr_left + 1
# define bins left of TMD
bins_left = list(range(-n_aa_surr_left,0))
# use linspace to make x number of bins for the TMD, evenly spaced between 0 and 1
bins_TM = list(np.linspace(-0.01,1.01,n_bins_in_TM_region))
# define bins to the right of the TMD
bins_right = list(range(2,n_aa_surr_right))
# combine left, TM and right bins together
newbins = bins_left + bins_TM + bins_right
# use the lower value of each bin as the label
newbins_lower_bin_edge = newbins[:-1]
# create a new dataframe to hold frequencies and bins
newbins_df = pd.DataFrame(newbins_lower_bin_edge, columns=["lower"])
# use numpy histogram to assign data into bins
freq, bins = np.histogram(pos_arr, bins=newbins)
# add frequency counts to dataframe
newbins_df["freq"] = freq
"""
The edges have artificially high number of gaps due to the alignment settings!
The settings need to be adjusted to lower gap penalties,
and the edge 2 residues need to be excluded from the dataset
"""
# the number of amino acids on the left, excluding the edges is the n_aa_before_tmd- n_excluded
aa_left_excl_edge = n_aa_before_tmd - n_aa_excluded_from_edge
# the number of amino acids on the right, excluding the edges is the n_aa_after_tmd - n_excluded
aa_right_excl_edge = n_aa_after_tmd - n_aa_excluded_from_edge
# the position at which the TMD starts should be the number of bins left of the TMD
start_TMD_in_index = n_aa_surr_left
# if there 8 lower+upper values gives 7 bins, so need to minus 1 to match
end_TMD_in_index = start_TMD_in_index + n_bins_in_TM_region - 1
# define where the data excluding the edges starts
start_data_excl_edge = start_TMD_in_index - aa_left_excl_edge
# define where the data excluding the edges ends
end_data_excl_edge = end_TMD_in_index + aa_right_excl_edge

########################################################################################
#                                                                                      #
#                              Create labels for x-axis                                #
#                                                                                      #
########################################################################################
# labels for the N-term remain unchanged from "pos"
newbins_df.loc[:start_TMD_in_index - 1, "pos"] = [str("%i"%x) for x in newbins_df.loc[:start_TMD_in_index - 1, "lower"]]
# labels for the TMD are simply "T"
newbins_df.loc[start_TMD_in_index:end_TMD_in_index, "pos"] = "T"
# labels for the C-term remain unchanged from "pos"
newbins_df.loc[end_TMD_in_index:, "pos"] = newbins_df.loc[end_TMD_in_index:, "lower"].astype(int)
# the frequency excluding the edge datapoints
newbins_df["freq_excl_edge"] = newbins_df[start_data_excl_edge:end_data_excl_edge].freq
# mean frequency is used to draw lines indicating the N-term, TMD and C-term regions
mean_freq = newbins_df.freq_excl_edge.mean()
# make line for N-term slightly higher, and line for Cterm slightly lower
newbins_df.loc[:start_TMD_in_index, "Nterm"] = mean_freq * 1.1
newbins_df.loc[end_TMD_in_index:, "Cterm"] = mean_freq * 0.9
newbins_df.loc[start_TMD_in_index : end_TMD_in_index, "TM"] = mean_freq


"""
Have a look at your dataframe, to see how it works
 	lower 	freq 	pos 	freq_ecl_edge
29 	-8.00 	14 	-8 	14.00 	26.49
30 	-7.00 	11 	-7 	11.00 	26.49
31 	-6.00 	18 	-6 	18.00 	26.49
32 	-5.00 	25 	-5 	25.00 	26.49
33 	-4.00 	20 	-4 	20.00 	26.49
34 	-3.00 	18 	-3 	18.00 	26.49
35 	-2.00 	17 	-2 	17.00 	26.49
36 	-1.00 	15 	-1 	15.00 	26.49
37 	-0.01 	4 	T 	4.00 	26.49
38 	0.04 	3 	T 	3.00 	nan
39 	0.09 	6 	T 	6.00 	nan
40 	0.14 	2 	T 	2.00 	nan
41 	0.18 	4 	T 	4.00 	nan
42 	0.23 	7 	T 	7.00 	nan
43 	0.28 	3 	T 	3.00 	nan
44 	0.33 	5 	T 	5.00 	nan
45 	0.38 	6 	T 	6.00 	nan
46 	0.43 	3 	T 	3.00 	nan
47 	0.48 	9 	T 	9.00 	nan
48 	0.52 	2 	T 	2.00 	nan
49 	0.57 	3 	T 	3.00 	nan
50 	0.62 	3 	T 	3.00 	nan
51 	0.67 	2 	T 	2.00 	nan
52 	0.72 	4 	T 	4.00 	nan
53 	0.77 	2 	T 	2.00 	nan
54 	0.82 	2 	T 	2.00 	nan
55 	0.86 	3 	T 	3.00 	nan
56 	0.91 	2 	T 	2.00 	nan
57 	0.96 	5 	T 	5.00 	nan
58 	1.01 	0 	1 	0.00 	nan
59 	2.00 	11 	2 	11.00 	nan
60 	3.00 	9 	3 	9.00 	nan
61 	4.00 	13 	4 	13.00 	nan
62 	5.00 	20 	5 	20.00 	nan
63 	6.00 	23 	6 	23.00 	nan
64 	7.00 	13 	7 	13.00 	nan
"""
# create a view with only the plotted data, excluding the edge positions
newbins_df_excl = newbins_df[newbins_df["freq_excl_edge"].notnull()]

# change plot defaults
plt.style.use('ggplot')
plt.rcParams["savefig.dpi"] = 240
plt.rcParams.update({'font.size': 8})

fig, ax = plt.subplots()
# plot frequency as line chart. The x-positions are the LOWER value for each bin.
ax.plot(newbins_df_excl.index, newbins_df_excl.freq)
# plot some lines to indicate regions
ax.plot(newbins_df_excl.index, newbins_df_excl.TM)
ax.plot(newbins_df_excl.index, newbins_df_excl.Nterm)
ax.plot(newbins_df_excl.index, newbins_df_excl.Cterm)
ax.legend()
ax.set_ylim(0)
ax.set_xlabel("position")
ax.set_ylabel("frequency")

fig_out = gap_all_pos_path[:-7] + "FigX07_fastagap_hist.png"
fig.savefig(fig_out)
