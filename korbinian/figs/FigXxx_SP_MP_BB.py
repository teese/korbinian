from scipy.stats import ttest_ind
import ast
import csv
import itertools
import korbinian
import korbinian.utils as utils
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.plotly as py
import plotly.graph_objs as go
import plotly

# % matplotlib inline

# plotly.offline.init_notebook_mode()

plt.rcParams["savefig.dpi"] = 240

backgroundcolour = '0.95'
# plt.style.use('ggplot')
# set default font size for plot
fontsize = 12
datapointsize = 8
alpha = 0.1

# mp_smallest_bin = 0.555
# mp_largest_bin = 1.455
# mp_number_of_bins = 31
# mp_final_highest_bin = 3

# linspace_binlist = np.linspace(mp_smallest_bin,
#                               mp_largest_bin,
#                               mp_number_of_bins)

# add 30 as the last bin, to make sure 100% of the data is added to the histogram, including major outliers
# binlist = np.append(linspace_binlist, mp_final_highest_bin)


fig, ax = plt.subplots()

data_to_plot = []

SP_csv = r"C:\Users\coffee oder tee\Dropbox\Undergraduate\Bachelor_Thesis\my_documents\plot_folder\original_csv\SP_MP_BB_40-100_excl_sp\List01_cr_summary.csv"
df_SP = pd.read_csv(SP_csv, sep=",", quoting=csv.QUOTE_NONNUMERIC, index_col=0)

MP_csv_file = r"C:\Users\coffee oder tee\Dropbox\Undergraduate\Bachelor_Thesis\my_documents\plot_folder\original_csv\SP_MP_BB_40-100_excl_sp\List02_cr_summary.csv"
df_MP = pd.read_csv(MP_csv_file, sep=",", quoting=csv.QUOTE_NONNUMERIC, index_col=0)

BB_input_csv = r"C:\Users\coffee oder tee\Dropbox\Undergraduate\Bachelor_Thesis\my_documents\plot_folder\original_csv\SP_MP_BB_40-100_excl_sp\List03_cr_summary.csv"
df_BB = pd.read_csv(BB_input_csv, sep=",", quoting=csv.QUOTE_NONNUMERIC, index_col=0)

#####################################################################
#                                                                   #
#                              BB                                   #
#                                                                   #
#####################################################################

# iterate through the proteins that have a list of TMDs
#for acc in df_BB.loc[df_BB['list_of_TMDs'].notnull()].index:
#    dict_AAIMON_ratio_mean_BB = {}
    #    print(df_BB.loc[acc, 'list_of_TMDs'])
#    for TMD in ast.literal_eval(df_BB.loc[acc, 'list_of_TMDs']):
#        dict_AAIMON_ratio_mean_BB[TMD] = df_BB.loc[acc, '%s_AAIMON_ratio_mean' % TMD]
#    df_BB.loc[acc, 'AAIMON_ratio_mean_all_TMDs'] = np.mean(list(dict_AAIMON_ratio_mean_BB.values()))
#    print(dict_AAIMON_ratio_mean_BB)

TM01_greater_5_homol = df_BB['TM01_AAIMON_n_homol'].loc[df_BB['TM01_AAIMON_n_homol'] > 5]
mean_AAIMON_filt_BB = df_BB['AAIMON_ratio_mean_all_TMDs'].loc[TM01_greater_5_homol.index]
print("mean_AAIMON_filt_mean_BB", mean_AAIMON_filt_BB.mean())

# create numpy array of membranous over nonmembranous conservation ratios (identity)
hist_data_AAIMON_mean = np.array(mean_AAIMON_filt_BB)
list_BB = list(hist_data_AAIMON_mean)
data_to_plot.append(list_BB)

# use numpy to create a histogram
# freq_counts_I, bin_array_I = np.histogram(hist_data_AAIMON_mean, bins=binlist)
# assuming all of the bins are exactly the same size, make the width of the column equal to XX% (e.g. 95%) of each bin
# col_width = float('%0.3f' % (0.95 * (bin_array_I[1] - bin_array_I[0])))
# when align='center', the central point of the bar in the x-axis is simply the middle of the bins ((bin_0-bin_1)/2, etc)
# centre_of_bar_in_x_axis = (bin_array_I[:-2] + bin_array_I[1:-1]) / 2
# add the final bin, which is physically located just after the last regular bin but represents all higher values
# bar_width = centre_of_bar_in_x_axis[3] - centre_of_bar_in_x_axis[2]
# centre_of_bar_in_x_axis = np.append(centre_of_bar_in_x_axis, centre_of_bar_in_x_axis[-1] + bar_width)
# linecontainer_AAIMON_mean = ax.plot(centre_of_bar_in_x_axis, freq_counts_I, color="#4B088A", alpha=0.5)


# ax.set_xlabel('average conservation ratio (membranous over nonmembranous)', fontsize=fontsize)
# move the x-axis label closer to the x-axis
# ax.xaxis.set_label_coords(0.45, -0.085)
# pylab.rcParams['figure.figsize'] = (50.0, 40.0)
# pylab.rcParams['figure.figsize'] = (20.0, 16.0)
# plt.show()
# xlim_min = 0.5
# take x-axis max from settings
# xlim_max = 1.75
# set x-axis min
# ax.set_xlim(xlim_min, xlim_max)
# set x-axis ticks
# use the slide selection to select every second item in the list as an xtick(axis label)
# ax.set_xticks([float('%0.1f' % c) for c in centre_of_bar_in_x_axis[::3]])
# ax.set_ylabel('freq', rotation='vertical', fontsize=fontsize)
# change axis font size
# ax.tick_params(labelsize=fontsize)
# create legend?#http://stackoverflow.com/questions/9834452/how-do-i-make-a-single-legend-for-many-subplots-with-matplotlib
# legend_obj = ax.legend(['BB'], loc='upper right', fontsize=fontsize)
# improve ggplot style for a canvas (fig) with 4 figures (plots)
# utils.improve_ggplot_for_4_plots(axarr, row_nr, col_nr, backgroundcolour, legend_obj)
# FROM RIMMA SCRIPT
# ax.yaxis.grid(True, zorder=0, linestyle=":", color="grey")
# for tic in ax.xaxis.get_major_ticks():
#    tic.tick1On = False
# for tic in ax.yaxis.get_major_ticks():
#    tic.tick1On = False
# ax.spines['top'].set_visible(False)
# ax.spines['right'].set_visible(False)


#####################################################################
#                                                                   #
#                              MP                                   #
#                                                                   #
#####################################################################

# fig, ax = plt.subplots()

# iterate through the proteins that have a list of TMDs
#for acc in df_MP.loc[df_MP['list_of_TMDs'].notnull()].index:
#    dict_AAIMON_ratio_mean_MP = {}
#    for TMD in ast.literal_eval(df_MP.loc[acc, 'list_of_TMDs']):
#        dict_AAIMON_ratio_mean_MP[TMD] = df_MP.loc[acc, '%s_AAIMON_ratio_mean' % TMD]
#    df_MP.loc[acc, 'AAIMON_ratio_mean_all_TMDs'] = np.mean(list(dict_AAIMON_ratio_mean_MP.values()))

homol_cutoff = 4500
df_MP_less_4500 = pd.DataFrame()
df_MP_less_4500['TM01_AAIMON_n_homol'] = df_MP['TM01_AAIMON_n_homol'].loc[df_MP['TM01_AAIMON_n_homol'] < homol_cutoff]
print('number_of_MP_proteins_excluded = ', len(df_MP.index) - len(df_MP_less_4500.index))

TM01_greater_50_homol = df_MP_less_4500['TM01_AAIMON_n_homol'].loc[df_MP['TM01_AAIMON_n_homol'] > 50]
mean_AAIMON_filt_MP = df_MP['AAIMON_ratio_mean_all_TMDs'].loc[TM01_greater_50_homol.index]
print("mean_AAIMON_filt_mean_MP", mean_AAIMON_filt_MP.mean())

# create numpy array of membranous over nonmembranous conservation ratios (identity)
hist_data_AAIMON_mean = np.array(mean_AAIMON_filt_MP)
list_MP = list(hist_data_AAIMON_mean)
data_to_plot.append(list_MP)

#####################################################################
#                                                                   #
#                              SP                                   #
#                                                                   #
#####################################################################

# fig, ax = plt.subplots()

# iterate through the proteins that have a list of TMDs
#for acc in df_SP.loc[df_SP['list_of_TMDs'].notnull()].index:
#    dict_AAIMON_ratio_mean_SP = {}
#    for TMD in ast.literal_eval(df_SP.loc[acc, 'list_of_TMDs']):
#        dict_AAIMON_ratio_mean_SP[TMD] = df_SP.loc[acc, '%s_AAIMON_ratio_mean' % TMD]
#    df_SP.loc[acc, 'AAIMON_ratio_mean_all_TMDs'] = np.mean(list(dict_AAIMON_ratio_mean_SP.values()))
homol_cutoff = 4500
df_SP_less_4500 = pd.DataFrame()
df_SP_less_4500['TM01_AAIMON_n_homol'] = df_SP['TM01_AAIMON_n_homol'].loc[df_SP['TM01_AAIMON_n_homol'] < homol_cutoff]
print('number_of_SP_proteins_excluded = ', len(df_SP.index) - len(df_SP_less_4500.index))

TM01_greater_50_homol = df_SP_less_4500['TM01_AAIMON_n_homol'].loc[df_SP['TM01_AAIMON_n_homol'] > 50]
mean_AAIMON_filt_SP = df_SP['AAIMON_ratio_mean_all_TMDs'].loc[TM01_greater_50_homol.index]
#print(mean_AAIMON_filt_SP)
print("mean_AAIMON_filt_mean_SP", mean_AAIMON_filt_SP.mean())

# create numpy array of membranous over nonmembranous conservation ratios (identity)
hist_data_AAIMON_mean = np.array(mean_AAIMON_filt_SP)
list_SP = list(hist_data_AAIMON_mean)
data_to_plot.append(list_SP)

# transfer data-list to a dataframe
dataframe = pd.DataFrame(data_to_plot, index=['BB', 'MP', 'SP'])
dataframe = dataframe.transpose()

# create histogram
trace1 = go.Histogram(
    x=data_to_plot[0],
    histnorm='count',
    name='BB',
    autobinx=False,
    xbins=dict(
        start=0.6,
        end=1.5,
        size=0.015
    ),
    marker=dict(
        color='blue',
        line=dict(
            color='grey',
            width=0.5
        )
    ),
    opacity=0.5
)
trace2 = go.Histogram(
    x=data_to_plot[1],
    histnorm='count',
    name='MP',
    autobinx=False,
    xbins=dict(
        start=0.6,
        end=1.5,
        size=0.015
    ),
    marker=dict(
        color='red',
        line=dict(
            color='grey',
            width=0.5
        )
    ),
    opacity=0.3
)
trace3 = go.Histogram(
    x=data_to_plot[2],
    histnorm='count',
    name='SP',
    autobinx=False,
    xbins=dict(
        start=0.6,
        end=1.5,
        size=0.015
    ),
    marker=dict(
        color='gray',
        line=dict(
            color='grey',
            width=0.5
        )
    ),
    opacity=0.4
)

data = [trace1, trace2, trace3]

layout = go.Layout(
    title='SP_MP_BB',
    font=dict(
        size=20
    ),
    xaxis=dict(
        title='AAIMON ratio mean',
        titlefont=dict(
            size=20
        ),
    ),
    yaxis=dict(
        title='freq',
        titlefont=dict(
            size=20
        ),
    ),
    barmode='overlay',
    bargap=0,
    bargroupgap=0.2
)
fig = go.Figure(data=data, layout=layout)
# py.iplot(fig)

#py.image.save_as(fig, filename='SP_MP_BB_40-100_perc.pdf')


#plotly.offline.iplot({
#    "data": data,
#    "layout": layout
#})