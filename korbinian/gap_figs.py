import ast
import csv
import itertools
import korbinian.utils as utils
import matplotlib.pyplot as plt
import numpy as np
import pickle
import sys

def create_graph_of_gap_density(pathdict, s, logging):
    logging.info('~~~~~~~~~~~~starting creating graphs of gap density~~~~~~~~~~~~')

    num_of_bins_in_tmd_region = s["num_of_bins_in_tmd_region"]

    with open(pathdict["gap_data_pickle"], "rb") as pkl_file:
        gap_data = pickle.load(pkl_file)

    data_names = ["flipped", "not_flipped", "hist_data_juxta_intracellular", "hist_data_juxta_extracellular", "min_value", "list_of_positionfrequency_extra", "list_of_positionfrequency_intra", "n_TMDs_in_all_proteins"]

    sys.stdout.write(type(gap_data))
    sys.stdout.write(len(gap_data))
    for i in range(len(gap_data)):
        list_or_int_in_gap_data = gap_data[i]
        name = data_names[i]
        if isinstance(list_or_int_in_gap_data, int):
            to_print = list_or_int_in_gap_data
            length = 1
        elif isinstance(list_or_int_in_gap_data, list) or isinstance(list_or_int_in_gap_data, np.ndarray):
            to_print = list_or_int_in_gap_data[0:10]
            length = len(list_or_int_in_gap_data)

        sys.stdout.write("name = {}, dtype = {}, len = {}, value (or value[0:10]) = {}".format(name, type(list_or_int_in_gap_data), length, to_print))

    flipped, not_flipped, hist_data_juxta_intracellular, hist_data_juxta_extracellular, min_value, list_of_positionfrequency_extra, list_of_positionfrequency_intra, n_TMDs_in_all_proteins = gap_data

    # the n_TMDs_in_all_proteins is way too large, normalisation is screwed up somehow


    ######################################################################################################################
    #                                                                                                                    #
    #                            Plot figures with the collected and normalised gap data                                 #
    #                                                                                                                    #
    ######################################################################################################################

    ######################################################################################################################
    #                                                                                                                    #
    #                                         fig_a,  fig_a : Test Figure ?                                              #
    #                                        based on flipped and not_flipped                                            #
    #  currently not working, as flipped and not_flipped are too small, missing append from a loop in gap.py somewhere   #
    #                                                                                                                    #
    ######################################################################################################################


    fig_a, ax_a = plt.subplots()

    tmd_gaps = flipped + not_flipped

    n, bins, patches = plt.hist(tmd_gaps, bins=num_of_bins_in_tmd_region, range=(0, num_of_bins_in_tmd_region),
                                normed=True, histtype='bar', color="darkblue", zorder=3, edgecolor="white")

    ax_a.set_ylabel("Propability for gapposition", fontsize=8, color="black")
    ax_a.set_xlabel("Depth of Residue", fontsize=8, color="black")
    # ax_a.set_title("Gap distribution along TMD", fontsize=16,color="black")

    ax_a.set_xlim(0, num_of_bins_in_tmd_region)
    ax_a.set_xticks(np.linspace(0.8, 19.2, 2))
    ax_a.set_xticklabels(("intracellular", "extracellular"), fontsize=8, color="black")
    ax_a.tick_params(reset=True, labelsize=8, color="black")

    for tic in ax_a.xaxis.get_major_ticks():
        tic.tick1On = tic.tick2On = False
    ax_a.grid(False, 'minor', color='0.99', linestyle='-', linewidth=0.7)
    # fig_a.patch.set_visible(False)
    ax_a.spines['top'].set_visible(False)
    ax_a.spines['right'].set_visible(False)
    ax_a.yaxis.tick_left()
    ax_a.patch.set_facecolor('white')
    ax_a.yaxis.grid(True, zorder=0, color="grey", linestyle=":")
    # save fig to the summaries folder
    fig_a.savefig(pathdict["gap_density_testfig_png"], format='png', dpi=200)

    ######################################################################################################################
    #                                                                                                                    #
    #               fig_b,  ax_b : Full figure showing +/- 30 residues, and number of datapoints                         #
    #                                                                                                                    #
    ######################################################################################################################

    fig_b, ax_b = plt.subplots()

    fontsize = 8
    ### Intracellular
    # NO IDEA WHY bins=hist_data_juxta_intracellular.max()
    """
    C:\Anaconda3\lib\site-packages\numpy\lib\function_base.py:564: VisibleDeprecationWarning: using a non-integer number instead of an integer will result in an error in the future
      n = np.zeros(bins, ntype)
    Traceback (most recent call last):
      File "D:/Schweris/Projects/Programming/Python/scripts/Pycharm_Projects/korbinian/korbinian/run_korbinian.py", line 151, in <module>
        korbinian.gap_figs.create_graph_of_gap_density(pathdict, s, logging)
      File "D:\Schweris\Projects\Programming\Python\scripts\Pycharm_Projects\korbinian\korbinian\gap_figs.py", line 293, in create_graph_of_gap_density
        freq_counts_I, bin_array_I = np.histogram(hist_data_juxta_intracellular, bins=hist_data_juxta_intracellular.max())
      File "C:\Anaconda3\lib\site-packages\numpy\lib\function_base.py", line 603, in histogram
        increment = (tmp_a_data >= bin_edges[indices + 1]) & (indices != bins - 1)
    IndexError: index 19 is out of bounds for axis 1 with size 19
"""
    # get the max value from the hist_data_juxta arrays, as an int
    hist_data_juxta_intracellular_ceil = int(np.ceil(hist_data_juxta_intracellular.max()))
    hist_data_juxta_extracellular_ceil = int(np.ceil(hist_data_juxta_extracellular.max()))

    # create frequency counts and array using np.histogram for the intracellular side
    # note that the bins=hist_data_juxta_intracellular_ceil is just giving the NUMBER of bins.
    freq_counts_I, bin_array_I = np.histogram(hist_data_juxta_intracellular, bins=hist_data_juxta_intracellular_ceil)

    sys.stdout.write("freq_counts_I :")
    sys.stdout.write(len(freq_counts_I), freq_counts_I[0:10])

    # In order to plot a histogram based on the central datapoint (e.g. as a line), the centre of each bin must be calculated

    centre_of_bar_in_x_axis_I = np.negative((bin_array_I[:-2] + bin_array_I[1:-1]) / 2)
    # when align='center', the central point of the bar in the x-axis is simply the middle of the bins ((bin_0-bin_1)/2, etc)
    bar_width_I = centre_of_bar_in_x_axis_I[3] - centre_of_bar_in_x_axis_I[2]
    # add the final bin, which is physically located just after the last regular bin but represents all higher values
    centre_of_bar_in_x_axis_I = np.append(centre_of_bar_in_x_axis_I, centre_of_bar_in_x_axis_I[-1] + bar_width_I)

    # create a list of normalised frequency counts, to plot on the graph
    # v = each frequency, normalised to the number of datapoints
    v_list_I = []
    for n in range(hist_data_juxta_intracellular_ceil):
        f = freq_counts_I[n]
        p = positionfreq_in_list_intra(list_of_positionfrequency_intra, n)
        v = f/p
        v_list_I.append(v)

    #ax_b.bar(left=centre_of_bar_in_x_axis_I, height=[(freq_counts_I.tolist()[n] / (positionfreq_in_list_intra(n))) for n in range(0, n_TMDs_max)], width=0.4, color="mediumblue", linewidth=0, zorder=3)  # edgecolor='black',
    # plot as a bar chart
    ax_b.bar(left=centre_of_bar_in_x_axis_I, height=v_list_I, width=0.4, color="mediumblue", linewidth=0, zorder=3)

    ######### TMD
    # data for tmd --> middle
    hist_data_tmds = np.array(tmd_gaps)
    # create frequency counts and array using np.histogram for the TM region
    # note that the bins=hist_data_juxta_intracellular_ceil is just giving the NUMBER of bins.
    freq_counts_TM, bin_array_TM = np.histogram(hist_data_tmds, bins=20, range=(0, 20))
    centre_of_bar_in_x_axis_TM = (bin_array_TM[:-2] + bin_array_TM[1:-1]) / 2
    bar_width_TM = centre_of_bar_in_x_axis_TM[3] - centre_of_bar_in_x_axis_TM[2]
    centre_of_bar_in_x_axis_TM = np.append(centre_of_bar_in_x_axis_TM, centre_of_bar_in_x_axis_TM[-1] + bar_width_TM)

    rimma_orig_height = np.array([n / n_TMDs_in_all_proteins for n in freq_counts_TM.tolist()])
    sys.stdout.write("rimma_orig_height", rimma_orig_height)
    rimma_orig_height_by_1000 = rimma_orig_height * 1000
    # plot the TMD data
    ax_b.bar(left=centre_of_bar_in_x_axis_TM,height=rimma_orig_height_by_1000, align='center', width=0.5, color="blue", linewidth=0, zorder=3)  # edgecolor='black',

    # freq counts and bin array for the extracellular side
    freq_counts_E, bin_array_E = np.histogram(hist_data_juxta_extracellular, bins=hist_data_juxta_extracellular_ceil)
    sys.stdout.write("freq_counts_E :")
    sys.stdout.write(len(freq_counts_E), freq_counts_E[0:10])
    # NO IDEA WHY bins=hist_data_juxta_intracellular.max()
    #freq_counts_E, bin_array_E = np.histogram(hist_data_juxta_extracellular)
    #####Extracellular
    v_list_E = []
    for n in range(hist_data_juxta_extracellular_ceil):
        f = freq_counts_E[n]
        p = positionfreq_in_list_extra(list_of_positionfrequency_extra, n)
        v = f/p
        v_list_E.append(v)

    centre_of_bar_in_x_axis_E = ((bin_array_E[:-2] + bin_array_E[1:-1]) / 2)
    bar_width_E = centre_of_bar_in_x_axis_E[3] - centre_of_bar_in_x_axis_E[2]
    centre_of_bar_in_x_axis_E = np.append(centre_of_bar_in_x_axis_E, centre_of_bar_in_x_axis_E[-1] + bar_width_E)

    # converted to v_list_E with the list of normalised values
    #ax_b.bar(left=centre_of_bar_in_x_axis_E + 19.5,height=[((freq_counts_E.tolist()[n]) / (positionfreq_in_list_extra(n))) for n in range(0, n_TMDs_max)], width=0.4, color="mediumblue", linewidth=0, zorder=3)  # edgecolor='black',
    ax_b.bar(left=centre_of_bar_in_x_axis_E + 19.5, height=v_list_E, width=0.4, color="mediumblue", linewidth=0, zorder=3)  # edgecolor='black',
    # ax_b.bar(left=centre_of_bar_in_x_axis_E+9.5, height=[(freq_counts_E.tolist()[n]/frequency_of_position_extracellular(n)) for n in range (0,207)], align = 'center', width=0.4, color="mediumblue" ,linewidth=0)  # edgecolor='black',
    # ax_b.bar(left=centre_of_bar_in_x_axis_E+9.5, height=[(freq_counts_E.tolist()[n]/frequency_of_position_extracellular(n)) for n in range (0,207)], align = 'center', width=0.4, color="mediumblue" ,linewidth=0)  # edgecolor='black',

    ##### Style
    ax_b.set_xlim(-30, 40)
    ax_b.set_ylim(0, 8)
    ax_b.set_xlabel("residue position relative to TM domain", fontsize=fontsize)
    ax_b.set_ylabel("gap propensity", color="mediumblue", fontsize=fontsize)

    ax_b.set_xticks([-30, -25, -20, -15, -10, -5, -1, 10, 21, 25, 30, 35, 40, 45, 50])
    labels = ["-30", "-25", "-20", "-15", "-10", "-5", "-1", "TM domain", "1", "5", "10", "15", "20", "25", "30"]
    ax_b.set_xticklabels(labels, color="mediumblue")

    ax_b.tick_params(reset=True, labelsize=8, color="mediumblue")

    for tic in ax_b.xaxis.get_major_ticks():
        tic.tick1On = tic.tick2On = False

    ax_b.grid(False)

    # fig_b.patch.set_visible(False)
    # ax_b.spines['top'].set_visible(False)
    # ax_b.spines['right'].set_visible(False)
    # ax_b.spines['left'].set_visible(True)
    # ax_b.spines['right'].set_visible(False)
    ax_b.yaxis.tick_left()

    ax_b.patch.set_facecolor('white')
    ax_b.yaxis.grid(True, zorder=0, linestyle=":", color="g")

    #ax_b.add_patch(patches.Rectangle((0, 0), 20, 10, alpha=0.1, linewidth=None, facecolor="lightseagreen"))

    ax_b.annotate('Intracellular', xy=(0, 0), xytext=(-20, 1.7), alpha=0.5, fontsize=fontsize)
    ax_b.annotate('Extracellular', xy=(0, 0), xytext=(30, 1.7), alpha=0.5, fontsize=fontsize)
    ax_b.annotate('TM', xy=(10, 0), xytext=(8.6, 1.3), alpha=0.5, fontsize=fontsize)
    ax_b.annotate('helix or sheet', xy=(10, 0), xytext=(8.2, 1.2), alpha=0.5, fontsize=fontsize)

    ax_b.spines['left'].set_color('mediumblue')
    # ax_b.spines['left'].set_visible(True)
    # ax_b.spines['bottom'].set_color('black')
    #### Second y- axis
    # ax_b.spines[]

    ax2 = ax_b.twinx()
    ax2.spines['right'].set_color("black")
    ax2.grid(False)

    # ax2.plot([n + 20 for n in range(0, 30)], [positionfreq_in_list_extra(n) for n in range(0, 30)], "black",
    #             linewidth=1.2)
    # # ax2.plot()
    #
    # ax2.plot([-n for n in range(0, 30)], [positionfreq_in_list_intra(n) for n in range(0, 30)], "black",
    #             linewidth=1.2)
    # ax2.plot()

    # ax2.plot(x_axis,y_axis,"black",linewidth=1.2)

    ax2.set_ylabel("Frequency of considered position in dataset", color="black", fontsize=fontsize)

    ax2.yaxis.label.set_size(fontsize)
    ax2.tick_params(axis='y', colors='black', labelsize=fontsize)
    ax2.spines["left"].set_color("mediumblue")

    ax_b.set_xlim(-30, 50)
    ax_b.set_ylim(0, 2)
    ax_b.set_yticks = ([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax_b.set_yticklabels(["0.0", "0.2", "0.4", "0.6", "0.8", "1.0"], color="mediumblue")
    # save fig to the summaries folder
    fig_b.savefig(pathdict["gap_density_fig_path"], format='png', dpi=200)

def positionfreq_in_list_extra(list_of_positionfrequency_extra, position):
    return (len([n for n in list_of_positionfrequency_extra if n >= position]) * 2)

def positionfreq_in_list_intra(list_of_positionfrequency_intra, position):
    return (len([n for n in list_of_positionfrequency_intra if n >= position]) * 2)

def juxta_function_orig(df, TMD):
    """
    Parameters
    ----------
    df : pd.Dataframe
        Dataframe for Sequences (df). This is the dataframe containing the full homologue sequences from the BLAST-like data analysis.
    TMD : str
        String denoting transmembrane domain number (e.g. "TM01")

    Returns
    -------

    """
    next_TM = "TM{:02d}".format(int(TMD[2:]) + 1)
    prev_TM = "TM{:02d}".format(int(TMD[2:])-1)

    df['start_juxta_after_%s'%TMD] = np.where(utils.isNaN(df['{}_start'.format(next_TM)]) == True, np.nan, df['%s_end'%TMD])

    df['end_juxta_before_%s'%TMD] = np.where(df["%s_start"%TMD] != 0, df["%s_start"%TMD], np.nan)

    df['end_juxta_after_%s'%TMD] = df["%s_end"%TMD] + ((df["{}_start".format(next_TM)]-df["%s_end"%TMD])/2).apply(lambda x :int(x) if not np.isnan(x) else np.nan)

    df['start_juxta_before_%s'%TMD] = np.where(df["end_juxta_after_{}".format(prev_TM)] == df['end_juxta_before_%s'%TMD], df["end_juxta_after_{}".format(prev_TM)],df["end_juxta_after_{}".format(prev_TM)])
    return df