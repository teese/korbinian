import ast
import csv
import itertools
import korbinian.mtutils as utils
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd

def create_graph_of_gap_density(pathdict, set_, logging):
    logging.info('~~~~~~~~~~~~starting creating graphs of gap density~~~~~~~~~~~~')

    # test if the dataframe has already been created, otherwise re-open from uniprot csv file
    if os.path.isfile(pathdict["dfout10_uniprot_gaps"]):
        df = pd.read_csv(pathdict["dfout10_uniprot_gaps"], sep=",", quoting=csv.QUOTE_NONNUMERIC, index_col=[0])
        logging.info('df loaded from %s' % pathdict["dfout10_uniprot_gaps"])
    else:
        raise FileNotFoundError(
            'No gap analysis has been done yet. %s is not found. Please run calculate calculate_gap_densities' % pathdict[
                "dfout10_uniprot_gaps"])

    num_of_bins_in_tmd_region = set_["num_of_bins_in_tmd_region"]
    # find the maximum number of TMDs amongst the proteins
    n_TMDs_max = df["number_of_TMDs"].max()

    flipped = []
    not_flipped = []

    for acc in df.index:
        print(acc)

        for num_TMD in range(1, int(df.loc[acc, "number_of_TMDs"]) + 1, 2):

            if not utils.isNaN(df.loc[acc, "TM%.2d_occuring_gaps" % num_TMD]):
                for n in ast.literal_eval(df.loc[acc, "TM%.2d_occuring_gaps" % num_TMD]):
                    # print (n)
                    # print ((df.loc[acc,"TM%.2d_len"%num_TMD]-1) )

                    if df.loc[acc, "n_term_ec"] == False:
                        not_flipped.append((n / (len(df.loc[acc, "TM%.2d_seq" % num_TMD]) - 1)) * num_of_bins_in_tmd_region)
                    if df.loc[acc, "n_term_ec"] == True:
                        flipped.append(num_of_bins_in_tmd_region - (
                        (n / (len(df.loc[acc, "TM%.2d_seq" % num_TMD]) - 1)) * num_of_bins_in_tmd_region))
                        # print ((n/(df.loc[acc,"TM%.2d_len"%num_TMD]-1))*10)

        for num_TMD in range(2, int(df.loc[acc, "number_of_TMDs"]) + 1, 2):
            # print(num_TMD)
            if not utils.isNaN(df.loc[acc, "TM%.2d_occuring_gaps" % num_TMD]):
                for n in ast.literal_eval(df.loc[acc, "TM%.2d_occuring_gaps" % num_TMD]):
                    # print(n)
                    # print(df.loc[acc,"TM%.2d_len"%num_TMD])

                    # m = np.prod(df.loc[acc,"TM%.2d_len"%num_TMD]*float(n))

                    if df.loc[acc, "n_term_ec"] == False:
                        flipped.append(num_of_bins_in_tmd_region - (
                        (n / (len(df.loc[acc, "TM%.2d_seq" % num_TMD]) - 1)) * num_of_bins_in_tmd_region))

                    if df.loc[acc, "n_term_ec"] == True:
                        not_flipped.append((n / (len(df.loc[acc, "TM%.2d_seq" % num_TMD]) - 1)) * num_of_bins_in_tmd_region)

        fig, axarr = plt.subplots()

        tmd_gaps = flipped + not_flipped

        n, bins, patches = plt.hist(tmd_gaps, bins=num_of_bins_in_tmd_region, range=(0, num_of_bins_in_tmd_region),
                                    normed=True, histtype='bar', color="darkblue", zorder=3, edgecolor="white")

        axarr.set_ylabel("Propability for gapposition", fontsize=8, color="black")
        axarr.set_xlabel("Depth of Residue", fontsize=8, color="black")
        # axarr.set_title("Gap distribution along TMD", fontsize=16,color="black")

        axarr.set_xlim(0, num_of_bins_in_tmd_region)
        axarr.set_xticks(np.linspace(0.8, 19.2, 2))
        axarr.set_xticklabels(("intracellular", "extracellular"), fontsize=8, color="black")
        axarr.tick_params(reset=True, labelsize=8, color="black")

        for tic in axarr.xaxis.get_major_ticks():
            tic.tick1On = tic.tick2On = False
        axarr.grid(False, 'minor', color='0.99', linestyle='-', linewidth=0.7)
        # fig.patch.set_visible(False)

        axarr.spines['top'].set_visible(False)
        axarr.spines['right'].set_visible(False)
        axarr.yaxis.tick_left()
        axarr.patch.set_facecolor('white')
        axarr.yaxis.grid(True, zorder=0, color="grey", linestyle=":")
        # fig.savefig("/nas/teeselab/students/rimma/omp/summaries/test.png", format='png', dpi=200)
        fig.savefig(pathdict["dfout11_gap_test_out_png"], format='png', dpi=200)

        TMD_range = range(1, n_TMDs_max + 1)
        TMD_range_plus_1 = range(1, n_TMDs_max + 2)
        TMD_range_2nd = range(1, n_TMDs_max + 1, 2)

        # nested_list_of_gaps_intracellular = [ast.literal_eval(m) for n in range (1,25) for m in df['juxta_TM%.2d_intracellular_possible_gappositions'%n].dropna().tolist()]
        nested_list_of_gaps_intracellular = [ast.literal_eval(m) for n in TMD_range for m in
                                             df['juxta_TM%.2d_intracellular_possible_gappositions' % n].dropna().tolist()]

        hist_data_juxta_intracellular = np.array(list(itertools.chain(*nested_list_of_gaps_intracellular)))

        min_value = int(abs(hist_data_juxta_intracellular.min()))
        # y-axis_intracell = [(freq_counts_I.tolist()[::-1][n]/frequency_of_position_intracellular(n)) for n in range (1,int(min_value))]

        # data for extracellular part --> right
        nested_list_of_gaps_extracellular = [ast.literal_eval(m) for n in TMD_range for m in
                                             df['juxta_TM%.2d_extracellular_possible_gappositions' % n].dropna().tolist()]
        hist_data_juxta_extracellular = np.array(list(itertools.chain(*nested_list_of_gaps_extracellular)))

        # data for tmd --> middle

        hist_data_tmds = np.array(tmd_gaps)
        # times 2, because TMDs in gap and in query are considered! --> double amount
        total_amount_of_TMDs_in_protein = df.loc[
                                              df.gaps_analysed == True, "number_of_TMDs"].sum() * 2
        # total_amount_of_TMDs_in_protein = len([TMD for acc in df.index for TMD in ast.literal_eval(df.loc[acc,"list_of_TMDs"])if (df.loc[acc,"gaps_analysed"]==True)and(utils.isNaN(df.loc[acc,"list_of_TMDs"]))])*2

        list_of_positionfrequency_extra = []
        list_of_positionfrequency_intra = []

        for acc in df.index:
            if df.loc[acc, "gaps_analysed"] == True:
                print(acc)
                if df.loc[acc, "n_term_ec"] == True:
                    for n in TMD_range_plus_1:
                        list_of_positionfrequency_extra.append(df.loc[acc, 'len_juxta_before_TM%.2d' % n].tolist())
                        list_of_positionfrequency_intra.append(df.loc[acc, 'len_juxta_after_TM%.2d' % n].tolist())
                    for n in TMD_range_2nd:
                        list_of_positionfrequency_extra.append(df.loc[acc, 'len_juxta_after_TM%.2d' % n].tolist())
                        list_of_positionfrequency_intra.append(df.loc[acc, 'len_juxta_before_TM%.2d' % n].tolist())

                utils.aaa(df)
                if df.loc[acc, "n_term_ec"] == False:
                    for n in TMD_range_plus_1:
                        list_of_positionfrequency_extra.append(df.loc[acc, 'len_juxta_after_TM%.2d' % n].tolist())
                        list_of_positionfrequency_intra.append(df.loc[acc, 'len_juxta_before_TM%.2d' % n].tolist())
                    for n in TMD_range_2nd:
                        list_of_positionfrequency_extra.append(df.loc[acc, 'len_juxta_before_TM%.2d' % n].tolist())
                        list_of_positionfrequency_intra.append(df.loc[acc, 'len_juxta_after_TM%.2d' % n].tolist())


        def positionfreq_in_list_extra(position):
            return (len([n for n in list_of_positionfrequency_extra if n >= position]) * 2)


        def positionfreq_in_list_intra(position):
            return (len([n for n in list_of_positionfrequency_intra if n >= position]) * 2)


        fig, axarr = plt.subplots()

        fontsize = 8

        ### Intracellular

        freq_counts_I, bin_array_I = np.histogram(hist_data_juxta_intracellular, bins=hist_data_juxta_intracellular.max())
        centre_of_bar_in_x_axis_I = -((bin_array_I[:-2] + bin_array_I[1:-1]) / 2)

        bar_width_I = centre_of_bar_in_x_axis_I[3] - centre_of_bar_in_x_axis_I[2]

        centre_of_bar_in_x_axis_I = np.append(centre_of_bar_in_x_axis_I, centre_of_bar_in_x_axis_I[-1] + bar_width_I)

        axarr.bar(left=centre_of_bar_in_x_axis_I,
                  height=[((freq_counts_I.tolist()[n]) / (positionfreq_in_list_intra(n))) for n in range(0, 32)], width=0.4,
                  color="mediumblue", linewidth=0, zorder=3)  # edgecolor='black',

        ######### TMD

        hist_data_tmds = np.array(tmd_gaps)
        freq_counts_II, bin_array_II = np.histogram(hist_data_tmds, bins=20, range=(0, 20))

        centre_of_bar_in_x_axis_II = (bin_array_II[:-2] + bin_array_II[1:-1]) / 2

        bar_width_II = centre_of_bar_in_x_axis_II[3] - centre_of_bar_in_x_axis_II[2]

        centre_of_bar_in_x_axis_II = np.append(centre_of_bar_in_x_axis_II, centre_of_bar_in_x_axis_II[-1] + bar_width_II)

        axarr.bar(left=centre_of_bar_in_x_axis_II,
                  height=[n / total_amount_of_TMDs_in_protein for n in freq_counts_II.tolist()], align='center', width=0.5,
                  color="blue", linewidth=0, zorder=3)  # edgecolor='black',

        #####Extracellular

        freq_counts_III, bin_array_III = np.histogram(hist_data_juxta_extracellular,
                                                      bins=hist_data_juxta_extracellular.max())

        centre_of_bar_in_x_axis_III = ((bin_array_III[:-2] + bin_array_III[1:-1]) / 2)

        bar_width_III = centre_of_bar_in_x_axis_III[3] - centre_of_bar_in_x_axis_III[2]

        centre_of_bar_in_x_axis_III = np.append(centre_of_bar_in_x_axis_III,
                                                centre_of_bar_in_x_axis_III[-1] + bar_width_III)

        axarr.bar(left=centre_of_bar_in_x_axis_III + 19.5,
                  height=[((freq_counts_III.tolist()[n]) / (positionfreq_in_list_extra(n))) for n in range(0, 32)],
                  width=0.4, color="mediumblue", linewidth=0, zorder=3)  # edgecolor='black',

        # axarr.bar(left=centre_of_bar_in_x_axis_III+9.5, height=[(freq_counts_III.tolist()[n]/frequency_of_position_extracellular(n)) for n in range (0,207)], align = 'center', width=0.4, color="mediumblue" ,linewidth=0)  # edgecolor='black',


        # axarr.bar(left=centre_of_bar_in_x_axis_III+9.5, height=[(freq_counts_III.tolist()[n]/frequency_of_position_extracellular(n)) for n in range (0,207)], align = 'center', width=0.4, color="mediumblue" ,linewidth=0)  # edgecolor='black',

        ##### Style
        axarr.set_xlim(-30, 40)
        axarr.set_ylim(0, 8)
        axarr.set_xlabel("residue position relative to TM domain", fontsize=fontsize)
        axarr.set_ylabel("gap propensity", color="mediumblue", fontsize=fontsize)

        axarr.set_xticks([-30, -25, -20, -15, -10, -5, -1, 10, 21, 25, 30, 35, 40, 45, 50])
        labels = ["-30", "-25", "-20", "-15", "-10", "-5", "-1", "TM domain", "1", "5", "10", "15", "20", "25", "30"]
        axarr.set_xticklabels(labels, color="mediumblue")

        axarr.tick_params(reset=True, labelsize=8, color="mediumblue")

        for tic in axarr.xaxis.get_major_ticks():
            tic.tick1On = tic.tick2On = False

        axarr.grid(False)

        # fig.patch.set_visible(False)

        # axarr.spines['top'].set_visible(False)
        # axarr.spines['right'].set_visible(False)
        # axarr.spines['left'].set_visible(True)
        # axarr.spines['right'].set_visible(False)

        axarr.yaxis.tick_left()

        axarr.patch.set_facecolor('white')
        axarr.yaxis.grid(True, zorder=0, linestyle=":", color="g")

        axarr.add_patch(patches.Rectangle((0, 0), 20, 10, alpha=0.1, linewidth=None, facecolor="lightseagreen"))

        axarr.annotate('Intracellular', xy=(0, 0), xytext=(-20, 1.7), alpha=0.5, fontsize=fontsize)
        axarr.annotate('Extracellular', xy=(0, 0), xytext=(30, 1.7), alpha=0.5, fontsize=fontsize)
        axarr.annotate('TM', xy=(10, 0), xytext=(8.6, 1.3), alpha=0.5, fontsize=fontsize)
        axarr.annotate('helix', xy=(10, 0), xytext=(8.2, 1.2), alpha=0.5, fontsize=fontsize)

        axarr.spines['left'].set_color('mediumblue')
        # axarr.spines['left'].set_visible(True)
        # axarr.spines['bottom'].set_color('black')
        #### Second y- axis
        # axarr.spines[]

        axarr2 = axarr.twinx()
        axarr2.spines['right'].set_color("black")
        axarr2.grid(False)

        axarr2.plot([n + 20 for n in range(0, 30)], [positionfreq_in_list_extra(n) for n in range(0, 30)], "black",
                    linewidth=1.2)
        # axarr2.plot()

        axarr2.plot([-n for n in range(0, 30)], [positionfreq_in_list_intra(n) for n in range(0, 30)], "black",
                    linewidth=1.2)
        # axarr2.plot()


        # axarr2.plot(x_axis,y_axis,"black",linewidth=1.2)

        axarr2.set_ylabel("Frequency of considered position in dataset", color="black", fontsize=fontsize)

        axarr2.yaxis.label.set_size(fontsize)
        axarr2.tick_params(axis='y', colors='black', labelsize=fontsize)
        axarr2.spines["left"].set_color("mediumblue")

        axarr.set_xlim(-30, 50)
        axarr.set_ylim(0, 2)
        axarr.set_yticks = ([0, 0.2, 0.4, 0.6, 0.8, 1.0])
        axarr.set_yticklabels(["0.0", "0.2", "0.4", "0.6", "0.8", "1.0"], color="mediumblue")

        # rcParams['savefig.dpi'] = 200
        #fig.savefig("/nas/teeselab/students/rimma/omp/summaries/Test234.png", format='png', dpi=200)
        fig.savefig(pathdict["create_graph_of_gap_density_png"], format='png', dpi=200)
