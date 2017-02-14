import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import ast
import sys
import csv
import datetime
import korbinian.utils as utils


def compare_lists (s):
    sys.stdout.write("\n\n~~~~~~~~~~~~         starting compare_lists           ~~~~~~~~~~~~\n\n")
    create_legend = True
    save_png = s['save_png']
    save_pdf = s['save_pdf']
    color_list = ['#949494', '#EE762C', '#005C96', '#A1B11A', '#9ECEEC', '#0076B8', '#454545']
    protein_lists = ast.literal_eval(s['protein_lists'])
    # initialise dataframe that hols all variables for every protein list
    dfv = pd.DataFrame(index=protein_lists)
    # get current date and time for folder name, join elements in protein_lists to string
    str_protein_lists = '-'.join(map(str, sorted(protein_lists)))
    now = datetime.datetime.now()
    folder_name = '{}{:02}{:02}-{:02}{:02}_Lists-{}'.format(now.year, now.month, now.day, now.hour, now.minute, str_protein_lists)

    if s["testing_mode"]:
        # folder_name for testing
        folder_name = 'Folder_Test'

    base_filepath = os.path.join(s["data_dir"], "compare_lists", folder_name)

    sys.stdout.write('data is saved to "{}"\n'.format(base_filepath))

    # create folder in list summary directory to hold keyword data
    if not os.path.exists(base_filepath):
        os.makedirs(base_filepath)

    # import datasets, load data and hold pandas dataframes in dictionary
    df_dict = {}
    for n, prot_list in enumerate(protein_lists):
        dfv.loc[prot_list, 'list_description'] = s['list_description'][prot_list]
        dfv.loc[prot_list, 'base_filename_lists'] = os.path.join(s["data_dir"], "summaries", '%02d' % prot_list, 'List%02d.csv' % prot_list)
        dfv.loc[prot_list, 'base_filename_cr_summaries'] = os.path.join(s["data_dir"], "summaries", '%02d' % prot_list, 'List%02d_cr_summary.csv' % prot_list)
        dfv.loc[prot_list, 'compare_dir'] = base_filepath
        dfv.loc[prot_list, 'color'] = color_list[n]
        dfv.loc[prot_list, 'min_homol'] = s['min_homol'][prot_list]

        # read list summary.csv and cr_summary.csv from disk, join them to one big dataframe
        dfx = pd.read_csv(dfv.loc[prot_list, 'base_filename_lists'], sep=",", quoting=csv.QUOTE_NONNUMERIC, index_col=0)
        dfy = pd.read_csv(dfv.loc[prot_list, 'base_filename_cr_summaries'], sep=",", quoting=csv.QUOTE_NONNUMERIC, index_col=0)
        df_temp = pd.merge(dfy, dfx, left_index=True, right_index=True, suffixes=('_dfy', ''))
        # count proteins before and after dropping proteins with less than x homologues (x from settings file "lists"), drop proteins
        proteins_before_dropping = len(df_temp)
        df_temp = df_temp[np.isfinite(df_temp['AAIMON_mean_all_TMDs'])]
        df_temp = df_temp[df_temp.TM01_AAIMON_n_homol >= dfv.loc[prot_list, 'min_homol']]
        proteins_after_dropping = len(df_temp)
        sys.stdout.write('List{:02} - {}: number of dropped proteins {}\n'.format(prot_list, dfv.loc[prot_list,'list_description'], proteins_before_dropping - proteins_after_dropping))
        # add dataframe to a dictionary of dataframes
        df_dict[prot_list] = df_temp

    ###############################################################
    #                                                             #
    #                           plot data                         #
    #                                                             #
    ###############################################################

    # set up general plotting parameters
    plt.style.use('seaborn-whitegrid')
    alpha = 1
    fontsize = 16
    linewidth = 2

    Fig_Nr = 1
    title = 'histograms AAIMON and AAIMON_n'
    Fig_name = 'Fig01_Histograms_of_mean_AAIMON_and_AAIMON_n'
    binlist = np.linspace(0, 2, 61)
    fig, ax = plt.subplots()
    offset = len(protein_lists) - 1

    for prot_list in protein_lists:
        ###   non-normalised AAIMON   ###
        # create numpy array of membranous over nonmembranous conservation ratios (identity)
        hist_data = np.array(df_dict[prot_list]['AAIMON_mean_all_TMDs'])
        # use numpy to create a histogram
        freq_counts, bin_array = np.histogram(hist_data, bins=binlist)
        freq_counts_normalised = freq_counts / freq_counts.max() + offset
        # assuming all of the bins are exactly the same size, make the width of the column equal to XX% (e.g. 95%) of each bin
        col_width = float('%0.3f' % (0.95 * (bin_array[1] - bin_array[0])))
        # when align='center', the central point of the bar in the x-axis is simply the middle of the bins ((bin_0-bin_1)/2, etc)
        centre_of_bar_in_x_axis = (bin_array[:-2] + bin_array[1:-1]) / 2
        # add the final bin, which is physically located just after the last regular bin but represents all higher values
        bar_width = centre_of_bar_in_x_axis[3] - centre_of_bar_in_x_axis[2]
        centre_of_bar_in_x_axis = np.append(centre_of_bar_in_x_axis, centre_of_bar_in_x_axis[-1] + bar_width)
        linecontainer_AAIMON_mean = ax.plot(centre_of_bar_in_x_axis, freq_counts_normalised, color=dfv.loc[prot_list, 'color'],
                                            alpha=alpha, linewidth=linewidth,
                                            label=dfv.loc[prot_list, 'list_description'])

        ###   normalised AAIMON   ###
        # create numpy array of membranous over nonmembranous conservation ratios (identity)
        hist_data = np.array(df_dict[prot_list]['AAIMON_mean_all_TMDs_n'])
        # use numpy to create a histogram
        freq_counts, bin_array = np.histogram(hist_data, bins=binlist)
        freq_counts_normalised = freq_counts / freq_counts.max() + offset
        # assuming all of the bins are exactly the same size, make the width of the column equal to XX% (e.g. 95%) of each bin
        col_width = float('%0.3f' % (0.95 * (bin_array[1] - bin_array[0])))
        # when align='center', the central point of the bar in the x-axis is simply the middle of the bins ((bin_0-bin_1)/2, etc)
        centre_of_bar_in_x_axis = (bin_array[:-2] + bin_array[1:-1]) / 2
        # add the final bin, which is physically located just after the last regular bin but represents all higher values
        bar_width = centre_of_bar_in_x_axis[3] - centre_of_bar_in_x_axis[2]
        centre_of_bar_in_x_axis = np.append(centre_of_bar_in_x_axis, centre_of_bar_in_x_axis[-1] + bar_width)
        linecontainer_AAIMON_mean = ax.plot(centre_of_bar_in_x_axis, freq_counts_normalised, ':', color=dfv.loc[prot_list, 'color'],
                                            alpha=alpha,
                                            linewidth=linewidth)
        offset = offset - 1

    ###############################################################
    #                                                             #
    #                       set up plot style                     #
    #                                                             #
    ###############################################################

    ax.set_xlabel('AAIMON', fontsize=fontsize)
    # move the x-axis label closer to the x-axis
    ax.xaxis.set_label_coords(0.5, -0.085)
    # x and y axes min and max
    xlim_min = 0
    xlim_max = 2
    ax.set_xlim(xlim_min, xlim_max)
    ylim_min = -0.01
    ylim_max = len(protein_lists) + 0.01
    ax.set_ylim(ylim_min, ylim_max)
    # set y-axis grid lines without tick labels
    ax.get_yaxis().set_ticks(list(np.arange(0, ylim_max, 1)))
    ax.yaxis.set_ticklabels([])
    ax.set_ylabel('relative frequency', rotation='vertical', fontsize=fontsize)
    # change axis font size
    ax.tick_params(labelsize=fontsize)
    # create legend
    ax.xaxis.set_label_coords(0.5, -0.07)
    ax.yaxis.set_label_coords(-0.005, 0.5)
    plt .xticks(np.arange(xlim_min, xlim_max + 0.1, 0.2))

    if create_legend:
        ### create legend with additional 2 elements corresponding to AAIMON and AAIMON_n ###
        # Get artists and labels for legend and chose which ones to display
        handles, labels = ax.get_legend_handles_labels()
        display = (list(range(0, len(protein_lists) + 1, 1)))
        # Create custom artists
        AAIMON = plt.Line2D((0, 1), (0, 0), color='k', linewidth=linewidth)
        AAIMON_norm = plt.Line2D((0, 1), (0, 0), color='k', linestyle=':', linewidth=linewidth)
        # Create legend from custom artist/label lists
        ax.legend([handle for i, handle in enumerate(handles) if i in display] + [AAIMON, AAIMON_norm],
                  [label for i, label in enumerate(labels) if i in display] + ['AAIMON', 'AAIMON norm.'],
                  fontsize=fontsize - 3, frameon=True, bbox_to_anchor=(1.07, 1.12))
    else:
        # Create custom artists
        AAIMON = plt.Line2D((0, 1), (0, 0), color='k', linewidth=linewidth)
        AAIMON_norm = plt.Line2D((0, 1), (0, 0), color='k', linestyle=':', linewidth=linewidth)
        # Create legend from custom artist/label lists
        ax.legend([AAIMON, AAIMON_norm],
                  ['AAIMON', 'AAIMON norm.'],
                  fontsize=fontsize - 3, frameon=True, bbox_to_anchor=(1.07, 1.12))

    utils.save_figure(fig, Fig_name, base_filepath=base_filepath, save_png=save_png, save_pdf=save_pdf)


    Fig_Nr = 2
    title = 'histograms AAIMON_slope and AAIMON_n_slope'
    Fig_name = 'Fig02_Histograms_of_mean_AAIMON_slope_and_AAIMON_n_slope'
    binlist = np.linspace(-0.04, 0.04, 61)
    fig, ax = plt.subplots()
    offset = len(protein_lists) - 1

    for prot_list in protein_lists:
        ###   non-normalised AAIMON   ###
        # create numpy array of membranous over nonmembranous conservation ratios (identity)
        hist_data = np.array(df_dict[prot_list]['AAIMON_slope_mean_all_TMDs'])
        # use numpy to create a histogram
        freq_counts, bin_array = np.histogram(hist_data, bins=binlist)
        freq_counts_normalised = freq_counts / freq_counts.max() + offset
        # assuming all of the bins are exactly the same size, make the width of the column equal to XX% (e.g. 95%) of each bin
        col_width = float('%0.3f' % (0.95 * (bin_array[1] - bin_array[0])))
        # when align='center', the central point of the bar in the x-axis is simply the middle of the bins ((bin_0-bin_1)/2, etc)
        centre_of_bar_in_x_axis = (bin_array[:-2] + bin_array[1:-1]) / 2
        # add the final bin, which is physically located just after the last regular bin but represents all higher values
        bar_width = centre_of_bar_in_x_axis[3] - centre_of_bar_in_x_axis[2]
        centre_of_bar_in_x_axis = np.append(centre_of_bar_in_x_axis, centre_of_bar_in_x_axis[-1] + bar_width)
        linecontainer_AAIMON_mean = ax.plot(centre_of_bar_in_x_axis, freq_counts_normalised, color=dfv.loc[prot_list, 'color'],
                                            alpha=alpha, linewidth=linewidth,
                                            label=dfv.loc[prot_list, 'list_description'])

        ###   normalised AAIMON   ###
        # create numpy array of membranous over nonmembranous conservation ratios (identity)
        hist_data = np.array(df_dict[prot_list]['AAIMON_n_slope_mean_all_TMDs'])
        # use numpy to create a histogram
        freq_counts, bin_array = np.histogram(hist_data, bins=binlist)
        freq_counts_normalised = freq_counts / freq_counts.max() + offset
        # assuming all of the bins are exactly the same size, make the width of the column equal to XX% (e.g. 95%) of each bin
        col_width = float('%0.3f' % (0.95 * (bin_array[1] - bin_array[0])))
        # when align='center', the central point of the bar in the x-axis is simply the middle of the bins ((bin_0-bin_1)/2, etc)
        centre_of_bar_in_x_axis = (bin_array[:-2] + bin_array[1:-1]) / 2
        # add the final bin, which is physically located just after the last regular bin but represents all higher values
        bar_width = centre_of_bar_in_x_axis[3] - centre_of_bar_in_x_axis[2]
        centre_of_bar_in_x_axis = np.append(centre_of_bar_in_x_axis, centre_of_bar_in_x_axis[-1] + bar_width)
        linecontainer_AAIMON_mean = ax.plot(centre_of_bar_in_x_axis, freq_counts_normalised, ':', color=dfv.loc[prot_list, 'color'],
                                            alpha=alpha,
                                            linewidth=linewidth)
        offset = offset - 1

    ###############################################################
    #                                                             #
    #                       set up plot style                     #
    #                                                             #
    ###############################################################

    ax.set_xlabel('AAIMON_slope', fontsize=fontsize)
    # move the x-axis label closer to the x-axis
    ax.xaxis.set_label_coords(0.5, -0.085)
    # x and y axes min and max
    xlim_min = -0.03
    xlim_max = 0.03
    ax.set_xlim(xlim_min, xlim_max)
    ylim_min = -0.01
    ylim_max = len(protein_lists) + 0.01
    ax.set_ylim(ylim_min, ylim_max)
    # set y-axis grid lines without tick labels
    ax.get_yaxis().set_ticks(list(np.arange(0, ylim_max, 1)))
    ax.yaxis.set_ticklabels([])
    ax.set_ylabel('relative frequency', rotation='vertical', fontsize=fontsize)
    # change axis font size
    ax.tick_params(labelsize=fontsize)
    # create legend
    ax.xaxis.set_label_coords(0.5, -0.07)
    ax.yaxis.set_label_coords(-0.005, 0.5)

    if create_legend:
        ### create legend with additional 2 elements corresponding to AAIMON and AAIMON_n ###
        # Get artists and labels for legend and chose which ones to display
        handles, labels = ax.get_legend_handles_labels()
        display = (list(range(0, len(protein_lists) + 1, 1)))
        # Create custom artists
        AAIMON = plt.Line2D((0, 1), (0, 0), color='k', linewidth=linewidth)
        AAIMON_norm = plt.Line2D((0, 1), (0, 0), color='k', linestyle=':', linewidth=linewidth)
        # Create legend from custom artist/label lists
        ax.legend([handle for i, handle in enumerate(handles) if i in display] + [AAIMON, AAIMON_norm],
                  [label for i, label in enumerate(labels) if i in display] + ['AAIMON_slope', 'AAIMON_slope norm.'],
                  fontsize=fontsize - 3, frameon=True, bbox_to_anchor=(1.07, 1.12))
    else:
        # Create custom artists
        AAIMON = plt.Line2D((0, 1), (0, 0), color='k', linewidth=linewidth)
        AAIMON_norm = plt.Line2D((0, 1), (0, 0), color='k', linestyle=':', linewidth=linewidth)
        # Create legend from custom artist/label lists
        ax.legend([AAIMON, AAIMON_norm],
                  ['AAIMON_slope', 'AAIMON_slope norm.'],
                  fontsize=fontsize - 3, frameon=True, bbox_to_anchor=(1.07, 1.12))


    utils.save_figure(fig, Fig_name, base_filepath, save_png, save_pdf)


    Fig_Nr = 3
    title = 'seqlen'
    Fig_name = 'Fig03_compare_seqlen'
    linspace_binlist = np.linspace(0, 2500, 100)
    binlist = np.append(linspace_binlist, 10000)
    fig, ax = plt.subplots()
    offset = len(protein_lists) - 1

    for prot_list in protein_lists:
        ###   non-normalised AAIMON   ###
        # create numpy array of membranous over nonmembranous conservation ratios (identity)
        # hist_data = np.array(df_dict[prot_list]['nonTMD_SW_align_len_excl_gaps_mean'])
        hist_data = np.array(df_dict[prot_list]['seqlen'])
        # use numpy to create a histogram
        freq_counts, bin_array = np.histogram(hist_data, bins=binlist)
        freq_counts_normalised = freq_counts / freq_counts.max() + offset
        # assuming all of the bins are exactly the same size, make the width of the column equal to XX% (e.g. 95%) of each bin
        col_width = float('%0.3f' % (0.95 * (bin_array[1] - bin_array[0])))
        # when align='center', the central point of the bar in the x-axis is simply the middle of the bins ((bin_0-bin_1)/2, etc)
        centre_of_bar_in_x_axis = (bin_array[:-2] + bin_array[1:-1]) / 2
        # add the final bin, which is physically located just after the last regular bin but represents all higher values
        # bar_width = centre_of_bar_in_x_axis[3] - centre_of_bar_in_x_axis[2]
        centre_of_bar_in_x_axis = np.append(centre_of_bar_in_x_axis, centre_of_bar_in_x_axis[-1] + 500)
        linecontainer_AAIMON_mean = ax.plot(centre_of_bar_in_x_axis, freq_counts_normalised, color=dfv.loc[prot_list, 'color'],
                                            alpha=alpha, linewidth=linewidth,
                                            label=dfv.loc[prot_list, 'list_description'])

        offset = offset - 1

    ###############################################################
    #                                                             #
    #                       set up plot style                     #
    #                                                             #
    ###############################################################

    ax.set_xlabel('sequence length', fontsize=fontsize)
    # move the x-axis label closer to the x-axis
    ax.xaxis.set_label_coords(0.5, -0.085)
    # x and y axes min and max
    xlim_min = 0
    xlim_max = 3000
    ax.set_xlim(xlim_min, xlim_max)
    ylim_min = -0.01
    ylim_max = len(protein_lists) + 0.01
    ax.set_ylim(ylim_min, ylim_max)
    # set y-axis grid lines without tick labels
    ax.get_yaxis().set_ticks(list(np.arange(0, ylim_max, 1)))
    ax.yaxis.set_ticklabels([])
    ax.set_ylabel('relative frequency', rotation='vertical', fontsize=fontsize)
    # change axis font size
    ax.tick_params(labelsize=fontsize)
    ax.xaxis.set_label_coords(0.5, -0.07)
    ax.yaxis.set_label_coords(-0.005, 0.5)

    # set custom x-ticks
    list_xticks = list(np.arange(xlim_min, xlim_max + 50, 500))
    ax.xaxis.set_ticks(list_xticks)
    # modify last x-tick with string
    list_xticks[-1] = '>2500'
    ax.set_xticklabels(list_xticks)

    if create_legend:
        # Get artists and labels for legend and chose which ones to display
        handles, labels = ax.get_legend_handles_labels()
        display = (list(range(0, len(protein_lists) + 1, 1)))
        # Create legend
        ax.legend([handle for i, handle in enumerate(handles) if i in display],
                  [label for i, label in enumerate(labels) if i in display],
                  fontsize=fontsize - 3, frameon=True, bbox_to_anchor=(1.07, 1.12))

    utils.save_figure(fig, Fig_name, base_filepath, save_png, save_pdf)

    Fig_Nr = 4
    title = 'number of TMDs'
    Fig_name = 'Fig04_comparing_number_of_TMDs'
    binlist = np.linspace(0, 25, 100)
    binlist = np.append(binlist, [50])
    fig, ax = plt.subplots()
    offset = len(protein_lists) - 1

    for prot_list in protein_lists:
        ###   non-normalised AAIMON   ###
        # create numpy array of membranous over nonmembranous conservation ratios (identity)
        # hist_data = np.array(df_dict[prot_list]['nonTMD_SW_align_len_excl_gaps_mean'])
        hist_data = np.array(df_dict[prot_list]['number_of_TMDs'])
        # use numpy to create a histogram
        freq_counts, bin_array = np.histogram(hist_data, bins=binlist)
        freq_counts_normalised = freq_counts / freq_counts.max() + offset
        # assuming all of the bins are exactly the same size, make the width of the column equal to XX% (e.g. 95%) of each bin
        col_width = float('%0.3f' % (0.95 * (bin_array[1] - bin_array[0])))
        # when align='center', the central point of the bar in the x-axis is simply the middle of the bins ((bin_0-bin_1)/2, etc)
        centre_of_bar_in_x_axis = (bin_array[:-2] + bin_array[1:-1]) / 2
        # add the final bin, which is physically located just after the last regular bin but represents all higher values
        # bar_width = centre_of_bar_in_x_axis[3] - centre_of_bar_in_x_axis[2]
        centre_of_bar_in_x_axis = np.append(centre_of_bar_in_x_axis, centre_of_bar_in_x_axis[-1] + 3)
        linecontainer_AAIMON_mean = ax.plot(centre_of_bar_in_x_axis, freq_counts_normalised, color=dfv.loc[prot_list, 'color'],
                                            alpha=alpha, linewidth=linewidth,
                                            label=dfv.loc[prot_list, 'list_description'])
        offset = offset - 1

    ###############################################################
    #                                                             #
    #                       set up plot style                     #
    #                                                             #
    ###############################################################

    ax.set_xlabel('number of TMDs', fontsize=fontsize)
    # move the x-axis label closer to the x-axis
    ax.xaxis.set_label_coords(0.5, -0.085)
    # x and y axes min and max
    xlim_min = 0
    xlim_max = 30
    ax.set_xlim(xlim_min, xlim_max)
    ylim_min = -0.01
    ylim_max = len(protein_lists) + 0.01
    ax.set_ylim(ylim_min, ylim_max)
    # set y-axis grid lines without tick labels
    ax.get_yaxis().set_ticks(list(np.arange(0, ylim_max, 1)))
    ax.yaxis.set_ticklabels([])
    ax.set_ylabel('relative frequency', rotation='vertical', fontsize=fontsize)
    # change axis font size
    ax.tick_params(labelsize=fontsize)
    ax.xaxis.set_label_coords(0.5, -0.07)
    ax.yaxis.set_label_coords(-0.005, 0.5)

    list_xticks = list(np.arange(xlim_min, xlim_max, 4))
    ax.xaxis.set_ticks(list_xticks)
    list_xticks[7] = '>24'
    ax.set_xticklabels(list_xticks)

    if create_legend:
        # Get artists and labels for legend and chose which ones to display
        handles, labels = ax.get_legend_handles_labels()
        display = (list(range(0, len(protein_lists) + 1, 1)))
        # Create legend
        ax.legend([handle for i, handle in enumerate(handles) if i in display],
                  [label for i, label in enumerate(labels) if i in display],
                  fontsize=fontsize - 3, frameon=True, bbox_to_anchor=(1.07, 1.12))

    utils.save_figure(fig, Fig_name, base_filepath, save_png, save_pdf)

    Fig_Nr = 5
    title = 'observed changes mean'
    Fig_name = 'Fig05_comparing_observed_changes_mean'
    binlist = np.linspace(0, 100, 31)
    fig, ax = plt.subplots()
    offset = len(protein_lists) - 1

    for prot_list in protein_lists:
        ###   non-normalised AAIMON   ###
        # create numpy array of membranous over nonmembranous conservation ratios (identity)
        # hist_data = np.array(df_dict[prot_list]['nonTMD_SW_align_len_excl_gaps_mean'])
        hist_data = np.array(df_dict[prot_list]['obs_changes_mean'])
        # use numpy to create a histogram
        freq_counts, bin_array = np.histogram(hist_data, bins=binlist)
        freq_counts_normalised = freq_counts / freq_counts.max() + offset
        # assuming all of the bins are exactly the same size, make the width of the column equal to XX% (e.g. 95%) of each bin
        col_width = float('%0.3f' % (0.95 * (bin_array[1] - bin_array[0])))
        # when align='center', the central point of the bar in the x-axis is simply the middle of the bins ((bin_0-bin_1)/2, etc)
        centre_of_bar_in_x_axis = (bin_array[:-2] + bin_array[1:-1]) / 2
        # add the final bin, which is physically located just after the last regular bin but represents all higher values
        bar_width = centre_of_bar_in_x_axis[3] - centre_of_bar_in_x_axis[2]
        centre_of_bar_in_x_axis = np.append(centre_of_bar_in_x_axis, centre_of_bar_in_x_axis[-1] + bar_width)
        linecontainer_AAIMON_mean = ax.plot(centre_of_bar_in_x_axis, freq_counts_normalised, color=dfv.loc[prot_list, 'color'],
                                            alpha=alpha, linewidth=linewidth,
                                            label=dfv.loc[prot_list, 'list_description'])

        offset = offset - 1

    ###############################################################
    #                                                             #
    #                       set up plot style                     #
    #                                                             #
    ###############################################################

    ax.set_xlabel('average % observed changes in homologues', fontsize=fontsize)
    # move the x-axis label closer to the x-axis
    ax.xaxis.set_label_coords(0.5, -0.085)
    # x and y axes min and max
    xlim_min = 0
    xlim_max = 100
    ax.set_xlim(xlim_min, xlim_max)
    ylim_min = -0.01
    ylim_max = len(protein_lists) + 0.01
    ax.set_ylim(ylim_min, ylim_max)
    # set y-axis grid lines without tick labels
    ax.get_yaxis().set_ticks(list(np.arange(0, ylim_max, 1)))
    ax.yaxis.set_ticklabels([])
    ax.set_ylabel('relative frequency', rotation='vertical', fontsize=fontsize)
    # change axis font size
    ax.tick_params(labelsize=fontsize)
    ax.xaxis.set_label_coords(0.5, -0.07)
    ax.yaxis.set_label_coords(-0.005, 0.5)

    # ax.xaxis.get_majorticklocs()
    list_xticks = list(np.arange(xlim_min, xlim_max + 1, 10))
    ax.xaxis.set_ticks(list_xticks)
    # list_xticks[11]='1-5000'
    # ax.set_xticklabels(list_xticks)

    if create_legend:
        # Get artists and labels for legend and chose which ones to display
        handles, labels = ax.get_legend_handles_labels()
        display = (list(range(0, len(protein_lists) + 1, 1)))
        # Create legend
        ax.legend([handle for i, handle in enumerate(handles) if i in display],
                  [label for i, label in enumerate(labels) if i in display],
                  fontsize=fontsize - 3, frameon=True, bbox_to_anchor=(1.07, 1.12))

    utils.save_figure(fig, Fig_name, base_filepath, save_png, save_pdf)

    Fig_Nr = 6
    title = 'number of homologues'
    Fig_name = 'Fig06_comparison_number_of_homologues'
    binlist = np.linspace(0, 5000, 51)
    offset = len(protein_lists) - 1

    fig, (ax, ax1) = plt.subplots(1, 2, gridspec_kw={'width_ratios': [6, 1]})

    for n, prot_list in enumerate(protein_lists):
        ax1.bar(n + 1, sum(df_dict[prot_list]['TM01_AAIMON_n_homol'] / 100000), width=0.75, align='center', color=dfv.loc[prot_list, 'color'], edgecolor='')

        ###   non-normalised AAIMON   ###
        # create numpy array of membranous over nonmembranous conservation ratios (identity)
        # hist_data = np.array(df_dict[prot_list]['nonTMD_SW_align_len_excl_gaps_mean'])
        hist_data = np.array(df_dict[prot_list]['TM01_AAIMON_n_homol'])
        # use numpy to create a histogram
        freq_counts, bin_array = np.histogram(hist_data, bins=binlist)
        freq_counts_normalised = freq_counts / freq_counts.max() + offset
        # assuming all of the bins are exactly the same size, make the width of the column equal to XX% (e.g. 95%) of each bin
        col_width = float('%0.3f' % (0.95 * (bin_array[1] - bin_array[0])))
        # when align='center', the central point of the bar in the x-axis is simply the middle of the bins ((bin_0-bin_1)/2, etc)
        centre_of_bar_in_x_axis = (bin_array[:-2] + bin_array[1:-1]) / 2
        # add the final bin, which is physically located just after the last regular bin but represents all higher values
        bar_width = centre_of_bar_in_x_axis[3] - centre_of_bar_in_x_axis[2]
        centre_of_bar_in_x_axis = np.append(centre_of_bar_in_x_axis, centre_of_bar_in_x_axis[-1] + bar_width)
        # centre_of_bar_in_x_axis [-len(linspace_binlist_2):] = linspace_binlist_2 - 50
        linecontainer_AAIMON_mean = ax.plot(centre_of_bar_in_x_axis, freq_counts_normalised, color=dfv.loc[prot_list, 'color'],
                                            alpha=alpha, linewidth=linewidth,
                                            label=dfv.loc[prot_list, 'list_description'])

        offset = offset - 1

    ###############################################################
    #                                                             #
    #                       set up plot style                     #
    #                                                             #
    ###############################################################

    ax.set_xlabel('number of homologues $*10^3$', fontsize=fontsize)
    # move the x-axis label closer to the x-axis
    ax.xaxis.set_label_coords(0.5, -0.085)
    # x and y axes min and max
    xlim_min = 0
    xlim_max = 5000
    ax.set_xlim(xlim_min, xlim_max)
    ylim_min = -0.01
    ylim_max = len(protein_lists)
    ax.set_ylim(ylim_min, ylim_max)
    # set y-axis grid lines without tick labels
    ax.get_yaxis().set_ticks(list(np.arange(0, ylim_max, 1)))
    ax.yaxis.set_ticklabels([])
    ax.set_ylabel('relative frequency', rotation='vertical', fontsize=fontsize)
    # change axis font size
    ax.tick_params(labelsize=fontsize)
    ax.xaxis.set_label_coords(0.5, -0.07)
    ax.yaxis.set_label_coords(-0.005, 0.5)

    ax.xaxis.get_majorticklocs()
    # list_xticks = list(np.arange(xlim_min, xlim_max+1 , 500))
    ax.xaxis.set_ticks(list(np.arange(xlim_min, xlim_max + 1, 500)))
    # list_xticks[-len(linspace_binlist_2):] = range(1100, 2700, 100)
    ax.set_xticklabels(list(np.arange(xlim_min, xlim_max + 1, 500) / 1000))

    ax1.yaxis.tick_right()
    ax1.set_ylabel('number of homologues in dataset $*10^5$', fontsize=fontsize)
    ax1.tick_params(labelsize=fontsize)
    ax1.yaxis.set_label_position("right")
    ax1.xaxis.set_ticks(range(1, len(df_dict) + 1, 1))
    ax1.set_xticklabels([])
    plt.subplots_adjust(wspace=0.05, hspace=0)

    if create_legend:
        # Get artists and labels for legend and chose which ones to display
        handles, labels = ax.get_legend_handles_labels()
        display = (list(range(0, len(protein_lists) + 1, 1)))
        # Create legend
        ax.legend([handle for i, handle in enumerate(handles) if i in display],
                  [label for i, label in enumerate(labels) if i in display],
                  fontsize=fontsize - 3, frameon=True)

    utils.save_figure(fig, Fig_name, base_filepath, save_png, save_pdf)


    Fig_Nr = 7
    title = 'hist_lipo_mean_all_TMDs'
    Fig_name = 'Fig07_hist_lipo_mean_all_TMDs'
    min_ = -0.5
    max_ = 0.8
    binlist = np.linspace(min_, max_, 41)
    fig, ax = plt.subplots()
    offset = len(protein_lists) - 1

    for prot_list in protein_lists:
        ###   non-normalised AAIMON   ###
        # create numpy array of membranous over nonmembranous conservation ratios (identity)
        hist_data = np.array(df_dict[prot_list]['lipo_mean_all_TMDs'])
        # use numpy to create a histogram
        freq_counts, bin_array = np.histogram(hist_data, bins=binlist)
        freq_counts_normalised = freq_counts / freq_counts.max() + offset
        # assuming all of the bins are exactly the same size, make the width of the column equal to XX% (e.g. 95%) of each bin
        col_width = float('%0.3f' % (0.95 * (bin_array[1] - bin_array[0])))
        # when align='center', the central point of the bar in the x-axis is simply the middle of the bins ((bin_0-bin_1)/2, etc)
        centre_of_bar_in_x_axis = (bin_array[:-2] + bin_array[1:-1]) / 2
        # add the final bin, which is physically located just after the last regular bin but represents all higher values
        bar_width = centre_of_bar_in_x_axis[3] - centre_of_bar_in_x_axis[2]
        centre_of_bar_in_x_axis = np.append(centre_of_bar_in_x_axis, centre_of_bar_in_x_axis[-1] + bar_width)
        linecontainer_AAIMON_mean = ax.plot(centre_of_bar_in_x_axis, freq_counts_normalised, color=dfv.loc[prot_list, 'color'],
                                            alpha=alpha, linewidth=linewidth,
                                            label=dfv.loc[prot_list, 'list_description'])

        ###   normalised AAIMON   ###
        # create numpy array of membranous over nonmembranous conservation ratios (identity)
        hist_data = np.array(df_dict[prot_list]['TM01_lipo'])
        # use numpy to create a histogram
        freq_counts, bin_array = np.histogram(hist_data, bins=binlist)
        freq_counts_normalised = freq_counts / freq_counts.max() + offset
        # assuming all of the bins are exactly the same size, make the width of the column equal to XX% (e.g. 95%) of each bin
        col_width = float('%0.3f' % (0.95 * (bin_array[1] - bin_array[0])))
        # when align='center', the central point of the bar in the x-axis is simply the middle of the bins ((bin_0-bin_1)/2, etc)
        centre_of_bar_in_x_axis = (bin_array[:-2] + bin_array[1:-1]) / 2
        # add the final bin, which is physically located just after the last regular bin but represents all higher values
        bar_width = centre_of_bar_in_x_axis[3] - centre_of_bar_in_x_axis[2]
        centre_of_bar_in_x_axis = np.append(centre_of_bar_in_x_axis, centre_of_bar_in_x_axis[-1] + bar_width)
        linecontainer_AAIMON_mean = ax.plot(centre_of_bar_in_x_axis, freq_counts_normalised, ':', color=dfv.loc[prot_list, 'color'],
                                            alpha=alpha,
                                            linewidth=linewidth)


        offset = offset - 1

    ###############################################################
    #                                                             #
    #                       set up plot style                     #
    #                                                             #
    ###############################################################

    ax.set_xlabel('lipophilicity (Hessa scale)', fontsize=fontsize)
    # move the x-axis label closer to the x-axis
    ax.xaxis.set_label_coords(0.5, -0.085)
    # x and y axes min and max
    xlim_min = min_
    xlim_max = max_
    ax.set_xlim(xlim_min, xlim_max)
    ylim_min = -0.01
    ylim_max = len(protein_lists) + 0.01
    ax.set_ylim(ylim_min, ylim_max)
    # set y-axis grid lines without tick labels
    ax.get_yaxis().set_ticks(list(np.arange(0, ylim_max, 1)))
    ax.yaxis.set_ticklabels([])
    ax.set_ylabel('relative frequency', rotation='vertical', fontsize=fontsize)
    # change axis font size
    ax.tick_params(labelsize=fontsize)
    # create legend
    ax.xaxis.set_label_coords(0.5, -0.07)
    ax.yaxis.set_label_coords(-0.005, 0.5)
    plt.xticks(np.arange(xlim_min, xlim_max + 0.1, 0.2))

    # add annotations
    ax.annotate(s="more lipophilic", xy=(0, -0.08), fontsize=fontsize, xytext=None, xycoords='axes fraction')
    ax.annotate(s="less lipophilic", xy=(1.0, -0.08), fontsize=fontsize, xytext=None, horizontalalignment='right', xycoords='axes fraction')

    if create_legend:
        ### create legend with additional 2 elements corresponding to AAIMON and AAIMON_n ###
        # Get artists and labels for legend and chose which ones to display
        handles, labels = ax.get_legend_handles_labels()
        display = (list(range(0, len(protein_lists)+1, 1)))
        # Create custom artists
        mean_ = plt.Line2D((0, 1), (0, 0), color='k', linewidth=linewidth)
        TM01_ = plt.Line2D((0, 1), (0, 0), color='k', linestyle=':', linewidth=linewidth)
        # Create legend from custom artist/label lists
        ax.legend([handle for i, handle in enumerate(handles) if i in display] + [mean_, TM01_],
                  [label for i, label in enumerate(labels) if i in display] + ['mean all TMDs', 'TM01'],
                  fontsize=fontsize-3, frameon=True, bbox_to_anchor=(1.07, 1.12))
    else:
        # Create custom artists
        mean_ = plt.Line2D((0, 1), (0, 0), color='k', linewidth=linewidth)
        TM01_ = plt.Line2D((0, 1), (0, 0), color='k', linestyle=':', linewidth=linewidth)
        # Create legend from custom artist/label lists
        ax.legend([mean_, TM01_],['mean all TMDs', 'TM01'],
                  fontsize=fontsize-3, frameon=True, bbox_to_anchor=(1.07, 1.12))

    utils.save_figure(fig, Fig_name, base_filepath=base_filepath, save_png=save_png, save_pdf=save_pdf)


    Fig_Nr = 8
    title = 'perc of TMD region in protein'
    Fig_name = 'Fig08_perc_of_TMD_region_in_protein'
    min_ = 0
    max_ = 100
    binlist = np.linspace(min_, max_, 101)
    fig, ax = plt.subplots()
    offset = len(protein_lists) - 1

    for prot_list in protein_lists:
        ###   non-normalised AAIMON   ###
        # create numpy array of membranous over nonmembranous conservation ratios (identity)
        hist_data = np.array(df_dict[prot_list]['perc_TMD'])
        # use numpy to create a histogram
        freq_counts, bin_array = np.histogram(hist_data, bins=binlist)
        freq_counts_normalised = freq_counts / freq_counts.max() + offset
        # assuming all of the bins are exactly the same size, make the width of the column equal to XX% (e.g. 95%) of each bin
        col_width = float('%0.3f' % (0.95 * (bin_array[1] - bin_array[0])))
        # when align='center', the central point of the bar in the x-axis is simply the middle of the bins ((bin_0-bin_1)/2, etc)
        centre_of_bar_in_x_axis = (bin_array[:-2] + bin_array[1:-1]) / 2
        # add the final bin, which is physically located just after the last regular bin but represents all higher values
        bar_width = centre_of_bar_in_x_axis[3] - centre_of_bar_in_x_axis[2]
        centre_of_bar_in_x_axis = np.append(centre_of_bar_in_x_axis, centre_of_bar_in_x_axis[-1] + bar_width)
        linecontainer_AAIMON_mean = ax.plot(centre_of_bar_in_x_axis, freq_counts_normalised, color=dfv.loc[prot_list, 'color'],
                                            alpha=alpha, linewidth=linewidth,
                                            label=dfv.loc[prot_list, 'list_description'])

        offset = offset - 1

    ###############################################################
    #                                                             #
    #                       set up plot style                     #
    #                                                             #
    ###############################################################

    ax.set_xlabel('% residues within TM region', fontsize=fontsize)
    # move the x-axis label closer to the x-axis
    ax.xaxis.set_label_coords(0.5, -0.085)
    # x and y axes min and max
    xlim_min = min_
    xlim_max = max_
    ax.set_xlim(xlim_min, xlim_max)
    ylim_min = -0.01
    ylim_max = len(protein_lists) + 0.01
    ax.set_ylim(ylim_min, ylim_max)
    # set y-axis grid lines without tick labels
    ax.get_yaxis().set_ticks(list(np.arange(0, ylim_max, 1)))
    ax.yaxis.set_ticklabels([])
    ax.set_ylabel('relative frequency', rotation='vertical', fontsize=fontsize)
    # change axis font size
    ax.tick_params(labelsize=fontsize)
    # create legend
    ax.xaxis.set_label_coords(0.5, -0.07)
    ax.yaxis.set_label_coords(-0.005, 0.5)
    plt.xticks(np.arange(xlim_min, xlim_max + 0.1, 10))

    if create_legend:
        # Get artists and labels for legend and chose which ones to display
        handles, labels = ax.get_legend_handles_labels()
        display = (list(range(0, len(protein_lists)+1, 1)))
        # Create legend from custom artist/label lists
        ax.legend([handle for i, handle in enumerate(handles) if i in display],
                  [label for i, label in enumerate(labels) if i in display],
                  fontsize=fontsize-3, frameon=True, bbox_to_anchor=(1.07, 1.12))

    utils.save_figure(fig, Fig_name, base_filepath=base_filepath, save_png=save_png, save_pdf=save_pdf)


    #dfv.to_csv(os.path.join(base_filepath, 'Lists_%s_variables.csv'%str_protein_lists))

    sys.stdout.write("\n~~~~~~~~~~~~         compare_lists finished           ~~~~~~~~~~~~\n")



