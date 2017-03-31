import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import ast
import sys
import csv
import datetime
import korbinian.utils as utils
import korbinian
import scipy


def compare_lists (s):
    """ Create figures comparing protein lists.

    Compares a number of parameters. Some from UniProt. Others from the
    conservation ratio analysis. Analysis of up to 7 lists is possible. Only tested for 3.

    Histograms of AAIMON ratios
    Histograms of AAIMON slope
    Full sequence length
    Length of nonTMD
    Number of homologues
    etc.

    Saves png and/or pdf files for each image. The first 3 colours are optimised for printing.

    Creates a new folder each time, unless "testing_mode" is activated in the settings file.

    Files are saved under the Databases folder, in a (new) subfolder compare_lists.

    To run: replace the list number with "compare"
    Put a python stringlist of the desired list numbers in protein_lists in the excel settings file (e.g. "[1,2,3]")

    Parameters
    ----------
    s : dict
        settings dictionary

    """
    sys.stdout.write("\n\n~~~~~~~~~~~~         starting compare_lists           ~~~~~~~~~~~~\n\n")
    create_legend = True
    save_png = s['save_png']
    save_pdf = s['save_pdf']

    ''' problems with RGB values in color_list when comparing less than 3 lists, colours were converted to HTML format'''
    #colour_dict = utils.create_colour_lists()
    #start with TUM colours, distinguishable in B&W
    #color_list = ["0.5", colour_dict["TUM_oranges"]['TUM1'], colour_dict["TUM_colours"]['TUM5'], (0.843, 0.098, 0.1098)]
    # colours for List 1, 2 and 3
    color_list = ['#808080', '#D59460', '#005293']
    # add other HTML colours
    color_list = color_list + ['#A1B11A', '#9ECEEC', '#0076B8', '#454545']

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
        dfx = pd.read_csv(dfv.loc[prot_list, 'base_filename_lists'], sep=",", quoting=csv.QUOTE_NONNUMERIC, index_col=0, low_memory=False)
        dfy = pd.read_csv(dfv.loc[prot_list, 'base_filename_cr_summaries'], sep=",", quoting=csv.QUOTE_NONNUMERIC, index_col=0, low_memory=False)
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
    fontsize = 12
    linewidth = 2


    #--------------------------------------------------------------------------------------------------------------------------------#
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
                  fontsize=fontsize - 3, frameon=True)#, bbox_to_anchor=(1.07, 1.12))
    else:
        # Create custom artists
        AAIMON = plt.Line2D((0, 1), (0, 0), color='k', linewidth=linewidth)
        AAIMON_norm = plt.Line2D((0, 1), (0, 0), color='k', linestyle=':', linewidth=linewidth)
        # Create legend from custom artist/label lists
        ax.legend([AAIMON, AAIMON_norm],
                  ['AAIMON', 'AAIMON norm.'],
                  fontsize=fontsize - 3, frameon=True)#, bbox_to_anchor=(1.07, 1.12))

    utils.save_figure(fig, Fig_name, base_filepath=base_filepath, save_png=save_png, save_pdf=save_pdf)

    #--------------------------------------------------------------------------------------------------------------------------------#
    Fig_Nr = 2
    title = 'histograms AAIMON_slope and AAIMON_n_slope'
    Fig_name = 'Fig02_Histograms_of_mean_AAIMON_slope_and_AAIMON_n_slope'
    binlist = np.linspace(-40, 40, 61)
    fig, ax = plt.subplots()
    offset = len(protein_lists) - 1

    for prot_list in protein_lists:
        ###   non-normalised AAIMON   ###
        # create numpy array of membranous over nonmembranous conservation ratios (identity)
        hist_data = np.array(df_dict[prot_list]['AAIMON_slope_mean_all_TMDs']*1000)
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
        hist_data = np.array(df_dict[prot_list]['AAIMON_n_slope_mean_all_TMDs']*1000)
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

    ax.set_xlabel(r'm$_{\rm TM/nonTM} *10^{\rm -3}$', fontsize=fontsize, labelpad=20)
    # move the x-axis label closer to the x-axis
    ax.xaxis.set_label_coords(0.5, -0.085)
    # x and y axes min and max
    xlim_min = -30
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
    # create legend
    ax.xaxis.set_label_coords(0.5, -0.07)
    ax.yaxis.set_label_coords(-0.005, 0.5)

    # add annotations
    ax.annotate(s="TM less conserved", xy=(0, -0.09), fontsize=fontsize, xytext=None, xycoords='axes fraction')
    ax.annotate(s="TM more conserved", xy=(1.0, -0.09), fontsize=fontsize, xytext=None, horizontalalignment='right', xycoords='axes fraction')

    create_legend_for_AAIMON_slope = False
    if create_legend or create_legend_for_AAIMON_slope:
        ### create legend with additional 2 elements corresponding to AAIMON and AAIMON_n ###
        # Get artists and labels for legend and chose which ones to display
        handles, labels = ax.get_legend_handles_labels()
        display = (list(range(0, len(protein_lists) + 1, 1)))
        # Create custom artists
        AAIMON = plt.Line2D((0, 1), (0, 0), color='k', linewidth=linewidth)
        AAIMON_norm = plt.Line2D((0, 1), (0, 0), color='k', linestyle=':', linewidth=linewidth)
        # Create legend from custom artist/label lists
        ax.legend([handle for i, handle in enumerate(handles) if i in display] + [AAIMON, AAIMON_norm],
                  [label for i, label in enumerate(labels) if i in display] + ['AAIMON slope', 'AAIMON slope norm.'],
                  fontsize=fontsize - 3, frameon=True, loc='upper right')#, bbox_to_anchor=(1.07, 1.12))
    else:
        # Create custom artists
        AAIMON = plt.Line2D((0, 1), (0, 0), color='k', linewidth=linewidth)
        AAIMON_norm = plt.Line2D((0, 1), (0, 0), color='k', linestyle=':', linewidth=linewidth)
        # Create legend from custom artist/label lists
        ax.legend([AAIMON, AAIMON_norm],
                  ['AAIMON slope', 'AAIMON slope norm.'],
                  fontsize=fontsize - 3, frameon=True, loc='upper right')#, bbox_to_anchor=(1.07, 1.12))
    #plt.gcf().subplots_adjust(bottom=0.15)
    utils.save_figure(fig, Fig_name, base_filepath, save_png, save_pdf)


    #--------------------------------------------------------------------------------------------------------------------------------#
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

    #--------------------------------------------------------------------------------------------------------------------------------#
    Fig_Nr = 4
    title = 'number of TMDs'
    Fig_name = 'Fig04_compare_number_of_TMDs'
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

    #--------------------------------------------------------------------------------------------------------------------------------#
    Fig_Nr = 5
    title = 'observed changes mean'
    Fig_name = 'Fig05_compare_observed_changes_mean'
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

    #--------------------------------------------------------------------------------------------------------------------------------#
    Fig_Nr = 6
    title = 'number of homologues'
    Fig_name = 'Fig06_comparison_number_of_homologues'
    #binlist = np.linspace(0, 5000, 101)
    offset = len(protein_lists) - 1

    # alternative binlist
    binlist_1 = np.linspace(0, 1000, 41)
    binlist_2 = np.linspace(1000, 5000, 81)
    binlist = np.append(binlist_1, binlist_2[1:])

    fig, (ax, ax1) = plt.subplots(1, 2, gridspec_kw={'width_ratios': [6, 1]})

    for n, prot_list in enumerate(protein_lists):
        n_of_homologues = sum(df_dict[prot_list]['TM01_AAIMON_n_homol'])
        sys.stdout.write('number of homologues {}: {}\n'.format(dfv.loc[prot_list, 'list_description'], n_of_homologues))
        ax1.bar(n + 1, n_of_homologues / 100000, width=0.75, align='center', color=dfv.loc[prot_list, 'color'], edgecolor='')

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
        centre_of_bar_in_x_axis[-80:] = np.linspace(10, 800, 80) + 1000
        linecontainer_AAIMON_mean = ax.plot(centre_of_bar_in_x_axis, freq_counts_normalised, color=dfv.loc[prot_list, 'color'],
                                            alpha=alpha, linewidth=linewidth,
                                            label=dfv.loc[prot_list, 'list_description'])

        offset = offset - 1

    ###############################################################
    #                                                             #
    #                       set up plot style                     #
    #                                                             #
    ###############################################################

    ax.set_xlabel('number of homologues $*10^3$ (unequal bins)', fontsize=fontsize)
    # move the x-axis label closer to the x-axis
    ax.xaxis.set_label_coords(0.5, -0.085)
    # x and y axes min and max
    xlim_min = 0
    xlim_max = 1800
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
    list_xticks = np.arange(xlim_min, xlim_max + 1, 200)
    list_xticks[-4:] = range(2000, 5001, 1000)
    ax.set_xticklabels((list_xticks) / 1000)

    ax1.yaxis.tick_right()
    ax1.set_ylabel('number of homologues in dataset $*10^5$', fontsize=fontsize)
    ax1.tick_params(labelsize=fontsize)
    ax1.yaxis.set_label_position("right")
    # ax1.xaxis.set_ticks(range(1, len(df_dict) + 1, 1))
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


    #--------------------------------------------------------------------------------------------------------------------------------#
    Fig_Nr = 7
    title = 'hist_lipo_mean_all_TMDs'
    Fig_name = 'Fig07_hist_lipo_mean_all_TMDs'
    min_ = -0.5
    max_ = 0.8
    binlist = np.linspace(min_, max_, 26) #41
    fig, ax = plt.subplots()
    offset = len(protein_lists) - 1

    for prot_list in protein_lists:
        ###   lipo mean all TMDs excl. TM01   ###
        # create numpy array of membranous over nonmembranous conservation ratios (identity)
        hist_data = np.array(df_dict[prot_list]['lipo_mean_excl_TM01'].dropna())
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

        ###   TM01 lipo   ###
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

        # ###   last TM lipo   ###
        # # create numpy array of membranous over nonmembranous conservation ratios (identity)
        # hist_data = np.array(df_dict[prot_list]['lipo_last_TMD'].dropna())
        # # use numpy to create a histogram
        # freq_counts, bin_array = np.histogram(hist_data, bins=binlist)
        # freq_counts_normalised = freq_counts / freq_counts.max() + offset
        # # assuming all of the bins are exactly the same size, make the width of the column equal to XX% (e.g. 95%) of each bin
        # col_width = float('%0.3f' % (0.95 * (bin_array[1] - bin_array[0])))
        # # when align='center', the central point of the bar in the x-axis is simply the middle of the bins ((bin_0-bin_1)/2, etc)
        # centre_of_bar_in_x_axis = (bin_array[:-2] + bin_array[1:-1]) / 2
        # # add the final bin, which is physically located just after the last regular bin but represents all higher values
        # bar_width = centre_of_bar_in_x_axis[3] - centre_of_bar_in_x_axis[2]
        # centre_of_bar_in_x_axis = np.append(centre_of_bar_in_x_axis, centre_of_bar_in_x_axis[-1] + bar_width)
        # linecontainer_AAIMON_mean = ax.plot(centre_of_bar_in_x_axis, freq_counts_normalised, '-.', color=dfv.loc[prot_list, 'color'],
        #                                     alpha=alpha,
        #                                     linewidth=linewidth)

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
    ax.annotate(s="more lipophilic", xy=(0, -0.1), fontsize=fontsize, xytext=None, xycoords='axes fraction')
    ax.annotate(s="less lipophilic", xy=(1.0, -0.1), fontsize=fontsize, xytext=None, horizontalalignment='right', xycoords='axes fraction')

    if create_legend:
        ### create legend with additional 2 elements corresponding to AAIMON and AAIMON_n ###
        # Get artists and labels for legend and chose which ones to display
        handles, labels = ax.get_legend_handles_labels()
        display = (list(range(0, len(protein_lists)+1, 1)))
        # Create custom artists
        excl_TM01 = plt.Line2D((0, 1), (0, 0), color='k', linewidth=linewidth)
        TM01 = plt.Line2D((0, 1), (0, 0), color='k', linestyle=':', linewidth=linewidth)
        # last_TMD = plt.Line2D((0, 1), (0, 0), color='k', linestyle='-.', linewidth=linewidth)
        # Create legend from custom artist/label lists
        ax.legend([handle for i, handle in enumerate(handles) if i in display] + [excl_TM01, TM01],
                  [label for i, label in enumerate(labels) if i in display] + ['mean excl. TM01', 'TM01'],
                  fontsize=fontsize-3, frameon=True, loc='upper right')#bbox_to_anchor=(1.07, 1.12))
    else:
        # Create custom artists
        excl_TM01 = plt.Line2D((0, 1), (0, 0), color='k', linewidth=linewidth)
        TM01 = plt.Line2D((0, 1), (0, 0), color='k', linestyle=':', linewidth=linewidth)
        # last_TMD = plt.Line2D((0, 1), (0, 0), color='k', linestyle='-.', linewidth=linewidth)
        # Create legend from custom artist/label lists
        ax.legend([excl_TM01, TM01],['mean excl. TM01', 'TM01'],
                  fontsize=fontsize-3, frameon=True, loc='upper right')#bbox_to_anchor=(1.07, 1.12))

    utils.save_figure(fig, Fig_name, base_filepath=base_filepath, save_png=save_png, save_pdf=save_pdf)


    #--------------------------------------------------------------------------------------------------------------------------------#
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

    #--------------------------------------------------------------------------------------------------------------------------------#
    Fig_Nr = 9
    title = 'length of TMD'
    Fig_name = 'Fig09_compare_length_of_TMD'
    binlist = np.linspace(0, 40, 161)
    fig, ax = plt.subplots()
    offset = len(protein_lists) - 1

    detailled_data = True

    # sys.stdout.write('\ncalculating mean length of TMD regions: - add this stuff to prepare_protein_list ?!?')
    # for prot_list in df_dict.keys():
    #     sys.stdout.write('\nprotein list: {}\n'.format(prot_list))
    #     for n, acc in enumerate(df_dict[prot_list].index):
    #         if n % 100 == 0:
    #             sys.stdout.write('. '), sys.stdout.flush()
    #         list_of_TMDs = ast.literal_eval(df_dict[prot_list].loc[acc, 'list_of_TMDs'])
    #         list_TMD_lenghts = []
    #         for TMD in list_of_TMDs:
    #             list_TMD_lenghts.append(len(df_dict[prot_list].loc[acc, '%s_seq' % TMD]))
    #         df_dict[prot_list].loc[acc, 'TMD_len_mean'] = np.mean(list_TMD_lenghts)

    for prot_list in protein_lists:
        # collect length data for all TMDs individually
        data = df_dict[prot_list]
        list_len_TMDs = []
        for n in range(1, data.number_of_TMDs_excl_SP.max().astype(int) + 1):
            TMD = 'TM%02d' % n
            list_len_TMDs.extend(data['%s_seqlen' % TMD].dropna())
        hist_data = np.array(list_len_TMDs)
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

    ax.set_xlabel('length of TMD regions', fontsize=fontsize)
    # move the x-axis label closer to the x-axis
    ax.xaxis.set_label_coords(0.5, -0.085)
    # x and y axes min and max
    xlim_min = 0
    xlim_max = 40
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
    ax.xaxis.set_label_coords(0.5, -0.08)
    ax.yaxis.set_label_coords(-0.005, 0.5)

    ax.xaxis.get_majorticklocs()
    list_xticks = np.arange(xlim_min, xlim_max + 1, 5)
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

    utils.save_figure(fig, Fig_name, base_filepath=base_filepath, save_png=save_png, save_pdf=save_pdf)

    #--------------------------------------------------------------------------------------------------------------------------------#
    Fig_Nr = 10
    title = 'norm. factors'
    Fig_name = 'Fig10_norm_factor_pairwise'
    fig, ax = plt.subplots()
    # define distance before randTM to stop drawing curves
    norm_fig_fontsize = 14


    for prot_list in protein_lists:
        ########################################################################################
        #                                                                                      #
        #                      get random identity from csv file                               #
        #                                                                                      #
        ########################################################################################

        rand_TM_path = os.path.join(s["data_dir"], "summaries", "{ln:02d}/List{ln:02d}_rand/List{ln:02d}_rand_TM.csv".format(ln=prot_list))
        rand_nonTM_path = os.path.join(s["data_dir"], "summaries", "{ln:02d}/List{ln:02d}_rand/List{ln:02d}_rand_nonTM.csv".format(ln=prot_list))
        serTM = pd.Series.from_csv(rand_TM_path, sep="\t")
        rand_TM = serTM["random_sequence_identity_output"]
        sernonTM = pd.Series.from_csv(rand_nonTM_path, sep="\t")
        rand_nonTM = sernonTM["random_sequence_identity_output"]

        # calculate the average percentage of residues within the TM region
        fraction_TM_residues = df_dict[prot_list]['perc_TMD'].mean() / 100
        sys.stdout.write("prot_list : {:02d}, fraction_TM_residues : {:0.03f}".format(prot_list, fraction_TM_residues))

        #create a new dataframe to hold the data for percentage identity and the normalisation factors
        df_norm = pd.DataFrame()
        # create datapoints for line, percent amino acid substitutions
        df_norm["perc_aa_sub"] = np.linspace(0,1.2,1000)
        # get the percentage identity
        df_norm["perc_aa_ident"] = 1 - df_norm["perc_aa_sub"]
        # calculate the norm factors
        df_norm["norm"] = df_norm["perc_aa_ident"].apply(korbinian.cons_ratio.norm.calc_AAIMON_aa_prop_norm_factor, args=(rand_TM, rand_nonTM, fraction_TM_residues))
        # get rid of any very high or low values
        df_norm = df_norm.loc[df_norm["norm"] > -10]
        df_norm = df_norm.loc[df_norm["norm"] < 10]

        # HYPERBOLIC GRAPH NEEDS TO BE SPLIT INTO THE TWO SIDES TO AVOID UGLY CROSSOVER EFFECT

        # get the indices for the highest and lowest norm factors
        index_max = df_norm["norm"].idxmax()
        index_min = df_norm["norm"].idxmin()

        label = dfv.loc[prot_list, 'list_description']

        # plot that section of the graph, from beginning to where values are high
        # convert to array to avoid auto legend)
        ax.plot(df_norm["perc_aa_sub"][:index_max]*100, np.array(df_norm["norm"][:index_max]), color=dfv.loc[prot_list, 'color'],
                alpha=alpha, linewidth=linewidth,
                label=None)
        # plot that section of the graph, from most negative value until the end
        ax.plot(df_norm["perc_aa_sub"][index_min:]*100, np.array(df_norm["norm"][index_min:]), color=dfv.loc[prot_list, 'color'],
                alpha=alpha, linewidth=linewidth, label=label)

    ###############################################################
    #                                                             #
    #                       set up plot style                     #
    #                                                             #
    ###############################################################

    ax.set_xlabel('% AA substitutions', fontsize=norm_fig_fontsize)
    # x and y axes min and max
    xlim_min = 0
    xlim_max = 120
    ax.set_xlim(xlim_min, xlim_max)
    ylim_min = 0
    ylim_max = 3
    ax.set_ylim(ylim_min, ylim_max)
    ax.set_ylabel('normalisation factor', rotation='vertical', fontsize=norm_fig_fontsize)
    # change axis font size
    ax.tick_params(labelsize=norm_fig_fontsize)
    ax.xaxis.set_label_coords(0.5, -0.08)

    # Get artists and labels for legend and chose which ones to display
    handles, labels = ax.get_legend_handles_labels()
    display = (list(range(0, len(protein_lists) + 1, 1)))
    # Create legend
    ax.legend(fontsize=norm_fig_fontsize, frameon=True, loc='upper left')

    utils.save_figure(fig, Fig_name, base_filepath=base_filepath, save_png=save_png, save_pdf=save_pdf)

    #--------------------------------------------------------------------------------------------------------------------------------#
    Fig_Nr = 11
    title = 'compare_random_identity'
    Fig_name = 'Fig11_random_ident'
    fig, ax = plt.subplots()

    df_rand = pd.DataFrame()
    df_aa_prop_TM = pd.DataFrame()
    df_aa_prop_nonTM = pd.DataFrame()

    for prot_list in protein_lists:
        rand_TM_path = os.path.join(s["data_dir"], "summaries", "{ln:02d}/List{ln:02d}_rand/List{ln:02d}_rand_TM.csv".format(ln=prot_list))
        rand_nonTM_path = os.path.join(s["data_dir"], "summaries", "{ln:02d}/List{ln:02d}_rand/List{ln:02d}_rand_nonTM.csv".format(ln=prot_list))

        # grab the random identity and aa propensity series for the TM
        ser_rand_TM = pd.Series.from_csv(rand_TM_path, sep="\t")
        # the first line is the random identity. Extract and delete.
        df_rand.loc[prot_list, "rand_TM"] = ser_rand_TM["random_sequence_identity_output"]
        del ser_rand_TM["random_sequence_identity_output"]
        # the rest of the index is A_input, etc. Delete the "_input" and the aa propensity table is ready
        ser_rand_TM.index = ser_rand_TM.index.str[:1]
        #df_rand.loc[prot_list, "TM_aa_prop"] = str(ser_rand_TM.to_dict())
        df_aa_prop_TM[prot_list] = ser_rand_TM

        # grab the random identity and aa propensity series for the nonTM
        ser_rand_nonTM = pd.Series.from_csv(rand_nonTM_path, sep="\t")
        # the first line is the random identity. Extract and delete.
        df_rand.loc[prot_list, "rand_nonTM"] = ser_rand_nonTM["random_sequence_identity_output"]
        del ser_rand_nonTM["random_sequence_identity_output"]
        # the rest of the index is A_input, etc. Delete the "_input" and the aa propensity table is ready
        ser_rand_nonTM.index = ser_rand_nonTM.index.str[:1]
        #df_rand.loc[prot_list, "nonTM_aa_prop"] = str(ser_rand_TM.to_dict())
        df_aa_prop_nonTM[prot_list] = ser_rand_nonTM

    # create list of list_numbers for x-axis
    x = df_rand.index.astype(str)
    # width of bar
    width = 0.3
    # indices (position of bar) are at 0.1, 1.1, 2.1
    ind = np.arange(len(x)) + 0.1
    # create bar for the rand_TM
    ax.bar(ind, df_rand.rand_TM, width, color=color_list[:len(x)], alpha=1, label="TM")#color=color_list[2]
    # create bar for the rand_nonTM
    ax.bar(ind + width, df_rand.rand_nonTM, width, color=color_list[:len(x)], alpha=0.5, label="nonTM")
    # put labels in between bars
    ax.set_xticks(ind + width)
    # extract list names
    x_labels = dfv.loc[protein_lists, 'list_description']
    ax.set_xticklabels(x_labels)
    ax.legend()
    ax.set_ylabel("random identity")
    ax.set_xlabel("list number")
    utils.save_figure(fig, Fig_name, base_filepath=base_filepath, save_png=save_png, save_pdf=save_pdf)

    #--------------------------------------------------------------------------------------------------------------------------------#
    Fig_Nr = 12
    title = 'compare_aa_propensity'
    Fig_name = 'Fig12_aa_propensity_TM'
    fig, ax = plt.subplots()
    n = len(protein_lists)
    # width of bar
    width = 0.8 / n
    # indices (position of bar)
    ind = np.arange(df_aa_prop_TM.shape[0]) + 0.1

    # iterate through the protein lists
    for m, prot_list in enumerate(df_aa_prop_TM.columns):
        #indices are 1.1, 2.1, etc until 20.1 for all amino acids
        ind_for_list = ind + width * m
        # create bar for that protein list
        ax.bar(ind_for_list, df_aa_prop_TM.loc[:,prot_list], width, color=color_list[m], label=dfv.loc[prot_list, 'list_description'])#color=color_list[2]

    # xticks are located in the middle of all the bars (depending on number of lists)
    xtick_locations = ind + n/2*width
    ax.set_xticks(xtick_locations)
    aa_list = df_aa_prop_TM.index
    ax.set_xticklabels(aa_list)
    ax.set_ylabel("amino acid propensity")
    ax.legend()
    utils.save_figure(fig, Fig_name, base_filepath=base_filepath, save_png=save_png, save_pdf=save_pdf)

    #--------------------------------------------------------------------------------------------------------------------------------#
    Fig_Nr = 13
    title = 'compare_aa_propensity'
    Fig_name = 'Fig13_aa_propensity_nonTM'
    fig, ax = plt.subplots()
    n = len(protein_lists)
    # width of bar
    width = 0.8 / n
    # indices (position of bar)
    ind = np.arange(df_aa_prop_nonTM.shape[0]) + 0.1

    aa_list = df_aa_prop_nonTM.index

    for m, prot_list in enumerate(df_aa_prop_nonTM.columns):
        ind_for_list = ind + width * m
        ax.bar(ind_for_list, df_aa_prop_nonTM.loc[:,prot_list], width, color=color_list[m], label=dfv.loc[prot_list, 'list_description'])#color=color_list[2]

    xtick_locations = ind + n/2*width
    ax.set_xticks(xtick_locations)
    ax.set_xticklabels(aa_list)
    ax.set_ylabel("amino acid propensity")
    ax.legend()
    utils.save_figure(fig, Fig_name, base_filepath=base_filepath, save_png=save_png, save_pdf=save_pdf)
    #--------------------------------------------------------------------------------------------------------------------------------#

    Fig_Nr = 14
    title = 'hist_lipo_TM01_vs_last_TM'
    Fig_name = 'Fig14_hist_lipo_TM01_vs_last_TM'
    min_ = -0.5
    max_ = 0.8
    binlist = np.linspace(min_, max_, 41)
    fig, ax = plt.subplots()
    offset = len(protein_lists) - 1

    for prot_list in protein_lists:
        ###   lipo TM01   ###
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
        linecontainer_AAIMON_mean = ax.plot(centre_of_bar_in_x_axis, freq_counts_normalised, color=dfv.loc[prot_list, 'color'],
                                            alpha=alpha, linewidth=linewidth,
                                            label=dfv.loc[prot_list, 'list_description'])

        ###   TM01 lipo   ###
        # create numpy array of membranous over nonmembranous conservation ratios (identity)
        hist_data = np.array(df_dict[prot_list]['lipo_last_TMD'])
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
        display = (list(range(0, len(protein_lists) + 1, 1)))
        # Create custom artists
        TM01 = plt.Line2D((0, 1), (0, 0), color='k', linewidth=linewidth)
        last_TM = plt.Line2D((0, 1), (0, 0), color='k', linestyle=':', linewidth=linewidth)
        # Create legend from custom artist/label lists
        ax.legend([handle for i, handle in enumerate(handles) if i in display] + [TM01, last_TM],
                  [label for i, label in enumerate(labels) if i in display] + ['TM01', 'last TM'],
                  fontsize=fontsize - 3, frameon=True, loc='upper right')  # bbox_to_anchor=(1.07, 1.12))
    else:
        # Create custom artists
        TM01 = plt.Line2D((0, 1), (0, 0), color='k', linewidth=linewidth)
        last_TM = plt.Line2D((0, 1), (0, 0), color='k', linestyle=':', linewidth=linewidth)
        # Create legend from custom artist/label lists
        ax.legend([TM01, last_TM], ['TM01', 'last TM'],
                  fontsize=fontsize - 3, frameon=True, loc='upper right')  # bbox_to_anchor=(1.07, 1.12))

    utils.save_figure(fig, Fig_name, base_filepath=base_filepath, save_png=save_png, save_pdf=save_pdf)

    # --------------------------------------------------------------------------------------------------------------------------------#
    Fig_Nr = 15
    title = 'histograms AAIMON_slope TM01 vs last TMD'
    Fig_name = 'Fig15_Histograms_AAIMON_slope_TM01_vs_lastTM'
    binlist = np.linspace(-40, 40, 61)
    fig, ax = plt.subplots()
    offset = len(protein_lists) - 1

    for prot_list in protein_lists:
        ###### for backwards compatibility ##### can be removed if all data is re-processed after march 5 2017
        if not 'AAIMON_slope_last_TMD' in df_dict[prot_list].columns:
            sys.stdout.write('List{:02d}: AAIMON_slope_last_TMD not in dataframe -> older version of data, re-run "gather_AAIMON_ratios"; adding data for figure'.format(prot_list))
            for n, acc in enumerate(df_dict[prot_list].index):
                if n % 200 == 0:
                    sys.stdout.write(". "), sys.stdout.flush()
                last_TMD = df_dict[prot_list].loc[acc, 'last_TMD']
                df_dict[prot_list].loc[acc, 'AAIMON_slope_last_TMD'] = df_dict[prot_list].loc[acc, '%s_AAIMON_slope' % last_TMD]


        ###   AAIMON_slope TM01   ###
        # create numpy array of membranous over nonmembranous conservation ratios (identity)
        hist_data = (df_dict[prot_list]['TM01_AAIMON_slope'] * 1000).dropna()
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

        ###   AAIMON_slope last TMD   ###
        # create numpy array of membranous over nonmembranous conservation ratios (identity)
        hist_data = (df_dict[prot_list]['AAIMON_slope_last_TMD'] * 1000).dropna()
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

    ax.set_xlabel(r'm$_{\rm TM/nonTM} *10^{\rm -3}$', fontsize=fontsize, labelpad=20)
    # move the x-axis label closer to the x-axis
    ax.xaxis.set_label_coords(0.5, -0.085)
    # x and y axes min and max
    xlim_min = -30
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
    # create legend
    ax.xaxis.set_label_coords(0.5, -0.07)
    ax.yaxis.set_label_coords(-0.005, 0.5)

    # add annotations
    ax.annotate(s="TM less conserved", xy=(0, -0.09), fontsize=fontsize, xytext=None, xycoords='axes fraction')
    ax.annotate(s="TM more conserved", xy=(1.0, -0.09), fontsize=fontsize, xytext=None, horizontalalignment='right', xycoords='axes fraction')

    create_legend_for_AAIMON_slope = False
    if create_legend or create_legend_for_AAIMON_slope:
        ### create legend with additional 2 elements corresponding to AAIMON and AAIMON_n ###
        # Get artists and labels for legend and chose which ones to display
        handles, labels = ax.get_legend_handles_labels()
        display = (list(range(0, len(protein_lists) + 1, 1)))
        # Create custom artists
        TM01 = plt.Line2D((0, 1), (0, 0), color='k', linewidth=linewidth)
        last_TM = plt.Line2D((0, 1), (0, 0), color='k', linestyle=':', linewidth=linewidth)
        # Create legend from custom artist/label lists
        ax.legend([handle for i, handle in enumerate(handles) if i in display] + [TM01, last_TM],
                  [label for i, label in enumerate(labels) if i in display] + ['TM01', 'last TM'],
                  fontsize=fontsize - 3, frameon=True, loc='upper right')  # , bbox_to_anchor=(1.07, 1.12))
    else:
        # Create custom artists
        TM01 = plt.Line2D((0, 1), (0, 0), color='k', linewidth=linewidth)
        last_TM = plt.Line2D((0, 1), (0, 0), color='k', linestyle=':', linewidth=linewidth)
        # Create legend from custom artist/label lists
        ax.legend([TM01, last_TM],
                  ['TM01', 'last TM'],
                  fontsize=fontsize - 3, frameon=True, loc='upper right')  # , bbox_to_anchor=(1.07, 1.12))
    # plt.gcf().subplots_adjust(bottom=0.15)
    utils.save_figure(fig, Fig_name, base_filepath, save_png, save_pdf)

    # --------------------------------------------------------------------------------------------------------------------------------#
    Fig_Nr = 16
    title = '% AA identity TMD vs nonTMD'
    Fig_name = 'Fig16_perc_AA_identity_TMD_vs_nonTMD'

    fig, axes = plt.subplots(ncols=len(df_dict), nrows=1, sharex=True, sharey=True, figsize=(len(df_dict) * 2.5, 5))
    vmax = 10
    for n, (ax, prot_list) in enumerate(zip(axes.flat, protein_lists)):

        ax.plot([0, 100], [0, 100], color='#CCCCCC', linewidth=1)
        # ax.set_aspect('equal')
        ax.grid(False, which='both')
        ax.tick_params(axis='both', which='major', length=3, width=1, color='#CCCCCC')
        ax.set(adjustable='box-forced', aspect='equal')
        ax.set_xticks([0, 25, 50, 75, 100])
        ax.set_yticks([0, 25, 50, 75, 100])
        # histogram definition
        # data range
        xyrange = [[0, 100], [0, 100]]
        # number of bins
        bins = [50, 50]
        # density threshold
        thresh = 1

        # data definition
        data_TMD = df_dict[prot_list].TMD_perc_identity_mean_all_TMDs * 100
        data_nonTMD = df_dict[prot_list].nonTMD_perc_ident_mean * 100
        xdat, ydat = data_TMD, data_nonTMD

        # histogram the data
        hh, locx, locy = scipy.histogram2d(xdat, ydat, range=xyrange, bins=bins)

        # fill the areas with low density by NaNs
        hh[hh < thresh] = np.nan

        # if there are just a few datapoints, multiply count to get visible values
        if len(xdat) < 50:
            hh = hh * 5
            ax.annotate('count in bin $*5$', fontsize=5, xy=(69, 1))

        # plot the data with imshow
        im = ax.imshow(np.flipud(hh.T), cmap='Oranges', extent=np.array(xyrange).flatten(),
                       interpolation='none', origin='upper', aspect='equal', vmin=1, vmax=vmax)

        ax.annotate(dfv.loc[prot_list, 'list_description'], xy=(3, 90), color=dfv.loc[prot_list, 'color'], fontsize=10)
        ax.tick_params(labelsize=10)
        ax.annotate('TMD less conserved', xy=(1, 7), rotation=90, fontsize=8, ha='left', va='bottom')
        ax.annotate('TMD more conserved', xy=(7, 1), rotation=0, fontsize=8, ha='left', va='bottom')

        if n == 0:
            ax.set_ylabel('% AA identity nonTMD')

    # get colorbar from latest imshow element (color scale should be the same for all subplots)
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.81, 0.35, 0.04 / len(df_dict), 0.3])
    fig.colorbar(im, cax=cbar_ax)
    labels = cbar_ax.get_ymajorticklabels()
    labels[-1] = '>{}'.format(vmax)
    cbar_ax.set_yticklabels(labels)

    plt.subplots_adjust(wspace=0.15, hspace=0.05)
    fig.text(0.465, 0.27, '% AA identity TMD', ha='center', fontsize=10)

    utils.save_figure(fig, Fig_name, base_filepath, save_png, save_pdf)

    #---------------------------------------------------------------------------------------------------------------------------------#
    Fig_Nr = 17
    title = 'number of TMDs in GPCRs'
    Fig_name = 'Fig17_number_of_TMDs_in_GPCRs'
    if len(df_dict) > 4:
        sys.stdout.write('cannot create plot, too many lists')
    else:
        fig, axes = plt.subplots(2, 2)

        for ax in axes.flat:
            ax.set_xlim(0, 15)
            ax.set_ylim(0, 1.15)
            ax.set_xlabel('number of TMDs', fontsize=fontsize)
            ax.set_ylabel('freq', fontsize=fontsize)
            ax.tick_params(labelsize=fontsize)
            ax.set_xticks(range(0, 16, 2))
            ax.set_yticklabels([])

        for ax, prot_list in zip(axes.flat, protein_lists):
            if "uniprot_KW" in df_dict[prot_list]:
                # check if protein is a GPCR
                list_GPCR_KW = ['G-protein coupled receptor']
                df_dict[prot_list]['GPCR'] = df_dict[prot_list]['uniprot_KW'].apply(korbinian.cons_ratio.keywords.KW_list_contains_any_desired_KW, args=(list_GPCR_KW,))
                df_GPCR = df_dict[prot_list][df_dict[prot_list].GPCR == True]
                # check if df_GPCR contains elements
                if df_GPCR.shape[0] == 0:
                    continue

                hist_data = df_GPCR.number_of_TMDs_excl_SP.value_counts().values
                hist_data_norm = hist_data / hist_data.max()
                bins = df_GPCR.number_of_TMDs_excl_SP.value_counts().index

                ax.bar(bins, hist_data_norm, align='center', color=dfv.loc[prot_list, 'color'])
                ax.annotate(dfv.loc[prot_list, 'list_description'], xy=(0.5, 1.05), color=dfv.loc[prot_list, 'color'], fontsize=fontsize - 2)
                ax.annotate('only GPCRs\nn = {}'.format(len(df_GPCR)), xy=(0.5, 0.80), fontsize=fontsize - 2, alpha=0.5)

        plt.tight_layout()
        utils.save_figure(fig, Fig_name, base_filepath, save_png, save_pdf)

    #dfv.to_csv(os.path.join(base_filepath, 'Lists_%s_variables.csv'%str_protein_lists))
    sys.stdout.write("\n~~~~~~~~~~~~         compare_lists finished           ~~~~~~~~~~~~\n")