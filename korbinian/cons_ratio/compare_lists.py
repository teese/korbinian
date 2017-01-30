import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import ast
import sys
import csv
import datetime
import korbinian.utils as utils


def compare_lists (s,):
    sys.stdout.write("\n~~~~~~~~~~~~         starting compare_lists           ~~~~~~~~~~~~\n")
    protein_lists = ast.literal_eval(s['protein_lists'])
    protein_list_names = ast.literal_eval(s['protein_list_names'])
    df_variables = pd.DataFrame(index=protein_lists)
    df_variables['protein_list_names'] = protein_list_names
    # get current date and time for folder name, join elements in protein_lists to string
    str_protein_lists = '-'.join(map(str, sorted(protein_lists)))
    now = datetime.datetime.now()
    folder_name = '{}{:02}{}-{}-{}_Lists-{}'.format(now.year, now.month, now.day, now.hour, now.minute, str_protein_lists)

    path_to_files = os.path.join(s["data_dir"], "compare_lists", folder_name)

    # create folder in list summary directory to hold keyword data
    if not os.path.exists(path_to_files):
        os.makedirs(path_to_files)

    # import datasets, load data and hold pandas dataframes in dictionary
    df_dict = {}
    for n, prot_list in enumerate(protein_lists):
        df_variables.loc[prot_list, 'base_filename_summaries'] = os.path.join(s["data_dir"], "summaries", '%02d' % prot_list, 'List%02d_summary.csv' % prot_list)
        df_variables.loc[prot_list, 'base_filename_cr_summaries'] = os.path.join(s["data_dir"], "summaries", '%02d' % prot_list, 'List%02d_cr_summary.csv' % prot_list)
        df_variables.loc[prot_list, 'compare_dir'] = os.path.join(path_to_files, )

        dfx = pd.read_csv(df_variables.loc[prot_list, 'base_filename_summaries'], sep=",", quoting=csv.QUOTE_NONNUMERIC, index_col=0)
        dfy = pd.read_csv(df_variables.loc[prot_list, 'base_filename_cr_summaries'], sep=",", quoting=csv.QUOTE_NONNUMERIC, index_col=0)
        df_temp = pd.merge(dfy, dfx, left_index=True, right_index=True, suffixes=('_dfy', ''))
        proteins_before_dropping = len(df_temp)
        df_temp = df_temp[np.isfinite(df_temp['AAIMON_mean_all_TMDs'])]
        df_temp = df_temp[df_temp.TM01_AAIMON_n_homol > 4]
        proteins_after_dropping = len(df_temp)
        sys.stdout.write('List{:02} number dropped proteins {}\n'.format(prot_list, proteins_before_dropping - proteins_after_dropping))
        df_dict[prot_list] = df_temp


    ###############################################################
    #                                                             #
    #                           plot data                         #
    #                                                             #
    ###############################################################

    binlist = np.linspace(0.4, 2, 23)
    plt.style.use('seaborn-whitegrid')
    fig, ax = plt.subplots()
    alpha = 1
    fontsize = 10
    linewidth = 2
    color_list = ['#949494', '#EE762C', '#005C96']
    offset = len(protein_lists) - 1

    for n, element in enumerate(df_dict.keys()):
        ###   non-normalised AAIMON   ###
        # create numpy array of membranous over nonmembranous conservation ratios (identity)
        hist_data_SP = np.array(df_dict[element]['AAIMON_mean_all_TMDs'])
        # use numpy to create a histogram
        freq_counts, bin_array = np.histogram(hist_data_SP, bins=binlist)
        freq_counts_normalised = freq_counts / freq_counts.max() + offset
        # assuming all of the bins are exactly the same size, make the width of the column equal to XX% (e.g. 95%) of each bin
        col_width = float('%0.3f' % (0.95 * (bin_array[1] - bin_array[0])))
        # when align='center', the central point of the bar in the x-axis is simply the middle of the bins ((bin_0-bin_1)/2, etc)
        centre_of_bar_in_x_axis = (bin_array[:-2] + bin_array[1:-1]) / 2
        # add the final bin, which is physically located just after the last regular bin but represents all higher values
        bar_width = centre_of_bar_in_x_axis[3] - centre_of_bar_in_x_axis[2]
        centre_of_bar_in_x_axis = np.append(centre_of_bar_in_x_axis, centre_of_bar_in_x_axis[-1] + bar_width)
        linecontainer_AAIMON_mean = ax.plot(centre_of_bar_in_x_axis, freq_counts_normalised, color=color_list[n],
                                            alpha=alpha, linewidth=linewidth,
                                            label=df_variables.loc[element, 'protein_list_names'])

        ###   normalised AAIMON   ###
        # create numpy array of membranous over nonmembranous conservation ratios (identity)
        hist_data_SP = np.array(df_dict[element]['AAIMON_mean_all_TMDs_n'])
        # use numpy to create a histogram
        freq_counts, bin_array = np.histogram(hist_data_SP, bins=binlist)
        freq_counts_normalised = freq_counts / freq_counts.max() + offset
        # assuming all of the bins are exactly the same size, make the width of the column equal to XX% (e.g. 95%) of each bin
        col_width = float('%0.3f' % (0.95 * (bin_array[1] - bin_array[0])))
        # when align='center', the central point of the bar in the x-axis is simply the middle of the bins ((bin_0-bin_1)/2, etc)
        centre_of_bar_in_x_axis = (bin_array[:-2] + bin_array[1:-1]) / 2
        # add the final bin, which is physically located just after the last regular bin but represents all higher values
        bar_width = centre_of_bar_in_x_axis[3] - centre_of_bar_in_x_axis[2]
        centre_of_bar_in_x_axis = np.append(centre_of_bar_in_x_axis, centre_of_bar_in_x_axis[-1] + bar_width)
        linecontainer_AAIMON_mean = ax.plot(centre_of_bar_in_x_axis, freq_counts_normalised, ':', color=color_list[n],
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
    ax.xaxis.set_label_coords(0.45, -0.085)
    # x and y axes min and max
    xlim_min = 0.45
    xlim_max = 1.95
    ax.set_xlim(xlim_min, xlim_max)
    ylim_min = -0.01
    ylim_max = len(protein_lists) + 0.1
    ax.set_ylim(ylim_min, ylim_max)
    # set x-axis ticks
    # use the slide selection to select every second item in the list as an xtick(axis label)
    # ax.set_xticks([float('%0.1f' % c) for c in centre_of_bar_in_x_axis[::3]])
    ax.get_yaxis().set_ticks([])
    ax.set_ylabel('relative freqency', rotation='vertical', fontsize=fontsize)
    # change axis font size
    ax.tick_params(labelsize=fontsize)
    # create legend
    ax.xaxis.set_label_coords(0.5, -0.07)
    ax.yaxis.set_label_coords(-0.005, 0.5)

    # Get artists and labels for legend and chose which ones to display
    handles, labels = ax.get_legend_handles_labels()
    display = (0, 1, 2)

    # Create custom artists
    AAIMON = plt.Line2D((0, 1), (0, 0), color='k', linewidth=linewidth)
    AAIMON_norm = plt.Line2D((0, 1), (0, 0), color='k', linestyle=':', linewidth=linewidth)

    # Create legend from custom artist/label lists
    ax.legend([handle for i, handle in enumerate(handles) if i in display] + [AAIMON, AAIMON_norm],
              [label for i, label in enumerate(labels) if i in display] + ['AAIMON', 'AAIMON norm.'],
              fontsize=fontsize, frameon=True)

    utils.save_figure(fig, Fig_name='comparison', base_filepath=path_to_files, save_png=True, save_pdf=False)


    #df_variables.to_csv(os.path.join(path_to_files, 'Lists_%s_variables.csv'%str_protein_lists))



