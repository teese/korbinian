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
import pickle
import scipy
import seaborn as sns
import scipy.stats as stats
import time
import zipfile
# import debugging tools
from korbinian.utils import pr, pc, pn, aaa


def compare_lists (s, df_lists_tab):
    """ Create figures comparing protein lists.

    Compares a number of parameters. Some from UniProt. Others from the
    conservation ratio analysis. Analysis of up to 7 lists is possible. Only tested for 3.

    Histograms of AAIMON ratios
    Histograms of AAIMON slope
    Full sequence length
    Length of EM
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

    df_lists_tab : pd.DataFrame
        Full dataframe of the lists tab of the settings file

    DataFrames
    ----------
    dfv : dataframe for paths and settings
        index : [1,2,44] (list numbers)
        columns : ['list_description', 'list_csv', 'cr_summary_csv', 'compare_dir', 'color'] 	etc

    df_merged : list01.csv and list01_cr_summary.csv merged together
                formerly df_temp


    """
    sys.stdout.write("\n\n~~~~~~~~~~~~         starting compare_lists           ~~~~~~~~~~~~\n\n")
    # set up general plotting parameters
    plt.style.use('seaborn-whitegrid')
    plt.rcParams["legend.framealpha"] = 1
    plt.rcParams['errorbar.capsize'] = 3
    alpha = 1
    fontsize = 12
    linewidth = 2
    anno_fontsize = fontsize - 3
    xyc = "axes fraction"

    # set numpy division error to "ignore"
    np.seterr(divide='ignore', invalid='ignore')
    create_legend = True
    save_png = s['save_png']
    save_pdf = s['save_pdf']
    n_bins_cons = 60
    n_bins_lipo = 25

    # get a color list (HTML works best)
    #color_list = utils.create_colour_lists()['HTML_list01']
    color_list = utils.create_colour_lists()['BGO_arr']
    color_list_dark = utils.create_colour_lists()['BGO_dark_arr']

    protein_lists = ast.literal_eval(s['protein_lists'])
    # number of protein lists
    n_prot_lists = len(protein_lists)

    df_lists_tab = df_lists_tab.reindex(index=protein_lists)
    df_lists_tab["color"] = color_list[0:n_prot_lists]
    df_lists_tab["color_dark"] = color_list_dark[0:n_prot_lists]

    # initialise dataframe that hols all variables for every protein list
    dfv = pd.DataFrame(index=protein_lists)
    # get current date and time for folder name, join elements in protein_lists to string
    str_protein_lists = '-'.join(map(str, sorted(protein_lists)))

    if s["compare_lists_create_new_unique_folder"]:
        # create a new folder with a datestamp, e.g. 20170727-1117_Lists-1-2-3
        now = datetime.datetime.now()
        folder_name = '{}{:02}{:02}-{:02}{:02}_Lists-{}'.format(now.year, now.month, now.day, now.hour, now.minute, str_protein_lists)
    else:

        folder_name = 'Lists-{}'.format(str_protein_lists)

    if s["compare_lists_testing_mode"]:
        # folder_name for testing
        folder_name = 'Folder_Test'

    if "compare_lists_run_only_this_fig" in s and isinstance(s["compare_lists_run_only_this_fig"], int):
        # run only this fig
        list_figs_to_run = [s["compare_lists_run_only_this_fig"]]
    else:
        #else run all possible figures, up to 100
        list_figs_to_run = range(100)

    base_filepath = os.path.join(s["data_dir"], "compare_lists", folder_name)

    sys.stdout.write('data is saved to "{}"\n'.format(base_filepath))

    # create folder in list summary directory to hold keyword data
    if not os.path.exists(base_filepath):
        os.makedirs(base_filepath)

    for n, prot_list in enumerate(protein_lists):
        # add the list description (also from the df_lists_tab)
        # this is legacy code : they could all be grabbed directly from df_lists_tab
        dfv.loc[prot_list, 'list_csv'] = os.path.join(s["data_dir"], "summaries", '%02d' % prot_list, 'List%02d.csv' % prot_list)
        dfv.loc[prot_list, 'cr_summary_csv'] = os.path.join(s["data_dir"], "summaries", '%02d' % prot_list, 'List%02d_cr_summary.csv' % prot_list)
        dfv.loc[prot_list, 'compare_dir'] = base_filepath

    df_dict_pickle_path = os.path.join(base_filepath, "df_dict.pickle")
    if os.path.isfile(df_dict_pickle_path):
        yesterday = time.time() - 60 * 60 * 24
        file_mtime = os.path.getmtime(df_dict_pickle_path)
        file_is_old = file_mtime < yesterday
    else:
        file_is_old = False

    if s["regenerate_df_dict_with_data_for_compare_lists"] or file_is_old:
        # import datasets, load data and hold pandas dataframes in dictionary
        df_dict = {}
        for n, prot_list in enumerate(protein_lists):
            # # add the list description (also from the df_lists_tab)
            # # this is legacy code : they could all be grabbed directly from df_lists_tab
            # dfv.loc[prot_list, 'list_description'] = df_lists_tab.loc[prot_list, 'list_description']
            # dfv.loc[prot_list, 'list_csv'] = os.path.join(s["data_dir"], "summaries", '%02d' % prot_list, 'List%02d.csv' % prot_list)
            # dfv.loc[prot_list, 'cr_summary_csv'] = os.path.join(s["data_dir"], "summaries", '%02d' % prot_list, 'List%02d_cr_summary.csv' % prot_list)
            # dfv.loc[prot_list, 'compare_dir'] = base_filepath
            # dfv.loc[prot_list, 'color'] = color_list[n]

            # read list summary.csv and cr_summary.csv from disk, join them to one big dataframe
            df_list = pd.read_csv(dfv.loc[prot_list, 'list_csv'], sep=",", quoting=csv.QUOTE_NONNUMERIC, index_col=0, low_memory=False)
            df_cr_summ = pd.read_csv(dfv.loc[prot_list, 'cr_summary_csv'], sep=",", quoting=csv.QUOTE_NONNUMERIC, index_col=0, low_memory=False)
            df_merged = pd.merge(df_list, df_cr_summ, left_index=True, right_index=True, suffixes=('_dfy', ''))
            # count proteins before and after dropping proteins with less than x homologues (x from settings file "lists"), drop proteins
            proteins_before_dropping = len(df_merged)
            # only keep proteins where an AAIMON_mean_all_TM_res has been calculated
            df_merged = df_merged[np.isfinite(df_merged['AAIMON_mean_all_TM_res'])]
            # only keep proteins where there is a minimum number of homologues
            df_merged = df_merged[df_merged.AAIMON_n_homol >= df_lists_tab.loc[prot_list, 'min_homol']]
            proteins_after_dropping = len(df_merged)
            sys.stdout.write('List{:02} - {}: {} proteins ({}/{} dropped)\n'.format(prot_list, df_lists_tab.loc[prot_list,'list_description'], proteins_after_dropping,
                                                                                    proteins_before_dropping - proteins_after_dropping, proteins_before_dropping))
            # make list of TMDs a python list
            df_merged.list_of_TMDs = df_merged.list_of_TMDs.apply(lambda x: ast.literal_eval(x))

            # SHIFTED TO prepare_protein_list
            # # check for uniprot keywords, make them a python list and search for GPCRs to exclude them later
            # if 'uniprot_KW' in df_merged.columns:
            #     df_merged['uniprot_KW'] = df_merged['uniprot_KW'].apply(lambda x: ast.literal_eval(x))
            #     df_merged['GPCR'] = df_merged['uniprot_KW'].apply(KW_list_contains_any_desired_KW, args=(['G-protein coupled receptor'],))
            # else:
            #     df_merged['GPCR'] = False

            # SHIFTED TO GATHER
            # if not 'AAIMON_slope_central_TMDs' in df_merged.columns:
            #     for acc in df_merged.index:
            #         if df_merged.loc[acc, 'number_of_TMDs'] >= 3:
            #             list_of_central_TMDs = df_merged.loc[acc, 'list_of_TMDs'][1:-1]
            #             list_mean_slope_central_TMDs = []
            #             for TMD in list_of_central_TMDs:
            #                 list_mean_slope_central_TMDs.append(pd.to_numeric(df_merged.loc[acc, '%s_AAIMON_slope' % TMD]))
            #             df_merged.loc[acc, 'AAIMON_slope_central_TMDs'] = np.mean(list_mean_slope_central_TMDs)
            # add dataframe to a dictionary of dataframes
            df_dict[prot_list] = df_merged

        # save the large df_dict as a pickle in the output folder
        with open(df_dict_pickle_path, "wb") as pkl:
            pickle.dump(df_dict, pkl, protocol=pickle.HIGHEST_PROTOCOL)

    else:
        list_str = "-".join([str(i) for i in sorted(protein_lists)])
        df_dict_pickle = r"D:\Databases\compare_lists\Lists-{}\df_dict.pickle".format(list_str)
        with open(df_dict_pickle, "rb") as pkl:
            df_dict = pickle.load(pkl)

    list_descriptions = df_lists_tab.loc[protein_lists, 'list_description']

    # check if structure data is plotted -> changed bins due to reduced number of proteins
    reduce_bins = False
    # iterate through "singlepass", "multipass", etc
    for element in df_lists_tab.loc[protein_lists, "list_description"]:
        #if 'STR' in element[:3]:
        if "crystal" in element:
            reduce_bins = True
    if reduce_bins == True:
        n_bins_cons = 20
        n_bins_lipo = 15
        sys.stdout.write('reduced number of bins due to structural data usage: n_bins_cons = {}; n_bins_lipo = {}\n\n'.format(n_bins_cons,n_bins_lipo))

    ###############################################################
    #                                                             #
    #                           plot data                         #
    #                                                             #
    ###############################################################
    # y-axis subfigure height (ysh)
    # note that to have a position of 80% of each subplot, the formula for y is
    # y = 0.8*ysh+ysh*offset
    # this gives the values between 0 and 1 used for the "axes fraction" xycoords of matplotlib
    ysh = 1 / n_prot_lists

    # for plots with first-central-last, set the minimum number of TMDs in the protein
    min_n_TMDs_first_last = 5

    # setup a dictionary with the offsets
    # take the range of the lists, minus one, and reverse (e.g for 3 lists, gives [2,1,0]
    offset_list = list(np.array(range(1, len(protein_lists) + 1)) - 1)[::-1]
    # add to a dict, so that it can be accessed using the protein list as a key
    offset_dict = {}
    for n, protein_list in enumerate(protein_lists):
        offset_dict[protein_list] = offset_list[n]

    #--------------------------------------------------------------------------------------------------------------------------------#
    Fig_Nr = 1
    if Fig_Nr in list_figs_to_run:
        title = 'histograms AAIMON and AAIMON_n'
        for suffix in ["", "_and_AAIMON_n"]:

            Fig_name = 'Fig01_Histograms_of_mean_AAIMON{}'.format(suffix)
            binlist = np.linspace(0, 2, n_bins_cons+1)
            fig, ax = plt.subplots()
            offset = len(protein_lists) - 1

            for prot_list in protein_lists:
                color = df_lists_tab.loc[prot_list, 'color']
                ###   non-normalised AAIMON   ###
                # create numpy array of membranous over nonmembranous conservation ratios (identity)
                hist_data = np.array(df_dict[prot_list]['AAIMON_mean_all_TM_res'])
                AAIMON_mean = hist_data.mean()
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
                linecontainer_AAIMON_mean = ax.plot(centre_of_bar_in_x_axis, freq_counts_normalised, color=color,
                                                    alpha=alpha, linewidth=linewidth,
                                                    label=df_lists_tab.loc[prot_list, 'list_description'])

                if suffix == "_and_AAIMON_n":
                    ###   normalised AAIMON   ###
                    # create numpy array of membranous over nonmembranous conservation ratios (identity)
                    hist_data = np.array(df_dict[prot_list]['AAIMON_n_mean_all_TM_res'])
                    AAIMON_n_mean = hist_data.mean()
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
                    linecontainer_AAIMON_mean = ax.plot(centre_of_bar_in_x_axis, freq_counts_normalised, ':', color=color,
                                                        alpha=alpha,
                                                        linewidth=linewidth)

                    ###############################################################
                    #                                                             #
                    #        annotate the mean values on the plots                #
                    #                                                             #
                    ###############################################################
                    ax.annotate("mean:", [0.005, 0.9*ysh+ysh*offset], fontsize=anno_fontsize, fontproperties="monospace", color=color, xycoords=xyc, weight="semibold")
                    ax.annotate("original{: >7.02f}".format(AAIMON_mean), [0.01, 0.8*ysh+ysh*offset], fontsize=anno_fontsize, fontproperties="monospace", color=color, xycoords=xyc, weight="semibold")
                    ax.annotate("normalised{: >5.02f}".format(AAIMON_n_mean), [0.01, 0.7*ysh+ysh*offset], fontsize=anno_fontsize, fontproperties="monospace", color=color, xycoords=xyc, weight="semibold")
                else:
                    ax.annotate("mean: {:.02f}".format(AAIMON_mean), [0.005, 0.9 * ysh + ysh * offset], fontsize=anno_fontsize, fontproperties="monospace", color=color, xycoords=xyc, weight="semibold")

                offset = offset - 1

            ###############################################################
            #                                                             #
            #                       set up plot style                     #
            #                                                             #
            ###############################################################

            ax.set_xlabel('TM/EM conservation', fontsize=fontsize)
            # move the x-axis label closer to the x-axis
            ax.xaxis.set_label_coords(0.5, -0.085)
            # x and y axes min and max
            xlim_min = 0.4
            xlim_max = 1.9
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
            ax.annotate(s="TM less conserved", xy=(0, -0.09), fontsize=fontsize, xytext=None, xycoords='axes fraction')
            ax.annotate(s="TM more conserved", xy=(1.0, -0.09), fontsize=fontsize, xytext=None, horizontalalignment='right', xycoords='axes fraction')

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
                          [label for i, label in enumerate(labels) if i in display] + ['non-normalised', 'normalised'],
                          fontsize=fontsize - 3, frameon=True, loc='upper right')#, bbox_to_anchor=(1.07, 1.12))
            else:
                # Create custom artists
                AAIMON = plt.Line2D((0, 1), (0, 0), color='k', linewidth=linewidth)
                AAIMON_norm = plt.Line2D((0, 1), (0, 0), color='k', linestyle=':', linewidth=linewidth)
                # Create legend from custom artist/label lists
                ax.legend([AAIMON, AAIMON_norm],
                          ['non-normalised', 'normalised'],
                          fontsize=fontsize - 3, frameon=True, loc='upper right')#, bbox_to_anchor=(1.07, 1.12))

            utils.save_figure(fig, Fig_name, base_filepath=base_filepath, save_png=save_png, save_pdf=save_pdf)

        # --------------------------------------------------------------------------------------------------------------------------------#
        title = 'hist AAIMON_single_yaxis'
        Fig_name = 'Fig01b_hist AAIMON_single_yaxis'
        binlist = np.linspace(0, 2, n_bins_cons+1)
        fig, ax = plt.subplots()

        for prot_list in protein_lists:
            offset = offset_dict[prot_list]
            label = df_lists_tab.loc[prot_list, 'list_description']
            color = df_lists_tab.loc[prot_list, 'color']
            ###   non-normalised AAIMON   ###
            # create numpy array of membranous over nonmembranous conservation ratios (identity)
            hist_data = np.array(df_dict[prot_list]['AAIMON_mean_all_TM_res'])
            AAIMON_slope_mean = hist_data.mean()
            # use numpy to create a histogram
            freq_counts, bin_array = np.histogram(hist_data, bins=binlist)
            freq_counts_normalised = freq_counts / freq_counts.max()
            # assuming all of the bins are exactly the same size, make the width of the column equal to XX% (e.g. 95%) of each bin
            col_width = float('%0.3f' % (0.95 * (bin_array[1] - bin_array[0])))
            # when align='center', the central point of the bar in the x-axis is simply the middle of the bins ((bin_0-bin_1)/2, etc)
            centre_of_bar_in_x_axis = (bin_array[:-2] + bin_array[1:-1]) / 2
            # add the final bin, which is physically located just after the last regular bin but represents all higher values
            bar_width = centre_of_bar_in_x_axis[3] - centre_of_bar_in_x_axis[2]
            centre_of_bar_in_x_axis = np.append(centre_of_bar_in_x_axis, centre_of_bar_in_x_axis[-1] + bar_width)
            linecontainer_AAIMON_mean = ax.plot(centre_of_bar_in_x_axis, freq_counts_normalised, color=color,
                                                alpha=alpha, linewidth=linewidth,
                                                label=label)


            ###############################################################
            #                                                             #
            #        annotate the mean values on the plots                #
            #                                                             #
            ###############################################################
            ax.annotate("mean:", [0.01, 0.9], fontsize=fontsize, fontproperties="monospace", color=color, xycoords=xyc, weight="semibold")
            ax.annotate("{:<10} {: >6.02f}".format(label, AAIMON_slope_mean), [0.02, 0.85 - 0.05*offset], fontsize=fontsize, fontproperties="monospace", color=color, xycoords=xyc, weight="semibold")


        ###############################################################
        #                                                             #
        #                       set up plot style                     #
        #                                                             #
        ###############################################################
        ax.set_xlabel('TM/EM conservation ratio', fontsize=fontsize)
        # move the x-axis label closer to the x-axis
        ax.xaxis.set_label_coords(0.5, -0.085)
        # x and y axes min and max
        xlim_min = 0.4
        xlim_max = 1.9
        ax.set_xlim(xlim_min, xlim_max)
        ylim_min = -0.01
        ylim_max = 1
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
        ax.annotate(s="TM less conserved", xy=(0, -0.09), fontsize=fontsize, xytext=None, xycoords='axes fraction')
        ax.annotate(s="TM more conserved", xy=(1.0, -0.09), fontsize=fontsize, xytext=None, horizontalalignment='right', xycoords='axes fraction')

        ax.legend(frameon=True, loc='upper right', fontsize=fontsize)

        utils.save_figure(fig, Fig_name, base_filepath, save_png, save_pdf)




    #--------------------------------------------------------------------------------------------------------------------------------#
    Fig_Nr = 2
    if Fig_Nr in list_figs_to_run:
        title = 'histograms AAIMON_slope and AAIMON_n_slope'
        for suffix in ["", "_and_AAIMON_n_slope"]:
            Fig_name = 'Fig02_Histograms_of_mean_AAIMON_slope{}'.format(suffix)
            binlist = np.linspace(-40, 40, n_bins_cons+1)
            fig, ax = plt.subplots()

            for prot_list in protein_lists:
                offset = offset_dict[prot_list]
                color = df_lists_tab.loc[prot_list, 'color']
                ###   non-normalised AAIMON   ###
                # create numpy array of membranous over nonmembranous conservation ratios (identity)
                hist_data = np.array(df_dict[prot_list]['AAIMON_slope_all_TM_res']*1000)
                AAIMON_slope_mean = hist_data.mean()
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
                linecontainer_AAIMON_mean = ax.plot(centre_of_bar_in_x_axis, freq_counts_normalised, color=color,
                                                    alpha=alpha, linewidth=linewidth,
                                                    label=df_lists_tab.loc[prot_list, 'list_description'])
                if suffix == "_and_AAIMON_n_slope":
                    ###   normalised AAIMON   ###
                    # create numpy array of membranous over nonmembranous conservation ratios (identity)
                    hist_data = np.array(df_dict[prot_list]['AAIMON_n_slope_all_TM_res']*1000)
                    AAIMON_slope_mean_n = hist_data.mean()
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
                    linecontainer_AAIMON_mean = ax.plot(centre_of_bar_in_x_axis, freq_counts_normalised, ':', color=color,
                                                        alpha=alpha,
                                                        linewidth=linewidth)

                    ###############################################################
                    #                                                             #
                    #        annotate the mean values on the plots                #
                    #                                                             #
                    ###############################################################
                    ax.annotate("mean:", [0.01, 0.9*ysh+ysh*offset], fontsize=anno_fontsize, fontproperties="monospace", color=color, xycoords=xyc, weight="semibold")
                    ax.annotate("original {: >8.02f}".format(AAIMON_slope_mean), [0.02, 0.8*ysh+ysh*offset], fontsize=anno_fontsize, fontproperties="monospace", color=color, xycoords=xyc, weight="semibold")
                    ax.annotate("normalised {: >6.02f}".format(AAIMON_slope_mean_n), [0.02, 0.7*ysh+ysh*offset], fontsize=anno_fontsize, fontproperties="monospace", color=color, xycoords=xyc, weight="semibold")
                else:
                    ax.annotate("mean: {:0.2f}".format(AAIMON_slope_mean), [0.01, 0.9 * ysh + ysh * offset], fontsize=anno_fontsize, fontproperties="monospace", color=color, xycoords=xyc, weight="semibold")

            ###############################################################
            #                                                             #
            #                       set up plot style                     #
            #                                                             #
            ###############################################################

            ax.set_xlabel(r'm$_{\rm TM/EM} *10^{\rm -3}$', fontsize=fontsize, labelpad=20)
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
                          [label for i, label in enumerate(labels) if i in display] + ['non-normalised', 'normalised'],
                          fontsize=fontsize - 3, frameon=True, loc='upper right')#, bbox_to_anchor=(1.07, 1.12))
            else:
                # Create custom artists
                AAIMON = plt.Line2D((0, 1), (0, 0), color='k', linewidth=linewidth)
                AAIMON_norm = plt.Line2D((0, 1), (0, 0), color='k', linestyle=':', linewidth=linewidth)
                # Create legend from custom artist/label lists
                ax.legend([AAIMON, AAIMON_norm],
                          ['non-normalised', 'normalised'],
                          fontsize=fontsize - 3, frameon=True, loc='upper right')#, bbox_to_anchor=(1.07, 1.12))
            #plt.gcf().subplots_adjust(bottom=0.15)
            utils.save_figure(fig, Fig_name, base_filepath, save_png, save_pdf)

        # --------------------------------------------------------------------------------------------------------------------------------#
        Fig_Nr = "2b"
        title = 'hist AAIMON_slope_single_yaxis'
        Fig_name = 'Fig02b_hist AAIMON_slope_single_yaxis'
        binlist = np.linspace(-40, 40, n_bins_cons + 1)
        fig, ax = plt.subplots()

        for prot_list in protein_lists:
            offset = offset_dict[prot_list]
            label = df_lists_tab.loc[prot_list, 'list_description']
            color = df_lists_tab.loc[prot_list, 'color']
            ###   non-normalised AAIMON   ###
            # create numpy array of membranous over nonmembranous conservation ratios (identity)
            hist_data = np.array(df_dict[prot_list]['AAIMON_slope_all_TM_res'] * 1000)
            AAIMON_slope_mean = hist_data.mean()
            # use numpy to create a histogram
            freq_counts, bin_array = np.histogram(hist_data, bins=binlist)
            freq_counts_normalised = freq_counts / freq_counts.max()
            # assuming all of the bins are exactly the same size, make the width of the column equal to XX% (e.g. 95%) of each bin
            col_width = float('%0.3f' % (0.95 * (bin_array[1] - bin_array[0])))
            # when align='center', the central point of the bar in the x-axis is simply the middle of the bins ((bin_0-bin_1)/2, etc)
            centre_of_bar_in_x_axis = (bin_array[:-2] + bin_array[1:-1]) / 2
            # add the final bin, which is physically located just after the last regular bin but represents all higher values
            bar_width = centre_of_bar_in_x_axis[3] - centre_of_bar_in_x_axis[2]
            centre_of_bar_in_x_axis = np.append(centre_of_bar_in_x_axis, centre_of_bar_in_x_axis[-1] + bar_width)
            linecontainer_AAIMON_mean = ax.plot(centre_of_bar_in_x_axis, freq_counts_normalised, color=color,
                                                alpha=alpha, linewidth=linewidth,
                                                label=label)


            ###############################################################
            #                                                             #
            #        annotate the mean values on the plots                #
            #                                                             #
            ###############################################################
            ax.annotate("mean:", [0.01, 0.9], fontsize=fontsize, fontproperties="monospace", color=color, xycoords=xyc, weight="semibold")
            ax.annotate("{:<10} {: >6.02f}".format(label, AAIMON_slope_mean), [0.02, 0.85 - 0.05*offset], fontsize=fontsize, fontproperties="monospace", color=color, xycoords=xyc, weight="semibold")


        ###############################################################
        #                                                             #
        #                       set up plot style                     #
        #                                                             #
        ###############################################################

        ax.set_xlabel(r'm$_{\rm TM/EM} *10^{\rm -3}$', fontsize=fontsize, labelpad=20)
        # move the x-axis label closer to the x-axis
        ax.xaxis.set_label_coords(0.5, -0.085)
        # x and y axes min and max
        xlim_min = -30
        xlim_max = 30
        ax.set_xlim(xlim_min, xlim_max)
        ylim_min = -0.01
        ylim_max = 1
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

        ax.legend(frameon=True, loc='upper right', fontsize=fontsize)

        utils.save_figure(fig, Fig_name, base_filepath, save_png, save_pdf)

    #--------------------------------------------------------------------------------------------------------------------------------#
    Fig_Nr = 3
    if Fig_Nr in list_figs_to_run:
        title = 'seqlen'
        Fig_name = 'Fig03_compare_seqlen'
        linspace_binlist = np.linspace(0, 2500, 100)
        binlist = np.append(linspace_binlist, 10000)
        fig, ax = plt.subplots()
        offset = len(protein_lists) - 1

        for prot_list in protein_lists:
            color = df_lists_tab.loc[prot_list, 'color']
            ###   non-normalised AAIMON   ###
            # create numpy array of membranous over nonmembranous conservation ratios (identity)
            # hist_data = np.array(df_dict[prot_list]['nonTMD_SW_align_len_excl_gaps_mean'])
            hist_data = np.array(df_dict[prot_list]['seqlen'])
            mean_seqlen = hist_data.mean()
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
            linecontainer_AAIMON_mean = ax.plot(centre_of_bar_in_x_axis, freq_counts_normalised, color=color,
                                                alpha=alpha, linewidth=linewidth,
                                                label=df_lists_tab.loc[prot_list, 'list_description'])

            ###############################################################
            #                                                             #
            #        annotate the mean values on the plots                #
            #                                                             #
            ###############################################################
            ax.annotate("mean: {: >2.0f}".format(mean_seqlen), [0.86, 0.15*ysh+ysh*offset], fontsize=anno_fontsize, fontproperties="monospace", color=color, xycoords=xyc, weight="semibold")

            offset = offset - 1

        ###############################################################
        #                                                             #
        #                       set up plot style                     #
        #                                                             #
        ###############################################################

        ax.set_xlabel('sequence length (full protein)', fontsize=fontsize)
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
            # ax.legend([handle for i, handle in enumerate(handles) if i in display],
            #           [label for i, label in enumerate(labels) if i in display],
            #           fontsize=fontsize - 3, frameon=True, bbox_to_anchor=(1.07, 1.12))

            ax.legend([handle for i, handle in enumerate(handles) if i in display],
                      [label for i, label in enumerate(labels) if i in display],
                      fontsize=fontsize - 3, frameon=True, loc='upper right')

        utils.save_figure(fig, Fig_name, base_filepath, save_png, save_pdf)

    #--------------------------------------------------------------------------------------------------------------------------------#
    Fig_Nr = 4
    if Fig_Nr in list_figs_to_run:
        title = 'number of TMDs'
        Fig_name = 'Fig04_compare_number_of_TMDs'
        binlist = np.linspace(0, 25, 100)
        binlist = np.append(binlist, [50])
        fig, ax = plt.subplots()
        offset = len(protein_lists) - 1

        for prot_list in protein_lists:
            color = df_lists_tab.loc[prot_list, 'color']
            ###   non-normalised AAIMON   ###
            # create numpy array of membranous over nonmembranous conservation ratios (identity)
            # hist_data = np.array(df_dict[prot_list]['nonTMD_SW_align_len_excl_gaps_mean'])
            hist_data = np.array(df_dict[prot_list]['number_of_TMDs'])
            number_of_TMDs_mean = hist_data.mean()
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
            linecontainer_AAIMON_mean = ax.plot(centre_of_bar_in_x_axis, freq_counts_normalised, color=color,
                                                alpha=alpha, linewidth=linewidth,
                                                label=df_lists_tab.loc[prot_list, 'list_description'])

            ###############################################################
            #                                                             #
            #        annotate the mean values on the plots                #
            #                                                             #
            ###############################################################
            ax.annotate("mean:{: >4.01f}".format(number_of_TMDs_mean), [0.81, 0.1*ysh+ysh*offset], fontsize=anno_fontsize, fontproperties="monospace", color=color, xycoords=xyc, weight="semibold")

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
                      fontsize=fontsize - 3, frameon=True, loc='upper right')#bbox_to_anchor=(1.07, 1.12)

        utils.save_figure(fig, Fig_name, base_filepath, save_png, save_pdf)

    #--------------------------------------------------------------------------------------------------------------------------------#
    Fig_Nr = 5
    if Fig_Nr in list_figs_to_run:
        title = 'observed changes mean'
        Fig_name = 'Fig05_compare_observed_changes_mean'
        binlist = np.linspace(0, 100, 31)
        fig, ax = plt.subplots()
        offset = len(protein_lists) - 1

        for prot_list in protein_lists:
            color = df_lists_tab.loc[prot_list, 'color']
            ###   non-normalised AAIMON   ###
            # create numpy array of membranous over nonmembranous conservation ratios (identity)
            # hist_data = np.array(df_dict[prot_list]['nonTMD_SW_align_len_excl_gaps_mean'])
            hist_data = np.array(df_dict[prot_list]['obs_changes_mean'])
            evol_dist_mean = hist_data.mean()
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
            linecontainer_AAIMON_mean = ax.plot(centre_of_bar_in_x_axis, freq_counts_normalised, color=color,
                                                alpha=alpha, linewidth=linewidth,
                                                label=df_lists_tab.loc[prot_list, 'list_description'])

            ###############################################################
            #                                                             #
            #        annotate the mean values on the plots                #
            #                                                             #
            ###############################################################
            ax.annotate("mean: {: >7.02f}".format(evol_dist_mean), [0.82, 0.1*ysh+ysh*offset], fontsize=anno_fontsize, fontproperties="monospace", color=color, xycoords=xyc, weight="semibold")

            offset = offset - 1

        ###############################################################
        #                                                             #
        #                       set up plot style                     #
        #                                                             #
        ###############################################################

        ax.set_xlabel('evolutionary distance (% substitutions)', fontsize=fontsize)
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
                      fontsize=fontsize - 3, frameon=True, loc='upper right')#bbox_to_anchor=(1.07, 1.12)



            # # Create custom artists
            # AAIMON = plt.Line2D((0, 1), (0, 0), color='k', linewidth=linewidth)
            # AAIMON_norm = plt.Line2D((0, 1), (0, 0), color='k', linestyle=':', linewidth=linewidth)
            # # Create legend from custom artist/label lists
            # ax.legend([AAIMON, AAIMON_norm],
            #           ['non-normalised', 'normalised'],
            #           fontsize=fontsize - 3, frameon=True)#, bbox_to_anchor=(1.07, 1.12))

        utils.save_figure(fig, Fig_name, base_filepath, save_png, save_pdf)

    #--------------------------------------------------------------------------------------------------------------------------------#
    Fig_Nr = 6
    if Fig_Nr in list_figs_to_run:
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
            color = df_lists_tab.loc[prot_list, 'color']
            n_of_homologues = sum(df_dict[prot_list]['AAIMON_n_homol'])
            sys.stdout.write('number of homologues {}: {:.0f}\n'.format(df_lists_tab.loc[prot_list, 'list_description'], n_of_homologues))
            #ax1.bar(n + 1, n_of_homologues / 100000, width=0.75, align='center', color=df_lists_tab.loc[prot_list, 'color'], edgecolor='')
            # IMPORTANT: barchart plots n_of_homologues /1000 (x 10^-3)
            ax1.bar(n + 1, n_of_homologues / 1000, width=0.75, align='center', color=df_lists_tab.loc[prot_list, 'color'], edgecolor='')

            ###   non-normalised AAIMON   ###
            # create numpy array of membranous over nonmembranous conservation ratios (identity)
            # hist_data = np.array(df_dict[prot_list]['nonTMD_SW_align_len_excl_gaps_mean'])
            hist_data = np.array(df_dict[prot_list]['AAIMON_n_homol'])
            AAIMON_n_homol_mean = hist_data.mean()
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
            linecontainer_AAIMON_mean = ax.plot(centre_of_bar_in_x_axis, freq_counts_normalised, color=color,
                                                alpha=alpha, linewidth=linewidth,
                                                label=df_lists_tab.loc[prot_list, 'list_description'])

            ###############################################################
            #                                                             #
            #        annotate the mean values on the plots                #
            #                                                             #
            ###############################################################
            if n_prot_lists == 3:
                ax.annotate("mean: {: >7.02f}".format(AAIMON_n_homol_mean), [0.8, 0.1*ysh+ysh*offset], fontsize=anno_fontsize, fontproperties="monospace", color=color, xycoords=xyc, weight="semibold")
            else:
                ax.annotate("mean:{: >4.0f}".format(AAIMON_n_homol_mean), [0.8, 0.4 * ysh + ysh * offset], fontsize=anno_fontsize, fontproperties="monospace", color=color, xycoords=xyc, weight="semibold")

            offset = offset - 1

        ###############################################################
        #                                                             #
        #                       set up plot style                     #
        #                                                             #
        ###############################################################

        ax.set_xlabel('number of valid homologues (note: unequal bins)', fontsize=fontsize)# $*10^3$
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
        #ax.set_xticklabels((list_xticks) / 1000)
        ax.set_xticklabels(list_xticks)

        ax1.yaxis.tick_right()
        ax1.set_ylabel('number of valid homologues in entire dataset ($*10^3$)', fontsize=fontsize)
        #ax1.tick_params(axis='y', pad=-11)
        ax1.tick_params(labelsize=fontsize)
        ax1.yaxis.set_label_position("right")
        # ax1.xaxis.set_ticks(range(1, len(df_dict) + 1, 1))
        ax1.set_xticklabels([])
        #ax1.set_yticklabels(ax.get_yticks() / 1000)
        plt.subplots_adjust(wspace=0.05, hspace=0)

        ax1.grid(axis="x", b=False)
        ax.grid(axis="x", b=False)

        create_legend = False
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
    if Fig_Nr in list_figs_to_run:
        title = 'hist_lipo_mean_all_TMDs_mean'
        Fig_name = 'Fig07_hist_lipo_mean_all_TMDs_mean'
        # set min for bins and xlim
        min_ = -0.5
        # check overlap of lists with betabarrel
        BB_lists = [3,15,53]
        overlap_BB = set(BB_lists).intersection(set(protein_lists))
        contains_betabarrel = True if len(overlap_BB) > 0 else False

        # set max for bins and xlim
        if contains_betabarrel:
            max_ = 1.5
        else:
            max_ = 0.8

        binlist = np.linspace(min_, max_, n_bins_lipo+1) #41
        fig, ax = plt.subplots()
        offset = len(protein_lists) - 1

        for prot_list in protein_lists:
            color = df_lists_tab.loc[prot_list, 'color']
            is_multipass = True if df_lists_tab.loc[prot_list, "max_TMDs"] > 2 else False

            ###   TM01 lipo   ###
            # create numpy array of membranous over nonmembranous conservation ratios (identity)
            hist_data = np.array(df_dict[prot_list]['TM01_lipo'])
            lipo_TM01_mean = hist_data.mean()
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
            linecontainer_AAIMON_mean = ax.plot(centre_of_bar_in_x_axis, freq_counts_normalised, ':', color=color,
                                                alpha=alpha,
                                                linewidth=linewidth)

            if is_multipass:
                ###   lipo mean all TMDs excl. TM01   ###
                # create numpy array of membranous over nonmembranous conservation ratios (identity)
                hist_data = np.array(df_dict[prot_list]['lipo_mean_excl_TM01'].dropna())
                lipo_excl_TM01_mean = hist_data.mean()
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
                linecontainer_AAIMON_mean = ax.plot(centre_of_bar_in_x_axis, freq_counts_normalised, color=color,
                                                    alpha=alpha, linewidth=linewidth,
                                                    label=df_lists_tab.loc[prot_list, 'list_description'])


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
            # linecontainer_AAIMON_mean = ax.plot(centre_of_bar_in_x_axis, freq_counts_normalised, '--', color=dfv.loc[prot_list, 'color'],
            #                                     alpha=alpha,
            #                                     linewidth=linewidth)

            ###############################################################
            #                                                             #
            #        annotate the mean values on the plots                #
            #                                                             #
            ###############################################################
            ax.annotate("mean:", [0.01, 0.9*ysh+ysh*offset], fontsize=anno_fontsize, fontproperties="monospace", color=color, xycoords=xyc, weight="semibold")
            ax.annotate("TM01 {: >14.02f}".format(lipo_TM01_mean), [0.02, 0.8*ysh+ysh*offset], fontsize=anno_fontsize, fontproperties="monospace", color=color, xycoords=xyc, weight="semibold")
            if is_multipass:
                ax.annotate("excluding TM01 {: >3.02f}".format(lipo_excl_TM01_mean), [0.02, 0.7*ysh+ysh*offset], fontsize=anno_fontsize, fontproperties="monospace", color=color, xycoords=xyc, weight="semibold")

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
    if Fig_Nr in list_figs_to_run:
        title = 'perc of TMD region in protein'
        Fig_name = 'Fig08_perc_of_TMD_region_in_protein'
        min_ = 0
        max_ = 100
        binlist = np.linspace(min_, max_, 101)
        fig, ax = plt.subplots()
        offset = len(protein_lists) - 1

        for prot_list in protein_lists:
            color = df_lists_tab.loc[prot_list, 'color']
            ###   non-normalised AAIMON   ###
            # create numpy array of membranous over nonmembranous conservation ratios (identity)
            hist_data = np.array(df_dict[prot_list]['perc_TMD'])
            perc_TMD_mean = hist_data.mean()/100
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
            linecontainer_AAIMON_mean = ax.plot(centre_of_bar_in_x_axis, freq_counts_normalised, color=df_lists_tab.loc[prot_list, 'color'],
                                                alpha=alpha, linewidth=linewidth,
                                                label=df_lists_tab.loc[prot_list, 'list_description'])

            ###############################################################
            #                                                             #
            #        annotate the mean values on the plots                #
            #                                                             #
            ###############################################################
            ax.annotate("mean: {: >2.0%}".format(perc_TMD_mean), [0.84, 0.1*ysh+ysh*offset], fontsize=anno_fontsize, fontproperties="monospace", color=color, xycoords=xyc, weight="semibold")

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
                      fontsize=fontsize-3, frameon=True, loc = 'upper right')# bbox_to_anchor=(1.07, 1.12))


        utils.save_figure(fig, Fig_name, base_filepath=base_filepath, save_png=save_png, save_pdf=save_pdf)

    #--------------------------------------------------------------------------------------------------------------------------------#
    Fig_Nr = 9
    if Fig_Nr in list_figs_to_run:
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
            color = df_lists_tab.loc[prot_list, 'color']
            # collect length data for all TMDs individually
            data = df_dict[prot_list]
            list_len_TMDs = []
            for n in range(1, data.number_of_TMDs_excl_SP.max().astype(int) + 1):
                TMD = 'TM%02d' % n
                list_len_TMDs.extend(data['%s_seqlen' % TMD].dropna())
            hist_data = np.array(list_len_TMDs)
            len_TMDs_mean = hist_data.mean()
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
            linecontainer_AAIMON_mean = ax.plot(centre_of_bar_in_x_axis, freq_counts_normalised, color=color,
                                                alpha=alpha, linewidth=linewidth,
                                                label=df_lists_tab.loc[prot_list, 'list_description'])

            ###############################################################
            #                                                             #
            #        annotate the mean values on the plots                #
            #                                                             #
            ###############################################################
            ax.annotate("mean: {: >7.02f}".format(len_TMDs_mean), [0.8, 0.1*ysh+ysh*offset], fontsize=anno_fontsize, fontproperties="monospace", color=color, xycoords=xyc, weight="semibold")

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
                      fontsize=fontsize - 3, frameon=True, loc="upper right")#bbox_to_anchor=(1.07, 1.12)

        utils.save_figure(fig, Fig_name, base_filepath=base_filepath, save_png=save_png, save_pdf=save_pdf)

    #--------------------------------------------------------------------------------------------------------------------------------#
    Fig_Nr = 10
    if Fig_Nr in list_figs_to_run:
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
            #sys.stdout.write("\nprot_list : {:02d}, fraction_TM_residues : {:0.03f}".format(prot_list, fraction_TM_residues))

            #create a new dataframe to hold the data for percentage identity and the normalisation factors
            df_norm = pd.DataFrame()
            # create datapoints for line, percent amino acid substitutions
            df_norm["perc_aa_sub"] = np.linspace(0,1.2,1000)
            # get the percentage identity
            df_norm["perc_aa_ident"] = 1 - df_norm["perc_aa_sub"]
            # calculate the norm factors
            df_norm["norm"] = df_norm["perc_aa_ident"].apply(korbinian.cons_ratio.norm.calc_aa_prop_norm_factor, args=(rand_TM, rand_nonTM, fraction_TM_residues))
            # get rid of any very high or low values
            df_norm = df_norm.loc[df_norm["norm"] > -10]
            df_norm = df_norm.loc[df_norm["norm"] < 10]

            # HYPERBOLIC GRAPH NEEDS TO BE SPLIT INTO THE TWO SIDES TO AVOID UGLY CROSSOVER EFFECT

            # get the indices for the highest and lowest norm factors
            index_max = df_norm["norm"].idxmax()
            index_min = df_norm["norm"].idxmin()

            label = df_lists_tab.loc[prot_list, 'list_description']

            # plot that section of the graph, from beginning to where values are high
            # convert to array to avoid auto legend)
            ax.plot(df_norm["perc_aa_sub"][:index_max]*100, np.array(df_norm["norm"][:index_max]), color=df_lists_tab.loc[prot_list, 'color'],
                    alpha=alpha, linewidth=linewidth,
                    label=None)
            # plot that section of the graph, from most negative value until the end
            ax.plot(df_norm["perc_aa_sub"][index_min:]*100, np.array(df_norm["norm"][index_min:]), color=df_lists_tab.loc[prot_list, 'color'],
                    alpha=alpha, linewidth=linewidth, label=label)

        ###############################################################
        #                                                             #
        #                       set up plot style                     #
        #                                                             #
        ###############################################################

        ax.set_xlabel('evolutionary distance (% substitutions)', fontsize=norm_fig_fontsize)
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
    if Fig_Nr in list_figs_to_run:
        title = 'compare_random_identity'
        Fig_name = 'Fig11_random_ident'
        #fig, ax = plt.subplots()
        fig, (ax, ax2) = plt.subplots(1, 2, sharey=False)

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


        ########################################################################################
        #                                                                                      #
        #                       Plot %LIVFA on the left side of the figure                  #
        #                                                                                      #
        ########################################################################################


        hydrophob_res = "LIVFA"
        list_hydrophob = list(hydrophob_res)
        TM_perc_LIVFA_for_each_dataset_ser = df_aa_prop_TM.loc[list_hydrophob, :].sum()
        nonTM_perc_LIVFA_for_each_dataset_ser = df_aa_prop_nonTM.loc[list_hydrophob, :].sum()

        # create list of list_numbers for x-axis
        x = TM_perc_LIVFA_for_each_dataset_ser.index.astype(str)
        # width of bar
        width = 0.4
        # indices (position of bar) are at 0.1, 1.1, 2.1
        ind = np.arange(len(x)) + 0.1
        # create bar for the rand_TM
        ax.bar(ind, TM_perc_LIVFA_for_each_dataset_ser*100, width, color=color_list_dark[:len(x)], alpha=1, label="TM")#color=color_list[2]
        # create bar for the rand_nonTM
        ax.bar(ind + width, nonTM_perc_LIVFA_for_each_dataset_ser*100, width, color=color_list[:len(x)], alpha=1, label="nonTM")
        # put labels in between bars
        ax.set_xticks(ind + width)
        # extract list names
        x_labels = df_lists_tab.loc[protein_lists, 'list_description']
        ax.set_xticklabels(x_labels, rotation=45, ha="right")

        if n_prot_lists == 3:
            ax.set_xlim(0, ind[-1]+ width*2.4)
        else:
            ax.set_xlim(-0.4, ind[-1] + width * 2.4)
        ax.legend(frameon=True)
        # turn off vertical grid
        ax.grid(axis="x", b=False)
        ax.set_ylabel("percent lipophilic residues (LIVFA)")

        ########################################################################################
        #                                                                                      #
        #                       Plot random identity on the RHS of the figure                  #
        #                                                                                      #
        ########################################################################################


        # create list of list_numbers for x-axis
        x = df_rand.index.astype(str)
        # width of bar
        width = 0.4
        # indices (position of bar) are at 0.1, 1.1, 2.1
        ind = np.arange(len(x)) + 0.1
        # create bar for the rand_TM
        ax2.bar(ind, df_rand.rand_TM, width, color=color_list_dark[:len(x)], alpha=1, label="TM")#color=color_list[2]
        # create bar for the rand_nonTM
        ax2.bar(ind + width, df_rand.rand_nonTM, width, color=color_list[:len(x)], alpha=1, label="nonTM")
        # put labels in between bars
        ax2.set_xticks(ind + width)
        # extract list names
        x_labels = df_lists_tab.loc[protein_lists, 'list_description']
        ax2.set_xticklabels(x_labels, rotation=45, ha="right")
        #ax2.set_xticklabels(x_labels)
        # turn off vertical grid
        ax2.grid(axis="x", b=False)
        ax2.legend(frameon=True)
        ax2.set_ylabel("random identity")

        fig.tight_layout()

        utils.save_figure(fig, Fig_name, base_filepath=base_filepath, save_png=save_png, save_pdf=save_pdf)


    #--------------------------------------------------------------------------------------------------------------------------------#
    Fig_Nr = 12
    if Fig_Nr in list_figs_to_run:
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
            ax.bar(ind_for_list, df_aa_prop_TM.loc[:,prot_list], width, color=color_list[m], label=df_lists_tab.loc[prot_list, 'list_description'])#color=color_list[2]

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
    if Fig_Nr in list_figs_to_run:
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
            ax.bar(ind_for_list, df_aa_prop_nonTM.loc[:,prot_list], width, color=color_list[m], label=df_lists_tab.loc[prot_list, 'list_description'])#color=color_list[2]

        xtick_locations = ind + n/2*width
        ax.set_xticks(xtick_locations)
        ax.set_xticklabels(aa_list)
        ax.set_ylabel("amino acid propensity")
        ax.legend()
        utils.save_figure(fig, Fig_name, base_filepath=base_filepath, save_png=save_png, save_pdf=save_pdf)
    #--------------------------------------------------------------------------------------------------------------------------------#

    Fig_Nr = 14
    if Fig_Nr in list_figs_to_run:
        def hist_lipo_first_vs_central_vs_last_TM_excl_GPCRs(df_dict, Fig_name):

            title = 'hist_lipo_first_vs_central_vs_last_TM_excl_GPCRs'
            min_ = -0.5
            max_ = 1.5
            binlist = np.linspace(min_, max_, n_bins_lipo + 1)

            fig, ax = plt.subplots()

            # move limits up, so they can be used for the annotations
            # x and y axes min and max
            xlim_min = min_
            xlim_max = max_
            ax.set_xlim(xlim_min, xlim_max)
            ylim_min = -0.01
            ylim_max = len(protein_lists) + 0.01

            # protein_lists_mp = []
            # for prot_list in protein_lists:
            #     for element in ['singlepass', 'multipass', 'betabarrel', 'SP', 'MP']:
            #         if element in df_lists_tab.loc[prot_list, 'list_description']:
            #             protein_lists_mp.append(prot_list)
            #
            #offset = len(protein_lists_mp) - 1
            offset = len(protein_lists) - 1

            for prot_list in protein_lists:
                # create a filtered dataframe (e.g. without GPCRs)
                df_filt = df_dict[prot_list]
                if "GPCR" in df_filt.columns:
                    df_filt = df_filt[df_filt.GPCR == False]
                # if it is multipass, remove proteins with 1-4 TMDs
                if df_lists_tab.loc[prot_list, "max_TMDs"] > 2:
                    df_filt = df_filt[df_filt['number_of_TMDs'] >= min_n_TMDs_first_last]

                # add "excl. GPCRs" to multipass label
                list_description = df_lists_tab.loc[prot_list, 'list_description']
                # make a list of classes where the excl. GPCR should be added
                if list_description in ["multipass"]:
                    list_description = list_description + " (excl. GPCRs)"

                # CREATE a boolean for singlepass and multipass
                is_multipass = True if df_lists_tab.loc[prot_list, "max_TMDs"] > 2 else False

                color = df_lists_tab.loc[prot_list, 'color']

                ###   lipo TM01   ###
                # create numpy array of membranous over nonmembranous conservation ratios (identity)
                hist_data = np.array(df_filt['TM01_lipo'].dropna())
                # get the mean value for annotations
                mean_lipo_TM01 = hist_data.mean()
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
                linecontainer_AAIMON_mean = ax.plot(centre_of_bar_in_x_axis, freq_counts_normalised, '--', color=color,
                                                    alpha=alpha, linewidth=linewidth, label=list_description)

                if df_lists_tab.loc[prot_list, "max_TMDs"] > 2:

                    ###   central TMDs   ###
                    # create numpy array of membranous over nonmembranous conservation ratios (identity)
                    hist_data = np.array(df_filt['lipo_mean_central_TMDs'].dropna())
                    # get the mean value for annotations
                    mean_lipo_central = hist_data.mean()
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
                    linecontainer_AAIMON_mean = ax.plot(centre_of_bar_in_x_axis, freq_counts_normalised, color=color,
                                                        alpha=alpha, linewidth=linewidth)

                    ###   last TM lipo   ###
                    # create numpy array of membranous over nonmembranous conservation ratios (identity)
                    hist_data = np.array(df_filt['lipo_last_TMD'].dropna())

                    if hist_data.shape[0] != 0:
                        # get the mean value for annotations
                        mean_lipo_last = hist_data.mean()
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
                        linecontainer_AAIMON_mean = ax.plot(centre_of_bar_in_x_axis, freq_counts_normalised, ':', color=color,
                                                            alpha=alpha,
                                                            linewidth=linewidth)
                    else:
                        sys.stdout.write("{} hist_data is empty".format(prot_list))

                ###############################################################
                #                                                             #
                #        annotate the mean values on the plots                #
                #                                                             #
                ###############################################################
                ax.annotate("mean:", [0.005, 0.9*ysh+ysh*offset], fontsize=anno_fontsize, fontproperties="monospace", color=color, xycoords=xyc, weight="semibold")
                ax.annotate("first{: >7.02f}".format(mean_lipo_TM01), [0.01, 0.8*ysh+ysh*offset], fontsize=anno_fontsize, fontproperties="monospace", color=color, xycoords=xyc, weight="semibold")
                if is_multipass:
                    ax.annotate("central{: >5.02f}".format(mean_lipo_central), [0.01, 0.7*ysh+ysh*offset], fontsize=anno_fontsize, fontproperties="monospace", color=color, xycoords=xyc, weight="semibold")
                    ax.annotate("last{: >8.02f}".format(mean_lipo_last), [0.01, 0.6*ysh+ysh*offset], fontsize=anno_fontsize, fontproperties="monospace", color=color, xycoords=xyc, weight="semibold")

                offset = offset - 1

            ###############################################################
            #                                                             #
            #                       set up plot style                     #
            #                                                             #
            ###############################################################

            ax.set_xlabel('lipophilicity (Hessa scale)', fontsize=fontsize)
            # move the x-axis label closer to the x-axis
            ax.xaxis.set_label_coords(0.5, -0.085)
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
                TM01 = plt.Line2D((0, 1), (0, 0), color='k', linestyle='--', linewidth=linewidth)
                central_TMs = plt.Line2D((0, 1), (0, 0), color='k', linestyle='-', linewidth=linewidth)
                last_TM = plt.Line2D((0, 1), (0, 0), color='k', linestyle=':', linewidth=linewidth)
                # Create legend from custom artist/label lists
                ax.legend([handle for i, handle in enumerate(handles) if i in display] + [TM01, central_TMs, last_TM],
                          [label for i, label in enumerate(labels) if i in display] + ['first TM', 'central TMs', 'last TM'],
                          fontsize=fontsize - 3, frameon=True, loc='upper right', handlelength=4)  # bbox_to_anchor=(1.07, 1.12))
            else:
                # Create custom artists
                TM01 = plt.Line2D((0, 1), (0, 0), color='k', linewidth=linewidth)
                central_TMs = plt.Line2D((0, 1), (0, 0), color='k', linestyle='-', linewidth=linewidth)
                last_TM = plt.Line2D((0, 1), (0, 0), color='k', linestyle=':', linewidth=linewidth)
                # Create legend from custom artist/label lists
                ax.legend([TM01, central_TMs, last_TM], ['first TM', 'central TMs', 'last TM'],
                          fontsize=fontsize - 3, frameon=True, loc='upper right', handlelength=4)  # bbox_to_anchor=(1.07, 1.12))

            utils.save_figure(fig, Fig_name, base_filepath=base_filepath, save_png=save_png, save_pdf=save_pdf)


        df_dict_list_csvs = {}
        for n, prot_list in enumerate(protein_lists):
            df_list = pd.read_csv(dfv.loc[prot_list, 'list_csv'], sep=",", quoting=csv.QUOTE_NONNUMERIC, index_col=0, low_memory=False)
            df_dict_list_csvs[prot_list] = df_list

        df_dict_list = [df_dict, df_dict_list_csvs]
        Fig_names = ['Fig14_hist_lipo_first_vs_central_vs_last_TM_excl_GPCRs', 'Fig14b_hist_lipo_f_c_l_all_prot_in_list_csv']
        for i in range(len(df_dict_list)):
            dict_for_graph = df_dict_list[i]
            Fig_name = Fig_names[i]
            hist_lipo_first_vs_central_vs_last_TM_excl_GPCRs(dict_for_graph, Fig_name)

    # --------------------------------------------------------------------------------------------------------------------------------#
    Fig_Nr = 15
    if Fig_Nr in list_figs_to_run:
        title = 'hist_AAIMON_slope_first_vs_central_vs_last_TM_excl_GPCRs'
        Fig_name = 'Fig15_hist_AAIMON_slope_first_vs_central_vs_last_TM_excl_GPCRs'
        binlist = np.linspace(-40, 40, n_bins_cons + 1)
        fig, ax = plt.subplots()
        offset = len(protein_lists) - 1

        # set up the axis min and max here, as it is necessary for annotations
        # x and y axes min and max
        xlim_min = -30
        xlim_max = 30
        ax.set_xlim(xlim_min, xlim_max)
        ylim_min = -0.01
        ylim_max = len(protein_lists) + 0.01
        ax.set_ylim(ylim_min, ylim_max)

        for prot_list in protein_lists:
            # create a filtered dataframe (e.g. without GPCRs)
            df_filt = df_dict[prot_list]
            if "GPCR" in df_filt.columns:
                df_filt = df_filt[df_filt.GPCR == False]
            # if it is multipass, remove proteins with 1-4 TMDs
            if df_lists_tab.loc[prot_list, "max_TMDs"] > 2:
                df_filt = df_filt[df_filt['number_of_TMDs'] >= min_n_TMDs_first_last]

            # add "excl. GPCRs" to multipass label
            list_description = df_lists_tab.loc[prot_list, 'list_description']
            # make a list of classes where the excl. GPCR should be added
            if list_description in ["multipass"]:
                list_description = list_description + " (excl. GPCRs)"

            is_multipass = True if df_lists_tab.loc[prot_list, "max_TMDs"] > 2 else False

            color = df_lists_tab.loc[prot_list, 'color']

            ###   AAIMON_slope TM01   ###
            # create numpy array of membranous over nonmembranous conservation ratios (identity)
            hist_data = (df_filt['TM01_AAIMON_slope'] * 1000).dropna()
            # get mean of first TMD
            mean_AAIMON_slope_TM01 = hist_data.mean()
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
            linecontainer_AAIMON_mean = ax.plot(centre_of_bar_in_x_axis, freq_counts_normalised, '--', color=color,
                                                alpha=alpha, linewidth=linewidth,
                                                label=df_lists_tab.loc[prot_list, 'list_description'])

            # purely singlepass datasets will probably not have 'AAIMON_slope_central_TMDs' in the columns
            if is_multipass:
                ###   AAIMON_slope central TMs   ###
                # create numpy array of membranous over nonmembranous conservation ratios (identity)
                hist_data = (df_filt['AAIMON_slope_central_TMDs'] * 1000).dropna()
                # get mean of central TMs
                mean_AAIMON_slope_central_TMDs = hist_data.mean()
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
                linecontainer_AAIMON_mean = ax.plot(centre_of_bar_in_x_axis, freq_counts_normalised, color=color,
                                                    alpha=alpha, linewidth=linewidth)


                ###   AAIMON_slope last TMD   ###
                # create numpy array of membranous over nonmembranous conservation ratios (identity)
                hist_data = (df_filt['AAIMON_slope_last_TMD'] * 1000).dropna()
                # get mean of last TMD
                mean_AAIMON_slope_last_TMD = hist_data.mean()
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
                linecontainer_AAIMON_mean = ax.plot(centre_of_bar_in_x_axis, freq_counts_normalised, ':', color=color,
                                                    alpha=alpha,
                                                    linewidth=linewidth)


            ###############################################################
            #                                                             #
            #        annotate the mean values on the plots                #
            #                                                             #
            ###############################################################
            ax.annotate("mean:", [0.01, 0.9*ysh+ysh*offset], fontsize=anno_fontsize, fontproperties="monospace", color=color, xycoords=xyc, weight="semibold")
            ax.annotate("first {: >6.02f}".format(mean_AAIMON_slope_TM01), [0.02, 0.8*ysh+ysh*offset], fontsize=anno_fontsize, fontproperties="monospace", color=color, xycoords=xyc, weight="semibold")
            if is_multipass:
                ax.annotate("central {: >4.02f}".format(mean_AAIMON_slope_central_TMDs), [0.02, 0.7*ysh+ysh*offset], fontsize=anno_fontsize, fontproperties="monospace", color=color, xycoords=xyc, weight="semibold")
                ax.annotate("last {: >7.02f}".format(mean_AAIMON_slope_last_TMD), [0.02, 0.6*ysh+ysh*offset], fontsize=anno_fontsize, fontproperties="monospace", color=color, xycoords=xyc, weight="semibold")

            offset = offset - 1

        ###############################################################
        #                                                             #
        #                       set up plot style                     #
        #                                                             #
        ###############################################################

        ax.set_xlabel(r'm$_{\rm TM/EM} *10^{\rm -3}$', fontsize=fontsize, labelpad=20)
        # move the x-axis label closer to the x-axis
        ax.xaxis.set_label_coords(0.5, -0.085)
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
            TM01 = plt.Line2D((0, 1), (0, 0), color='k', linestyle='-.', linewidth=linewidth)
            central_TM = plt.Line2D((0, 1), (0, 0), color='k', linewidth=linewidth)
            last_TM = plt.Line2D((0, 1), (0, 0), color='k', linestyle=':', linewidth=linewidth)
            # Create legend from custom artist/label lists
            ax.legend([handle for i, handle in enumerate(handles) if i in display] + [TM01, central_TM, last_TM],
                      [label for i, label in enumerate(labels) if i in display] + ['first TM', 'central TMs', 'last TM'],
                      fontsize=fontsize - 3, frameon=True, loc='upper right', handlelength=4)  # , bbox_to_anchor=(1.07, 1.12))
        else:
            # Create custom artists
            TM01 = plt.Line2D((0, 1), (0, 0), color='k', linestyle='-.', linewidth=linewidth)
            central_TM = plt.Line2D((0, 1), (0, 0), color='k', linewidth=linewidth)
            last_TM = plt.Line2D((0, 1), (0, 0), color='k', linestyle=':', linewidth=linewidth)
            # Create legend from custom artist/label lists
            ax.legend([TM01, central_TM, last_TM],
                      ['first TM', 'central TMs', 'last TM'],
                      fontsize=fontsize - 3, frameon=True, loc='upper right', handlelength=4)  # , bbox_to_anchor=(1.07, 1.12))
        utils.save_figure(fig, Fig_name, base_filepath, save_png, save_pdf)

    # --------------------------------------------------------------------------------------------------------------------------------#
    Fig_Nr = 16
    if Fig_Nr in list_figs_to_run:
        title = '% AA identity TMD vs EM'
        Fig_name = 'Fig16_perc_AA_identity_TMD_vs_nonTMD'
        sys.stdout.write("Fig16 may be empty until data is rerun.")

        fig, axes = plt.subplots(ncols=len(df_dict), nrows=1, sharex=True, sharey=True, figsize=(len(df_dict) * 2.5, 5))
        vmax = 15
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

            # data definition
            #min_ = df_dict[prot_list].perc_ident_mean.min()
            #max_ = df_dict[prot_list].perc_ident_mean.max()
            #sys.stdout.write("Fig16 prot_list {} : min_ = {}, max_ = {}".format(prot_list, min_, max_))

            data_TMD = df_dict[prot_list].TMD_perc_ident_mean * 100
            data_nonTMD = df_dict[prot_list].nonTMD_perc_ident_mean * 100

            #sys.stdout.write("data_TMD[0:5] {}".format(data_TMD[0:5]))
            #sys.stdout.flush()
            # histogram the data
            hh, locx, locy = scipy.histogram2d(data_TMD, data_nonTMD, range=xyrange, bins=bins)

            # density threshold
            thresh = 1
            # fill the areas with low density by NaNs
            hh[hh < thresh] = np.nan

            # if there are just a few datapoints, multiply count to get visible values
            if len(data_TMD) < 50:
                hh = hh * 5
                ax.annotate('count in bin $*5$', fontsize=5, xy=(69, 1))

            # plot the data with imshow
            im = ax.imshow(np.flipud(hh.T), cmap='Oranges', extent=np.array(xyrange).flatten(),
                           interpolation='none', origin='upper', aspect='equal', vmin=1, vmax=vmax)

            ax.annotate(df_lists_tab.loc[prot_list, 'list_description'], xy=(3, 90), color=df_lists_tab.loc[prot_list, 'color'], fontsize=10)
            ax.tick_params(labelsize=10)
            ax.annotate('TMD less conserved', xy=(1, 7), rotation=90, fontsize=8, ha='left', va='bottom')
            ax.annotate('TMD more conserved', xy=(7, 1), rotation=0, fontsize=8, ha='left', va='bottom')

            if n == 0:
                ax.set_ylabel('mean % identity\nof EM region', fontsize=10)

        # get colorbar from latest imshow element (color scale should be the same for all subplots)
        fig.subplots_adjust(right=0.8)
        cbar_ax = fig.add_axes([0.81, 0.35, 0.04 / len(df_dict), 0.3])
        fig.colorbar(im, cax=cbar_ax)
        labels = cbar_ax.get_ymajorticklabels()
        labels[-1] = '>{}'.format(vmax)
        cbar_ax.set_yticklabels(labels)

        plt.subplots_adjust(wspace=0.15, hspace=0.05)
        fig.text(0.465, 0.27, 'mean % identity of TM region amongst homologues', ha='center', fontsize=10)

        utils.save_figure(fig, Fig_name, base_filepath, save_png, save_pdf)

    #---------------------------------------------------------------------------------------------------------------------------------#
    Fig_Nr = 17
    if Fig_Nr in list_figs_to_run:
        title = 'number of TMDs in GPCRs'
        Fig_name = 'Fig17_number_of_TMDs_in_GPCRs'
        if len(df_dict) > 4:
            sys.stdout.write('cannot create plot, too many lists')
        else:
            list_GPCR_KW = ['G-protein coupled receptor']
            for prot_list in protein_lists:
                if 'GPCR' not in df_dict[prot_list].columns:
                    df_dict[prot_list]['GPCR'] = df_dict[prot_list]['uniprot_KW'].apply(korbinian.cons_ratio.keywords.KW_list_contains_any_desired_KW, args=(list_GPCR_KW,))
            a_list_contains_GPCRs = False
            for prot_list in protein_lists:
                if True in df_dict[prot_list]['GPCR']:
                    a_list_contains_GPCRs = True
                    break
            if a_list_contains_GPCRs:
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
                        df_GPCR = df_dict[prot_list][df_dict[prot_list].GPCR == True]
                        # check if df_GPCR contains elements
                        if df_GPCR.shape[0] == 0:
                            continue

                        hist_data = df_GPCR.number_of_TMDs_excl_SP.value_counts().values
                        hist_data_norm = hist_data / hist_data.max()
                        bins = df_GPCR.number_of_TMDs_excl_SP.value_counts().index

                        ax.bar(bins, hist_data_norm, align='center', color=df_lists_tab.loc[prot_list, 'color'])
                        ax.annotate(df_lists_tab.loc[prot_list, 'list_description'], xy=(0.5, 1.05), color=df_lists_tab.loc[prot_list, 'color'], fontsize=fontsize - 2)
                        ax.annotate('only GPCRs\nn = {}'.format(len(df_GPCR)), xy=(0.5, 0.80), fontsize=fontsize - 2, alpha=0.5)

                plt.tight_layout()
                utils.save_figure(fig, Fig_name, base_filepath, save_png, save_pdf)
            else:
                sys.stdout.write("Fig17 skipped, no GPCRs in datasets.")
                sys.stdout.flush()

    # ---------------------------------------------------------------------------------------------------------------------------------#
    Fig_Nr = 18
    if Fig_Nr in list_figs_to_run:
        title = 'delta mTM/EM normalised vs non-normalised'
        Fig_name = 'Fig18_d_mTM_nonTM_norm_vs_nonnorm'

        # difference of AAIMON and AAIMON_n
        fig, ax = plt.subplots()

        mean_diff = []
        for n, prot_list in enumerate(protein_lists):
            slope = df_dict[prot_list].AAIMON_slope_all_TMDs_mean * 1000
            slope_n = df_dict[prot_list].AAIMON_n_slope_all_TMDs_mean * 1000
            diff = np.mean(abs(slope - slope_n))
            mean_diff.append(diff)

        ax.bar(range(1, len(protein_lists) + 1), mean_diff, edgecolor='k', color='#0076B8')
        ax.set_xticks(range(1, len(protein_lists) + 1))
        labels = [element.replace('_', '\n') for element in df_lists_tab.loc[protein_lists, "list_description"].tolist()]
        ax.set_xticklabels(labels)
        ax.set_ylabel(r'$\Delta$ m$_{\rm TM/EM} *10^{\rm -3}$ norm. vs. non-norm.', fontsize=fontsize)
        ax.tick_params(labelsize=fontsize)
        plt.tight_layout()

        utils.save_figure(fig, Fig_name, base_filepath, save_png, save_pdf)

    # ---------------------------------------------------------------------------------------------------------------------------------#
    Fig_Nr = 19
    if Fig_Nr in list_figs_to_run:
        title = 'heatmap deltas between datasets'
        Fig_name = 'Fig19_Heatmap_d_AAIMON_slope_between_datasets'

        # add mean values to dfv
        for prot_list in protein_lists:
            dfv.loc[prot_list, 'mean_AAIMON_slope'] = np.mean(df_dict[prot_list].AAIMON_slope_all_TMDs_mean * 1000)
        # initialise empty lists to hold data, annotations and labels
        heatmap_data = []
        heatmap_text = []
        labels = []
        for prot_list_1 in protein_lists:
            labels.append('List{:02d}'.format(prot_list_1))
            list_data = []
            list_text = []
            for prot_list_2 in protein_lists:
                data = (dfv.loc[prot_list_2, 'mean_AAIMON_slope'] - dfv.loc[prot_list_1, 'mean_AAIMON_slope'])
                list_text.append('{:.02f}'.format(data))
                list_data.append(data)
            heatmap_data.append(list_data)
            heatmap_text.append(list_text)
        heatmap_data = np.array(heatmap_data)
        heatmap_text = np.array(heatmap_text)

        fig, ax = plt.subplots()

        diagonal = np.zeros_like(heatmap_data, dtype=np.bool)
        np.fill_diagonal(diagonal, True)
        label = [element.replace('_', '\n') for element in df_lists_tab.loc[protein_lists, "list_description"].tolist()]
        ax = sns.heatmap(heatmap_data, annot=heatmap_text, fmt="", square=True, mask=diagonal, linewidths=2,
                         annot_kws={"size": fontsize + 2})
        ax.set_xticklabels(label, fontsize=fontsize)
        ax.set_yticklabels(reversed(label), rotation=0, fontsize=fontsize)
        # get colorbar axis for labelling
        cbar_ax = ax.collections[0].colorbar
        cbar_ax.set_label(r'$\Delta$ m$_{\rm TM/EM}*10^{\rm -3}$ between datasets', fontsize=fontsize)
        cbar_ax.ax.tick_params(labelsize=fontsize)
        plt.tight_layout()

        utils.save_figure(fig, Fig_name, base_filepath, save_png, save_pdf)

    # ---------------------------------------------------------------------------------------------------------------------------------#
    Fig_Nr = 20
    if Fig_Nr in list_figs_to_run:
        title = 'compare percentage of different keywords'
        Fig_name = 'Fig20_perc_of_different_keywords' # optimised for two lists containing keywords

        KW_to_compare = ['Cell membrane', 'Endoplasmic reticulum', 'Golgi apparatus', 'Enzyme', 'G-protein coupled receptor']
        if 'Enzyme' in KW_to_compare:
            compare_enzyme = True
            KW_to_compare.remove('Enzyme')
            list_enzyme, list_ignored, PFAM_dict = korbinian.cons_ratio.keywords.get_list_enzyme_KW_and_list_ignored_KW()
        else:
            compare_enzyme = False
        protein_lists_KW = []
        df = pd.DataFrame(index=dfv.index)
        #df['list_description'] = dfv['list_description']
        df['list_description'] = list_descriptions
        df['color'] = color_list[0:n_prot_lists]
        df['color dark'] = color_list_dark[0:n_prot_lists]

        for n, prot_list in enumerate(protein_lists):
            # skip lists that don't contain uniprot_KW
            if not 'uniprot_KW' in df_dict[prot_list].columns:
                continue
            protein_lists_KW.append(prot_list)
            # if the keyword list is a stringlist, change to python list
            if type(df_dict[prot_list].ix[0, 'uniprot_KW']) == str:
                df_dict[prot_list]['uniprot_KW'] = df_dict[prot_list]['uniprot_KW'].apply(lambda x: ast.literal_eval(x))
            # calculate percentage of dataset with KW
            df_dict[prot_list]['perc_of_dataset'] = 100 / len(df_dict[prot_list])
            for keyword in KW_to_compare:
                df_dict[prot_list][keyword] = df_dict[prot_list].uniprot_KW.apply(korbinian.cons_ratio.keywords.KW_list_contains_any_desired_KW, args=([keyword],))
                df.loc[prot_list, '_'.join(keyword.split(' '))] = sum(df_dict[prot_list][df_dict[prot_list][keyword] == True]['perc_of_dataset'])
            if compare_enzyme == True:
                df_dict[prot_list]['Enzyme'] = df_dict[prot_list].uniprot_KW.apply(korbinian.cons_ratio.keywords.KW_list_contains_any_desired_KW, args=(list_enzyme,))
                df.loc[prot_list, 'Enzyme'] = sum(df_dict[prot_list][df_dict[prot_list]['Enzyme'] == True]['perc_of_dataset'])

        df = df.set_index(['list_description']).dropna().T

        if "G-protein_coupled_receptor" in df.index:
            df.drop("G-protein_coupled_receptor", axis=0, inplace=True)
            sys.stdout.write("G-protein_coupled_receptor dropped from KW analysis")

        df.index = df.index.str.replace("_", "\n")

        if protein_lists_KW != []:
            manual = False
            if manual:
                fig, ax = plt.subplots()
                n = len(protein_lists_KW)
                # width of bar
                width = 0.7 / n
                # indices (position of bar)
                ind = np.arange(df.shape[0])
                for m, prot_list in enumerate(protein_lists_KW):
                    list_description = df_lists_tab.loc[prot_list, 'list_description']
                    ind_for_list = ind + width * m
                    ax.bar(ind_for_list, df[list_description], width=width, color=df_lists_tab.loc[prot_list, 'color'], label=df_lists_tab.loc[prot_list, 'list_description'],
                           edgecolor=df_lists_tab.loc[prot_list, 'color_dark'])
                xtick_locations = ind + width# / 2  # - n/2*width
                ax.set_xticks(xtick_locations)
                label = [element.replace('_', '\n') for element in df.index]
                ax.set_xticklabels(label, rotation=0, ha="center")
                ax.tick_params(labelsize=fontsize, pad=3)
                ax.set_ylabel('% of dataset', fontsize=fontsize)
                ax.legend(frameon=True, fontsize=fontsize)
                plt.tight_layout()
                utils.save_figure(fig, Fig_name, base_filepath, save_png, save_pdf)

            else:
                fig, ax = plt.subplots()
                df_to_plot = df.iloc[2:,:]
                df_to_plot.plot(kind="bar", ax=ax, color=df.loc["color", :], width=0.7)#, edgecolor=df.loc["color dark", :]
                ax.set_xticklabels(df_to_plot.index, rotation=0)
                ax.legend(frameon=True, fontsize=fontsize)
                ax.grid(axis="x", b=False)
                ax.set_ylabel('percentage of dataset', fontsize=fontsize)
                plt.tight_layout()
                utils.save_figure(fig, Fig_name, base_filepath, save_png, save_pdf)


        else:
            sys.stdout.write('no uniprot keywords found, {} cannot be processed!'. format(Fig_name))


    # ---------------------------------------------------------------------------------------------------------------------------------#
    Fig_Nr = 21
    if Fig_Nr in list_figs_to_run:
        title = 'barchart comparing mean slopes'
        Fig_name = 'Fig21_barchart_comparing_mean_slopes'

        fig, ax = plt.subplots(figsize=(3, 4.5))

        for n, prot_list in enumerate(protein_lists):
            # make sure that the first column cluster centers around 1
            n += 1
            # define data
            df = df_dict[prot_list]
            edgecolor = df_lists_tab.loc[prot_list, "color_dark"]
            color = df_lists_tab.loc[prot_list, "color"]


            # calculate mean, SEM and plot data for non-normalised slopes
            data_slope = df.AAIMON_slope_all_TM_res.dropna() * 1000
            mean_slope = np.mean(data_slope)
            sem_slope = stats.sem(data_slope)
            width = 0.4
            if n == 1:
                label = "original"
            else:
                label = None
            ax.bar(n - 0.2, mean_slope, color=color, width=width, align="center", label=label, edgecolor=edgecolor)#, yerr=std_slope
            ax.errorbar(n - 0.2, mean_slope, yerr=sem_slope, fmt="none", ecolor="k", ls="none", capthick=1, elinewidth=1, capsize=4, label=None)

            # calculate mean, SEM and plot data for normalised slopes
            data_slope_n = df.AAIMON_n_slope_all_TM_res.dropna() * 1000
            mean_slope_n = np.mean(data_slope_n)
            sem_slope_n = stats.sem(data_slope_n)
            #ax.bar(n + 0.2, mean_slope_n, color=color_list[prot_list - 1], width=0.4, yerr=sem_slope_n, alpha=0.4)
            if n == 1:
                label = "normalised"
            else:
                label = None
            # # for the normalised data, create a bar that is lighter than the original data
            # #lighter_colour = np.array(utils.HTMLColorToRGB(color_list[n]))/255 + 0.2
            # lighter_colour = color_list[n] + 0.2
            # # replace all values above 1 with 1
            # lighter_colour[lighter_colour > 1] = 1
            ax.bar(n + 0.2, mean_slope_n, color=color, width=width, align="center", label=label, hatch="////", edgecolor=edgecolor)  # , yerr=std_slope , alpha=0.4
            ax.errorbar(n + 0.2, mean_slope_n, yerr=sem_slope_n, fmt="none", ecolor="k", ls="none", capthick=1, elinewidth=1, capsize=4, label=None)


        # set x-axis limits
        ax.set_xlim(0.5, len(protein_lists) + 0.5)
        # set custom x-axis ticks and labels
        ax.set_xticks(range(1, len(protein_lists) + 1))
        ax.set_xticklabels(df_lists_tab.list_description, fontsize=fontsize, ha="right", rotation=45)
        # set y-axis label
        ax.set_ylabel(r'mean m$_{\rm TM/EM} *10^{\rm -3}$', fontsize=fontsize)
        # set y-limits - can be changed to fit purposes
        #ax.set_ylim(-1, 7)

        ax.tick_params(labelsize=fontsize-3, pad=2)
        #ax.legend(['non-normalised', 'normalised'], frameon=True, fontsize=fontsize - 3)#, loc='upper center')
        leg = ax.legend(frameon=True, fontsize=fontsize - 3)
        leg.legendHandles[0].set_color('0.5')
        leg.legendHandles[1].set_color('0.5')

        import matplotlib.patches as mpatches
        circ1 = mpatches.Patch(facecolor="0.5", label='original', edgecolor="k")
        circ2 = mpatches.Patch(facecolor="0.5", label='normalised', hatch="////", edgecolor="k")

        ax.legend(handles=[circ1, circ2], frameon=True)

        ax.yaxis.set_label_coords(-0.07, 0.5)
        ax.grid(False)
        plt.tight_layout()

        utils.save_figure(fig, Fig_name, base_filepath, save_png, save_pdf)



    # ---------------------------------------------------------------------------------------------------------------------------------#
    Fig_Nr = 22
    if Fig_Nr in list_figs_to_run:
        Fig_name = 'Fig22_linechart_comparing_AAIMON_over_evol_dist'

        data_dir = r"D:\Databases"
        fig, ax = plt.subplots(figsize=(6.84, 6.84))
        linestyles = ["-", ":", "--", "-.", "-"]*10

        list_max_evol_dist = []

        for n, prot_list in enumerate(protein_lists):

            base_filename_summaries = os.path.join(data_dir, "summaries", "{:02d}".format(prot_list), "List{:02d}".format(prot_list))
            in_zipfile = '%s_characterising_each_homol_TMD.zip' % base_filename_summaries
            # read data from disk
            if os.path.isfile(in_zipfile):
                with zipfile.ZipFile(in_zipfile, "r", zipfile.ZIP_DEFLATED) as openzip:
                    binned_data = pickle.load(openzip.open("binned_data_characterising_each_homol_TMD.pickle", "r"))

            x_perc_subst_data = binned_data[:, 0]
            y_AAIMON = binned_data[:, 1]

            # get number of valid proteins
            n_prot = df_dict[prot_list].shape[0]
            # add number of proteins to label
            label = "{}, n={}".format(df_lists_tab.loc[prot_list, 'list_description'], n_prot)

            ax.plot(x_perc_subst_data, y_AAIMON, label=label, linestyle=linestyles[n])

            # get the largest bin on the x-axis (first datapoint). This is used for xlim.
            list_max_evol_dist.append(x_perc_subst_data[0])


        # get the max value for the x-axis limit, based on last perc subst. bin in all datasets
        max_perc_ident_all_datasets = np.ceil(np.array(list_max_evol_dist).max())

        ax.legend(frameon=True, fontsize=fontsize)
        ax.set_xlabel("evolutionary distance (% substitutions)", fontsize=fontsize)
        ax.set_ylabel("TM/EM conservation ratio", fontsize=fontsize)
        ax.set_xlim(0, max_perc_ident_all_datasets)

        if protein_lists == (1, 12, 13, 14, 3):
            ax.set_ylim(0.9, 1.5)

        # add annotations
        ax.annotate(s="TM less\nconserved", xy=(-0.1, 0.1), fontsize=fontsize-2, xytext=None, xycoords='axes fraction', rotation=90)
        ax.annotate(s="TM more\nconserved", xy=(-0.1, 0.9), fontsize=fontsize-2, xytext=None, xycoords='axes fraction', rotation=90)

        #ax.autoscale_view(tight=True, scalex=True, scaley=True)
        utils.save_figure(fig, Fig_name, base_filepath, save_png, save_pdf)

    # ---------------------------------------------------------------------------------------------------------------------------------#
    def create_boxplot(col, ylabel, df_dict, protein_lists, list_descriptions, mult=1, add_xlabel=False, log=False, sym="o"):
        colors = korbinian.utils.create_colour_lists()['BGO_arr']
        darker_colors = [utils.darken_or_lighten(c, -0.33) for c in colors]
        # colors = [tuple(c) for c in colors]
        double_colors = []
        for i in darker_colors:
            double_colors.extend([i, i])

        nested_data = []
        for n in protein_lists:
            dft = df_dict[n]
            # multiply by a value?
            data = dft[col] * mult
            nested_data.append(data.tolist())
        plt.close("all")
        fig, ax = plt.subplots(figsize=(5, 5))
        meanpointprops = dict(marker='.', markerfacecolor='k', markeredgecolor="k", markersize=4)
        boxplotcontainer = ax.boxplot(nested_data, showmeans=True, meanprops=meanpointprops, patch_artist=True,
                                      notch=True, bootstrap=2000, widths=0.7)
        ax.set_ylabel(ylabel)
        if add_xlabel:
            #ax.set_xticklabels(df.name, rotation=90)
            ax.set_xticklabels(list_descriptions, rotation=90)

        # ax.grid(False)
        if log:
            ax.set_yscale('log')
        # fill with colors
        # colors = color_list*3

        for n, patch in enumerate(boxplotcontainer["boxes"]):
            # fc = darker_colors[n]
            # lighter_colour = np.array(korbinian.utils.HTMLColorToRGB(fc))/255 + 0.3

            # lighter_colour = fc + 0.3
            # replace all values above 1 with 1
            # lighter_colour[lighter_colour > 1] = 1
            # plt.setp(boxplotcontainer["facecolor"], color=fc)
            patch.set_color(darker_colors[n])
            patch.set_facecolor(colors[n])

        for element in ["whiskers", "caps"]:
            for n, patch in enumerate(boxplotcontainer[element]):
                fc = double_colors[n]
                patch.set_color(fc)

        ## change color and linewidth of the medians
        for n, median in enumerate(boxplotcontainer['medians']):
            median.set_color(darker_colors[n])
        for n, mean in enumerate(boxplotcontainer['means']):
            mean.set_markeredgecolor(darker_colors[n])
            mean.set_markerfacecolor(darker_colors[n])

        for n, flier in enumerate(boxplotcontainer['fliers']):
            flier.set(marker="o", alpha=0.8, markersize=3)
            flier.set_markeredgecolor(darker_colors[n])
            # flier.set_linewidth(0)
            # flier.set_fillstyle("full")
            # flier.set_markerfacecolor(colors[n])

            # change the style of fliers and their fill
            # for flier in boxplotcontainer['fliers']:
            #    flier.set(marker=sym, color='0.8', alpha=0.5, markerfacecolor='k', markersize=3)
            # plt.setp(boxplotcontainer["fliers"], color="red")

            #     for patch, color in zip(boxplotcontainer["boxes"], colors):
            #         patch.set_facecolor(color)
            #         patch.set_color(color)
        ax.grid(axis="x", b=False)
        fig.tight_layout()
        return fig, ax

    Fig_Nr = 23
    if Fig_Nr in list_figs_to_run:
        Fig_name = 'Fig23_boxplot_AAIMON'
        #df, protein_lists, col, ylabel, df_dict, list_descriptions,
        fig, ax = create_boxplot("AAIMON_mean_all_TM_res", r'TM/EM conservation ratio', df_dict, protein_lists, list_descriptions, add_xlabel=False)
        #fig.savefig("AAIMON_slope_all_TM_res.png")
        #fig.savefig("AAIMON_slope_all_TM_res.pdf")

        utils.save_figure(fig, Fig_name, base_filepath, save_png, save_pdf)

    Fig_Nr = 24
    if Fig_Nr in list_figs_to_run:
        Fig_name = 'Fig24_boxplot_AAIMON_slope'
        #df, protein_lists, col, ylabel, df_dict, list_descriptions,
        fig, ax = create_boxplot("AAIMON_slope_all_TM_res", r'm$_{\rm TM/EM} *10^{\rm -3}$', df_dict, protein_lists, list_descriptions, mult=1000, add_xlabel=True)
        #fig.savefig("AAIMON_slope_all_TM_res.png")
        #fig.savefig("AAIMON_slope_all_TM_res.pdf")

        utils.save_figure(fig, Fig_name, base_filepath, save_png, save_pdf)


    Fig_Nr = 25
    if Fig_Nr in list_figs_to_run:
        Fig_name = 'Fig25_boxplot_lipo_mean_all_TM_res'
        fig, ax = create_boxplot("lipo_mean_all_TM_res", "mean TM lipophilicity", df_dict, protein_lists, list_descriptions, add_xlabel=True)
        # add annotations
        fontsize = 8
        ax.annotate(s="more\nlipophilic", xy=(-0.14, 0.1), fontsize=fontsize, xytext=None, xycoords='axes fraction', rotation=90)
        ax.annotate(s="less\nlipophilic", xy=(-0.08, 0.9), horizontalalignment="right", fontsize=fontsize, xytext=None, xycoords='axes fraction', rotation=90)
        ax.set_ylabel("mean TM lipophilicity", labelpad=12)
        utils.save_figure(fig, Fig_name, base_filepath, save_png, save_pdf)

    Fig_Nr = 26
    if Fig_Nr in list_figs_to_run:
        Fig_name = 'Fig26_boxplot_seqlen'
        fig, ax = create_boxplot("seqlen", "sequence length", df_dict, protein_lists, list_descriptions, add_xlabel=True, log=True)
        utils.save_figure(fig, Fig_name, base_filepath, save_png, save_pdf)

    Fig_Nr = 27
    if Fig_Nr in list_figs_to_run:
        Fig_name = 'Fig27_boxplot_evol_dist'
        fig, ax = create_boxplot("obs_changes_mean", "evolutionary distance\n(% substitutions)", df_dict, protein_lists, list_descriptions, add_xlabel=True)
        utils.save_figure(fig, Fig_name, base_filepath, save_png, save_pdf)

    Fig_Nr = 28
    if Fig_Nr in list_figs_to_run:
        Fig_name = 'Fig28_boxplot_n_homol'
        fig, ax = create_boxplot("AAIMON_n_homol", "number of homologues", df_dict, protein_lists, list_descriptions, add_xlabel=True, log=True)
        #ax.set_yscale("log")
        utils.save_figure(fig, Fig_name, base_filepath, save_png, save_pdf)



    #dfv.to_csv(os.path.join(base_filepath, 'Lists_%s_variables.csv'%str_protein_lists))
    sys.stdout.write("\n~~~~~~~~~~~~         compare_lists finished           ~~~~~~~~~~~~\n")