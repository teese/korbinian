from scipy.stats import ttest_ind
import ast
import csv
import itertools
import korbinian.utils as utils
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
import zipfile
import os
import pickle

def save_figures_describing_proteins_in_list(pathdict, s, logging):
    logging.info("~~~~~~~~~~~~         starting run_save_figures_describing_proteins_in_list          ~~~~~~~~~~~~")
    # define save settings
    save_png = s["save_png"]
    save_pdf = s["save_pdf"]
    base_filepath = pathdict["single_list_fig_path"]
    # set resolution for plots in png format
    dpi = 300

    print ('Preparing data for plotting', flush=True)

    backgroundcolour = '0.95'
    plt.style.use('seaborn-whitegrid')

    # set default font size for plot
    fontsize = 8
    datapointsize = 8
    #alpha = 0.1

    '''Prepare data for the following plots'''

    # load cr_summary file
    dfc = pd.read_csv(pathdict["list_cr_summary_csv"], sep=",", quoting=csv.QUOTE_NONNUMERIC, index_col=0)
    # load summary file
    dfu = pd.read_csv(pathdict["list_summary_csv"], sep=",", quoting=csv.QUOTE_NONNUMERIC, index_col=0)
    # merge cr_summary and summary file, if columns are equal in both files, suffix _dfc will be added in cr_summary column names for backwards compatibility
    df = pd.merge(dfc, dfu, left_index=True, right_index=True, suffixes=('_dfc', ''))

    # create number of datapoint dependent alpha_dpd
    alpha_dpd = utils.calc_alpha_from_datapoints(df['AAIMON_mean_all_TMDs'])
    sys.stdout.write('opacity of datapoints: {a:.2f}\n'.format(a=alpha_dpd))
    # filter to remove proteins that have less than ~5 homologues
    # this is only important for the beta-barrel dataset, which has a lot of these proteins!
    min_n_homol = s["min_n_homol_for_figs"]
    n_prot_before_n_homol_cutoff = df.shape[0]
    df = df.loc[df['TM01_AAIMON_n_homol'] >= min_n_homol]
    n_prot_after_n_homol_cutoff = df.shape[0]
    n_removed = n_prot_before_n_homol_cutoff - n_prot_after_n_homol_cutoff
    # if any proteins have been removed, then print the exact number.
    if n_removed >= 1:
        sys.stdout.write("{}/{} proteins were removed, as they contained less than {} valid homologues. "
              "Final number of proteins = {}".format(n_removed, n_prot_before_n_homol_cutoff, min_n_homol, n_prot_after_n_homol_cutoff))

    # open list_summary_csv file
    #df_uniprot = pd.read_csv(pathdict["list_summary_csv"], sep=",", quoting=csv.QUOTE_NONNUMERIC, index_col=0)

    if 'uniprot_KW' in df.columns:
        # convert the keywords from a stringlist to a python list
        if isinstance(df['uniprot_KW'][0], str):
            df['uniprot_KW'] = df['uniprot_KW'].apply(lambda x: ast.literal_eval(x))
        # create a new column showing whether the protein is a GPCR
        df['G-protein_coupled_receptor'] = df['uniprot_KW'].apply(lambda x: 'G-protein coupled receptor' in x)
        df_GPCR = df.loc[df['G-protein_coupled_receptor'] == True]

        # bool to check if dataframe contains GPCRs
        if df['G-protein_coupled_receptor'].any():
            GPCR_in_df = True
        else:
            GPCR_in_df = False

    else:
        sys.stdout.write('No uniprot keywords available! cannot create figures 19-21 \n')

    # # save dataframe
    # df.to_csv(pathdict["base_filename_summaries"] + '_df_figs.csv', sep=",", quoting=csv.QUOTE_NONNUMERIC)

    # create binlist
    linspace_binlist = np.linspace(s["mp_smallest_bin"],
                                   s["mp_largest_bin"],
                                   s["mp_number_of_bins"])

    # add 30 as the last bin, to make sure 100% of the data is added to the histogram, including major outliers
    binlist = np.append(linspace_binlist, s["mp_final_highest_bin"])

    # THE MEAN CALCULATIONS ARE QUITE SLOW FOR THOUSANDS OF PROTEINS,
    # # iterate through the proteins that have a list of TMDs
    # for acc in df.loc[df['list_of_TMDs'].notnull()].loc[df['list_of_TMDs'] != 'nan'].index:
    #     dict_AAIMON_mean = {}
    #     dict_AAIMON_std = {}
    #     dict_AASMON_ratio_mean = {}
    #     dict_AASMON_ratio_std = {}
    #     for TMD in ast.literal_eval(df.loc[acc, 'list_of_TMDs']):
    #         dict_AAIMON_mean[TMD] = df.loc[acc, '%s_AAIMON_mean' % TMD]
    #         dict_AAIMON_std[TMD] = df.loc[acc, '%s_AAIMON_std' % TMD]
    #         dict_AASMON_ratio_mean[TMD] = df.loc[acc, '%s_AASMON_ratio_mean' % TMD]
    #         dict_AASMON_ratio_std[TMD] = df.loc[acc, '%s_AASMON_ratio_std' % TMD]
    #     df.loc[acc, 'AAIMON_mean_all_TMDs'] = np.mean(list(dict_AAIMON_mean.values()))
    #     df.loc[acc, 'AAIMON_std_all_TMDs'] = np.mean(list(dict_AAIMON_std.values()))
    #     df.loc[acc, 'AASMON_ratio_mean_all_TMDs'] = np.mean(list(dict_AASMON_ratio_mean.values()))
    #     df.loc[acc, 'AASMON_ratio_std_all_TMDs'] = np.mean(list(dict_AASMON_ratio_std.values()))

    # create list of colours to use in figures
    colour_lists = utils.create_colour_lists()
    tableau20 = colour_lists['tableau20']

    # create dataframe mean_AAIMON_each_TM
    df_mean_AAIMON_each_TM = pd.DataFrame()
    # add AAIMON each TMD to dataframe
    for acc in df.index:
        for TMD in ast.literal_eval(df.loc[acc, 'list_of_TMDs']):
            df_mean_AAIMON_each_TM.loc[acc, '{a}_AAIMON_mean'.format(a=TMD)] = df.loc[
                acc, '{b}_AAIMON_mean'.format(b=TMD)]

    # logging saved data types
    sys.stdout.write('Saving figures as: ')
    if s['save_pdf']:
        sys.stdout.write(' .pdf ')
    if s['save_png']:
        sys.stdout.write(' .png ')
    sys.stdout.write('\n')

    if s['Fig01_Histogram_of_mean_AAIMON_and_AASMON_ratios_SP_vs_MP']:
        Fig_Nr = 1
        title = 'Mean ratios'
        Fig_name = 'Fig01_Histogram_of_mean_AAIMON_and_AASMON_ratios_SP_vs_MP'
        # create a new figure
        fig, ax = plt.subplots()
        # create numpy array of membranous over nonmembranous conservation ratios (identity)
        hist_data_AAIMON_mean = np.array(df['AAIMON_mean_all_TMDs'].dropna())
        # use numpy to create a histogram
        freq_counts_I, bin_array_I = np.histogram(hist_data_AAIMON_mean, bins=binlist)
        # assuming all of the bins are exactly the same size, make the width of the column equal to 70% of each bin
        col_width = float('%0.3f' % (0.95 * (bin_array_I[1] - bin_array_I[0])))
        # when align='center', the central point of the bar in the x-axis is simply the middle of the bins ((bin_0-bin_1)/2, etc)
        centre_of_bar_in_x_axis = (bin_array_I[:-2] + bin_array_I[1:-1]) / 2
        # add the final bin, which is physically located just after the last regular bin but represents all higher values
        bar_width = centre_of_bar_in_x_axis[3] - centre_of_bar_in_x_axis[2]
        centre_of_bar_in_x_axis = np.append(centre_of_bar_in_x_axis, centre_of_bar_in_x_axis[-1] + bar_width)
        barcontainer_AAIMON_mean = ax.bar(left=centre_of_bar_in_x_axis, height=freq_counts_I,
                                          align='center', width=col_width, color="#0489B1",
                                          alpha=0.5)  # edgecolor='black',
        # create numpy array of normalised membranous over nonmembranous conservation ratios (identity)
        hist_data_AAIMON_n_mean = np.array(df['AAIMON_mean_all_TMDs_n'].dropna())
        # use numpy to create a histogram
        freq_counts_I, bin_array_I = np.histogram(hist_data_AAIMON_n_mean, bins=binlist)
        # assuming all of the bins are exactly the same size, make the width of the column equal to 70% of each bin
        col_width = float('%0.3f' % (0.95 * (bin_array_I[1] - bin_array_I[0])))
        # when align='center', the central point of the bar in the x-axis is simply the middle of the bins ((bin_0-bin_1)/2, etc)
        centre_of_bar_in_x_axis = (bin_array_I[:-2] + bin_array_I[1:-1]) / 2
        # add the final bin, which is physically located just after the last regular bin but represents all higher values
        bar_width = centre_of_bar_in_x_axis[3] - centre_of_bar_in_x_axis[2]
        centre_of_bar_in_x_axis = np.append(centre_of_bar_in_x_axis, centre_of_bar_in_x_axis[-1] + bar_width)
        barcontainer_AAIMON_mean = ax.bar(left=centre_of_bar_in_x_axis, height=freq_counts_I,
                                          align='center', width=col_width, color="#EE762C",
                                          alpha=0.5)
        # create numpy array of membranous over nonmembranous conservation ratios (identity + similarity)
        hist_data_AASMON_mean = np.array(df['AASMON_ratio_mean_all_TMDs'].dropna())
        # use numpy to create a histogram
        freq_counts_S, bin_array_S = np.histogram(hist_data_AASMON_mean, bins=binlist)
        # barcontainer_S = axarr[row_nr,col_nr].bar(left=centre_of_bar_in_x_axis, height=freq_counts_S, align='center', width=col_width, color="#0101DF", edgecolor="#0101DF", alpha = 0.5)
        # create a line graph rather than a bar graph for the AASMON (ident + similarity)
        linecontainer_AASMON_mean = ax.plot(centre_of_bar_in_x_axis, freq_counts_S, color="#0101DF",
                                            alpha=0.5)
        # other colours that are compatible with colourblind readers: #8A084B Dark red, #B45F04 deep orange, reddish purple #4B088A
        # http://html-color-codes.info/
        # label the x-axis for each plot, based on the TMD
        ax.set_xlabel('average conservation ratio (membranous over nonmembranous)', fontsize=fontsize)
        # move the x-axis label closer to the x-axis
        ax.xaxis.set_label_coords(0.45, -0.085)
        xlim_min = s["mp_xlim_min01"]
        # take x-axis max from settings
        xlim_max = s["mp_xlim_max01"]
        # set x-axis min
        ax.set_xlim(xlim_min, xlim_max)
        # set x-axis ticks
        # use the slide selection to select every second item in the list as an xtick(axis label)
        ax.set_xticks([float('%0.1f' % c) for c in centre_of_bar_in_x_axis[::3]])
        ax.set_ylabel('freq', rotation='vertical', fontsize=fontsize)
        # change axis font size
        ax.tick_params(labelsize=fontsize)
        # create legend?#http://stackoverflow.com/questions/9834452/how-do-i-make-a-single-legend-for-many-subplots-with-matplotlib
        legend_obj = ax.legend(['AASMON (identity + similarity)', 'AAIMON (identity)', 'AAIMON norm (identity)'], loc='upper right',
                               fontsize=fontsize)
        # add figure number to top left of subplot
        ax.annotate(s=str(Fig_Nr) + '.', xy=(0.04, 0.9), fontsize=fontsize, xytext=None, xycoords='axes fraction',
                    alpha=0.75)
        # add figure title to top left of subplot
        ax.annotate(s=title, xy=(0.1, 0.9), fontsize=fontsize, xytext=None, xycoords='axes fraction', alpha=0.75)
        # # FROM RIMMA SCRIPT - not necessary
        # ax.yaxis.grid(True, zorder=0, linestyle=":", color="grey")
        # ax
        # for tic in ax.xaxis.get_major_ticks():
        #     tic.tick1On = False
        # for tic in ax.yaxis.get_major_ticks():
        #     tic.tick1On = False
        # ax.spines['top'].set_visible(False)
        # ax.spines['right'].set_visible(False)
        # # END FROM RIMMA SCRIPT

        utils.save_figure(fig, Fig_name, base_filepath, save_png, save_pdf, dpi)

    if s['Fig02_Histogram_of_standard_deviations_for_AAIMON_and_AASMON']:
        Fig_Nr = 2
        title = 'Standard Deviaton, SP vs MP'
        Fig_name = 'Fig02_Histogram_of_standard_deviations_for_AAIMON_and_AASMON'
        # create a new figure
        fig, ax = plt.subplots()
        
        # create numpy array of membranous over nonmembranous conservation ratios (identity)
        hist_data_AAIMON_std = np.array(df['AAIMON_std_all_TMDs'].dropna())
        # use numpy to create a histogram
        number_of_bins = 50
        freq_counts_S, bin_array_S = np.histogram(hist_data_AAIMON_std, bins=number_of_bins)
        # assuming all of the bins are exactly the same size, make the width of the column equal to 70% of each bin
        col_width = float('%0.3f' % (0.95 * (bin_array_S[1] - bin_array_S[0])))
        # when align='center', the central point of the bar in the x-axis is simply the middle of the bins ((bin_0-bin_1)/2, etc)
        centre_of_bar_in_x_axis = (bin_array_S[:-1] + bin_array_S[1:]) / 2
        barcontainer_AAIMON_std = ax.bar(left=centre_of_bar_in_x_axis, height=freq_counts_S,
                                                            align='center', width=col_width, color="#0489B1",
                                                            alpha=0.5)  # edgecolor='black',
        # create numpy array of membranous over nonmembranous conservation ratios (identity + similarity)
        hist_data_AASMON_std = np.array(df['AASMON_ratio_std_all_TMDs'].dropna())
        # use numpy to create a histogram
        # N.B. use the bins from the previous plot
        freq_counts, bin_array = np.histogram(hist_data_AASMON_std, bins=bin_array_S)
        # barcontainer_S = axarr[row_nr,col_nr].bar(left=centre_of_bar_in_x_axis, height=freq_counts, align='center', width=col_width, color="#0101DF", edgecolor="#0101DF", alpha = 0.5)
        # create a line graph rather than a bar graph for the AASMON (ident + similarity)
        linecontainer_AASMON_std = ax.plot(centre_of_bar_in_x_axis, freq_counts, color="#0101DF",
                                                              alpha=0.5)
        # other colours that are compatible with colourblind readers: #8A084B Dark red, #B45F04 deep orange, reddish purple #4B088A
        # http://html-color-codes.info/
        # label the x-axis for each plot, based on the TMD
        ax.set_xlabel('average standard deviation', fontsize=fontsize)
        # move the x-axis label closer to the x-axis
        ax.xaxis.set_label_coords(0.45, -0.085)
        # xlim_min = s["mp_xlim_min01"]
        # take x-axis max from settings
        # xlim_max = s["mp_xlim_max01"]
        # set x-axis min
        # axarr[row_nr,col_nr].set_xlim(xlim_min,xlim_max)
        # set x-axis ticks
        # use the slide selection to select every second item in the list as an xtick(axis label)
        ax.set_xticks([float('%0.1f' % c) for c in centre_of_bar_in_x_axis[::3]])
        ax.set_ylabel('freq', rotation='vertical', fontsize=fontsize)
        # change axis font size
        ax.tick_params(labelsize=fontsize)
        # create legend?#http://stackoverflow.com/questions/9834452/how-do-i-make-a-single-legend-for-many-subplots-with-matplotlib
        legend_obj = ax.legend(['AASMON (identity + similarity)', 'AAIMON (identity)'],
                                                  loc='upper right',
                                                  fontsize=fontsize)
        # add figure number to top left of subplot
        ax.annotate(s=str(Fig_Nr) + '.', xy=(0.04, 0.9), fontsize=fontsize, xytext=None,
                                       xycoords='axes fraction',
                                       alpha=0.75)
        # add figure title to top left of subplot
        ax.annotate(s=title, xy=(0.1, 0.9), fontsize=fontsize, xytext=None, xycoords='axes fraction',
                                       alpha=0.75)

        utils.save_figure(fig, Fig_name, base_filepath, save_png, save_pdf, dpi)

    if s['Fig03_Scattergram_comparing_mean_AAIMON_and_AASMON']:
        Fig_Nr = 3
        title = 'AAIMON vs AASMON'
        Fig_name = 'Fig03_Scattergram_comparing_mean_AAIMON_and_AASMON'
        fig, ax = plt.subplots()
        # pylab.rcParams['figure.figsize'] = (50.0, 40.0)
        x = np.array(df['AAIMON_mean_all_TMDs'])
        y = np.array(df['AASMON_ratio_mean_all_TMDs'])
        scattercontainer_AAIMON_AASMON_mean = ax.scatter(x=x, y=y, color="#8A084B", alpha=alpha_dpd,
                                                                            s=datapointsize)
        # label the x-axis for each plot, based on the TMD
        ax.set_xlabel('AAIMON (aa identity)', fontsize=fontsize)
        # move the x-axis label closer to the x-axis
        ax.xaxis.set_label_coords(0.45, -0.085)
        ax.set_ylabel('AASMON_ratio (aa identity + similarity)', fontsize=fontsize)
        # change axis font size
        ax.tick_params(labelsize=fontsize)
        # create legend?#http://stackoverflow.com/questions/9834452/how-do-i-make-a-single-legend-for-many-subplots-with-matplotlib
        legend_obj = ax.legend(['mean'], loc='upper right', fontsize=fontsize)
        # add figure number to top left of subplot
        ax.annotate(s=str(Fig_Nr) + '.', xy=(0.04, 0.9), fontsize=fontsize, xytext=None,
                                       xycoords='axes fraction',
                                       alpha=0.75)
        # add figure title to top left of subplot
        ax.annotate(s=title, xy=(0.1, 0.9), fontsize=fontsize, xytext=None, xycoords='axes fraction',
                                       alpha=0.75)

        utils.save_figure(fig, Fig_name, base_filepath, save_png, save_pdf, dpi)

    if s['Fig04_Scattergram_comparing_standard_deviation_AAIMON_and_AASMON']:
        Fig_Nr = 4
        title = 'standard deviation AAIMON vs AASMON'
        Fig_name = 'Fig04_Scattergram_comparing_standard_deviation_AAIMON_and_AASMON'
        fig, ax = plt.subplots()

        x = np.array(df['AAIMON_std_all_TMDs'])
        y = np.array(df['AASMON_ratio_std_all_TMDs'])
        scattercontainer_AAIMON_AASMON_std = ax.scatter(x=x, y=y, color="#B45F04", alpha=alpha_dpd,
                                                                           s=datapointsize)
        # other colours that are compatible with colourblind readers: #8A084B Dark red, #B45F04 deep orange, reddish purple #4B088A
        # http://html-color-codes.info/
        # label the x-axis for each plot, based on the TMD
        ax.set_xlabel('AAIMON', fontsize=fontsize)
        # move the x-axis label closer to the x-axis
        ax.xaxis.set_label_coords(0.45, -0.085)
        ax.set_ylabel('AASMON_ratio', fontsize=fontsize)
        # set x-axis ticks
        # use the slide selection to select every second item in the list as an xtick(axis label)
        # axarr[row_nr,col_nr].set_xticks([float('%0.1f' % c) for c in centre_of_bar_in_x_axis[::3]])
        # axarr[row_nr,col_nr].set_ylabel('freq',rotation = 'vertical', fontsize = fontsize)
        # change axis font size
        ax.tick_params(labelsize=fontsize)
        # create legend?#http://stackoverflow.com/questions/9834452/how-do-i-make-a-single-legend-for-many-subplots-with-matplotlib
        legend_obj = ax.legend(['standard deviation'], loc='upper right', fontsize=fontsize)
        # add figure number to top left of subplot
        ax.annotate(s=str(Fig_Nr) + '.', xy=(0.04, 0.9), fontsize=fontsize, xytext=None,
                                       xycoords='axes fraction',
                                       alpha=0.75)
        # add figure title to top left of subplot
        ax.annotate(s=title, xy=(0.1, 0.9), fontsize=fontsize, xytext=None, xycoords='axes fraction',
                                       alpha=0.75)

        utils.save_figure(fig, Fig_name, base_filepath, save_png, save_pdf, dpi)

    if s['Fig05_Scattergram_comparing_number_of_TMDs_with_mean_AAIMON']:
        Fig_Nr = 5
        title = 'num_TMDs vs AAIMON'
        Fig_name = 'Fig05_Scattergram_comparing_number_of_TMDs_with_mean_AAIMON'
        fig, ax = plt.subplots()
        

        # for backwards compatibility, check for old name.
        # commented out by MO - assumed that it is not needed
        # if 'number_of_TMDs' in df.columns:
        #     x = np.array(df['number_of_TMDs'])
        # else:
        #     x = np.array(df['number_of_TMDs'])
        x = np.array(df['number_of_TMDs'])
        y = np.array(df['AAIMON_mean_all_TMDs'])
        scattercontainer_AAIMON_AASMON_std = ax.scatter(x=x, y=y, color="#0489B1", alpha=alpha_dpd,
                                                                           s=datapointsize)
        # other colours that are compatible with colourblind readers: #8A084B Dark red, #B45F04 deep orange, reddish purple #4B088A
        # http://html-color-codes.info/
        # label the x-axis for each plot, based on the TMD
        ax.set_xlabel('number of TMDs in protein', fontsize=fontsize)
        # move the x-axis label closer to the x-axis
        ax.xaxis.set_label_coords(0.45, -0.085)
        ax.set_ylabel('Average AAIMON ratio for all TMDs', fontsize=fontsize)
        # change axis font size
        ax.tick_params(labelsize=fontsize)
        # create legend?#http://stackoverflow.com/questions/9834452/how-do-i-make-a-single-legend-for-many-subplots-with-matplotlib
        legend_obj = ax.legend(['AAIMON_mean'], loc='upper right', fontsize=fontsize)
        # add background grid
        ax.grid(True, color='0.75', alpha=0.3)
        # add figure number to top left of subplot
        ax.annotate(s=str(Fig_Nr) + '.', xy=(0.04, 0.9), fontsize=fontsize, xytext=None,
                                       xycoords='axes fraction',
                                       alpha=0.75)
        # add figure title to top left of subplot
        ax.annotate(s=title, xy=(0.1, 0.9), fontsize=fontsize, xytext=None, xycoords='axes fraction',
                                       alpha=0.75)

        utils.save_figure(fig, Fig_name, base_filepath, save_png, save_pdf, dpi)

    if s['Fig06_Scattergram_comparing_seqlen_with_mean_AAIMON']:
        Fig_Nr = 6
        title = 'seqlen vs AAIMON'
        Fig_name = 'Fig06_Scattergram_comparing_seqlen_with_mean_AAIMON'
        fig, ax = plt.subplots()
        
        # pylab.rcParams['figure.figsize'] = (100.0, 80.0)
        x = np.array(df['seqlen'])
        y = np.array(df['AAIMON_mean_all_TMDs'])
        scattercontainer_AAIMON_AASMON_std = ax.scatter(x=x, y=y, color="#0489B1", alpha=alpha_dpd,
                                                                           s=datapointsize)
        # label the x-axis for each plot, based on the TMD
        ax.set_xlabel('Length of protein', fontsize=fontsize)
        # move the x-axis label closer to the x-axis
        ax.xaxis.set_label_coords(0.45, -0.085)
        ax.set_ylabel('AAIMON', fontsize=fontsize)
        # change axis font size
        ax.tick_params(labelsize=fontsize)
        # add figure number to top left of subplot
        ax.annotate(s=str(Fig_Nr) + '.', xy=(0.04, 0.9), fontsize=fontsize, xytext=None,
                                       xycoords='axes fraction',
                                       alpha=0.75)
        # add figure title to top left of subplot
        ax.annotate(s=title, xy=(0.1, 0.9), fontsize=fontsize, xytext=None, xycoords='axes fraction',
                                       alpha=0.75)
        utils.save_figure(fig, Fig_name, base_filepath, save_png, save_pdf, dpi)

    if s['Fig07_Scattergram_comparing_nonTMD_SW_align_len_mean_with_mean_AAIMON']:
        Fig_Nr = 7
        title = 'length nonTMD region'
        Fig_name = 'Fig07_Scattergram_comparing_nonTMD_SW_align_len_mean_with_mean_AAIMON'
        fig, ax = plt.subplots()
        
        # pylab.rcParams['figure.figsize'] = (100.0, 80.0)
        x = np.array(df['nonTMD_SW_align_len_mean'])
        y = np.array(df['AAIMON_mean_all_TMDs'])
        scattercontainer_AAIMON_AASMON_std = ax.scatter(x=x, y=y, color="#0489B1", alpha=alpha_dpd,
                                                                           s=datapointsize)
        # label the x-axis for each plot, based on the TMD
        ax.set_xlabel('Average length of nonTMD region in homologues', fontsize=fontsize)
        # move the x-axis label closer to the x-axis
        ax.xaxis.set_label_coords(0.45, -0.085)
        ax.set_ylabel('AAIMON', fontsize=fontsize)
        # change axis font size
        ax.tick_params(labelsize=fontsize)
        # add figure number to top left of subplot
        ax.annotate(s=str(Fig_Nr) + '.', xy=(0.04, 0.9), fontsize=fontsize, xytext=None,
                                       xycoords='axes fraction',
                                       alpha=0.75)
        # add figure title to top left of subplot
        ax.annotate(s=title, xy=(0.1, 0.9), fontsize=fontsize, xytext=None, xycoords='axes fraction',
                                       alpha=0.75)

        utils.save_figure(fig, Fig_name, base_filepath, save_png, save_pdf, dpi)

    if s['Fig08_Scattergram_comparing_total_number_of_simap_hits_with_mean_AAIMON']:
        '''
        note that the total hits comes from SIMAP, so this doesn't say anything about how much data is available for each protein
        '''
        Fig_Nr = 8
        title = 'number SIMAP hits'
        Fig_name = 'Fig08_Scattergram_comparing_total_number_of_simap_hits_with_mean_AAIMON'
        fig, ax = plt.subplots()

        # plot AAIMON
        x = np.array(df['TM01_AAIMON_n_homol']) # total_number_of_simap_hits can be replaced with TM01_AAIMON_n_homol
        y = np.array(df['AAIMON_mean_all_TMDs'])
        scattercontainer_AAIMON_AASMON_std = ax.scatter(x=x, y=y, color="#0489B1", alpha=alpha_dpd,
                                                                           s=datapointsize)

        # plot AAIMON normalised
        x = np.array(df['TM01_AAIMON_n_homol']) # total_number_of_simap_hits can be replaced with TM01_AAIMON_n_homol
        y = np.array(df['AAIMON_mean_all_TMDs_n'])
        scattercontainer_AAIMON_AASMON_std = ax.scatter(x=x, y=y, color="#EE762C", alpha=alpha_dpd,
                                                                           s=datapointsize)
        # label the x-axis for each plot, based on the TMD
        ax.set_xlabel('total number of homologues', fontsize=fontsize)
        # move the x-axis label closer to the x-axis
        ax.xaxis.set_label_coords(0.45, -0.085)
        ax.set_ylabel('AAIMON', fontsize=fontsize)

        # create legend
        legend_obj = ax.legend(['AAIMON', 'AAIMON norm'],
                               loc='upper right',
                               fontsize=fontsize)
        # change axis font size
        ax.tick_params(labelsize=fontsize)
        # add figure number to top left of subplot
        ax.annotate(s=str(Fig_Nr) + '.', xy=(0.04, 0.9), fontsize=fontsize, xytext=None,
                                       xycoords='axes fraction',
                                       alpha=0.75)
        # add figure title to top left of subplot
        ax.annotate(s=title, xy=(0.1, 0.9), fontsize=fontsize, xytext=None, xycoords='axes fraction',
                                       alpha=0.75)
        utils.save_figure(fig, Fig_name, base_filepath, save_png, save_pdf, dpi)

    if s['Fig09_Histogram_of_mean_AAIMON_ratios_for_each_TMD_separately']:
        Fig_Nr = 9
        title = 'Histogram of mean AAIMON ratios'
        Fig_name = 'Fig09_Histogram_of_mean_AAIMON_ratios_for_each_TMD_separately'
        fig, ax = plt.subplots()

        title = 'AAIMON each TMD separately'
        num_bins = 30
        # "#0489B1"
        alpha = 0.25
        col_width_value = 0.95
        ylabel = 'freq'
        xlabel = 'average conservation ratio (membranous over nonmembranous)'

        for n, TM in enumerate(df_mean_AAIMON_each_TM.columns):
            # define the colour for that TMD
            # if there are more TMDs than colours, simply start from the beginning of the list again
            if n < len(tableau20):
                color_num = n
            else:
                color_num = n - len(tableau20)
            color = tableau20[color_num]

            hist_data = np.array(df_mean_AAIMON_each_TM[TM].dropna())
            '''
            Calculated the bins for a histogram, even for highly non-normal data
            '''
            # calculate 5th percentile
            percentile_5 = np.percentile(hist_data, 5)  # hist_data.min()
            # calculate 9th percentile
            percentile_95 = np.percentile(hist_data, 95)  # hist_data.max()
            # calculate difference
            percentile_95_minus_5 = percentile_95 - percentile_5
            # create buffer for bins
            extra_xaxis_range = percentile_95_minus_5 / 4
            # lowest bin is the 5th percentile minus the buffer, except where that is below zero
            data_min = percentile_5 - extra_xaxis_range  # hist_data.min()
            # ata_min = 0 if data_max < 0 else data_max
            # highest bin is the 95th percentile
            data_max = percentile_95 + extra_xaxis_range  # hist_data.max()
            # create bins using the min and max
            binlist = np.linspace(data_min, data_max, num_bins)
            # use numpy to create a histogram
            freq_counts_I, bin_array_I = np.histogram(hist_data, bins=binlist)
            # assuming all of the bins are exactly the same size, make the width of the column equal to 70% of each bin
            col_width = float('%0.3f' % (col_width_value * (bin_array_I[1] - bin_array_I[0])))
            # when align='center', the central point of the bar in the x-axis is simply the middle of the bins ((bin_0-bin_1)/2, etc)
            centre_of_bar_in_x_axis = (bin_array_I[:-2] + bin_array_I[1:-1]) / 2
            # add the final bin, which is physically located just after the last regular bin but represents all higher values
            bar_width = centre_of_bar_in_x_axis[3] - centre_of_bar_in_x_axis[2]
            centre_of_bar_in_x_axis = np.append(centre_of_bar_in_x_axis, centre_of_bar_in_x_axis[-1] + bar_width)
            barcontainer = ax.bar(left=centre_of_bar_in_x_axis, height=freq_counts_I,
                                                     align='center', width=col_width, facecolor=color,
                                                     alpha=alpha, edgecolor='black',
                                                     linewidth=0.2)  # edgecolor='black',
            # barcontainer = ax.bar(left=centre_of_bar_in_x_axis, height=freq_counts_I,
            #                                                     align='center', width=col_width, facecolor=color,
            #                                                     alpha=alpha, edgecolor='black', linewidth=0.1)  # edgecolor='black',    #label the x-axis for each plot, based on the TMD
        ax.set_xlabel(xlabel, fontsize=fontsize)
        # move the x-axis label closer to the x-axis
        ax.xaxis.set_label_coords(0.45, -0.085)
        ax.set_ylabel(ylabel, rotation='vertical', fontsize=fontsize)
        # change axis font size
        ax.tick_params(labelsize=fontsize)
        # create legend?#http://stackoverflow.com/questions/9834452/how-do-i-make-a-single-legend-for-many-subplots-with-matplotlib ; commented out MO
        # legend_obj = ax.legend(legend, loc='upper right', fontsize=fontsize)
        # add title
        # ax.set_title(title,fontsize=fontsize)
        # add figure number to top left of subplot
        ax.annotate(s=str(Fig_Nr) + '.', xy=(0.04, 0.9), fontsize=fontsize, xytext=None,
                                       xycoords='axes fraction',
                                       alpha=0.75)
        # add figure title to top left of subplot
        ax.annotate(s=title, xy=(0.1, 0.9), fontsize=fontsize, xytext=None, xycoords='axes fraction',
                                       alpha=0.75)

        utils.save_figure(fig, Fig_name, base_filepath, save_png, save_pdf, dpi)

    if s['Fig10_Line_histogram_of_mean_AAIMON_ratios_for_each_TMD_separately']:
        Fig_Nr = 10
        title = 'Line histogram each TMD'
        Fig_name = 'Fig10_Line_histogram_of_mean_AAIMON_ratios_for_each_TMD_separately'
        fig, ax = plt.subplots()

        num_bins = 50
        # "#0489B1"
        alpha = 0.9
        col_width_value = 0.95
        ylabel = 'freq'
        xlabel = 'average conservation ratio (membranous over nonmembranous)'

        for n, TM in enumerate(df_mean_AAIMON_each_TM.columns):
            # define the colour for that TMD
            # if there are more TMDs than colours, simply start from the beginning of the list again
            if n < len(tableau20):
                color_num = n
            else:
                color_num = n - len(tableau20)
            color = tableau20[color_num]

            hist_data = np.array(df_mean_AAIMON_each_TM[TM].dropna())
            '''
            Calculated the bins for a histogram, even for highly non-normal data
            '''
            # calculate 5th percentile
            percentile_5 = np.percentile(hist_data, 5)  # hist_data.min()
            # calculate 9th percentile
            percentile_95 = np.percentile(hist_data, 95)  # hist_data.max()
            # calculate difference
            percentile_95_minus_5 = percentile_95 - percentile_5
            # create buffer for bins
            extra_xaxis_range = percentile_95_minus_5 / 4
            # lowest bin is the 5th percentile minus the buffer, except where that is below zero
            data_min = percentile_5 - extra_xaxis_range  # hist_data.min()
            # ata_min = 0 if data_max < 0 else data_max
            # highest bin is the 95th percentile
            data_max = percentile_95 + extra_xaxis_range  # hist_data.max()
            # create bins using the min and max
            binlist = np.linspace(data_min, data_max, num_bins)
            # use numpy to create a histogram
            freq_counts_I, bin_array_I = np.histogram(hist_data, bins=binlist)
            # assuming all of the bins are exactly the same size, make the width of the column equal to 70% of each bin
            col_width = float('%0.3f' % (col_width_value * (bin_array_I[1] - bin_array_I[0])))
            # when align='center', the central point of the bar in the x-axis is simply the middle of the bins ((bin_0-bin_1)/2, etc)
            centre_of_bar_in_x_axis = (bin_array_I[:-2] + bin_array_I[1:-1]) / 2
            # add the final bin, which is physically located just after the last regular bin but represents all higher values
            bar_width = centre_of_bar_in_x_axis[3] - centre_of_bar_in_x_axis[2]
            centre_of_bar_in_x_axis = np.append(centre_of_bar_in_x_axis, centre_of_bar_in_x_axis[-1] + bar_width)
            linecontainer = ax.plot(centre_of_bar_in_x_axis, freq_counts_I, alpha=alpha,
                                                       linewidth=0.3)  # edgecolor='black',
        # label the x-axis for each plot, based on the TMD
        ax.set_xlabel(xlabel, fontsize=fontsize)
        # move the x-axis label closer to the x-axis
        ax.xaxis.set_label_coords(0.45, -0.085)
        ax.set_ylabel(ylabel, rotation='vertical', fontsize=fontsize)
        # change axis font size
        ax.tick_params(labelsize=fontsize)
        # create legend?#http://stackoverflow.com/questions/9834452/how-do-i-make-a-single-legend-for-many-subplots-with-matplotlib; commented out MO
        # legend_obj = ax.legend(legend, loc='upper right', fontsize=fontsize)
        # add title
        # ax.set_title(title,fontsize=fontsize)
        # add figure number to top left of subplot
        ax.annotate(s=str(Fig_Nr) + '.', xy=(0.04, 0.9), fontsize=fontsize, xytext=None,
                                       xycoords='axes fraction',
                                       alpha=0.75)
        # add figure title to top left of subplot
        ax.annotate(s=title, xy=(0.1, 0.9), fontsize=fontsize, xytext=None, xycoords='axes fraction',
                                       alpha=0.75)

        utils.save_figure(fig, Fig_name, base_filepath, save_png, save_pdf, dpi)

    if s['Fig11_Line_histogram_of_mean_AAIMON_ratios_for_selected_TMDs,_highlighting_difference_for_TM07']:
        # these graphs are only applicable for multi-pass proteins. Use where at least 2 proteins have a 7th TMD
        if "TM07_AAIMON_mean" in df_mean_AAIMON_each_TM.columns:
            if df_mean_AAIMON_each_TM['TM07_AAIMON_mean'].dropna().shape[0] >= 2:
                dataset_contains_multipass_prots = True
            else:
                dataset_contains_multipass_prots = False
                sys.stdout.write('Dataset does not contain proteins with 7 TMDs; cannot create figure 11 \n')
        else:
            dataset_contains_multipass_prots = False
            sys.stdout.write('Dataset does not contain multipass proteins; figure 11 cannot be created \n')

        if dataset_contains_multipass_prots:
            Fig_Nr = 11
            title = 'Select TMDs, all data'
            Fig_name = 'Fig11_Line_histogram_of_mean_AAIMON_ratios_for_selected_TMDs'
            #cols_for_analysis = ['TM01', 'TM07', 'TM08', 'last_TM_AAIMON_mean']
            cols_for_analysis = ['TM01_AAIMON_mean', 'TM07_AAIMON_mean', 'TM08_AAIMON_mean', 'TM{last_TM:02d}_AAIMON_mean'.format(last_TM=len(df_mean_AAIMON_each_TM.columns))]
            fig, ax = plt.subplots()
            num_bins = 50
            # "#0489B1"
            alpha = 0.7
            col_width_value = 0.95
            ylabel = 'freq'
            xlabel = 'average conservation ratio (membranous over nonmembranous)'

            for n, TM in enumerate(cols_for_analysis):
                # define the colour for that TMD
                # if there are more TMDs than colours, simply start from the beginning of the list again
                if n < len(tableau20):
                    color_num = n
                else:
                    color_num = n - len(tableau20)
                color = tableau20[color_num]

                hist_data = np.array(df_mean_AAIMON_each_TM[TM].dropna())
                '''
                Calculated the bins for a histogram, even for highly non-normal data
                '''
                # calculate 5th percentile
                percentile_5 = np.percentile(hist_data, 5)  # hist_data.min()
                # calculate 9th percentile
                percentile_95 = np.percentile(hist_data, 95)  # hist_data.max()
                # calculate difference
                percentile_95_minus_5 = percentile_95 - percentile_5
                # create buffer for bins
                extra_xaxis_range = percentile_95_minus_5 / 4
                # lowest bin is the 5th percentile minus the buffer, except where that is below zero
                data_min = percentile_5 - extra_xaxis_range  # hist_data.min()
                # ata_min = 0 if data_max < 0 else data_max
                # highest bin is the 95th percentile
                data_max = percentile_95 + extra_xaxis_range  # hist_data.max()
                # create bins using the min and max
                binlist = np.linspace(data_min, data_max, num_bins)
                # use numpy to create a histogram
                freq_counts_I, bin_array_I = np.histogram(hist_data, bins=binlist)
                # assuming all of the bins are exactly the same size, make the width of the column equal to 70% of each bin
                col_width = float('%0.3f' % (col_width_value * (bin_array_I[1] - bin_array_I[0])))
                # when align='center', the central point of the bar in the x-axis is simply the middle of the bins ((bin_0-bin_1)/2, etc)
                centre_of_bar_in_x_axis = (bin_array_I[:-2] + bin_array_I[1:-1]) / 2
                # add the final bin, which is physically located just after the last regular bin but represents all higher values
                bar_width = centre_of_bar_in_x_axis[3] - centre_of_bar_in_x_axis[2]
                centre_of_bar_in_x_axis = np.append(centre_of_bar_in_x_axis, centre_of_bar_in_x_axis[-1] + bar_width)
                linecontainer = ax.plot(centre_of_bar_in_x_axis, freq_counts_I, alpha=alpha,
                                                           linewidth=0.3)  # edgecolor='black',
            # label the x-axis for each plot, based on the TMD
            ax.set_xlabel(xlabel, fontsize=fontsize)
            # move the x-axis label closer to the x-axis
            ax.xaxis.set_label_coords(0.45, -0.085)
            ax.set_ylabel(ylabel, rotation='vertical', fontsize=fontsize)
            # change axis font size
            ax.tick_params(labelsize=fontsize)
            # create legend?#http://stackoverflow.com/questions/9834452/how-do-i-make-a-single-legend-for-many-subplots-with-matplotlib
            legend_obj = ax.legend(cols_for_analysis, loc='upper right', fontsize=fontsize)
            # add title
            # ax.set_title(title,fontsize=fontsize)
            # add figure number to top left of subplot
            ax.annotate(s=str(Fig_Nr) + '.', xy=(0.04, 0.9), fontsize=fontsize, xytext=None,
                                           xycoords='axes fraction', alpha=0.75)
            # add figure title to top left of subplot
            ax.annotate(s=title, xy=(0.1, 0.9), fontsize=fontsize, xytext=None, xycoords='axes fraction',
                                           alpha=0.75)
            utils.save_figure(fig, Fig_name, base_filepath, save_png, save_pdf, dpi)

    if s['Fig12_TMD_1-5_only']:
        Fig_Nr = 12
        col_start = 0
        col_end = 5
        cols_to_analyse = df_mean_AAIMON_each_TM.columns[col_start:col_end]
        title = 'TMD 1 to 5, all data'
        Fig_name = 'Fig12_TMD_1-5_only'
        fig, ax = plt.subplots()
        num_bins = 30
        # "#0489B1"
        alpha = 0.9
        col_width_value = 0.95
        ylabel = 'freq'
        xlabel = 'average conservation ratio (membranous over nonmembranous)'

        for n, TM in enumerate(cols_to_analyse):
            # define the colour for that TMD
            # if there are more TMDs than colours, simply start from the beginning of the list again
            if n < len(tableau20):
                color_num = n
            else:
                color_num = n - len(tableau20)
            color = tableau20[color_num]

            hist_data = np.array(df_mean_AAIMON_each_TM[TM].dropna())
            '''
            Calculated the bins for a histogram, even for highly non-normal data
            '''
            # calculate 5th percentile
            percentile_5 = np.percentile(hist_data, 5)  # hist_data.min()
            # calculate 9th percentile
            percentile_95 = np.percentile(hist_data, 95)  # hist_data.max()
            # calculate difference
            percentile_95_minus_5 = percentile_95 - percentile_5
            # create buffer for bins
            extra_xaxis_range = percentile_95_minus_5 / 4
            # lowest bin is the 5th percentile minus the buffer, except where that is below zero
            data_min = percentile_5 - extra_xaxis_range  # hist_data.min()
            # ata_min = 0 if data_max < 0 else data_max
            # highest bin is the 95th percentile
            data_max = percentile_95 + extra_xaxis_range  # hist_data.max()
            # create bins using the min and max
            binlist = np.linspace(data_min, data_max, num_bins)
            # use numpy to create a histogram
            freq_counts_I, bin_array_I = np.histogram(hist_data, bins=binlist)
            # assuming all of the bins are exactly the same size, make the width of the column equal to 70% of each bin
            col_width = float('%0.3f' % (col_width_value * (bin_array_I[1] - bin_array_I[0])))
            # when align='center', the central point of the bar in the x-axis is simply the middle of the bins ((bin_0-bin_1)/2, etc)
            centre_of_bar_in_x_axis = (bin_array_I[:-2] + bin_array_I[1:-1]) / 2
            # add the final bin, which is physically located just after the last regular bin but represents all higher values
            bar_width = centre_of_bar_in_x_axis[3] - centre_of_bar_in_x_axis[2]
            centre_of_bar_in_x_axis = np.append(centre_of_bar_in_x_axis, centre_of_bar_in_x_axis[-1] + bar_width)
            linecontainer = ax.plot(centre_of_bar_in_x_axis, freq_counts_I, alpha=alpha,
                                                       linewidth=0.3)  # edgecolor='black',
        # label the x-axis for each plot, based on the TMD
        ax.set_xlabel(xlabel, fontsize=fontsize)
        # move the x-axis label closer to the x-axis
        ax.xaxis.set_label_coords(0.45, -0.085)
        ax.set_ylabel(ylabel, rotation='vertical', fontsize=fontsize)
        # change axis font size
        ax.tick_params(labelsize=fontsize)
        # create legend?#http://stackoverflow.com/questions/9834452/how-do-i-make-a-single-legend-for-many-subplots-with-matplotlib
        legend_obj = ax.legend(cols_to_analyse, loc='upper right', fontsize=fontsize)
        # add title
        # ax.set_title(title,fontsize=fontsize)
        # add figure number to top left of subplot
        ax.annotate(s=str(Fig_Nr) + '.', xy=(0.04, 0.9), fontsize=fontsize, xytext=None,
                                       xycoords='axes fraction', alpha=0.75)
        # add figure title to top left of subplot
        ax.annotate(s=title, xy=(0.1, 0.9), fontsize=fontsize, xytext=None, xycoords='axes fraction',
                                       alpha=0.75)
        utils.save_figure(fig, Fig_name, base_filepath, save_png, save_pdf, dpi)

    if s['Fig13_TMD_5-10_only']:
        Fig_Nr = 13
        title = 'TMD 5-10, all data'
        Fig_name = 'Fig13_TMD_5-10_only'
        col_start = 5
        col_end = 10
        # analyzing TM01 (as control?) and TM05-10
        cols_to_analyse = ['TM01_AAIMON_mean'] + list(df_mean_AAIMON_each_TM.columns[col_start:col_end])
        
        fig, ax = plt.subplots()

        num_bins = 30
        # "#0489B1"
        alpha = 0.9
        col_width_value = 0.95
        ylabel = 'freq'
        xlabel = 'average conservation ratio (membranous over nonmembranous)'

        for n, TM in enumerate(cols_to_analyse):
            # define the colour for that TMD
            # if there are more TMDs than colours, simply start from the beginning of the list again
            if n < len(tableau20):
                color_num = n
            else:
                color_num = n - len(tableau20)
            color = tableau20[color_num]

            hist_data = np.array(df_mean_AAIMON_each_TM[TM].dropna())
            '''
            Calculated the bins for a histogram, even for highly non-normal data
            '''
            # calculate 5th percentile
            percentile_5 = np.percentile(hist_data, 5)  # hist_data.min()
            # calculate 9th percentile
            percentile_95 = np.percentile(hist_data, 95)  # hist_data.max()
            # calculate difference
            percentile_95_minus_5 = percentile_95 - percentile_5
            # create buffer for bins
            extra_xaxis_range = percentile_95_minus_5 / 4
            # lowest bin is the 5th percentile minus the buffer, except where that is below zero
            data_min = percentile_5 - extra_xaxis_range  # hist_data.min()
            # ata_min = 0 if data_max < 0 else data_max
            # highest bin is the 95th percentile
            data_max = percentile_95 + extra_xaxis_range  # hist_data.max()
            # create bins using the min and max
            binlist = np.linspace(data_min, data_max, num_bins)
            # use numpy to create a histogram
            freq_counts_I, bin_array_I = np.histogram(hist_data, bins=binlist)
            # assuming all of the bins are exactly the same size, make the width of the column equal to 70% of each bin
            col_width = float('%0.3f' % (col_width_value * (bin_array_I[1] - bin_array_I[0])))
            # when align='center', the central point of the bar in the x-axis is simply the middle of the bins ((bin_0-bin_1)/2, etc)
            centre_of_bar_in_x_axis = (bin_array_I[:-2] + bin_array_I[1:-1]) / 2
            # add the final bin, which is physically located just after the last regular bin but represents all higher values
            bar_width = centre_of_bar_in_x_axis[3] - centre_of_bar_in_x_axis[2]
            centre_of_bar_in_x_axis = np.append(centre_of_bar_in_x_axis, centre_of_bar_in_x_axis[-1] + bar_width)
            linecontainer = ax.plot(centre_of_bar_in_x_axis, freq_counts_I, alpha=alpha,
                                                       linewidth=0.3)  # edgecolor='black',
        # label the x-axis for each plot, based on the TMD
        ax.set_xlabel(xlabel, fontsize=fontsize)
        # move the x-axis label closer to the x-axis
        ax.xaxis.set_label_coords(0.45, -0.085)
        ax.set_ylabel(ylabel, rotation='vertical', fontsize=fontsize)
        # change axis font size
        ax.tick_params(labelsize=fontsize)
        # create legend?#http://stackoverflow.com/questions/9834452/how-do-i-make-a-single-legend-for-many-subplots-with-matplotlib
        legend_obj = ax.legend(cols_to_analyse, loc='upper right', fontsize=fontsize)
        # add title
        # ax.set_title(title,fontsize=fontsize)
        # add figure number to top left of subplot
        ax.annotate(s=str(Fig_Nr) + '.', xy=(0.04, 0.9), fontsize=fontsize, xytext=None,
                                       xycoords='axes fraction', alpha=0.75)
        # add figure title to top left of subplot
        ax.annotate(s=title, xy=(0.1, 0.9), fontsize=fontsize, xytext=None, xycoords='axes fraction',
                                       alpha=0.75)
        utils.save_figure(fig, Fig_name, base_filepath, save_png, save_pdf, dpi)

    if s['Fig14_TMD_10-15_only']:
        Fig_Nr = 14
        title = 'TMD 10-15, all data'
        Fig_name = 'Fig14_TMD_10-15_only'
        col_start = 10
        col_end = 15
        cols_to_analyse = ['TM01_AAIMON_mean'] + list(df_mean_AAIMON_each_TM.columns[col_start:col_end])
        fig, ax = plt.subplots()
        num_bins = 10
        # "#0489B1"
        alpha = 0.9
        col_width_value = 0.95
        ylabel = 'freq'
        xlabel = 'average conservation ratio (membranous over nonmembranous)'

        for n, TM in enumerate(cols_to_analyse):
            # define the colour for that TMD
            # if there are more TMDs than colours, simply start from the beginning of the list again
            if n < len(tableau20):
                color_num = n
            else:
                color_num = n - len(tableau20)
            color = tableau20[color_num]

            hist_data = np.array(df_mean_AAIMON_each_TM[TM].dropna())
            '''
            Calculated the bins for a histogram, even for highly non-normal data
            '''
            # calculate 5th percentile
            percentile_5 = np.percentile(hist_data, 5)  # hist_data.min()
            # calculate 9th percentile
            percentile_95 = np.percentile(hist_data, 95)  # hist_data.max()
            # calculate difference
            percentile_95_minus_5 = percentile_95 - percentile_5
            # create buffer for bins
            extra_xaxis_range = percentile_95_minus_5 / 4
            # lowest bin is the 5th percentile minus the buffer, except where that is below zero
            data_min = percentile_5 - extra_xaxis_range  # hist_data.min()
            # ata_min = 0 if data_max < 0 else data_max
            # highest bin is the 95th percentile
            data_max = percentile_95 + extra_xaxis_range  # hist_data.max()
            # create bins using the min and max
            binlist = np.linspace(data_min, data_max, num_bins)
            # use numpy to create a histogram
            freq_counts_I, bin_array_I = np.histogram(hist_data, bins=binlist)
            # normalize the frequency counts
            freq_counts_normalised = freq_counts_I / freq_counts_I.max()

            # assuming all of the bins are exactly the same size, make the width of the column equal to 70% of each bin
            col_width = float('%0.3f' % (col_width_value * (bin_array_I[1] - bin_array_I[0])))
            # when align='center', the central point of the bar in the x-axis is simply the middle of the bins ((bin_0-bin_1)/2, etc)
            centre_of_bar_in_x_axis = (bin_array_I[:-2] + bin_array_I[1:-1]) / 2
            # add the final bin, which is physically located just after the last regular bin but represents all higher values
            bar_width = centre_of_bar_in_x_axis[3] - centre_of_bar_in_x_axis[2]
            centre_of_bar_in_x_axis = np.append(centre_of_bar_in_x_axis, centre_of_bar_in_x_axis[-1] + bar_width)
            linecontainer = ax.plot(centre_of_bar_in_x_axis, freq_counts_normalised, alpha=alpha,
                                                       linewidth=0.3)  # edgecolor='black',
        # label the x-axis for each plot, based on the TMD
        ax.set_xlabel(xlabel, fontsize=fontsize)
        # move the x-axis label closer to the x-axis
        ax.xaxis.set_label_coords(0.45, -0.085)
        ax.set_ylabel(ylabel, rotation='vertical', fontsize=fontsize)
        # change axis font size
        ax.tick_params(labelsize=fontsize)
        # create legend?#http://stackoverflow.com/questions/9834452/how-do-i-make-a-single-legend-for-many-subplots-with-matplotlib
        legend_obj = ax.legend(cols_to_analyse, loc='upper right', fontsize=fontsize)
        # add title
        # ax.set_title(title,fontsize=fontsize)
        # add figure number to top left of subplot
        ax.annotate(s=str(Fig_Nr) + '.', xy=(0.04, 0.9), fontsize=fontsize, xytext=None,
                                       xycoords='axes fraction', alpha=0.75)
        # add figure title to top left of subplot
        ax.annotate(s=title, xy=(0.1, 0.9), fontsize=fontsize, xytext=None, xycoords='axes fraction',
                                       alpha=0.75)
        utils.save_figure(fig, Fig_name, base_filepath, save_png, save_pdf, dpi)

    if s['Fig15_TMD_15-20_only']:
        Fig_Nr = 15
        title = 'TMD 15-20, all data'
        Fig_name = 'Fig15_TMD_15-20_only'
        col_start = 15
        col_end = 20
        cols_to_analyse = ['TM01_AAIMON_mean'] + list(df_mean_AAIMON_each_TM.columns[col_start:col_end])
        fig, ax = plt.subplots()
        num_bins = 10
        # "#0489B1"
        alpha = 0.9
        col_width_value = 0.95
        ylabel = 'freq'
        xlabel = 'average conservation ratio (membranous over nonmembranous)'

        for n, TM in enumerate(cols_to_analyse):
            # define the colour for that TMD
            # if there are more TMDs than colours, simply start from the beginning of the list again
            if n < len(tableau20):
                color_num = n
            else:
                color_num = n - len(tableau20)
            color = tableau20[color_num]

            hist_data = np.array(df_mean_AAIMON_each_TM[TM].dropna())
            '''
            Calculated the bins for a histogram, even for highly non-normal data
            '''
            # calculate 5th percentile
            percentile_5 = np.percentile(hist_data, 5)  # hist_data.min()
            # calculate 9th percentile
            percentile_95 = np.percentile(hist_data, 95)  # hist_data.max()
            # calculate difference
            percentile_95_minus_5 = percentile_95 - percentile_5
            # create buffer for bins
            extra_xaxis_range = percentile_95_minus_5 / 4
            # lowest bin is the 5th percentile minus the buffer, except where that is below zero
            data_min = percentile_5 - extra_xaxis_range  # hist_data.min()
            # ata_min = 0 if data_max < 0 else data_max
            # highest bin is the 95th percentile
            data_max = percentile_95 + extra_xaxis_range  # hist_data.max()
            # create bins using the min and max
            binlist = np.linspace(data_min, data_max, num_bins)
            # use numpy to create a histogram
            freq_counts_I, bin_array_I = np.histogram(hist_data, bins=binlist)
            # normalize the frequency counts
            freq_counts_normalised = freq_counts_I / freq_counts_I.max()
            # assuming all of the bins are exactly the same size, make the width of the column equal to 70% of each bin
            col_width = float('%0.3f' % (col_width_value * (bin_array_I[1] - bin_array_I[0])))
            # when align='center', the central point of the bar in the x-axis is simply the middle of the bins ((bin_0-bin_1)/2, etc)
            centre_of_bar_in_x_axis = (bin_array_I[:-2] + bin_array_I[1:-1]) / 2
            # add the final bin, which is physically located just after the last regular bin but represents all higher values
            bar_width = centre_of_bar_in_x_axis[3] - centre_of_bar_in_x_axis[2]
            centre_of_bar_in_x_axis = np.append(centre_of_bar_in_x_axis, centre_of_bar_in_x_axis[-1] + bar_width)
            linecontainer = ax.plot(centre_of_bar_in_x_axis, freq_counts_normalised, alpha=alpha,
                                                       linewidth=0.3)  # edgecolor='black',
        # label the x-axis for each plot, based on the TMD
        ax.set_xlabel(xlabel, fontsize=fontsize)
        # move the x-axis label closer to the x-axis
        ax.xaxis.set_label_coords(0.45, -0.085)
        ax.set_ylabel(ylabel, rotation='vertical', fontsize=fontsize)
        # change axis font size
        ax.tick_params(labelsize=fontsize)
        # create legend?#http://stackoverflow.com/questions/9834452/how-do-i-make-a-single-legend-for-many-subplots-with-matplotlib
        legend_obj = ax.legend(cols_to_analyse, loc='upper right', fontsize=fontsize)
        # add title
        # ax.set_title(title,fontsize=fontsize)
        # add figure number to top left of subplot
        ax.annotate(s=str(Fig_Nr) + '.', xy=(0.04, 0.9), fontsize=fontsize, xytext=None,
                                       xycoords='axes fraction', alpha=0.75)
        # add figure title to top left of subplot
        ax.annotate(s=title, xy=(0.1, 0.9), fontsize=fontsize, xytext=None, xycoords='axes fraction',
                                       alpha=0.75)
        utils.save_figure(fig, Fig_name, base_filepath, save_png, save_pdf, dpi)

    if s['Fig16_ONLY_proteins_with_7_TMDs']:
        Fig_Nr = 16
        title = 'ONLY prot with 7 TMDs'
        Fig_name = 'Fig16_ONLY_proteins_with_7_TMDs'
        # create a new dataframe containing data only for proteins with 7 TMDs
        df_seven = df.loc[df.number_of_TMDs == 7].copy()
        # logging.info('df_seven.shape: %s' % str(df_seven.shape))
        # if there are any proteins with 7 TM helices in the dataset (may not be true for beta-barrel datasets)
        if df_seven.shape[0] != 0:
            df_mean_AAIMON_each_TM_7TM, max_num_TMDs, legend = utils.create_df_with_mean_AAIMON_each_TM(df_seven)
            # cols_for_analysis = ['TM01', 'TM06', 'TM07']
            cols_for_analysis = ['TM01_AAIMON_mean', 'TM06_AAIMON_mean', 'TM07_AAIMON_mean']
            fig, ax = plt.subplots()

            num_bins = 30
            # "#0489B1"
            alpha = 0.7
            col_width_value = 0.95
            ylabel = 'freq'
            xlabel = 'average conservation ratio (membranous over nonmembranous)'
            # legend =

            for n, TM in enumerate(cols_for_analysis):
                # define the colour for that TMD
                # if there are more TMDs than colours, simply start from the beginning of the list again
                if n < len(tableau20):
                    color_num = n
                else:
                    color_num = n - len(tableau20)
                color = tableau20[color_num]

                hist_data = np.array(df_mean_AAIMON_each_TM_7TM[TM].dropna())
                '''
                Calculated the bins for a histogram, even for highly non-normal data
                '''
                # calculate 5th percentile
                percentile_5 = np.percentile(hist_data, 5)  # hist_data.min()
                # calculate 9th percentile
                percentile_95 = np.percentile(hist_data, 95)  # hist_data.max()
                # calculate difference
                percentile_95_minus_5 = percentile_95 - percentile_5
                # create buffer for bins
                extra_xaxis_range = percentile_95_minus_5 / 4
                # lowest bin is the 5th percentile minus the buffer, except where that is below zero
                data_min = percentile_5 - extra_xaxis_range  # hist_data.min()
                # ata_min = 0 if data_max < 0 else data_max
                # highest bin is the 95th percentile
                data_max = percentile_95 + extra_xaxis_range  # hist_data.max()
                # create bins using the min and max
                binlist = np.linspace(data_min, data_max, num_bins)
                # use numpy to create a histogram
                freq_counts_I, bin_array_I = np.histogram(hist_data, bins=binlist)
                # assuming all of the bins are exactly the same size, make the width of the column equal to 70% of each bin
                col_width = float('%0.3f' % (col_width_value * (bin_array_I[1] - bin_array_I[0])))
                # when align='center', the central point of the bar in the x-axis is simply the middle of the bins ((bin_0-bin_1)/2, etc)
                centre_of_bar_in_x_axis = (bin_array_I[:-2] + bin_array_I[1:-1]) / 2
                # add the final bin, which is physically located just after the last regular bin but represents all higher values
                bar_width = centre_of_bar_in_x_axis[3] - centre_of_bar_in_x_axis[2]
                centre_of_bar_in_x_axis = np.append(centre_of_bar_in_x_axis, centre_of_bar_in_x_axis[-1] + bar_width)
                linecontainer = ax.plot(centre_of_bar_in_x_axis, freq_counts_I, alpha=alpha,
                                                           linewidth=0.3)  # edgecolor='black',
            # label the x-axis for each plot, based on the TMD
            ax.set_xlabel(xlabel, fontsize=fontsize)
            # move the x-axis label closer to the x-axis
            ax.xaxis.set_label_coords(0.45, -0.085)
            ax.set_ylabel(ylabel, rotation='vertical', fontsize=fontsize)
            # change axis font size
            ax.tick_params(labelsize=fontsize)
            # create legend?#http://stackoverflow.com/questions/9834452/how-do-i-make-a-single-legend-for-many-subplots-with-matplotlib
            legend_obj = ax.legend(cols_for_analysis, loc='upper right', fontsize=fontsize)
            # add title
            # ax.set_title(title,fontsize=fontsize)
            # add figure number to top left of subplot
            ax.annotate(s=str(Fig_Nr) + '.', xy=(0.04, 0.9), fontsize=fontsize, xytext=None,
                                           xycoords='axes fraction', alpha=0.75)
            # add figure title to top left of subplot
            ax.annotate(s=title, xy=(0.1, 0.9), fontsize=fontsize, xytext=None, xycoords='axes fraction',
                                           alpha=0.75)

        utils.save_figure(fig, Fig_name, base_filepath, save_png, save_pdf, dpi)

    if s['Fig17_Less_than_12_TMDs_vs_at_least_12_TMDs']:
        Fig_Nr = 17
        title = '<12 TMDs vs >12 TMDs'
        Fig_name = 'Fig17_Less_than_12_TMDs_vs_at_least_12_TMDs'
        df_under_12 = df.loc[df.number_of_TMDs < 12]

        fig, ax = plt.subplots()
        
        # create numpy array of membranous over nonmembranous conservation ratios (identity)
        hist_data_AAIMON_mean = np.array(df_under_12['AAIMON_mean_all_TMDs'].dropna())
        # use numpy to create a histogram
        freq_counts_I, bin_array_I = np.histogram(hist_data_AAIMON_mean, bins=binlist)
        # assuming all of the bins are exactly the same size, make the width of the column equal to 70% of each bin
        col_width = float('%0.3f' % (0.95 * (bin_array_I[1] - bin_array_I[0])))
        # when align='center', the central point of the bar in the x-axis is simply the middle of the bins ((bin_0-bin_1)/2, etc)
        centre_of_bar_in_x_axis = (bin_array_I[:-2] + bin_array_I[1:-1]) / 2
        # add the final bin, which is physically located just after the last regular bin but represents all higher values
        bar_width = centre_of_bar_in_x_axis[3] - centre_of_bar_in_x_axis[2]
        centre_of_bar_in_x_axis = np.append(centre_of_bar_in_x_axis, centre_of_bar_in_x_axis[-1] + bar_width)
        barcontainer_AAIMON_mean = ax.bar(left=centre_of_bar_in_x_axis, height=freq_counts_I,
                                                             align='center', width=col_width, color="#0489B1",
                                                             alpha=0.5, linewidth=0.1)  # edgecolor='black',
        # other colours that are compatible with colourblind readers: #8A084B Dark red, #B45F04 deep orange, reddish purple #4B088A
        # http://html-color-codes.info/
        # label the x-axis for each plot, based on the TMD
        ax.set_xlabel('average conservation ratio (membranous over nonmembranous)', fontsize=fontsize)
        # move the x-axis label closer to the x-axis
        ax.xaxis.set_label_coords(0.45, -0.085)
        xlim_min = 0.8
        xlim_max = 1.5
        ax.set_xlim(xlim_min, xlim_max)
        # set x-axis ticks
        # use the slide selection to select every second item in the list as an xtick(axis label)
        ax.set_xticks([float('%0.1f' % c) for c in centre_of_bar_in_x_axis[::3]])
        ax.set_ylabel('freq', rotation='vertical', fontsize=fontsize)
        # change axis font size
        ax.tick_params(labelsize=fontsize)


        '''
        At least 12 TMDs
        '''
        df_at_least_12 = df.loc[df.number_of_TMDs >= 12]

        # create numpy array of membranous over nonmembranous conservation ratios (identity)
        hist_data_AAIMON_mean = np.array(df_at_least_12['AAIMON_mean_all_TMDs'].dropna())
        # use numpy to create a histogram
        freq_counts_I, bin_array_I = np.histogram(hist_data_AAIMON_mean, bins=binlist)
        # assuming all of the bins are exactly the same size, make the width of the column equal to 70% of each bin
        col_width = float('%0.3f' % (0.95 * (bin_array_I[1] - bin_array_I[0])))
        # when align='center', the central point of the bar in the x-axis is simply the middle of the bins ((bin_0-bin_1)/2, etc)
        centre_of_bar_in_x_axis = (bin_array_I[:-2] + bin_array_I[1:-1]) / 2
        # add the final bin, which is physically located just after the last regular bin but represents all higher values
        bar_width = centre_of_bar_in_x_axis[3] - centre_of_bar_in_x_axis[2]
        centre_of_bar_in_x_axis = np.append(centre_of_bar_in_x_axis, centre_of_bar_in_x_axis[-1] + bar_width)
        barcontainer_AAIMON_mean = ax.bar(left=centre_of_bar_in_x_axis, height=freq_counts_I * 5,
                                                             align='center', width=col_width, color='#B45F04',
                                                             alpha=0.5, linewidth=0.1)  # edgecolor='black',
        # other colours that are compatible with colourblind readers: #8A084B Dark red, #B45F04 deep orange, reddish purple #4B088A
        # http://html-color-codes.info/
        # label the x-axis for each plot, based on the TMD
        ax.set_xlabel('average conservation ratio (membranous over nonmembranous)', fontsize=fontsize)
        # move the x-axis label closer to the x-axis
        ax.xaxis.set_label_coords(0.45, -0.085)
        xlim_min = 0.8
        xlim_max = 1.5
        ax.set_xlim(xlim_min, xlim_max)
        # set x-axis ticks
        # use the slide selection to select every second item in the list as an xtick(axis label)
        ax.set_xticks([float('%0.1f' % c) for c in centre_of_bar_in_x_axis[::3]])
        ax.set_ylabel('freq', rotation='vertical', fontsize=fontsize)
        # change axis font size
        ax.tick_params(labelsize=fontsize)
        # create legend?#http://stackoverflow.com/questions/9834452/how-do-i-make-a-single-legend-for-many-subplots-with-matplotlib
        legend_obj = ax.legend(['AAIMON 01-11 TMDs', 'AAIMON >= 12 TMDs'], loc='upper right',
                                                  fontsize=fontsize)
        # add figure number to top left of subplot
        ax.annotate(s=str(Fig_Nr) + '.', xy=(0.04, 0.9), fontsize=fontsize, xytext=None,
                                       xycoords='axes fraction', alpha=0.75)
        # add figure title to top left of subplot
        ax.annotate(s=title, xy=(0.1, 0.9), fontsize=fontsize, xytext=None, xycoords='axes fraction',
                                       alpha=0.75)

        utils.save_figure(fig, Fig_name, base_filepath, save_png, save_pdf, dpi)
        
    if s['Fig18_Boxplot_of_all_TMDs']:
        Fig_Nr = 18
        title = 'Boxplot of all TMDs'
        Fig_name = 'Fig18_Boxplot_of_all_TMDs'
        fig, ax = plt.subplots()
        num_bins = 30
        # "#0489B1"
        alpha = 0.25
        col_width_value = 0.95
        ylabel = 'freq'
        xlabel = 'average conservation ratio (membranous over nonmembranous)'
        # legend =

        max_num_TMDs = df.number_of_TMDs.max()
        legend = []
        data_to_plot = []
        for i in range(1, max_num_TMDs.astype(np.int64) + 1):
            TM = 'TM%02d' % i
            hist_data_AAIMON_each_TM = df['TM%02d_AAIMON_mean' % i].dropna()
            if len(hist_data_AAIMON_each_TM) > 0:
                data_to_plot.append(hist_data_AAIMON_each_TM)
                legend.append(TM)

        meanpointprops = dict(marker='o', markerfacecolor='black', markersize=2)  # markeredgecolor='0.75',

        flierprops = dict(marker='o', markerfacecolor='green', markersize=12,
                          linestyle='none')
        boxplotcontainer = ax.boxplot(data_to_plot, sym='+', whis=1.5, showmeans=True,
                                                         meanprops=meanpointprops)
        ax.tick_params(labelsize=fontsize)
        for box in boxplotcontainer['boxes']:
            # change outline color
            box.set(color='black', linewidth=0.4)  # '7570b3'
            # change fill color
            # box.set( facecolor = '#1b9e77' )
            box.set_linewidth(0.4)

        ## change color and linewidth of the whiskers
        for whisker in boxplotcontainer['whiskers']:
            whisker.set(color='black', linewidth=0.4, dashes=(1, 1))

        ## change color and linewidth of the caps
        for cap in boxplotcontainer['caps']:
            cap.set(color='black', linewidth=0.4)

        ## change color and linewidth of the medians
        for median in boxplotcontainer['medians']:
            median.set(color='black', linewidth=0.4)

        # change the style of fliers and their fill
        for flier in boxplotcontainer['fliers']:
            flier.set(marker='o', color='0.8', alpha=0.1, markerfacecolor='0.3', markersize=3)

        ax.set_ylabel('AAIMON', rotation='vertical', fontsize=fontsize)

        ## Remove top axes and right axes ticks
        ax.get_xaxis().tick_bottom()
        ax.get_yaxis().tick_left()
        ## Custom x-axis labels
        ax.set_xticklabels(legend, rotation=45)
        # add figure number to top left of subplot
        ax.annotate(s=str(Fig_Nr) + '.', xy=(0.04, 0.9), fontsize=fontsize, xytext=None,
                                       xycoords='axes fraction', alpha=0.75)
        # add figure title to top left of subplot
        ax.annotate(s=title, xy=(0.1, 0.9), fontsize=fontsize, xytext=None, xycoords='axes fraction',
                                       alpha=0.75)

        utils.save_figure(fig, Fig_name, base_filepath, save_png, save_pdf, dpi)

    if 'uniprot_KW' in df.columns:

        if s['Fig19_Show_only_GPCRs_in_full_dataset']:
            if GPCR_in_df:
                Fig_Nr = 19
                title = 'only GPCR in uniprot KW, NORM'
                Fig_name = 'Fig19_Show_only_GPCRs_in_full_dataset'
                fig, ax = plt.subplots()
                # create numpy array of membranous over nonmembranous conservation ratios (identity)
                hist_data_AAIMON_mean = np.array(df_GPCR['AAIMON_mean_all_TMDs'].dropna())
                # use numpy to create a histogram
                freq_counts_I, bin_array_I = np.histogram(hist_data_AAIMON_mean, bins=binlist)
                # normalize the frequency counts
                freq_counts_normalised = freq_counts_I / freq_counts_I.max()
                # assuming all of the bins are exactly the same size, make the width of the column equal to 70% of each bin
                col_width = float('%0.3f' % (0.95 * (bin_array_I[1] - bin_array_I[0])))
                # when align='center', the central point of the bar in the x-axis is simply the middle of the bins ((bin_0-bin_1)/2, etc)
                centre_of_bar_in_x_axis = (bin_array_I[:-2] + bin_array_I[1:-1]) / 2
                # add the final bin, which is physically located just after the last regular bin but represents all higher values
                bar_width = centre_of_bar_in_x_axis[3] - centre_of_bar_in_x_axis[2]
                centre_of_bar_in_x_axis = np.append(centre_of_bar_in_x_axis, centre_of_bar_in_x_axis[-1] + bar_width)
                barcontainer_AAIMON_mean = ax.bar(left=centre_of_bar_in_x_axis,
                                                  height=freq_counts_normalised,
                                                  align='center', width=col_width, color="#0489B1",
                                                  alpha=0.5, linewidth=0.1)  # edgecolor='black',
                # other colours that are compatible with colourblind readers: #8A084B Dark red, #B45F04 deep orange, reddish purple #4B088A
                # http://html-color-codes.info/
                # label the x-axis for each plot, based on the TMD
                ax.set_xlabel('average conservation ratio (membranous over nonmembranous)',
                              fontsize=fontsize)
                # move the x-axis label closer to the x-axis
                ax.xaxis.set_label_coords(0.45, -0.085)
                xlim_min = 0.8
                xlim_max = 1.5
                ax.set_xlim(xlim_min, xlim_max)
                # set x-axis ticks
                # use the slide selection to select every second item in the list as an xtick(axis label)
                ax.set_xticks([float('%0.1f' % c) for c in centre_of_bar_in_x_axis[::3]])
                ax.set_ylabel('freq', rotation='vertical', fontsize=fontsize)
                # change axis font size
                ax.tick_params(labelsize=fontsize)

                '''
                NON-GPCRS
                '''
                df_nonGPCR = df.loc[df['G-protein_coupled_receptor'] == False]

                # create numpy array of membranous over nonmembranous conservation ratios (identity)
                hist_data_AAIMON_mean = np.array(df_nonGPCR['AAIMON_mean_all_TMDs'].dropna())
                # use numpy to create a histogram
                freq_counts_I, bin_array_I = np.histogram(hist_data_AAIMON_mean, bins=binlist)
                # normalize the frequency counts
                freq_counts_normalised = freq_counts_I / freq_counts_I.max()
                # assuming all of the bins are exactly the same size, make the width of the column equal to 70% of each bin
                col_width = float('%0.3f' % (0.95 * (bin_array_I[1] - bin_array_I[0])))
                # when align='center', the central point of the bar in the x-axis is simply the middle of the bins ((bin_0-bin_1)/2, etc)
                centre_of_bar_in_x_axis = (bin_array_I[:-2] + bin_array_I[1:-1]) / 2
                # add the final bin, which is physically located just after the last regular bin but represents all higher values
                bar_width = centre_of_bar_in_x_axis[3] - centre_of_bar_in_x_axis[2]
                centre_of_bar_in_x_axis = np.append(centre_of_bar_in_x_axis, centre_of_bar_in_x_axis[-1] + bar_width)
                barcontainer_AAIMON_mean = ax.bar(left=centre_of_bar_in_x_axis,
                                                  height=freq_counts_normalised,
                                                  align='center', width=col_width, color='#B45F04',
                                                  alpha=0.5, linewidth=0.1)  # edgecolor='black',
                # other colours that are compatible with colourblind readers: #8A084B Dark red, #B45F04 deep orange, reddish purple #4B088A
                # http://html-color-codes.info/
                # label the x-axis for each plot, based on the TMD
                ax.set_xlabel('average conservation ratio (membranous over nonmembranous)',
                              fontsize=fontsize)
                # move the x-axis label closer to the x-axis
                ax.xaxis.set_label_coords(0.45, -0.085)
                xlim_min = 0.8
                xlim_max = 1.5
                ax.set_xlim(xlim_min, xlim_max)
                # set x-axis ticks
                # use the slide selection to select every second item in the list as an xtick(axis label)
                ax.set_xticks([float('%0.1f' % c) for c in centre_of_bar_in_x_axis[::3]])
                ax.set_ylabel('freq', rotation='vertical', fontsize=fontsize)
                # change axis font size
                ax.tick_params(labelsize=fontsize)
                # create legend?#http://stackoverflow.com/questions/9834452/how-do-i-make-a-single-legend-for-many-subplots-with-matplotlib
                legend_obj = ax.legend(['AAIMON GPCRs', 'AAIMON non-GPCRs'], loc='upper right',
                                       fontsize=fontsize)
                # add figure number to top left of subplot
                ax.annotate(s=str(Fig_Nr) + '.', xy=(0.04, 0.9), fontsize=fontsize, xytext=None,
                            xycoords='axes fraction', alpha=0.75)
                # add figure title to top left of subplot
                ax.annotate(s=title, xy=(0.1, 0.9), fontsize=fontsize, xytext=None, xycoords='axes fraction',
                            alpha=0.75)

                utils.save_figure(fig, Fig_name, base_filepath, save_png, save_pdf, dpi)
            else:
                sys.stdout.write('Dataset does not contain GPCRs; cannot create figure 19 \n')

        if s['Fig20_Show_GPCRs_vs_full_dataset']:
            if GPCR_in_df:
                Fig_Nr = 20
                title = 'GPCR vs full dataset, NORM'
                Fig_name = 'Fig20_Show_GPCRs_vs_full_dataset'
                df_GPCR = df_GPCR  # see above, already defined
                fig, ax = plt.subplots()

                '''GPCR'''

                # create numpy array of membranous over nonmembranous conservation ratios (identity)
                hist_data_AAIMON_mean = np.array(df_GPCR['AAIMON_mean_all_TMDs'].dropna())
                # use numpy to create a histogram
                freq_counts_I, bin_array_I = np.histogram(hist_data_AAIMON_mean, bins=binlist)
                # normalize the frequency counts
                freq_counts_normalised = freq_counts_I / freq_counts_I.max()
                # assuming all of the bins are exactly the same size, make the width of the column equal to 70% of each bin
                col_width = float('%0.3f' % (0.95 * (bin_array_I[1] - bin_array_I[0])))
                # when align='center', the central point of the bar in the x-axis is simply the middle of the bins ((bin_0-bin_1)/2, etc)
                centre_of_bar_in_x_axis = (bin_array_I[:-2] + bin_array_I[1:-1]) / 2
                # add the final bin, which is physically located just after the last regular bin but represents all higher values
                bar_width = centre_of_bar_in_x_axis[3] - centre_of_bar_in_x_axis[2]
                centre_of_bar_in_x_axis = np.append(centre_of_bar_in_x_axis, centre_of_bar_in_x_axis[-1] + bar_width)
                barcontainer_AAIMON_mean = ax.bar(left=centre_of_bar_in_x_axis,
                                                  height=freq_counts_normalised,
                                                  align='center', width=col_width, color="#0489B1",
                                                  alpha=0.5, linewidth=0.1)  # edgecolor='black',
                # other colours that are compatible with colourblind readers: #8A084B Dark red, #B45F04 deep orange, reddish purple #4B088A
                # http://html-color-codes.info/
                # label the x-axis for each plot, based on the TMD
                ax.set_xlabel('average conservation ratio (membranous over nonmembranous)',
                              fontsize=fontsize)
                # move the x-axis label closer to the x-axis
                ax.xaxis.set_label_coords(0.45, -0.085)
                xlim_min = 0.8
                xlim_max = 1.5
                ax.set_xlim(xlim_min, xlim_max)
                # set x-axis ticks
                # use the slide selection to select every second item in the list as an xtick(axis label)
                ax.set_xticks([float('%0.1f' % c) for c in centre_of_bar_in_x_axis[::3]])
                ax.set_ylabel('freq', rotation='vertical', fontsize=fontsize)
                # change axis font size
                ax.tick_params(labelsize=fontsize)

                '''
                full dataset
                '''
                # df_nonGPCR = df.loc[df['Gprotein'] == False]

                # create numpy array of membranous over nonmembranous conservation ratios (identity)
                hist_data_AAIMON_mean = np.array(df['AAIMON_mean_all_TMDs'].dropna())
                # use numpy to create a histogram
                freq_counts_I, bin_array_I = np.histogram(hist_data_AAIMON_mean, bins=binlist)
                # normalize the frequency counts
                freq_counts_normalised = freq_counts_I / freq_counts_I.max()
                # assuming all of the bins are exactly the same size, make the width of the column equal to 70% of each bin
                col_width = float('%0.3f' % (0.95 * (bin_array_I[1] - bin_array_I[0])))
                # when align='center', the central point of the bar in the x-axis is simply the middle of the bins ((bin_0-bin_1)/2, etc)
                centre_of_bar_in_x_axis = (bin_array_I[:-2] + bin_array_I[1:-1]) / 2
                # add the final bin, which is physically located just after the last regular bin but represents all higher values
                bar_width = centre_of_bar_in_x_axis[3] - centre_of_bar_in_x_axis[2]
                centre_of_bar_in_x_axis = np.append(centre_of_bar_in_x_axis, centre_of_bar_in_x_axis[-1] + bar_width)
                barcontainer_AAIMON_mean = ax.bar(left=centre_of_bar_in_x_axis,
                                                  height=freq_counts_normalised,
                                                  align='center', width=col_width, color='#B45F04',
                                                  alpha=0.5, linewidth=0.1)  # edgecolor='black',
                # other colours that are compatible with colourblind readers: #8A084B Dark red, #B45F04 deep orange, reddish purple #4B088A
                # http://html-color-codes.info/
                # label the x-axis for each plot, based on the TMD
                ax.set_xlabel('average conservation ratio (membranous over nonmembranous)',
                              fontsize=fontsize)
                # move the x-axis label closer to the x-axis
                ax.xaxis.set_label_coords(0.45, -0.085)
                xlim_min = 0.8
                xlim_max = 1.5
                ax.set_xlim(xlim_min, xlim_max)
                # set x-axis ticks
                # use the slide selection to select every second item in the list as an xtick(axis label)
                ax.set_xticks([float('%0.1f' % c) for c in centre_of_bar_in_x_axis[::3]])
                ax.set_ylabel('freq', rotation='vertical', fontsize=fontsize)
                # change axis font size
                ax.tick_params(labelsize=fontsize)
                # create legend?#http://stackoverflow.com/questions/9834452/how-do-i-make-a-single-legend-for-many-subplots-with-matplotlib
                legend_obj = ax.legend(['AAIMON GPCR', 'AAIMON ALL'], loc='upper right',
                                       fontsize=fontsize)
                # add figure number to top left of subplot
                ax.annotate(s=str(Fig_Nr) + '.', xy=(0.04, 0.9), fontsize=fontsize, xytext=None,
                            xycoords='axes fraction', alpha=0.75)
                # add figure title to top left of subplot
                ax.annotate(s=title, xy=(0.1, 0.9), fontsize=fontsize, xytext=None, xycoords='axes fraction',
                            alpha=0.75)

                utils.save_figure(fig, Fig_name, base_filepath, save_png, save_pdf, dpi)

            else:
                sys.stdout.write('Dataset does not contain GPCRs; cannot create figure 20 \n')

        if s['Fig21_Boxplot_only_GPCRs']:
            if GPCR_in_df:
                Fig_Nr = 21
                title = 'Only GPCRs, boxplot for each TMD'
                Fig_name = 'Fig21_Boxplot_only_GPCRs'
                fig, ax = plt.subplots()

                num_bins = 30
                # "#0489B1"
                alpha = 0.25
                col_width_value = 0.95
                ylabel = 'freq'
                xlabel = 'average conservation ratio (membranous over nonmembranous)'
                # legend =

                max_num_TMDs = df_GPCR.number_of_TMDs.max()
                legend = []
                data_to_plot = []
                for i in range(1, max_num_TMDs.astype(np.int64) + 1):
                    TM = 'TM%02d' % i
                    hist_data_AAIMON_each_TM = df_GPCR['TM%02d_AAIMON_mean' % i].dropna()
                    if len(hist_data_AAIMON_each_TM) > 0:
                        data_to_plot.append(hist_data_AAIMON_each_TM)
                        legend.append(TM)

                meanpointprops = dict(marker='o', markerfacecolor='black', markersize=3)  # markeredgecolor='0.75',

                # flierprops = dict(marker='o', color = 'black', markerfacecolor='black', markersize=1)

                # flierprops = dict(marker='o',color='0.1', alpha=0.1)
                flierprops = dict(marker='o', markerfacecolor='green', markersize=12,
                                  linestyle='none')
                boxplotcontainer = ax.boxplot(data_to_plot, sym='+', whis=1.5, showmeans=True,
                                              meanprops=meanpointprops)
                ax.tick_params(labelsize=fontsize)
                for box in boxplotcontainer['boxes']:
                    # change outline color
                    box.set(color='black', linewidth=0.4)  # '7570b3'
                    # change fill color
                    # box.set( facecolor = '#1b9e77' )
                    box.set_linewidth(0.4)

                ## change color and linewidth of the whiskers
                for whisker in boxplotcontainer['whiskers']:
                    whisker.set(color='black', linewidth=0.4, dashes=(1, 1))

                ## change color and linewidth of the caps
                for cap in boxplotcontainer['caps']:
                    cap.set(color='black', linewidth=0.4)

                ## change color and linewidth of the medians
                for median in boxplotcontainer['medians']:
                    median.set(color='black', linewidth=0.4)

                # change the style of fliers and their fill
                for flier in boxplotcontainer['fliers']:
                    flier.set(marker='o', color='0.8', alpha=0.1, markerfacecolor='0.3', markersize=3)

                ## Remove top axes and right axes ticks
                ax.get_xaxis().tick_bottom()
                ax.get_yaxis().tick_left()
                ## Custom x-axis labels
                ax.set_xticklabels(legend, rotation=45)
                # add figure number to top left of subplot
                ax.annotate(s=str(Fig_Nr) + '.', xy=(0.04, 0.9), fontsize=fontsize, xytext=None,
                            xycoords='axes fraction', alpha=0.75)
                # add figure title to top left of subplot
                ax.annotate(s=title, xy=(0.1, 0.9), fontsize=fontsize, xytext=None, xycoords='axes fraction',
                            alpha=0.75)

                utils.save_figure(fig, Fig_name, base_filepath, save_png, save_pdf, dpi)
            else:
                sys.stdout.write('Dataset does not contain GPCRs; cannot create figure 21 \n')

    if s['Fig22_Boxplot_comparing_number_of_TMDs_with_mean_AAIMON']:
        Fig_Nr = 22
        title = 'num_TMDs vs AAIMON'
        Fig_name = 'Fig22_Boxplot_comparing_number_of_TMDs_with_mean_AAIMON'
        fig, ax = plt.subplots()

        alpha = 0.25
        col_width_value = 0.95
        max_num_TMDs = df.number_of_TMDs.max()
        legend = []
        data_to_plot = []
        # iterate through df and get all AAIMONs with specified number of TMD
        for i in range(1, max_num_TMDs.astype(np.int64) + 1):
            hist_data = []
            for acc in df.loc[df['list_of_TMDs'].notnull()].loc[df['list_of_TMDs'] != 'nan'].index:
                if df.loc[acc, 'number_of_TMDs'] == i:
                    hist_data.append(df.loc[acc, 'AAIMON_mean_all_TMDs'])
            data_to_plot.append(hist_data)
            legend.append(i)
        meanpointprops = dict(marker='o', markerfacecolor='black', markersize=3)  # markeredgecolor='0.75',
        flierprops = dict(marker='o', markerfacecolor='green', markersize=12,
                          linestyle='none')
        boxplotcontainer = ax.boxplot(data_to_plot, sym='+', whis=1.5, showmeans=True,
                                      meanprops=meanpointprops)
        ax.tick_params(labelsize=fontsize)
        for box in boxplotcontainer['boxes']:
            # change outline color
            box.set(color='black', linewidth=0.4)  # '7570b3'
            # change fill color
            # box.set( facecolor = '#1b9e77' )
            box.set_linewidth(0.4)

        ## change color and linewidth of the whiskers
        for whisker in boxplotcontainer['whiskers']:
            whisker.set(color='black', linewidth=0.4, dashes=(1, 1))

        ## change color and linewidth of the caps
        for cap in boxplotcontainer['caps']:
            cap.set(color='black', linewidth=0.4)

        ## change color and linewidth of the medians
        for median in boxplotcontainer['medians']:
            median.set(color='black', linewidth=0.4)

        # change the style of fliers and their fill
        for flier in boxplotcontainer['fliers']:
            flier.set(marker='o', color='0.8', alpha=0.1, markerfacecolor='0.3', markersize=3)

        ax.set_xlabel('number of TMDs in protein', fontsize=fontsize)
        # move the x-axis label closer to the x-axis
        ax.xaxis.set_label_coords(0.45, -0.085)
        ax.set_ylabel('Average AAIMON ratio for all TMDs', fontsize=fontsize)
        # Remove top axes and right axes ticks
        ax.get_xaxis().tick_bottom()
        ax.get_yaxis().tick_left()
        # Custom x-axis labels
        ax.set_xticklabels(legend)
        # add figure number to top left of subplot
        ax.annotate(s=str(Fig_Nr) + '.', xy=(0.04, 0.9), fontsize=fontsize, xytext=None,
                    xycoords='axes fraction', alpha=0.75)
        # add figure title to top left of subplot
        ax.annotate(s=title, xy=(0.1, 0.9), fontsize=fontsize, xytext=None, xycoords='axes fraction',
                    alpha=0.75)

        utils.save_figure(fig, Fig_name, base_filepath, save_png, save_pdf, dpi)

    if s['Fig23_Boxplot_comparing_seqlen_with_mean_AAIMON']:
        Fig_Nr = 23
        title = 'num_TMDs vs seqlen'
        Fig_name = 'Fig23_Boxplot_comparing_seqlen_with_mean_AAIMON'
        fig, ax = plt.subplots()

        num_bins = 30
        # "#0489B1"
        alpha = 0.25
        col_width_value = 0.95

        max_seqlen = df.seqlen.max()
        size_of_bin = 500  # amino acids
        legend = []
        data_to_plot = []

        # iterate through df and get all AAIMONs with specified number of TMD

        for i in range(1, max_seqlen.astype(np.int64) + 1000, size_of_bin):
            hist_data = []
            for acc in df.loc[df['list_of_TMDs'].notnull()].loc[df['list_of_TMDs'] != 'nan'].index:
                if df.loc[acc, 'seqlen'] > i and df.loc[acc, 'seqlen'] < i + size_of_bin:
                    hist_data.append(df.loc[acc, 'AAIMON_mean_all_TMDs'])
            data_to_plot.append(hist_data)
            legend.append(i + (size_of_bin / 2 - 1))

        meanpointprops = dict(marker='o', markerfacecolor='black', markersize=3)  # markeredgecolor='0.75',

        # flierprops = dict(marker='o', color = 'black', markerfacecolor='black', markersize=1)

        # flierprops = dict(marker='o',color='0.1', alpha=0.1)
        flierprops = dict(marker='o', markerfacecolor='green', markersize=12,
                          linestyle='none')
        boxplotcontainer = ax.boxplot(data_to_plot, sym='+', whis=1.5, showmeans=True,
                                      meanprops=meanpointprops)
        ax.tick_params(labelsize=fontsize)
        for box in boxplotcontainer['boxes']:
            # change outline color
            box.set(color='black', linewidth=0.4)  # '7570b3'
            # change fill color
            # box.set( facecolor = '#1b9e77' )
            box.set_linewidth(0.4)

        # change color and linewidth of the whiskers
        for whisker in boxplotcontainer['whiskers']:
            whisker.set(color='black', linewidth=0.4, dashes=(1, 1))

        # change color and linewidth of the caps
        for cap in boxplotcontainer['caps']:
            cap.set(color='black', linewidth=0.4)

        # change color and linewidth of the medians
        for median in boxplotcontainer['medians']:
            median.set(color='black', linewidth=0.4)

        # change the style of fliers and their fill
        for flier in boxplotcontainer['fliers']:
            flier.set(marker='o', color='0.8', alpha=0.1, markerfacecolor='0.3', markersize=3)

        ax.set_xlabel('Length of protein in bins of 500 amino acids', fontsize=fontsize)
        # move the x-axis label closer to the x-axis
        ax.xaxis.set_label_coords(0.45, -0.085)
        ax.set_ylabel('Average AAIMON ratio for all TMDs', fontsize=fontsize)
        ## Remove top axes and right axes ticks
        ax.get_xaxis().tick_bottom()
        ax.get_yaxis().tick_left()
        ## Custom x-axis labels
        ax.set_xticklabels(legend, rotation=25)
        # add figure number to top left of subplot
        ax.annotate(s=str(Fig_Nr) + '.', xy=(0.04, 0.9), fontsize=fontsize, xytext=None,
                    xycoords='axes fraction', alpha=0.75)
        # add figure title to top left of subplot
        ax.annotate(s=title, xy=(0.1, 0.9), fontsize=fontsize, xytext=None, xycoords='axes fraction',
                    alpha=0.75)

        utils.save_figure(fig, Fig_name, base_filepath, save_png, save_pdf, dpi)

    if s['Fig24_Scatterplot_AAIMON_n_vs._AAIMON_n_slope']:
        Fig_Nr = 24
        title = 'AAIMON_n vs. AAIMON_n_slope'
        Fig_name = 'Fig24_Scatterplot_AAIMON_n_vs._AAIMON_n_slope'
        fig, ax = plt.subplots()

        x = df['AAIMON_mean_all_TMDs_n']
        y = df['AAIMON_n_slope_mean_all_TMDs']

        if len(x) > 5:
            # calculate linear regression for fitted line
            linear_regression = np.polyfit(x, y, 1)
            fit_fn = np.poly1d(linear_regression)
            fitted_data_x = fit_fn(x)
            ax.plot(x, fitted_data_x, alpha=0.75, color='k')
            ax.annotate(s='y = {a:.5f}x + {b:.5f}'.format(a=linear_regression[0], b=linear_regression[1]), xy=(0.85, 0.95),
                        fontsize=fontsize-2, xytext=None, xycoords='axes fraction', alpha=0.75)
        else:
            logging.info("The dataset has less than 5 proteins. Lines of best fit will not be calculated.")

        ax.scatter(x, y, alpha=alpha_dpd, s=datapointsize)
        ax.set_ylabel('AAIMON_n_slope', rotation='vertical', fontsize=fontsize)
        ax.set_xlabel('AAIMON_n', fontsize=fontsize)
        ax.annotate(s=str(Fig_Nr) + '.', xy=(0.04, 0.9), fontsize=fontsize, xytext=None,
                    xycoords='axes fraction', alpha=0.75)
        # add figure title to top left of subplot
        ax.annotate(s=title, xy=(0.1, 0.9), fontsize=fontsize, xytext=None, xycoords='axes fraction',
                    alpha=0.75)

        # change axis font size
        ax.tick_params(labelsize=fontsize)

        utils.save_figure(fig, Fig_name, base_filepath, save_png, save_pdf, dpi)


    if s['Fig25_Scatterplot_AAIMON_n_vs_obs_changes_mean']:
        Fig_Nr = 25
        title = 'AAIMON_n  vs. obs_changes_mean'
        Fig_name = 'Fig25_Scatterplot_AAIMON_n_vs_obs_changes_mean'
        fig, ax = plt.subplots()

        ax.scatter(df['obs_changes_mean'], df['AAIMON_mean_all_TMDs_n'], color='b', alpha=alpha_dpd, s=datapointsize)
        ax.set_ylabel('AAIMON_n', rotation='vertical', fontsize=fontsize)
        ax.set_xlabel('obs_changes_mean', fontsize=fontsize)

        #ax2 = ax.twinx()
        #ax2.scatter(df['obs_changes_mean'], df['AAIMON_mean_all_TMDs_n'], color="r", alpha=alpha_dpd, s=datapointsize)
        #ax2.set_ylabel('AAIMON_n', rotation='vertical', color='r', fontsize=fontsize)

        ax.set_xlim(0, 60)
        ax.set_ylim(0.2, 1.8)
        #ax2.set_ylim(0.4, 1.6)

        ax.annotate(s=str(Fig_Nr) + '.', xy=(0.04, 0.9), fontsize=fontsize, xytext=None,
                    xycoords='axes fraction', alpha=0.75)
        # add figure title to top left of subplot
        ax.annotate(s=title, xy=(0.1, 0.9), fontsize=fontsize, xytext=None, xycoords='axes fraction',
                    alpha=0.75)

        # change axis font size
        ax.tick_params(labelsize=fontsize)
        #ax2.tick_params(labelsize=fontsize)

        utils.save_figure(fig, Fig_name, base_filepath, save_png, save_pdf, dpi)

    if s['Fig26_Scatterplot_AAIMON_n_slope_vs_obs_changes_mean']:
        Fig_Nr = 26
        title = 'AAIMON_n_slope  vs. obs_changes_mean'
        Fig_name = 'Fig26_Scatterplot_AAIMON_n_slope_vs_obs_changes_mean'
        fig, ax = plt.subplots()

        ax.scatter(df['obs_changes_mean'], df['AAIMON_n_slope_mean_all_TMDs'], color="r", alpha=alpha_dpd, s=datapointsize)
        ax.set_ylabel('AAIMON_n_slope', rotation='vertical', fontsize=fontsize)
        ax.set_xlabel('obs_changes_mean', fontsize=fontsize)


        ax.set_xlim(0, 60)
        ax.set_ylim(-0.020, 0.020)
        # ax2.set_ylim(0.4, 1.6)

        ax.annotate(s=str(Fig_Nr) + '.', xy=(0.04, 0.9), fontsize=fontsize, xytext=None,
                    xycoords='axes fraction', alpha=0.75)
        # add figure title to top left of subplot
        ax.annotate(s=title, xy=(0.1, 0.9), fontsize=fontsize, xytext=None, xycoords='axes fraction',
                    alpha=0.75)

        # change axis font size
        ax.tick_params(labelsize=fontsize)
        # ax2.tick_params(labelsize=fontsize)

        utils.save_figure(fig, Fig_name, base_filepath, save_png, save_pdf, dpi)

    if s['Fig27_Scatterplot_perc_identity_nonTMD_vs_TMD']:
        Fig_Nr = 27
        title = 'perc identity all TMDs  vs. perc identity nonTMD'
        Fig_name = 'Fig27_Scatterplot_perc_identity_nonTMD_vs_TMD'
        fig, ax = plt.subplots()

        x = df['TMD_perc_identity_mean_all_TMDs'] * 100
        y = df['nonTMD_perc_ident_mean'] * 100

        if len(x) > 5:
            linear_regression = np.polyfit(x, y, 1)
            fit_fn = np.poly1d(linear_regression)
            fitted_data_x = fit_fn(x)
            ax.plot(x, fitted_data_x, alpha=0.75, color='k')
            ax.annotate(s='y = {a:.5f}x + {b:.5f}'.format(a=linear_regression[0], b=linear_regression[1]),
                        xy=(0.85, 0.95), fontsize=fontsize - 2, xytext=None, xycoords='axes fraction',alpha=0.75)

        ax.scatter(x, y, s=datapointsize, alpha=alpha_dpd, color='r')

        ax.set_xlabel('TMD_perc_identity_all_TMDs', fontsize=fontsize)
        ax.set_ylabel('nonTMD_perc_ident_mean', rotation='vertical', fontsize=fontsize)
        ax.tick_params(labelsize=fontsize)
        ax.set_xlim(40, 100)
        ax.set_ylim(40, 100)

        ax.annotate(s=str(Fig_Nr) + '.', xy=(0.04, 0.9), fontsize=fontsize, xytext=None,
                    xycoords='axes fraction', alpha=0.75)
        # add figure title to top left of subplot
        ax.annotate(s=title, xy=(0.1, 0.9), fontsize=fontsize, xytext=None, xycoords='axes fraction',
                    alpha=0.75)

        # change axis font size
        ax.tick_params(labelsize=fontsize)

        utils.save_figure(fig, Fig_name, base_filepath, save_png, save_pdf, dpi)


    if s['Fig98_Scatterplot_AAIMON_vs_perc_ident_all_homol_all_proteins']:
        Fig_Nr = 98
        Fig_name = 'Fig98_Scatterplot_AAIMON_vs_perc_ident_all_homol_all_proteins'

        # read data from disk
        in_zipfile = pathdict["save_df_characterising_each_homol_TMD"]
        if os.path.isfile(in_zipfile):
            with zipfile.ZipFile(in_zipfile, "r", zipfile.ZIP_DEFLATED) as openzip:
                data = pickle.load(openzip.open("data_characterising_each_homol_TMD.pickle", "r"))
                binned_data = pickle.load(openzip.open("binned_data_characterising_each_homol_TMD.pickle", "r"))
        else:
            raise FileNotFoundError("{} not found".format(in_zipfile))

        fontsize = 14
        datapointsize = 1
        alpha = 0.05
        linewidth = 2
        color_nonnorm = "#EE762C"
        color_norm = "#0076B8"
        fig, ax = plt.subplots()

        # set color of axis label to black
        ax.tick_params(axis='x', colors='black')
        ax.tick_params(axis='y', colors='black')
        ax.yaxis.label.set_color('black')
        ax.xaxis.label.set_color('black')

        # pylab.rcParams['figure.figsize'] = (50.0, 40.0)
        x = data[:, 0]  # FASTA_gapped_identity
        y = data[:, 1]  # AAIMON for each TMD
        ax.scatter(x=x, y=y, color=color_nonnorm, alpha=alpha, s=datapointsize)  # color="#003366" is TUM-blue
        plt.ylim(ymin=0, ymax=3)
        plt.xlim(xmin=0, xmax=60)
        # label the x-axis for each plot, based on the TMD
        ax.set_xlabel('% observed changes', fontsize=fontsize)
        # move the x-axis label closer to the x-axis
        # ax.xaxis.set_label_coords(0.45, -0.085)
        ax.set_ylabel('AAIMON ratio', fontsize=fontsize)
        # change axis font size
        ax.tick_params(labelsize=fontsize)
        x_line = binned_data[:, 0]
        y_line = binned_data[:, 1]
        plt.plot(x_line, y_line, linewidth=linewidth, color=color_nonnorm)  # plot linegraph

        # plot normalised data
        x = data[:, 0]  # FASTA_gapped_identity
        y = data[:, 2]  # AAIMON_n for each TMD
        ax.scatter(x=x, y=y, color=color_norm, alpha=alpha, s=datapointsize)  # color="#FF6633" is TUM-orange
        x_line = binned_data[:, 0]
        y_line = binned_data[:, 2]
        plt.plot(x_line, y_line, linewidth=linewidth, color=color_norm)  # plot linegraph

        # remove unwanted axis ticks on top and right
        ax.get_xaxis().tick_bottom()
        ax.get_yaxis().tick_left()

        #ax.xaxis.set_label_coords(0.45, -0.1)  # move x-axis label towards graph
        #ax.yaxis.set_label_coords(-0.09, 0.45)  # move y-axis label towards graph

        ax.legend(['AAIMON', 'AAIMON norm.'], loc='upper right', fontsize=fontsize)  # create legend

        utils.save_figure(fig, Fig_name, base_filepath, save_png=True, save_pdf=False, dpi=150, close=False)
        utils.save_figure(fig, Fig_name + "_highres", base_filepath, save_png=True, save_pdf=False, dpi=600)

    if s['Fig99_Linegraph_CI_95_AAIMON_vs_perc_ident_all_homol_all_proteins']:
        Fig_Nr = 99
        Fig_name = 'Fig99_Linegraph_CI_95_AAIMON_vs_perc_ident_all_homol_all_proteins'

        # read data from disk
        in_zipfile = pathdict["save_df_characterising_each_homol_TMD"]
        if os.path.isfile(in_zipfile):
            with zipfile.ZipFile(in_zipfile, "r", zipfile.ZIP_DEFLATED) as openzip:
                binned_data = pickle.load(openzip.open("binned_data_characterising_each_homol_TMD.pickle", "r"))
        else:
            raise FileNotFoundError("{} not found".format(in_zipfile))

        fontsize = 14
        datapointsize = 1
        alpha = 0.05
        alpha_line = 1
        linewidth = 1
        color_nonnorm = "#EE762C"
        color_norm = "#0076B8"
        fig, ax = plt.subplots()

        # set color of axis label to black
        ax.tick_params(axis='x', colors='black')
        ax.tick_params(axis='y', colors='black')
        ax.yaxis.label.set_color('black')
        ax.xaxis.label.set_color('black')

        # pylab.rcParams['figure.figsize'] = (50.0, 40.0)
        # x = data[:,0] # FASTA_gapped_identity
        # y = data[:,1] # AAIMON for each TMD
        # ax.scatter(x=x, y=y, color=color_nonnorm, alpha=alpha, s=datapointsize) # color="#003366" is TUM-blue
        plt.ylim(ymin=0.5, ymax=2)
        plt.xlim(xmin=0, xmax=60)
        # label the x-axis for each plot, based on the TMD
        ax.set_xlabel('% observed changes', fontsize=fontsize)
        ax.set_ylabel('AAIMON ratio', fontsize=fontsize)
        # change axis font size
        ax.tick_params(labelsize=fontsize)

        # plot AAIMON
        plt.plot(binned_data[:, 0], binned_data[:, 1], linewidth=linewidth, color=color_nonnorm,
                 alpha=alpha_line)  # plot linegraph
        # plot AAIMON_n
        plt.plot(binned_data[:, 0], binned_data[:, 2], linewidth=linewidth, color=color_norm,
                 alpha=alpha_line)  # plot linegraph

        # plot 95% confidence intervals
        # AAIMON
        plt.plot(binned_data[:, 0], binned_data[:, 3], linewidth=linewidth, color=color_nonnorm, alpha=alpha_line - 0.5)
        plt.plot(binned_data[:, 0], binned_data[:, 4], linewidth=linewidth, color=color_nonnorm, alpha=alpha_line - 0.5)
        # AAIMON_n
        plt.plot(binned_data[:, 0], binned_data[:, 5], linewidth=linewidth, color=color_norm, alpha=alpha_line - 0.5)
        plt.plot(binned_data[:, 0], binned_data[:, 6], linewidth=linewidth, color=color_norm, alpha=alpha_line - 0.5)

        # remove unwanted axis ticks on top and right
        ax.get_xaxis().tick_bottom()
        ax.get_yaxis().tick_left()

        # ax.xaxis.set_label_coords(0.45, -0.1)  # move x-axis label towards graph

        ax.legend(['AAIMON', 'AAIMON norm.'], loc='upper right', fontsize=fontsize)  # create legend

        utils.save_figure(fig, Fig_name, base_filepath, save_png, save_pdf, dpi)

    logging.info("~~~~~~~~~~~~        run_save_figures_describing_proteins_in_list is finished        ~~~~~~~~~~~~")