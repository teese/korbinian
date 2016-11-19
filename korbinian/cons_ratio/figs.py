from scipy.stats import ttest_ind
import ast
import csv
import itertools
import korbinian.utils as utils
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys

def save_figures_describing_proteins_in_list(pathdict, s, logging):
    logging.info("~~~~~~~~~~~~         starting run_save_figures_describing_proteins_in_list          ~~~~~~~~~~~~")

    backgroundcolour = '0.95'
    plt.style.use('ggplot')

    # set default font size for plot
    fontsize = 8
    datapointsize = 8
    alpha = 0.1

    # set resolution for plots in png format
    dpi = 300

    '''Prepare data for the following plots'''

    # open cons_ratio summary file
    df = pd.read_csv(pathdict["list_cr_summary_csv"], sep=",", quoting=csv.QUOTE_NONNUMERIC, index_col=0)

    # create binlist
    linspace_binlist = np.linspace(s["mp_smallest_bin"],
                                   s["mp_largest_bin"],
                                   s["mp_number_of_bins"])

    # add 30 as the last bin, to make sure 100% of the data is added to the histogram, including major outliers
    binlist = np.append(linspace_binlist, s["mp_final_highest_bin"])

    # iterate through the proteins that have a list of TMDs
    for acc in df.loc[df['list_of_TMDs'].notnull()].loc[df['list_of_TMDs'] != 'nan'].index:
        dict_AAIMON_ratio_mean = {}
        dict_AAIMON_ratio_std = {}
        dict_AASMON_ratio_mean = {}
        dict_AASMON_ratio_std = {}
        for TMD in ast.literal_eval(df.loc[acc, 'list_of_TMDs']):
            dict_AAIMON_ratio_mean[TMD] = df.loc[acc, '%s_AAIMON_ratio_mean' % TMD]
            dict_AAIMON_ratio_std[TMD] = df.loc[acc, '%s_AAIMON_ratio_std' % TMD]
            dict_AASMON_ratio_mean[TMD] = df.loc[acc, '%s_AASMON_ratio_mean' % TMD]
            dict_AASMON_ratio_std[TMD] = df.loc[acc, '%s_AASMON_ratio_std' % TMD]
        df.loc[acc, 'AAIMON_ratio_mean_all_TMDs'] = np.mean(list(dict_AAIMON_ratio_mean.values()))
        df.loc[acc, 'AAIMON_ratio_std_all_TMDs'] = np.mean(list(dict_AAIMON_ratio_std.values()))
        df.loc[acc, 'AASMON_ratio_mean_all_TMDs'] = np.mean(list(dict_AASMON_ratio_mean.values()))
        df.loc[acc, 'AASMON_ratio_std_all_TMDs'] = np.mean(list(dict_AASMON_ratio_std.values()))

    # logging saved data types
    sys.stdout.write('saving figures as: ')
    if s['save_fig_to_pdf']:
        sys.stdout.write(' .pdf ')
    if s['save_fig_to_png']:
        sys.stdout.write(' .png ')
    sys.stdout.write('\n')

    # create list of colours to use in figures
    colour_lists = utils.create_colour_lists()
    tableau20 = colour_lists['tableau20']

    if s['Fig01_Histogram_of_mean_AAIMON_and_AASMON_ratios_SP_vs_MP']:
        Fig_Nr = 1

        sys.stdout.write('Figure Processed: Fig01_Histogram_of_mean_AAIMON_and_AASMON_ratios_SP_vs_MP \n')
        title = 'Mean ratios'

        # for the first figure, create empty objects that will be overwritten into matplotlib objects.
        #fig, ax = 'empty', 'objects'
        # create a new figure
        fig, ax = plt.subplots()
        # create numpy array of membranous over nonmembranous conservation ratios (identity)
        hist_data_AAIMON_mean = np.array(df['AAIMON_ratio_mean_all_TMDs'].dropna())
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
        legend_obj = ax.legend(['AASMON (identity + similarity)', 'AAIMON (identity)'], loc='upper right',
                               fontsize=fontsize)
        # add figure number to top left of subplot
        ax.annotate(s=str(Fig_Nr) + '.', xy=(0.04, 0.9), fontsize=8, xytext=None, xycoords='axes fraction',
                    alpha=0.75)
        # add figure title to top left of subplot
        ax.annotate(s=title, xy=(0.1, 0.9), fontsize=5, xytext=None, xycoords='axes fraction', alpha=0.75)
        # improve ggplot style for a canvas (fig) with 4 figures (plots)
        # FROM RIMMA SCRIPT
        ax.yaxis.grid(True, zorder=0, linestyle=":", color="grey")
        ax
        for tic in ax.xaxis.get_major_ticks():
            tic.tick1On = False
        for tic in ax.yaxis.get_major_ticks():
            tic.tick1On = False
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        # END FROM RIMMA SCRIPT

        utils.save_figure(s, fig, Fig_Nr, base_filepath=pathdict["figures_describing_proteins_in_list"], dpi=dpi)

    if s['Fig02_Histogram_of_standard_deviations_for_AAIMON_and_AASMON']:
        Fig_Nr = 2
        sys.stdout.write('Figure Processed: Fig02_Histogram_of_standard_deviations_for_AAIMON_and_AASMON \n')
        title = 'Standard Deviaton, SP vs MP'
        # create a new figure
        fig, ax = plt.subplots()
        
        # create numpy array of membranous over nonmembranous conservation ratios (identity)
        hist_data_AAIMON_std = np.array(df['AAIMON_ratio_std_all_TMDs'].dropna())
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
        ax.annotate(s=str(Fig_Nr) + '.', xy=(0.04, 0.9), fontsize=5, xytext=None,
                                       xycoords='axes fraction',
                                       alpha=0.75)
        # add figure title to top left of subplot
        ax.annotate(s=title, xy=(0.1, 0.9), fontsize=5, xytext=None, xycoords='axes fraction',
                                       alpha=0.75)

        utils.save_figure(s, fig, Fig_Nr, base_filepath=pathdict["figures_describing_proteins_in_list"], dpi=dpi)


    if s['Fig03_Scattergram_comparing_mean_AAIMON_and_AASMON']:
        Fig_Nr = 3
        sys.stdout.write('Figure Processed: Fig03_Scattergram_comparing_mean_AAIMON_and_AASMON \n')
        title = 'AAIMON vs AASMON'

        fig, ax = plt.subplots()
        # pylab.rcParams['figure.figsize'] = (50.0, 40.0)
        x = np.array(df['AAIMON_ratio_mean_all_TMDs'])
        y = np.array(df['AASMON_ratio_mean_all_TMDs'])
        scattercontainer_AAIMON_AASMON_mean = ax.scatter(x=x, y=y, color="#8A084B", alpha=alpha,
                                                                            s=datapointsize)
        # label the x-axis for each plot, based on the TMD
        ax.set_xlabel('AAIMON_ratio (aa identity)', fontsize=fontsize)
        # move the x-axis label closer to the x-axis
        ax.xaxis.set_label_coords(0.45, -0.085)
        ax.set_ylabel('AASMON_ratio (aa identity + similarity)', fontsize=fontsize)
        # change axis font size
        ax.tick_params(labelsize=fontsize)
        # create legend?#http://stackoverflow.com/questions/9834452/how-do-i-make-a-single-legend-for-many-subplots-with-matplotlib
        legend_obj = ax.legend(['mean'], loc='upper right', fontsize=fontsize)
        # add figure number to top left of subplot
        ax.annotate(s=str(Fig_Nr) + '.', xy=(0.04, 0.9), fontsize=5, xytext=None,
                                       xycoords='axes fraction',
                                       alpha=0.75)
        # add figure title to top left of subplot
        ax.annotate(s=title, xy=(0.1, 0.9), fontsize=5, xytext=None, xycoords='axes fraction',
                                       alpha=0.75)

        utils.save_figure(s, fig, Fig_Nr, base_filepath=pathdict["figures_describing_proteins_in_list"], dpi=dpi)

    if s['Fig04_Scattergram_comparing_standard_deviation_AAIMON_and_AASMON']:
        Fig_Nr = 4
        sys.stdout.write('Figure Processed: Fig04_Scattergram_comparing_standard_deviation_AAIMON_and_AASMON \n')
        title = 'standard deviation AAIMON vs AASMON'
        fig, ax = plt.subplots()

        x = np.array(df['AAIMON_ratio_std_all_TMDs'])
        y = np.array(df['AASMON_ratio_std_all_TMDs'])
        scattercontainer_AAIMON_AASMON_std = ax.scatter(x=x, y=y, color="#B45F04", alpha=alpha,
                                                                           s=datapointsize)
        # other colours that are compatible with colourblind readers: #8A084B Dark red, #B45F04 deep orange, reddish purple #4B088A
        # http://html-color-codes.info/
        # label the x-axis for each plot, based on the TMD
        ax.set_xlabel('AAIMON_ratio', fontsize=fontsize)
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
        ax.annotate(s=str(Fig_Nr) + '.', xy=(0.04, 0.9), fontsize=5, xytext=None,
                                       xycoords='axes fraction',
                                       alpha=0.75)
        # add figure title to top left of subplot
        ax.annotate(s=title, xy=(0.1, 0.9), fontsize=5, xytext=None, xycoords='axes fraction',
                                       alpha=0.75)

        utils.save_figure(s, fig, Fig_Nr, base_filepath=pathdict["figures_describing_proteins_in_list"], dpi=dpi)

    if s['Fig05_Scattergram_comparing_number_of_TMDs_with_mean_AAIMON']:
        Fig_Nr = 5
        sys.stdout.write('Figure Processed: Fig05_Scattergram_comparing_number_of_TMDs_with_mean_AAIMON \n')
        title = 'num_TMDs vs AAIMON'
        fig, ax = plt.subplots()
        

        # for backwards compatibility, check for old name.
        # commented out by MO - assumed that it is not needed
        # if 'number_of_TMDs' in df.columns:
        #     x = np.array(df['number_of_TMDs'])
        # else:
        #     x = np.array(df['number_of_TMDs'])
        x = np.array(df['number_of_TMDs'])
        y = np.array(df['AAIMON_ratio_mean_all_TMDs'])
        scattercontainer_AAIMON_AASMON_std = ax.scatter(x=x, y=y, color="#0489B1", alpha=alpha,
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
        legend_obj = ax.legend(['AAIMON_ratio_mean'], loc='upper right', fontsize=fontsize)
        # add background grid
        ax.grid(True, color='0.75', alpha=0.3)
        # add figure number to top left of subplot
        ax.annotate(s=str(Fig_Nr) + '.', xy=(0.04, 0.9), fontsize=5, xytext=None,
                                       xycoords='axes fraction',
                                       alpha=0.75)
        # add figure title to top left of subplot
        ax.annotate(s=title, xy=(0.1, 0.9), fontsize=5, xytext=None, xycoords='axes fraction',
                                       alpha=0.75)

        utils.save_figure(s, fig, Fig_Nr, base_filepath=pathdict["figures_describing_proteins_in_list"], dpi=dpi)

    if s['Fig06_Scattergram_comparing_seqlen_with_mean_AAIMON']:
        Fig_Nr = 6
        sys.stdout.write('Figure Processed: Fig06_Scattergram_comparing_seqlen_with_mean_AAIMON \n')
        title = 'seqlen vs AAIMON'
        fig, ax = plt.subplots()
        
        # pylab.rcParams['figure.figsize'] = (100.0, 80.0)
        x = np.array(df['seqlen'])
        y = np.array(df['AAIMON_ratio_mean_all_TMDs'])
        scattercontainer_AAIMON_AASMON_std = ax.scatter(x=x, y=y, color="#0489B1", alpha=alpha,
                                                                           s=datapointsize)
        # label the x-axis for each plot, based on the TMD
        ax.set_xlabel('Length of protein', fontsize=fontsize)
        # move the x-axis label closer to the x-axis
        ax.xaxis.set_label_coords(0.45, -0.085)
        ax.set_ylabel('AAIMON_ratio', fontsize=fontsize)
        # change axis font size
        ax.tick_params(labelsize=fontsize)
        # add figure number to top left of subplot
        ax.annotate(s=str(Fig_Nr) + '.', xy=(0.04, 0.9), fontsize=5, xytext=None,
                                       xycoords='axes fraction',
                                       alpha=0.75)
        # add figure title to top left of subplot
        ax.annotate(s=title, xy=(0.1, 0.9), fontsize=5, xytext=None, xycoords='axes fraction',
                                       alpha=0.75)

        utils.save_figure(s, fig, Fig_Nr, base_filepath=pathdict["figures_describing_proteins_in_list"], dpi=dpi)

    if s['Fig07_Scattergram_comparing_nonTMD_SW_align_len_mean_with_mean_AAIMON']:
        Fig_Nr = 7
        sys.stdout.write('Figure Processed: Fig07_Scattergram_comparing_nonTMD_SW_align_len_mean_with_mean_AAIMON \n')
        title = 'length nonTMD region'
        fig, ax = plt.subplots()
        
        # pylab.rcParams['figure.figsize'] = (100.0, 80.0)
        x = np.array(df['nonTMD_SW_align_len_mean'])
        y = np.array(df['AAIMON_ratio_mean_all_TMDs'])
        scattercontainer_AAIMON_AASMON_std = ax.scatter(x=x, y=y, color="#0489B1", alpha=alpha,
                                                                           s=datapointsize)
        # label the x-axis for each plot, based on the TMD
        ax.set_xlabel('Average length of nonTMD region in homologues', fontsize=fontsize)
        # move the x-axis label closer to the x-axis
        ax.xaxis.set_label_coords(0.45, -0.085)
        ax.set_ylabel('AAIMON_ratio', fontsize=fontsize)
        # change axis font size
        ax.tick_params(labelsize=fontsize)
        # add figure number to top left of subplot
        ax.annotate(s=str(Fig_Nr) + '.', xy=(0.04, 0.9), fontsize=5, xytext=None,
                                       xycoords='axes fraction',
                                       alpha=0.75)
        # add figure title to top left of subplot
        ax.annotate(s=title, xy=(0.1, 0.9), fontsize=5, xytext=None, xycoords='axes fraction',
                                       alpha=0.75)

        utils.save_figure(s, fig, Fig_Nr, base_filepath=pathdict["figures_describing_proteins_in_list"], dpi=dpi)

    if s['Fig08_Scattergram_comparing_total_number_of_simap_hits_with_mean_AAIMON']:
        '''
        note that the total hits comes from SIMAP, so this doesn't say anything about how much data is available for each protein
        '''
        Fig_Nr = 8
        sys.stdout.write('Figure Processed: Fig08_Scattergram_comparing_total_number_of_simap_hits_with_mean_AAIMON \n')
        title = 'number SIMAP hits'
        fig, ax = plt.subplots()
        
        # pylab.rcParams['figure.figsize'] = (100.0, 80.0)
        x = np.array(df['total_number_of_simap_hits'])
        y = np.array(df['AAIMON_ratio_mean_all_TMDs'])
        scattercontainer_AAIMON_AASMON_std = ax.scatter(x=x, y=y, color="#0489B1", alpha=alpha,
                                                                           s=datapointsize)
        # label the x-axis for each plot, based on the TMD
        ax.set_xlabel('total number of homologues', fontsize=fontsize)
        # move the x-axis label closer to the x-axis
        ax.xaxis.set_label_coords(0.45, -0.085)
        ax.set_ylabel('AAIMON_ratio', fontsize=fontsize)
        # change axis font size
        ax.tick_params(labelsize=fontsize)
        # add figure number to top left of subplot
        ax.annotate(s=str(Fig_Nr) + '.', xy=(0.04, 0.9), fontsize=5, xytext=None,
                                       xycoords='axes fraction',
                                       alpha=0.75)
        # add figure title to top left of subplot
        ax.annotate(s=title, xy=(0.1, 0.9), fontsize=5, xytext=None, xycoords='axes fraction',
                                       alpha=0.75)
        utils.save_figure(s, fig, Fig_Nr, base_filepath=pathdict["figures_describing_proteins_in_list"], dpi=dpi)

    if s['Fig09_Histogram_of_mean_AAIMON_ratios_for_each_TMD_separately']:
        Fig_Nr = 9
        sys.stdout.write('Figure Processed: Fig09_Histogram_of_mean_AAIMON_ratios_for_each_TMD_separately \n')
        title = 'Histogram of mean AAIMON ratios'
        fig, ax = plt.subplots()

        title = 'AAIMON each TMD separately'
        num_bins = 30
        # "#0489B1"
        alpha = 0.25
        col_width_value = 0.95
        ylabel = 'freq'
        xlabel = 'average conservation ratio (membranous over nonmembranous)'

        # create dataframe mean_AAIMON_each_TM
        df_mean_AAIMON_each_TM = pd.DataFrame()
        # add AAIMON each TMD to dataframe
        for acc in df.index:
            for TMD in ast.literal_eval(df.loc[acc, 'list_of_TMDs']):
                df_mean_AAIMON_each_TM.loc[acc, '{a}_AAIMON_ratio_mean'.format(a=TMD)] = df.loc[
                    acc, '{b}_AAIMON_ratio_mean'.format(b=TMD)]

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
        ax.annotate(s=str(Fig_Nr) + '.', xy=(0.04, 0.9), fontsize=5, xytext=None,
                                       xycoords='axes fraction',
                                       alpha=0.75)
        # add figure title to top left of subplot
        ax.annotate(s=title, xy=(0.1, 0.9), fontsize=5, xytext=None, xycoords='axes fraction',
                                       alpha=0.75)

        utils.save_figure(s, fig, Fig_Nr, base_filepath=pathdict["figures_describing_proteins_in_list"], dpi=dpi)

    if s['Fig10_Line_histogram_of_mean_AAIMON_ratios_for_each_TMD_separately']:
        Fig_Nr = 10
        sys.stdout.write('Figure Processed: Fig10_Line_histogram_of_mean_AAIMON_ratios_for_each_TMD_separately \n')
        title = 'Line histogram each TMD'
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
        ax.annotate(s=str(Fig_Nr) + '.', xy=(0.04, 0.9), fontsize=5, xytext=None,
                                       xycoords='axes fraction',
                                       alpha=0.75)
        # add figure title to top left of subplot
        ax.annotate(s=title, xy=(0.1, 0.9), fontsize=5, xytext=None, xycoords='axes fraction',
                                       alpha=0.75)

        utils.save_figure(s, fig, Fig_Nr, base_filepath=pathdict["figures_describing_proteins_in_list"], dpi=dpi)


    logging.info("~~~~~~~~~~~~        run_save_figures_describing_proteins_in_list is finished        ~~~~~~~~~~~~")