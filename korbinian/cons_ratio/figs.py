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
    backgroundcolour = '0.95'
    plt.style.use('ggplot')
    # add non-functional object to aid document navigation in some IDEs (e.g. Spyder)
    fig_title = ''

    '''
    Prepare subplots and default fontsizes etc
    '''
    # set default font size for plot
    fontsize = 4
    datapointsize = 2
    alpha = 0.1
    # use a dictionary to organise the saving of multiple plots in multiple figures, with a certain number of plots per figure
    n_plots_per_fig = 4
    nrows_in_each_fig = 2
    ncols_in_each_fig = 2
    dict_organising_subplots = utils.create_dict_organising_subplots(n_plots_per_fig=n_plots_per_fig,
                                                                     n_rows=nrows_in_each_fig,
                                                                     n_cols=ncols_in_each_fig)
    '''
    Prepare data for following figures
    '''
    df = pd.read_csv(pathdict["list_cr_summary_csv"], sep = ",", quoting = csv.QUOTE_NONNUMERIC, index_col = 0)
    # if "uniprot_acc" in df.columns:
    #     df.set_index("uniprot_acc", drop=False, inplace=True)
    # else:
    #     df["uniprot_acc"] = df.index

    # iterate over the datafra    # filter to remove sequences where no TMDs are found (will contain either np.nan, or 'nan')
    # df = df.loc[df['list_of_TMDs'].notnull()]
    # df = df.loc[df['list_of_TMDs'] != 'nan']me. Note that acc = uniprot accession here.
    linspace_binlist = np.linspace(s["mp_smallest_bin"],
                                   s["mp_largest_bin"],
                                   s["mp_number_of_bins"])

    # add 30 as the last bin, to make sure 100% of the data is added to the histogram, including major outliers
    binlist = np.append(linspace_binlist, s["mp_final_highest_bin"])

    '''
    The beta-barrel dataset contained a lot of proteins with an average AAIMON of 1.000000. This can only mean that there are not enough homologues.
    The singlepass dataset contained only 2 proteins, the alpha-helicas multipass only 8 with 1.000000.
    All these are excluded from the dataset (and all following graphs). Note that it would be better to exactly count the valid homologues, rather than exclude them afterwards like this.
    '''
    # iterate through the proteins that have a list of TMDs
    for acc in df.loc[df['list_of_TMDs'].notnull()].loc[df['list_of_TMDs'] != 'nan'].index:
        dict_AAIMON_ratio_mean = {}
        dict_AAIMON_ratio_std = {}
        dict_AASMON_ratio_mean = {}
        dict_AASMON_ratio_std = {}
        for TMD in ast.literal_eval(df.loc[acc, 'list_of_TMDs']):
            dict_AAIMON_ratio_mean[TMD] = df.loc[acc, '%s_AAIMON_ratio_mean'%TMD]
            dict_AAIMON_ratio_std[TMD] = df.loc[acc, '%s_AAIMON_ratio_std'%TMD]
            dict_AASMON_ratio_mean[TMD] = df.loc[acc, '%s_AASMON_ratio_mean'%TMD]
            dict_AASMON_ratio_std[TMD] = df.loc[acc, '%s_AASMON_ratio_std'%TMD]
        df.loc[acc, 'AAIMON_ratio_mean_all_TMDs'] = np.mean(list(dict_AAIMON_ratio_mean.values()))
        df.loc[acc, 'AAIMON_ratio_std_all_TMDs'] = np.mean(list(dict_AAIMON_ratio_std.values()))
        df.loc[acc, 'AASMON_ratio_mean_all_TMDs'] = np.mean(list(dict_AASMON_ratio_mean.values()))
        df.loc[acc, 'AASMON_ratio_std_all_TMDs'] = np.mean(list(dict_AASMON_ratio_std.values()))

    # count numbers of each AAIMON ratio
    vc_AAIMON = df.AAIMON_ratio_mean_all_TMDs.value_counts()
    # replace df with a filtered dataframe, with all rows excluded where AAIMON_ratio_mean_all_TMDs is 1.000000
    if 1.000000 in vc_AAIMON:
        num_proteins_with_AAIMON_of_ONE = vc_AAIMON[1.000000]
        total_num_prot_with_data = len(df['AAIMON_ratio_mean_all_TMDs'].dropna())
        logging.info('num_proteins_with_AAIMON_of_ONE in orig dataframe : %i from %i total' % (
        num_proteins_with_AAIMON_of_ONE, total_num_prot_with_data))
        # replace the original dataframe with the dataframe lacking data containing a TMD/rest ratio of exactly 1.
        df = df.loc[df['AAIMON_ratio_mean_all_TMDs'] != 1.000000]

    '''
    Fig01: Histogram of mean AAIMON and AASMON ratios, SP vs MP
    '''
    Fig_Nr = 1
    # add non-functional "if" line to aid document navigation in some IDEs (e.g. Spyder)
    if fig_title == "Fig01: Histogram of mean AAIMON and AASMON ratios, SP vs MP":
        pass
    sys.stdout.write('Figures Processed' + str(Fig_Nr) + ', ')
    title = 'Mean ratios'
    newfig, savefig, fig_nr, plot_nr_in_fig, row_nr, col_nr = dict_organising_subplots[Fig_Nr]
    # for the first figure, create empty objects that will be overwritten into matplotlib objects.
    fig, axarr = 'empty', 'objects'
    # create a new figure
    fig, axarr = utils.create_new_fig_if_necessary(newfig, fig, axarr, nrows_in_each_fig, ncols_in_each_fig, dpi=300)
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
    barcontainer_AAIMON_mean = axarr[row_nr, col_nr].bar(left=centre_of_bar_in_x_axis, height=freq_counts_I,
                                                         align='center', width=col_width, color="#0489B1",
                                                         alpha=0.5)  # edgecolor='black',
    # create numpy array of membranous over nonmembranous conservation ratios (identity + similarity)
    hist_data_AASMON_mean = np.array(df['AASMON_ratio_mean_all_TMDs'].dropna())
    # use numpy to create a histogram
    freq_counts_S, bin_array_S = np.histogram(hist_data_AASMON_mean, bins=binlist)
    # barcontainer_S = axarr[row_nr,col_nr].bar(left=centre_of_bar_in_x_axis, height=freq_counts_S, align='center', width=col_width, color="#0101DF", edgecolor="#0101DF", alpha = 0.5)
    # create a line graph rather than a bar graph for the AASMON (ident + similarity)
    linecontainer_AASMON_mean = axarr[row_nr, col_nr].plot(centre_of_bar_in_x_axis, freq_counts_S, color="#0101DF",
                                                           alpha=0.5)
    # other colours that are compatible with colourblind readers: #8A084B Dark red, #B45F04 deep orange, reddish purple #4B088A
    # http://html-color-codes.info/
    # label the x-axis for each plot, based on the TMD
    axarr[row_nr, col_nr].set_xlabel('average conservation ratio (membranous over nonmembranous)', fontsize=fontsize)
    # move the x-axis label closer to the x-axis
    axarr[row_nr, col_nr].xaxis.set_label_coords(0.45, -0.085)
    # pylab.rcParams['figure.figsize'] = (50.0, 40.0)
    # pylab.rcParams['figure.figsize'] = (20.0, 16.0)
    # plt.show()
    xlim_min = s["mp_xlim_min01"]
    # take x-axis max from settings
    xlim_max = s["mp_xlim_max01"]
    # set x-axis min
    axarr[row_nr, col_nr].set_xlim(xlim_min, xlim_max)
    # set x-axis ticks
    # use the slide selection to select every second item in the list as an xtick(axis label)
    axarr[row_nr, col_nr].set_xticks([float('%0.1f' % c) for c in centre_of_bar_in_x_axis[::3]])
    axarr[row_nr, col_nr].set_ylabel('freq', rotation='vertical', fontsize=fontsize)
    # change axis font size
    axarr[row_nr, col_nr].tick_params(labelsize=fontsize)
    # create legend?#http://stackoverflow.com/questions/9834452/how-do-i-make-a-single-legend-for-many-subplots-with-matplotlib
    legend_obj = axarr[row_nr, col_nr].legend(['AASMON (identity + similarity)', 'AAIMON (identity)'], loc='upper right',
                                              fontsize=fontsize)
    # add figure number to top left of subplot
    axarr[row_nr, col_nr].annotate(s=str(Fig_Nr) + '.', xy=(0.04, 0.9), fontsize=5, xytext=None, xycoords='axes fraction',
                                   alpha=0.75)
    # add figure title to top left of subplot
    axarr[row_nr, col_nr].annotate(s=title, xy=(0.1, 0.9), fontsize=5, xytext=None, xycoords='axes fraction', alpha=0.75)
    # improve ggplot style for a canvas (fig) with 4 figures (plots)
    utils.improve_ggplot_for_4_plots(axarr, row_nr, col_nr, backgroundcolour, legend_obj)
    # FROM RIMMA SCRIPT
    axarr[row_nr, col_nr].yaxis.grid(True, zorder=0, linestyle=":", color="grey")
    axarr[row_nr, col_nr]
    for tic in axarr[row_nr, col_nr].xaxis.get_major_ticks():
        tic.tick1On = False
    for tic in axarr[row_nr, col_nr].yaxis.get_major_ticks():
        tic.tick1On = False
    axarr[row_nr, col_nr].spines['top'].set_visible(False)
    axarr[row_nr, col_nr].spines['right'].set_visible(False)
    # END FROM RIMMA SCRIPT

    '''
    Fig02: Histogram of standard deviations for AAIMON and AASMON (std among homologues for each protein)
    '''
    Fig_Nr = 2
    # add non-functional "if" line to aid document navigation in some IDEs (e.g. Spyder)
    if fig_title == "Fig02: Histogram of standard deviations for AAIMON and AASMON":
        pass
    sys.stdout.write(str(Fig_Nr) + ', ')
    title = 'Standard Deviaton, SP vs MP'
    newfig, savefig, fig_nr, plot_nr_in_fig, row_nr, col_nr = dict_organising_subplots[Fig_Nr]
    # if a new figure should be created (either because the orig is full, or the last TMD is analysed)
    fig, axarr = utils.create_new_fig_if_necessary(newfig, fig, axarr, nrows_in_each_fig, ncols_in_each_fig, dpi=300)
    # create numpy array of membranous over nonmembranous conservation ratios (identity)
    hist_data_AAIMON_std = np.array(df['AAIMON_ratio_std_all_TMDs'].dropna())
    # use numpy to create a histogram
    number_of_bins = 50
    freq_counts_S, bin_array_S = np.histogram(hist_data_AAIMON_std, bins=number_of_bins)
    # assuming all of the bins are exactly the same size, make the width of the column equal to 70% of each bin
    col_width = float('%0.3f' % (0.95 * (bin_array_S[1] - bin_array_S[0])))
    # when align='center', the central point of the bar in the x-axis is simply the middle of the bins ((bin_0-bin_1)/2, etc)
    centre_of_bar_in_x_axis = (bin_array_S[:-1] + bin_array_S[1:]) / 2
    barcontainer_AAIMON_std = axarr[row_nr, col_nr].bar(left=centre_of_bar_in_x_axis, height=freq_counts_S,
                                                        align='center', width=col_width, color="#0489B1",
                                                        alpha=0.5)  # edgecolor='black',
    # create numpy array of membranous over nonmembranous conservation ratios (identity + similarity)
    hist_data_AASMON_std = np.array(df['AASMON_ratio_std_all_TMDs'].dropna())
    # use numpy to create a histogram
    # N.B. use the bins from the previous plot
    freq_counts, bin_array = np.histogram(hist_data_AASMON_std, bins=bin_array_S)
    # barcontainer_S = axarr[row_nr,col_nr].bar(left=centre_of_bar_in_x_axis, height=freq_counts, align='center', width=col_width, color="#0101DF", edgecolor="#0101DF", alpha = 0.5)
    # create a line graph rather than a bar graph for the AASMON (ident + similarity)
    linecontainer_AASMON_std = axarr[row_nr, col_nr].plot(centre_of_bar_in_x_axis, freq_counts, color="#0101DF",
                                                          alpha=0.5)
    # other colours that are compatible with colourblind readers: #8A084B Dark red, #B45F04 deep orange, reddish purple #4B088A
    # http://html-color-codes.info/
    # label the x-axis for each plot, based on the TMD
    axarr[row_nr, col_nr].set_xlabel('average standard deviation', fontsize=fontsize)
    # move the x-axis label closer to the x-axis
    axarr[row_nr, col_nr].xaxis.set_label_coords(0.45, -0.085)
    # xlim_min = s["mp_xlim_min01"]
    # take x-axis max from settings
    # xlim_max = s["mp_xlim_max01"]
    # set x-axis min
    # axarr[row_nr,col_nr].set_xlim(xlim_min,xlim_max)
    # set x-axis ticks
    # use the slide selection to select every second item in the list as an xtick(axis label)
    axarr[row_nr, col_nr].set_xticks([float('%0.1f' % c) for c in centre_of_bar_in_x_axis[::3]])
    axarr[row_nr, col_nr].set_ylabel('freq', rotation='vertical', fontsize=fontsize)
    # change axis font size
    axarr[row_nr, col_nr].tick_params(labelsize=fontsize)
    # create legend?#http://stackoverflow.com/questions/9834452/how-do-i-make-a-single-legend-for-many-subplots-with-matplotlib
    legend_obj = axarr[row_nr, col_nr].legend(['AASMON (identity + similarity)', 'AAIMON (identity)'], loc='upper right',
                                              fontsize=fontsize)
    # add figure number to top left of subplot
    axarr[row_nr, col_nr].annotate(s=str(Fig_Nr) + '.', xy=(0.04, 0.9), fontsize=5, xytext=None, xycoords='axes fraction',
                                   alpha=0.75)
    # add figure title to top left of subplot
    axarr[row_nr, col_nr].annotate(s=title, xy=(0.1, 0.9), fontsize=5, xytext=None, xycoords='axes fraction', alpha=0.75)
    # improve ggplot style for a canvas (fig) with 4 figures (plots)
    utils.improve_ggplot_for_4_plots(axarr, row_nr, col_nr, backgroundcolour, legend_obj)

    '''
    Fig03: Scattergram comparing mean AAIMON and AASMON
    '''
    Fig_Nr = 3
    # add non-functional "if" line to aid document navigation in some IDEs (e.g. Spyder)
    if fig_title == "Fig03: Scattergram comparing mean AAIMON and AASMON":
        pass
    sys.stdout.write(str(Fig_Nr) + ', ')
    title = 'AAIMON vs AASMON'
    newfig, savefig, fig_nr, plot_nr_in_fig, row_nr, col_nr = dict_organising_subplots[Fig_Nr]
    # if a new figure should be created (either because the orig is full, or the last TMD is analysed)
    fig, axarr = utils.create_new_fig_if_necessary(newfig, fig, axarr, nrows_in_each_fig, ncols_in_each_fig, dpi=300)
    # pylab.rcParams['figure.figsize'] = (50.0, 40.0)
    x = np.array(df['AAIMON_ratio_mean_all_TMDs'])
    y = np.array(df['AASMON_ratio_mean_all_TMDs'])
    scattercontainer_AAIMON_AASMON_mean = axarr[row_nr, col_nr].scatter(x=x, y=y, color="#8A084B", alpha=alpha,
                                                                        s=datapointsize)
    # label the x-axis for each plot, based on the TMD
    axarr[row_nr, col_nr].set_xlabel('AAIMON_ratio (aa identity)', fontsize=fontsize)
    # move the x-axis label closer to the x-axis
    axarr[row_nr, col_nr].xaxis.set_label_coords(0.45, -0.085)
    axarr[row_nr, col_nr].set_ylabel('AASMON_ratio (aa identity + similarity)', fontsize=fontsize)
    # change axis font size
    axarr[row_nr, col_nr].tick_params(labelsize=fontsize)
    # create legend?#http://stackoverflow.com/questions/9834452/how-do-i-make-a-single-legend-for-many-subplots-with-matplotlib
    legend_obj = axarr[row_nr, col_nr].legend(['mean'], loc='upper right', fontsize=fontsize)
    # add figure number to top left of subplot
    axarr[row_nr, col_nr].annotate(s=str(Fig_Nr) + '.', xy=(0.04, 0.9), fontsize=5, xytext=None, xycoords='axes fraction',
                                   alpha=0.75)
    # add figure title to top left of subplot
    axarr[row_nr, col_nr].annotate(s=title, xy=(0.1, 0.9), fontsize=5, xytext=None, xycoords='axes fraction', alpha=0.75)
    # improve ggplot style for a canvas (fig) with 4 figures (plots)
    utils.improve_ggplot_for_4_plots(axarr, row_nr, col_nr, backgroundcolour, legend_obj)

    '''
    Fig04: Scattergram comparing standard deviation AAIMON and AASMON
    '''
    Fig_Nr = 4
    # add non-functional "if" line to aid document navigation in some IDEs (e.g. Spyder)
    if fig_title == "Fig04: Scattergram comparing standard deviation AAIMON and AASMON":
        pass
    title = 'standard deviation AAIMON vs AASMON'
    newfig, savefig, fig_nr, plot_nr_in_fig, row_nr, col_nr = dict_organising_subplots[Fig_Nr]
    # if a new figure should be created (either because the orig is full, or the last TMD is analysed)
    fig, axarr = utils.create_new_fig_if_necessary(newfig, fig, axarr, nrows_in_each_fig, ncols_in_each_fig, dpi=300)
    x = np.array(df['AAIMON_ratio_std_all_TMDs'])
    y = np.array(df['AASMON_ratio_std_all_TMDs'])
    scattercontainer_AAIMON_AASMON_std = axarr[row_nr, col_nr].scatter(x=x, y=y, color="#B45F04", alpha=alpha,
                                                                       s=datapointsize)
    # other colours that are compatible with colourblind readers: #8A084B Dark red, #B45F04 deep orange, reddish purple #4B088A
    # http://html-color-codes.info/
    # label the x-axis for each plot, based on the TMD
    axarr[row_nr, col_nr].set_xlabel('AAIMON_ratio', fontsize=fontsize)
    # move the x-axis label closer to the x-axis
    axarr[row_nr, col_nr].xaxis.set_label_coords(0.45, -0.085)
    axarr[row_nr, col_nr].set_ylabel('AASMON_ratio', fontsize=fontsize)
    # set x-axis ticks
    # use the slide selection to select every second item in the list as an xtick(axis label)
    # axarr[row_nr,col_nr].set_xticks([float('%0.1f' % c) for c in centre_of_bar_in_x_axis[::3]])
    # axarr[row_nr,col_nr].set_ylabel('freq',rotation = 'vertical', fontsize = fontsize)
    # change axis font size
    axarr[row_nr, col_nr].tick_params(labelsize=fontsize)
    # create legend?#http://stackoverflow.com/questions/9834452/how-do-i-make-a-single-legend-for-many-subplots-with-matplotlib
    legend_obj = axarr[row_nr, col_nr].legend(['standard deviation'], loc='upper right', fontsize=fontsize)
    # add figure number to top left of subplot
    axarr[row_nr, col_nr].annotate(s=str(Fig_Nr) + '.', xy=(0.04, 0.9), fontsize=5, xytext=None, xycoords='axes fraction',
                                   alpha=0.75)
    # add figure title to top left of subplot
    axarr[row_nr, col_nr].annotate(s=title, xy=(0.1, 0.9), fontsize=5, xytext=None, xycoords='axes fraction', alpha=0.75)
    # improve ggplot style for a canvas (fig) with 4 figures (plots)
    utils.improve_ggplot_for_4_plots(axarr, row_nr, col_nr, backgroundcolour, legend_obj)
    # save the figure if necessary (i.e., if the maximum number of plots per figure has been obtained)
    utils.savefig_if_necessary(savefig, fig, fig_nr, base_filepath=pathdict["base_filename_summaries"])

    '''
    Fig05: Scattergram comparing number_of_TMDs with mean AAIMON
    '''
    Fig_Nr = 5
    sys.stdout.write(str(Fig_Nr) + ', ')
    title = 'num_TMDs vs AAIMON'
    newfig, savefig, fig_nr, plot_nr_in_fig, row_nr, col_nr = dict_organising_subplots[Fig_Nr]
    # if a new figure should be created (either because the orig is full, or the last TMD is analysed)
    fig, axarr = utils.create_new_fig_if_necessary(newfig, fig, axarr, nrows_in_each_fig, ncols_in_each_fig, dpi=300)
    # for backwards compatibility, check for old name.
    # commented out by MO - assumed that it is not needed
    # if 'number_of_TMDs' in df.columns:
    #     x = np.array(df['number_of_TMDs'])
    # else:
    #     x = np.array(df['number_of_TMDs'])
    x = np.array(df['number_of_TMDs'])
    y = np.array(df['AAIMON_ratio_mean_all_TMDs'])
    scattercontainer_AAIMON_AASMON_std = axarr[row_nr, col_nr].scatter(x=x, y=y, color="#0489B1", alpha=alpha,
                                                                       s=datapointsize)
    # other colours that are compatible with colourblind readers: #8A084B Dark red, #B45F04 deep orange, reddish purple #4B088A
    # http://html-color-codes.info/
    # label the x-axis for each plot, based on the TMD
    axarr[row_nr, col_nr].set_xlabel('number of TMDs in protein', fontsize=fontsize)
    # move the x-axis label closer to the x-axis
    axarr[row_nr, col_nr].xaxis.set_label_coords(0.45, -0.085)
    axarr[row_nr, col_nr].set_ylabel('Average AAIMON ratio for all TMDs', fontsize=fontsize)
    # change axis font size
    axarr[row_nr, col_nr].tick_params(labelsize=fontsize)
    # create legend?#http://stackoverflow.com/questions/9834452/how-do-i-make-a-single-legend-for-many-subplots-with-matplotlib
    legend_obj = axarr[row_nr, col_nr].legend(['AAIMON_ratio_mean'], loc='upper right', fontsize=fontsize)
    # add background grid
    axarr[row_nr, col_nr].grid(True, color='0.75', alpha=0.3)
    # add figure number to top left of subplot
    axarr[row_nr, col_nr].annotate(s=str(Fig_Nr) + '.', xy=(0.04, 0.9), fontsize=5, xytext=None, xycoords='axes fraction',
                                   alpha=0.75)
    # add figure title to top left of subplot
    axarr[row_nr, col_nr].annotate(s=title, xy=(0.1, 0.9), fontsize=5, xytext=None, xycoords='axes fraction', alpha=0.75)
    # improve ggplot style for a canvas (fig) with 4 figures (plots)
    utils.improve_ggplot_for_4_plots(axarr, row_nr, col_nr, backgroundcolour, legend_obj)

    '''
    Fig06: Scattergram comparing seqlen with mean AAIMON
    '''
    Fig_Nr = 6
    # add non-functional "if" line to aid document navigation in some IDEs (e.g. Spyder)
    if fig_title == "Fig06: Scattergram comparing seqlen with mean AAIMON":
        pass
    sys.stdout.write(str(Fig_Nr) + ', ')
    title = 'seqlen vs AAIMON'
    newfig, savefig, fig_nr, plot_nr_in_fig, row_nr, col_nr = dict_organising_subplots[Fig_Nr]
    # if a new figure should be created (either because the orig is full, or the last TMD is analysed)
    fig, axarr = utils.create_new_fig_if_necessary(newfig, fig, axarr, nrows_in_each_fig, ncols_in_each_fig, dpi=300)
    # pylab.rcParams['figure.figsize'] = (100.0, 80.0)
    x = np.array(df['seqlen'])
    y = np.array(df['AAIMON_ratio_mean_all_TMDs'])
    scattercontainer_AAIMON_AASMON_std = axarr[row_nr, col_nr].scatter(x=x, y=y, color="#0489B1", alpha=alpha,
                                                                       s=datapointsize)
    # label the x-axis for each plot, based on the TMD
    axarr[row_nr, col_nr].set_xlabel('Length of protein', fontsize=fontsize)
    # move the x-axis label closer to the x-axis
    axarr[row_nr, col_nr].xaxis.set_label_coords(0.45, -0.085)
    axarr[row_nr, col_nr].set_ylabel('AAIMON_ratio', fontsize=fontsize)
    # change axis font size
    axarr[row_nr, col_nr].tick_params(labelsize=fontsize)
    # add figure number to top left of subplot
    axarr[row_nr, col_nr].annotate(s=str(Fig_Nr) + '.', xy=(0.04, 0.9), fontsize=5, xytext=None, xycoords='axes fraction',
                                   alpha=0.75)
    # add figure title to top left of subplot
    axarr[row_nr, col_nr].annotate(s=title, xy=(0.1, 0.9), fontsize=5, xytext=None, xycoords='axes fraction', alpha=0.75)
    # improve ggplot style for a canvas (fig) with 4 figures (plots)
    utils.improve_ggplot_for_4_plots(axarr, row_nr, col_nr, backgroundcolour, legend_obj)

    '''
    Fig07: Scattergram comparing nonTMD_SW_align_len_mean with mean AAIMON
    '''
    Fig_Nr = 7
    sys.stdout.write(str(Fig_Nr) + ', ')
    title = 'length nonTMD region'
    newfig, savefig, fig_nr, plot_nr_in_fig, row_nr, col_nr = dict_organising_subplots[Fig_Nr]
    # if a new figure should be created (either because the orig is full, or the last TMD is analysed)
    fig, axarr = utils.create_new_fig_if_necessary(newfig, fig, axarr, nrows_in_each_fig, ncols_in_each_fig, dpi=300)
    # pylab.rcParams['figure.figsize'] = (100.0, 80.0)
    x = np.array(df['nonTMD_SW_align_len_mean'])
    y = np.array(df['AAIMON_ratio_mean_all_TMDs'])
    scattercontainer_AAIMON_AASMON_std = axarr[row_nr, col_nr].scatter(x=x, y=y, color="#0489B1", alpha=alpha,
                                                                       s=datapointsize)
    # label the x-axis for each plot, based on the TMD
    axarr[row_nr, col_nr].set_xlabel('Average length of nonTMD region in homologues', fontsize=fontsize)
    # move the x-axis label closer to the x-axis
    axarr[row_nr, col_nr].xaxis.set_label_coords(0.45, -0.085)
    axarr[row_nr, col_nr].set_ylabel('AAIMON_ratio', fontsize=fontsize)
    # change axis font size
    axarr[row_nr, col_nr].tick_params(labelsize=fontsize)
    # add figure number to top left of subplot
    axarr[row_nr, col_nr].annotate(s=str(Fig_Nr) + '.', xy=(0.04, 0.9), fontsize=5, xytext=None, xycoords='axes fraction',
                                   alpha=0.75)
    # add figure title to top left of subplot
    axarr[row_nr, col_nr].annotate(s=title, xy=(0.1, 0.9), fontsize=5, xytext=None, xycoords='axes fraction', alpha=0.75)
    # improve ggplot style for a canvas (fig) with 4 figures (plots)
    utils.improve_ggplot_for_4_plots(axarr, row_nr, col_nr, backgroundcolour, legend_obj)

    '''
    Fig08: Scattergram comparing total_number_of_simap_hits with mean AAIMON
    note that the total hits comes from SIMAP, so this doesn't say anything about how much data is available for each protein
    '''
    Fig_Nr = 8
    # add non-functional "if" line to aid document navigation in some IDEs (e.g. Spyder)
    if fig_title == "Fig08: Scattergram comparing total_number_of_simap_hits with mean AAIMON":
        pass
    sys.stdout.write(str(Fig_Nr) + ', ')
    title = 'number SIMAP hits'
    newfig, savefig, fig_nr, plot_nr_in_fig, row_nr, col_nr = dict_organising_subplots[Fig_Nr]
    # if a new figure should be created (either because the orig is full, or the last TMD is analysed)
    fig, axarr = utils.create_new_fig_if_necessary(newfig, fig, axarr, nrows_in_each_fig, ncols_in_each_fig, dpi=300)
    # pylab.rcParams['figure.figsize'] = (100.0, 80.0)
    # commented out as data for total_number_of_simap_hits is currently not available
    # x = np.array(df['total_number_of_simap_hits'])
    # y = np.array(df['AAIMON_ratio_mean_all_TMDs'])
    scattercontainer_AAIMON_AASMON_std = axarr[row_nr, col_nr].scatter(x=x, y=y, color="#0489B1", alpha=alpha,
                                                                       s=datapointsize)
    # label the x-axis for each plot, based on the TMD
    axarr[row_nr, col_nr].set_xlabel('total number of homologues', fontsize=fontsize)
    # move the x-axis label closer to the x-axis
    axarr[row_nr, col_nr].xaxis.set_label_coords(0.45, -0.085)
    axarr[row_nr, col_nr].set_ylabel('AAIMON_ratio', fontsize=fontsize)
    # change axis font size
    axarr[row_nr, col_nr].tick_params(labelsize=fontsize)
    # add figure number to top left of subplot
    axarr[row_nr, col_nr].annotate(s=str(Fig_Nr) + '.', xy=(0.04, 0.9), fontsize=5, xytext=None, xycoords='axes fraction',
                                   alpha=0.75)
    # add figure title to top left of subplot
    axarr[row_nr, col_nr].annotate(s=title, xy=(0.1, 0.9), fontsize=5, xytext=None, xycoords='axes fraction', alpha=0.75)
    # improve ggplot style for a canvas (fig) with 4 figures (plots)
    utils.improve_ggplot_for_4_plots(axarr, row_nr, col_nr, backgroundcolour, legend_obj)
    # save the figure if necessary (i.e., if the maximum number of plots per figure has been obtained)
    utils.savefig_if_necessary(savefig, fig, fig_nr, base_filepath=pathdict["base_filename_summaries"])

    '''
    SHOW THE DATA FOR EACH TMD IN THE DATASET
    '''
    # create list of colours to use in figures
    colour_lists = utils.create_colour_lists()
    tableau20 = colour_lists['tableau20']

    '''
    Fig09: Histogram of mean AAIMON ratios for each TMD separately
    '''
    Fig_Nr = 9
    sys.stdout.write(str(Fig_Nr) + ', ')
    title = 'Histogram of mean AAIMON ratios'

    if 'number_of_TMDs' not in df.columns:
        df['number_of_TMDs'] = df.list_of_TMDs.apply(lambda list_tmds: len(eval(list_tmds)))
    df_mean_AAIMON_each_TM, max_num_TMDs, legend = utils.create_df_with_mean_AAIMON_each_TM(df)

    newfig, savefig, fig_nr, plot_nr_in_fig, row_nr, col_nr = dict_organising_subplots[Fig_Nr]
    # if a new figure should be created (either because the orig is full, or the last TMD is analysed)
    fig, axarr = utils.create_new_fig_if_necessary(newfig, fig, axarr, nrows_in_each_fig, ncols_in_each_fig, dpi=300)

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
        barcontainer = axarr[row_nr, col_nr].bar(left=centre_of_bar_in_x_axis, height=freq_counts_I,
                                                 align='center', width=col_width, facecolor=color,
                                                 alpha=alpha, edgecolor='black', linewidth=0.2)  # edgecolor='black',
        # barcontainer = axarr[row_nr, col_nr].bar(left=centre_of_bar_in_x_axis, height=freq_counts_I,
        #                                                     align='center', width=col_width, facecolor=color,
        #                                                     alpha=alpha, edgecolor='black', linewidth=0.1)  # edgecolor='black',    #label the x-axis for each plot, based on the TMD
    axarr[row_nr, col_nr].set_xlabel(xlabel, fontsize=fontsize)
    # move the x-axis label closer to the x-axis
    axarr[row_nr, col_nr].xaxis.set_label_coords(0.45, -0.085)
    axarr[row_nr, col_nr].set_ylabel(ylabel, rotation='vertical', fontsize=fontsize)
    # change axis font size
    axarr[row_nr, col_nr].tick_params(labelsize=fontsize)
    # create legend?#http://stackoverflow.com/questions/9834452/how-do-i-make-a-single-legend-for-many-subplots-with-matplotlib
    legend_obj = axarr[row_nr, col_nr].legend(legend, loc='upper right', fontsize=fontsize)
    # add title
    # axarr[row_nr, col_nr].set_title(title,fontsize=fontsize)
    # add figure number to top left of subplot
    axarr[row_nr, col_nr].annotate(s=str(Fig_Nr) + '.', xy=(0.04, 0.9), fontsize=5, xytext=None, xycoords='axes fraction',
                                   alpha=0.75)
    # add figure title to top left of subplot
    axarr[row_nr, col_nr].annotate(s=title, xy=(0.1, 0.9), fontsize=5, xytext=None, xycoords='axes fraction', alpha=0.75)
    # improve ggplot style for a canvas (fig) with 4 figures (plots)
    utils.improve_ggplot_for_4_plots(axarr, row_nr, col_nr, backgroundcolour, legend_obj)

    '''
    Fig10: Line histogram of mean AAIMON ratios for each TMD separately
    '''
    Fig_Nr = 10
    # add non-functional "if" line to aid document navigation in some IDEs (e.g. Spyder)
    if fig_title == "Fig10: Line histogram of mean AAIMON ratios for each TMD separately":
        pass
    sys.stdout.write(str(Fig_Nr) + ', ')
    title = 'Line histogram each TMD'
    newfig, savefig, fig_nr, plot_nr_in_fig, row_nr, col_nr = dict_organising_subplots[Fig_Nr]
    # if a new figure should be created (either because the orig is full, or the last TMD is analysed)
    fig, axarr = utils.create_new_fig_if_necessary(newfig, fig, axarr, nrows_in_each_fig, ncols_in_each_fig, dpi=300)
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
        linecontainer = axarr[row_nr, col_nr].plot(centre_of_bar_in_x_axis, freq_counts_I, alpha=alpha,
                                                   linewidth=0.3)  # edgecolor='black',
    # label the x-axis for each plot, based on the TMD
    axarr[row_nr, col_nr].set_xlabel(xlabel, fontsize=fontsize)
    # move the x-axis label closer to the x-axis
    axarr[row_nr, col_nr].xaxis.set_label_coords(0.45, -0.085)
    axarr[row_nr, col_nr].set_ylabel(ylabel, rotation='vertical', fontsize=fontsize)
    # change axis font size
    axarr[row_nr, col_nr].tick_params(labelsize=fontsize)
    # create legend?#http://stackoverflow.com/questions/9834452/how-do-i-make-a-single-legend-for-many-subplots-with-matplotlib
    legend_obj = axarr[row_nr, col_nr].legend(legend, loc='upper right', fontsize=fontsize)
    # add title
    # axarr[row_nr, col_nr].set_title(title,fontsize=fontsize)
    # add figure number to top left of subplot
    axarr[row_nr, col_nr].annotate(s=str(Fig_Nr) + '.', xy=(0.04, 0.9), fontsize=5, xytext=None, xycoords='axes fraction',
                                   alpha=0.75)
    # add figure title to top left of subplot
    axarr[row_nr, col_nr].annotate(s=title, xy=(0.1, 0.9), fontsize=5, xytext=None, xycoords='axes fraction', alpha=0.75)
    # improve ggplot style for a canvas (fig) with 4 figures (plots)
    utils.improve_ggplot_for_4_plots(axarr, row_nr, col_nr, backgroundcolour, legend_obj)

    # these graphs are only applicable for multi-pass proteins. Use where at least 2 proteins have a 7th TMD
    if "TM07" in df_mean_AAIMON_each_TM.columns:
        if df_mean_AAIMON_each_TM['TM07'].dropna().shape[0] >= 2:
            dataset_contains_multipass_prots = True
        else:
            dataset_contains_multipass_prots = False
    else:
        dataset_contains_multipass_prots = False
    if dataset_contains_multipass_prots:
        '''
        Fig11: Line histogram of mean AAIMON ratios for selected TMDs, highlighting difference for TM07
        '''
        Fig_Nr = 11
        sys.stdout.write(str(Fig_Nr) + ', ')
        title = 'Select TMDs, all data'
        cols_for_analysis = ['TM01', 'TM07', 'TM08', 'last_TM_AAIMON_ratio_mean']
        newfig, savefig, fig_nr, plot_nr_in_fig, row_nr, col_nr = dict_organising_subplots[Fig_Nr]
        # if a new figure should be created (either because the orig is full, or the last TMD is analysed)
        fig, axarr = utils.create_new_fig_if_necessary(newfig, fig, axarr, nrows_in_each_fig, ncols_in_each_fig, dpi=300)
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
            linecontainer = axarr[row_nr, col_nr].plot(centre_of_bar_in_x_axis, freq_counts_I, alpha=alpha,
                                                       linewidth=0.3)  # edgecolor='black',
        # label the x-axis for each plot, based on the TMD
        axarr[row_nr, col_nr].set_xlabel(xlabel, fontsize=fontsize)
        # move the x-axis label closer to the x-axis
        axarr[row_nr, col_nr].xaxis.set_label_coords(0.45, -0.085)
        axarr[row_nr, col_nr].set_ylabel(ylabel, rotation='vertical', fontsize=fontsize)
        # change axis font size
        axarr[row_nr, col_nr].tick_params(labelsize=fontsize)
        # create legend?#http://stackoverflow.com/questions/9834452/how-do-i-make-a-single-legend-for-many-subplots-with-matplotlib
        legend_obj = axarr[row_nr, col_nr].legend(cols_for_analysis, loc='upper right', fontsize=fontsize)
        # add title
        # axarr[row_nr, col_nr].set_title(title,fontsize=fontsize)
        # add figure number to top left of subplot
        axarr[row_nr, col_nr].annotate(s=str(Fig_Nr) + '.', xy=(0.04, 0.9), fontsize=5, xytext=None,
                                       xycoords='axes fraction', alpha=0.75)
        # add figure title to top left of subplot
        axarr[row_nr, col_nr].annotate(s=title, xy=(0.1, 0.9), fontsize=5, xytext=None, xycoords='axes fraction',
                                       alpha=0.75)
        # improve ggplot style for a canvas (fig) with 4 figures (plots)
        utils.improve_ggplot_for_4_plots(axarr, row_nr, col_nr, backgroundcolour, legend_obj)

        '''
        Fig12: TMD 1-5 only
        '''
        Fig_Nr = 12
        # add non-functional "if" line to aid document navigation in some IDEs (e.g. Spyder)
        if fig_title == "Fig12: TMD 1-5 only":
            pass
        sys.stdout.write(str(Fig_Nr) + ', ')
        col_start = 0
        col_end = 5
        cols_to_analyse = df_mean_AAIMON_each_TM.columns[col_start:col_end]
        title = 'TMD 1 to 5, all data'
        newfig, savefig, fig_nr, plot_nr_in_fig, row_nr, col_nr = dict_organising_subplots[Fig_Nr]
        # if a new figure should be created (either because the orig is full, or the last TMD is analysed)
        fig, axarr = utils.create_new_fig_if_necessary(newfig, fig, axarr, nrows_in_each_fig, ncols_in_each_fig, dpi=300)
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
            linecontainer = axarr[row_nr, col_nr].plot(centre_of_bar_in_x_axis, freq_counts_I, alpha=alpha,
                                                       linewidth=0.3)  # edgecolor='black',
        # label the x-axis for each plot, based on the TMD
        axarr[row_nr, col_nr].set_xlabel(xlabel, fontsize=fontsize)
        # move the x-axis label closer to the x-axis
        axarr[row_nr, col_nr].xaxis.set_label_coords(0.45, -0.085)
        axarr[row_nr, col_nr].set_ylabel(ylabel, rotation='vertical', fontsize=fontsize)
        # change axis font size
        axarr[row_nr, col_nr].tick_params(labelsize=fontsize)
        # create legend?#http://stackoverflow.com/questions/9834452/how-do-i-make-a-single-legend-for-many-subplots-with-matplotlib
        legend_obj = axarr[row_nr, col_nr].legend(cols_to_analyse, loc='upper right', fontsize=fontsize)
        # add title
        # axarr[row_nr, col_nr].set_title(title,fontsize=fontsize)
        # add figure number to top left of subplot
        axarr[row_nr, col_nr].annotate(s=str(Fig_Nr) + '.', xy=(0.04, 0.9), fontsize=5, xytext=None,
                                       xycoords='axes fraction', alpha=0.75)
        # add figure title to top left of subplot
        axarr[row_nr, col_nr].annotate(s=title, xy=(0.1, 0.9), fontsize=5, xytext=None, xycoords='axes fraction',
                                       alpha=0.75)
        # improve ggplot style for a canvas (fig) with 4 figures (plots)
        utils.improve_ggplot_for_4_plots(axarr, row_nr, col_nr, backgroundcolour, legend_obj)

        # save the figure as it is
        savefig = True
        # save the figure if necessary (i.e., if the maximum number of plots per figure has been obtained)
        utils.savefig_if_necessary(savefig, fig, fig_nr, base_filepath=pathdict["base_filename_summaries"])

        '''
        Fig13: TMD 5-10 only
        '''
        Fig_Nr = 13
        title = 'TMD 5-10, all data'
        sys.stdout.write(str(Fig_Nr) + ', ')
        col_start = 5
        col_end = 10
        cols_to_analyse = ['TM01'] + list(df_mean_AAIMON_each_TM.columns[col_start:col_end])

        newfig, savefig, fig_nr, plot_nr_in_fig, row_nr, col_nr = dict_organising_subplots[Fig_Nr]
        # if a new figure should be created (either because the orig is full, or the last TMD is analysed)
        fig, axarr = utils.create_new_fig_if_necessary(newfig, fig, axarr, nrows_in_each_fig, ncols_in_each_fig, dpi=300)
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
            linecontainer = axarr[row_nr, col_nr].plot(centre_of_bar_in_x_axis, freq_counts_I, alpha=alpha,
                                                       linewidth=0.3)  # edgecolor='black',
        # label the x-axis for each plot, based on the TMD
        axarr[row_nr, col_nr].set_xlabel(xlabel, fontsize=fontsize)
        # move the x-axis label closer to the x-axis
        axarr[row_nr, col_nr].xaxis.set_label_coords(0.45, -0.085)
        axarr[row_nr, col_nr].set_ylabel(ylabel, rotation='vertical', fontsize=fontsize)
        # change axis font size
        axarr[row_nr, col_nr].tick_params(labelsize=fontsize)
        # create legend?#http://stackoverflow.com/questions/9834452/how-do-i-make-a-single-legend-for-many-subplots-with-matplotlib
        legend_obj = axarr[row_nr, col_nr].legend(cols_to_analyse, loc='upper right', fontsize=fontsize)
        # add title
        # axarr[row_nr, col_nr].set_title(title,fontsize=fontsize)
        # add figure number to top left of subplot
        axarr[row_nr, col_nr].annotate(s=str(Fig_Nr) + '.', xy=(0.04, 0.9), fontsize=5, xytext=None,
                                       xycoords='axes fraction', alpha=0.75)
        # add figure title to top left of subplot
        axarr[row_nr, col_nr].annotate(s=title, xy=(0.1, 0.9), fontsize=5, xytext=None, xycoords='axes fraction',
                                       alpha=0.75)
        # improve ggplot style for a canvas (fig) with 4 figures (plots)
        utils.improve_ggplot_for_4_plots(axarr, row_nr, col_nr, backgroundcolour, legend_obj)

        '''
        Fig14: TMD 10-15 only
        '''
        Fig_Nr = 14
        # add non-functional "if" line to aid document navigation in some IDEs (e.g. Spyder)
        if fig_title == "Fig14: TMD 10-15 only":
            pass
        title = 'TMD 10-15, all data'
        sys.stdout.write(str(Fig_Nr) + ', ')
        col_start = 10
        col_end = 15
        cols_to_analyse = ['TM01'] + list(df_mean_AAIMON_each_TM.columns[col_start:col_end])
        newfig, savefig, fig_nr, plot_nr_in_fig, row_nr, col_nr = dict_organising_subplots[Fig_Nr]
        # if a new figure should be created (either because the orig is full, or the last TMD is analysed)
        fig, axarr = utils.create_new_fig_if_necessary(newfig, fig, axarr, nrows_in_each_fig, ncols_in_each_fig, dpi=300)
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
            linecontainer = axarr[row_nr, col_nr].plot(centre_of_bar_in_x_axis, freq_counts_normalised, alpha=alpha,
                                                       linewidth=0.3)  # edgecolor='black',
        # label the x-axis for each plot, based on the TMD
        axarr[row_nr, col_nr].set_xlabel(xlabel, fontsize=fontsize)
        # move the x-axis label closer to the x-axis
        axarr[row_nr, col_nr].xaxis.set_label_coords(0.45, -0.085)
        axarr[row_nr, col_nr].set_ylabel(ylabel, rotation='vertical', fontsize=fontsize)
        # change axis font size
        axarr[row_nr, col_nr].tick_params(labelsize=fontsize)
        # create legend?#http://stackoverflow.com/questions/9834452/how-do-i-make-a-single-legend-for-many-subplots-with-matplotlib
        legend_obj = axarr[row_nr, col_nr].legend(cols_to_analyse, loc='upper right', fontsize=fontsize)
        # add title
        # axarr[row_nr, col_nr].set_title(title,fontsize=fontsize)
        # add figure number to top left of subplot
        axarr[row_nr, col_nr].annotate(s=str(Fig_Nr) + '.', xy=(0.04, 0.9), fontsize=5, xytext=None,
                                       xycoords='axes fraction', alpha=0.75)
        # add figure title to top left of subplot
        axarr[row_nr, col_nr].annotate(s=title, xy=(0.1, 0.9), fontsize=5, xytext=None, xycoords='axes fraction',
                                       alpha=0.75)
        # improve ggplot style for a canvas (fig) with 4 figures (plots)
        utils.improve_ggplot_for_4_plots(axarr, row_nr, col_nr, backgroundcolour, legend_obj)

        '''
        Fig15: TMD 15-20 only
        '''
        Fig_Nr = 15
        title = 'TMD 15-20, all data'
        sys.stdout.write(str(Fig_Nr) + ', ')
        col_start = 15
        col_end = 20
        cols_to_analyse = ['TM01'] + list(df_mean_AAIMON_each_TM.columns[col_start:col_end])
        newfig, savefig, fig_nr, plot_nr_in_fig, row_nr, col_nr = dict_organising_subplots[Fig_Nr]
        # if a new figure should be created (either because the orig is full, or the last TMD is analysed)
        fig, axarr = utils.create_new_fig_if_necessary(newfig, fig, axarr, nrows_in_each_fig, ncols_in_each_fig, dpi=300)
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
            linecontainer = axarr[row_nr, col_nr].plot(centre_of_bar_in_x_axis, freq_counts_normalised, alpha=alpha,
                                                       linewidth=0.3)  # edgecolor='black',
        # label the x-axis for each plot, based on the TMD
        axarr[row_nr, col_nr].set_xlabel(xlabel, fontsize=fontsize)
        # move the x-axis label closer to the x-axis
        axarr[row_nr, col_nr].xaxis.set_label_coords(0.45, -0.085)
        axarr[row_nr, col_nr].set_ylabel(ylabel, rotation='vertical', fontsize=fontsize)
        # change axis font size
        axarr[row_nr, col_nr].tick_params(labelsize=fontsize)
        # create legend?#http://stackoverflow.com/questions/9834452/how-do-i-make-a-single-legend-for-many-subplots-with-matplotlib
        legend_obj = axarr[row_nr, col_nr].legend(cols_to_analyse, loc='upper right', fontsize=fontsize)
        # add title
        # axarr[row_nr, col_nr].set_title(title,fontsize=fontsize)
        # add figure number to top left of subplot
        axarr[row_nr, col_nr].annotate(s=str(Fig_Nr) + '.', xy=(0.04, 0.9), fontsize=5, xytext=None,
                                       xycoords='axes fraction', alpha=0.75)
        # add figure title to top left of subplot
        axarr[row_nr, col_nr].annotate(s=title, xy=(0.1, 0.9), fontsize=5, xytext=None, xycoords='axes fraction',
                                       alpha=0.75)
        # improve ggplot style for a canvas (fig) with 4 figures (plots)
        utils.improve_ggplot_for_4_plots(axarr, row_nr, col_nr, backgroundcolour, legend_obj)

        '''
        Fig16: ONLY proteins with 7 TMDs
        '''
        Fig_Nr = 16
        # add non-functional "if" line to aid document navigation in some IDEs (e.g. Spyder)
        if fig_title == "Fig16: ONLY proteins with 7 TMDs":
            pass
        title = 'ONLY prot with 7 TMDs'
        sys.stdout.write(str(Fig_Nr) + ', ')
        # create a new dataframe containing data only for proteins with 7 TMDs
        df_seven = df.loc[df.number_of_TMDs == 7].copy()
        logging.info('df_seven.shape: %s' % str(df_seven.shape))

        # if there are any proteins with 7 TM helices in the dataset (may not be true for beta-barrel datasets)
        if df_seven.shape[0] != 0:
            df_mean_AAIMON_each_TM_7TM, max_num_TMDs, legend = utils.create_df_with_mean_AAIMON_each_TM(df_seven)
            cols_for_analysis = ['TM01', 'TM06', 'TM07']
            newfig, savefig, fig_nr, plot_nr_in_fig, row_nr, col_nr = dict_organising_subplots[Fig_Nr]
            # if a new figure should be created (either because the orig is full, or the last TMD is analysed)
            fig, axarr = utils.create_new_fig_if_necessary(newfig, fig, axarr, nrows_in_each_fig, ncols_in_each_fig,
                                                           dpi=300)

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
                linecontainer = axarr[row_nr, col_nr].plot(centre_of_bar_in_x_axis, freq_counts_I, alpha=alpha,
                                                           linewidth=0.3)  # edgecolor='black',
            # label the x-axis for each plot, based on the TMD
            axarr[row_nr, col_nr].set_xlabel(xlabel, fontsize=fontsize)
            # move the x-axis label closer to the x-axis
            axarr[row_nr, col_nr].xaxis.set_label_coords(0.45, -0.085)
            axarr[row_nr, col_nr].set_ylabel(ylabel, rotation='vertical', fontsize=fontsize)
            # change axis font size
            axarr[row_nr, col_nr].tick_params(labelsize=fontsize)
            # create legend?#http://stackoverflow.com/questions/9834452/how-do-i-make-a-single-legend-for-many-subplots-with-matplotlib
            legend_obj = axarr[row_nr, col_nr].legend(cols_for_analysis, loc='upper right', fontsize=fontsize)
            # add title
            # axarr[row_nr, col_nr].set_title(title,fontsize=fontsize)
            # add figure number to top left of subplot
            axarr[row_nr, col_nr].annotate(s=str(Fig_Nr) + '.', xy=(0.04, 0.9), fontsize=5, xytext=None,
                                           xycoords='axes fraction', alpha=0.75)
            # add figure title to top left of subplot
            axarr[row_nr, col_nr].annotate(s=title, xy=(0.1, 0.9), fontsize=5, xytext=None, xycoords='axes fraction',
                                           alpha=0.75)
            # improve ggplot style for a canvas (fig) with 4 figures (plots)
            utils.improve_ggplot_for_4_plots(axarr, row_nr, col_nr, backgroundcolour, legend_obj)

        # save the figure as it is
        savefig = True
        # save the figure if necessary (i.e., if the maximum number of plots per figure has been obtained)
        utils.savefig_if_necessary(savefig, fig, fig_nr, base_filepath=pathdict["base_filename_summaries"])

        '''
        Fig17: Less than 12 TMDs vs at least 12 TMDs
        '''
        Fig_Nr = 17
        title = '<12 TMDs vs >12 TMDs'
        sys.stdout.write(str(Fig_Nr) + ', ')
        df_under_12 = df.loc[df.number_of_TMDs < 12]

        newfig, savefig, fig_nr, plot_nr_in_fig, row_nr, col_nr = dict_organising_subplots[Fig_Nr]
        # if a new figure should be created (either because the orig is full, or the last TMD is analysed)
        fig, axarr = utils.create_new_fig_if_necessary(newfig, fig, axarr, nrows_in_each_fig, ncols_in_each_fig, dpi=300)
        # create numpy array of membranous over nonmembranous conservation ratios (identity)
        hist_data_AAIMON_mean = np.array(df_under_12['AAIMON_ratio_mean_all_TMDs'].dropna())
        # use numpy to create a histogram
        freq_counts_I, bin_array_I = np.histogram(hist_data_AAIMON_mean, bins=binlist)
        # assuming all of the bins are exactly the same size, make the width of the column equal to 70% of each bin
        col_width = float('%0.3f' % (0.95 * (bin_array_I[1] - bin_array_I[0])))
        # when align='center', the central point of the bar in the x-axis is simply the middle of the bins ((bin_0-bin_1)/2, etc)
        centre_of_bar_in_x_axis = (bin_array_I[:-2] + bin_array_I[1:-1]) / 2
        # add the final bin, which is physically located just after the last regular bin but represents all higher values
        bar_width = centre_of_bar_in_x_axis[3] - centre_of_bar_in_x_axis[2]
        centre_of_bar_in_x_axis = np.append(centre_of_bar_in_x_axis, centre_of_bar_in_x_axis[-1] + bar_width)
        barcontainer_AAIMON_mean = axarr[row_nr, col_nr].bar(left=centre_of_bar_in_x_axis, height=freq_counts_I,
                                                             align='center', width=col_width, color="#0489B1",
                                                             alpha=0.5, linewidth=0.1)  # edgecolor='black',
        # other colours that are compatible with colourblind readers: #8A084B Dark red, #B45F04 deep orange, reddish purple #4B088A
        # http://html-color-codes.info/
        # label the x-axis for each plot, based on the TMD
        axarr[row_nr, col_nr].set_xlabel('average conservation ratio (membranous over nonmembranous)', fontsize=fontsize)
        # move the x-axis label closer to the x-axis
        axarr[row_nr, col_nr].xaxis.set_label_coords(0.45, -0.085)
        xlim_min = 0.8
        xlim_max = 1.5
        axarr[row_nr, col_nr].set_xlim(xlim_min, xlim_max)
        # set x-axis ticks
        # use the slide selection to select every second item in the list as an xtick(axis label)
        axarr[row_nr, col_nr].set_xticks([float('%0.1f' % c) for c in centre_of_bar_in_x_axis[::3]])
        axarr[row_nr, col_nr].set_ylabel('freq', rotation='vertical', fontsize=fontsize)
        # change axis font size
        axarr[row_nr, col_nr].tick_params(labelsize=fontsize)
        # improve ggplot style for a canvas (fig) with 4 figures (plots)
        utils.improve_ggplot_for_4_plots(axarr, row_nr, col_nr, backgroundcolour, legend_obj)

        '''
        At least 12 TMDs
        '''
        df_at_least_12 = df.loc[df.number_of_TMDs >= 12]

        # create numpy array of membranous over nonmembranous conservation ratios (identity)
        hist_data_AAIMON_mean = np.array(df_at_least_12['AAIMON_ratio_mean_all_TMDs'].dropna())
        # use numpy to create a histogram
        freq_counts_I, bin_array_I = np.histogram(hist_data_AAIMON_mean, bins=binlist)
        # assuming all of the bins are exactly the same size, make the width of the column equal to 70% of each bin
        col_width = float('%0.3f' % (0.95 * (bin_array_I[1] - bin_array_I[0])))
        # when align='center', the central point of the bar in the x-axis is simply the middle of the bins ((bin_0-bin_1)/2, etc)
        centre_of_bar_in_x_axis = (bin_array_I[:-2] + bin_array_I[1:-1]) / 2
        # add the final bin, which is physically located just after the last regular bin but represents all higher values
        bar_width = centre_of_bar_in_x_axis[3] - centre_of_bar_in_x_axis[2]
        centre_of_bar_in_x_axis = np.append(centre_of_bar_in_x_axis, centre_of_bar_in_x_axis[-1] + bar_width)
        barcontainer_AAIMON_mean = axarr[row_nr, col_nr].bar(left=centre_of_bar_in_x_axis, height=freq_counts_I * 5,
                                                             align='center', width=col_width, color='#B45F04',
                                                             alpha=0.5, linewidth=0.1)  # edgecolor='black',
        # other colours that are compatible with colourblind readers: #8A084B Dark red, #B45F04 deep orange, reddish purple #4B088A
        # http://html-color-codes.info/
        # label the x-axis for each plot, based on the TMD
        axarr[row_nr, col_nr].set_xlabel('average conservation ratio (membranous over nonmembranous)', fontsize=fontsize)
        # move the x-axis label closer to the x-axis
        axarr[row_nr, col_nr].xaxis.set_label_coords(0.45, -0.085)
        xlim_min = 0.8
        xlim_max = 1.5
        axarr[row_nr, col_nr].set_xlim(xlim_min, xlim_max)
        # set x-axis ticks
        # use the slide selection to select every second item in the list as an xtick(axis label)
        axarr[row_nr, col_nr].set_xticks([float('%0.1f' % c) for c in centre_of_bar_in_x_axis[::3]])
        axarr[row_nr, col_nr].set_ylabel('freq', rotation='vertical', fontsize=fontsize)
        # change axis font size
        axarr[row_nr, col_nr].tick_params(labelsize=fontsize)
        # create legend?#http://stackoverflow.com/questions/9834452/how-do-i-make-a-single-legend-for-many-subplots-with-matplotlib
        legend_obj = axarr[row_nr, col_nr].legend(['AAIMON 01-11 TMDs', 'AAIMON >= 12 TMDs'], loc='upper right',
                                                  fontsize=fontsize)
        # add figure number to top left of subplot
        axarr[row_nr, col_nr].annotate(s=str(Fig_Nr) + '.', xy=(0.04, 0.9), fontsize=5, xytext=None,
                                       xycoords='axes fraction', alpha=0.75)
        # add figure title to top left of subplot
        axarr[row_nr, col_nr].annotate(s=title, xy=(0.1, 0.9), fontsize=5, xytext=None, xycoords='axes fraction',
                                       alpha=0.75)
        # improve ggplot style for a canvas (fig) with 4 figures (plots)
        utils.improve_ggplot_for_4_plots(axarr, row_nr, col_nr, backgroundcolour, legend_obj)

        '''
        Fig20: Boxplot of all TMDs
        '''
        Fig_Nr = 18
        title = 'Boxplot of all TMDs'
        # add non-functional "if" line to aid document navigation in some IDEs (e.g. Spyder)
        if fig_title == "Fig20: Boxplot of all TMDs":
            pass
        sys.stdout.write(str(Fig_Nr) + ', ')
        newfig, savefig, fig_nr, plot_nr_in_fig, row_nr, col_nr = dict_organising_subplots[Fig_Nr]
        # if a new figure should be created (either because the orig is full, or the last TMD is analysed)
        fig, axarr = utils.create_new_fig_if_necessary(newfig, fig, axarr, nrows_in_each_fig, ncols_in_each_fig, dpi=300)

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
        for i in range(1, max_num_TMDs + 1):
            TM = 'TM%02d' % i
            hist_data_AAIMON_each_TM = df['TM%02d_AAIMON_ratio_mean' % i].dropna()
            if len(hist_data_AAIMON_each_TM) > 0:
                data_to_plot.append(hist_data_AAIMON_each_TM)
                legend.append(TM)

        meanpointprops = dict(marker='o', markerfacecolor='black', markersize=2)  # markeredgecolor='0.75',

        flierprops = dict(marker='o', markerfacecolor='green', markersize=12,
                          linestyle='none')
        boxplotcontainer = axarr[row_nr, col_nr].boxplot(data_to_plot, sym='+', whis=1.5, showmeans=True,
                                                         meanprops=meanpointprops)
        axarr[row_nr, col_nr].tick_params(labelsize=fontsize)
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
        axarr[row_nr, col_nr].get_xaxis().tick_bottom()
        axarr[row_nr, col_nr].get_yaxis().tick_left()
        ## Custom x-axis labels
        axarr[row_nr, col_nr].set_xticklabels(legend, rotation=45)
        # add figure number to top left of subplot
        axarr[row_nr, col_nr].annotate(s=str(Fig_Nr) + '.', xy=(0.04, 0.9), fontsize=5, xytext=None,
                                       xycoords='axes fraction', alpha=0.75)
        # add figure title to top left of subplot
        axarr[row_nr, col_nr].annotate(s=title, xy=(0.1, 0.9), fontsize=5, xytext=None, xycoords='axes fraction',
                                       alpha=0.75)
        # improve ggplot style for a canvas (fig) with 4 figures (plots)
        utils.improve_ggplot_for_4_plots(axarr, row_nr, col_nr, backgroundcolour, legend_obj)

        # conduct a keyword analysis for sequences that come from uniprot, or have functional keywords annotated
        if 'uniprot_KW' in df.columns:
            '''
            Fig18: SHOW ONLY GPCRS IN FULL DATASET
            '''

            Fig_Nr = 19
            title = 'only GPCR in uniprot KW, NORM'
            # add non-functional "if" line to aid document navigation in some IDEs (e.g. Spyder)
            if fig_title == "Fig18: SHOW ONLY GPCRS IN FULL DATASET":
                pass
            sys.stdout.write(str(Fig_Nr) + ', ')
            # if it hasn't been done already, convert the keywords from a stringlist to a python list
            if isinstance(df['uniprot_KW'][0], str):
                df['uniprot_KW'] = df['uniprot_KW'].apply(lambda x: ast.literal_eval(x))
            # create a new column showing whether the protein is a GPCR
            df['G-protein coupled receptor'] = df['uniprot_KW'].apply(lambda x: 'G-protein coupled receptor' in x)
            df_GPCR = df.loc[df['G-protein coupled receptor'] == True]

            newfig, savefig, fig_nr, plot_nr_in_fig, row_nr, col_nr = dict_organising_subplots[Fig_Nr]
            # if a new figure should be created (either because the orig is full, or the last TMD is analysed)
            fig, axarr = utils.create_new_fig_if_necessary(newfig, fig, axarr, nrows_in_each_fig, ncols_in_each_fig,
                                                           dpi=300)
            # create numpy array of membranous over nonmembranous conservation ratios (identity)
            hist_data_AAIMON_mean = np.array(df_GPCR['AAIMON_ratio_mean_all_TMDs'].dropna())
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
            barcontainer_AAIMON_mean = axarr[row_nr, col_nr].bar(left=centre_of_bar_in_x_axis,
                                                                 height=freq_counts_normalised,
                                                                 align='center', width=col_width, color="#0489B1",
                                                                 alpha=0.5, linewidth=0.1)  # edgecolor='black',
            # other colours that are compatible with colourblind readers: #8A084B Dark red, #B45F04 deep orange, reddish purple #4B088A
            # http://html-color-codes.info/
            # label the x-axis for each plot, based on the TMD
            axarr[row_nr, col_nr].set_xlabel('average conservation ratio (membranous over nonmembranous)',
                                             fontsize=fontsize)
            # move the x-axis label closer to the x-axis
            axarr[row_nr, col_nr].xaxis.set_label_coords(0.45, -0.085)
            xlim_min = 0.8
            xlim_max = 1.5
            axarr[row_nr, col_nr].set_xlim(xlim_min, xlim_max)
            # set x-axis ticks
            # use the slide selection to select every second item in the list as an xtick(axis label)
            axarr[row_nr, col_nr].set_xticks([float('%0.1f' % c) for c in centre_of_bar_in_x_axis[::3]])
            axarr[row_nr, col_nr].set_ylabel('freq', rotation='vertical', fontsize=fontsize)
            # change axis font size
            axarr[row_nr, col_nr].tick_params(labelsize=fontsize)

            '''
            NON-GPCRS
            '''
            df_nonGPCR = df.loc[df['G-protein coupled receptor'] == False]

            # create numpy array of membranous over nonmembranous conservation ratios (identity)
            hist_data_AAIMON_mean = np.array(df_nonGPCR['AAIMON_ratio_mean_all_TMDs'].dropna())
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
            barcontainer_AAIMON_mean = axarr[row_nr, col_nr].bar(left=centre_of_bar_in_x_axis,
                                                                 height=freq_counts_normalised,
                                                                 align='center', width=col_width, color='#B45F04',
                                                                 alpha=0.5, linewidth=0.1)  # edgecolor='black',
            # other colours that are compatible with colourblind readers: #8A084B Dark red, #B45F04 deep orange, reddish purple #4B088A
            # http://html-color-codes.info/
            # label the x-axis for each plot, based on the TMD
            axarr[row_nr, col_nr].set_xlabel('average conservation ratio (membranous over nonmembranous)',
                                             fontsize=fontsize)
            # move the x-axis label closer to the x-axis
            axarr[row_nr, col_nr].xaxis.set_label_coords(0.45, -0.085)
            xlim_min = 0.8
            xlim_max = 1.5
            axarr[row_nr, col_nr].set_xlim(xlim_min, xlim_max)
            # set x-axis ticks
            # use the slide selection to select every second item in the list as an xtick(axis label)
            axarr[row_nr, col_nr].set_xticks([float('%0.1f' % c) for c in centre_of_bar_in_x_axis[::3]])
            axarr[row_nr, col_nr].set_ylabel('freq', rotation='vertical', fontsize=fontsize)
            # change axis font size
            axarr[row_nr, col_nr].tick_params(labelsize=fontsize)
            # create legend?#http://stackoverflow.com/questions/9834452/how-do-i-make-a-single-legend-for-many-subplots-with-matplotlib
            legend_obj = axarr[row_nr, col_nr].legend(['AAIMON GPCRs', 'AAIMON non-GPCRs'], loc='upper right',
                                                      fontsize=fontsize)
            # add figure number to top left of subplot
            axarr[row_nr, col_nr].annotate(s=str(Fig_Nr) + '.', xy=(0.04, 0.9), fontsize=5, xytext=None,
                                           xycoords='axes fraction', alpha=0.75)
            # add figure title to top left of subplot
            axarr[row_nr, col_nr].annotate(s=title, xy=(0.1, 0.9), fontsize=5, xytext=None, xycoords='axes fraction',
                                           alpha=0.75)
            # improve ggplot style for a canvas (fig) with 4 figures (plots)
            utils.improve_ggplot_for_4_plots(axarr, row_nr, col_nr, backgroundcolour, legend_obj)

            '''
            Fig19: SHOW GPCRS VS FULL DATASET
            '''
            Fig_Nr = 20
            title = 'GPCR vs full dataset, NORM'
            sys.stdout.write(str(Fig_Nr) + ', ')
            df_GPCR = df_GPCR  # see above, already defined

            '''GPCR'''
            newfig, savefig, fig_nr, plot_nr_in_fig, row_nr, col_nr = dict_organising_subplots[Fig_Nr]
            # if a new figure should be created (either because the orig is full, or the last TMD is analysed)
            fig, axarr = utils.create_new_fig_if_necessary(newfig, fig, axarr, nrows_in_each_fig, ncols_in_each_fig,
                                                           dpi=300)
            # create numpy array of membranous over nonmembranous conservation ratios (identity)
            hist_data_AAIMON_mean = np.array(df_GPCR['AAIMON_ratio_mean_all_TMDs'].dropna())
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
            barcontainer_AAIMON_mean = axarr[row_nr, col_nr].bar(left=centre_of_bar_in_x_axis,
                                                                 height=freq_counts_normalised,
                                                                 align='center', width=col_width, color="#0489B1",
                                                                 alpha=0.5, linewidth=0.1)  # edgecolor='black',
            # other colours that are compatible with colourblind readers: #8A084B Dark red, #B45F04 deep orange, reddish purple #4B088A
            # http://html-color-codes.info/
            # label the x-axis for each plot, based on the TMD
            axarr[row_nr, col_nr].set_xlabel('average conservation ratio (membranous over nonmembranous)',
                                             fontsize=fontsize)
            # move the x-axis label closer to the x-axis
            axarr[row_nr, col_nr].xaxis.set_label_coords(0.45, -0.085)
            xlim_min = 0.8
            xlim_max = 1.5
            axarr[row_nr, col_nr].set_xlim(xlim_min, xlim_max)
            # set x-axis ticks
            # use the slide selection to select every second item in the list as an xtick(axis label)
            axarr[row_nr, col_nr].set_xticks([float('%0.1f' % c) for c in centre_of_bar_in_x_axis[::3]])
            axarr[row_nr, col_nr].set_ylabel('freq', rotation='vertical', fontsize=fontsize)
            # change axis font size
            axarr[row_nr, col_nr].tick_params(labelsize=fontsize)

            '''
            full dataset
            '''
            # df_nonGPCR = df.loc[df['Gprotein'] == False]

            # create numpy array of membranous over nonmembranous conservation ratios (identity)
            hist_data_AAIMON_mean = np.array(df['AAIMON_ratio_mean_all_TMDs'].dropna())
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
            barcontainer_AAIMON_mean = axarr[row_nr, col_nr].bar(left=centre_of_bar_in_x_axis,
                                                                 height=freq_counts_normalised,
                                                                 align='center', width=col_width, color='#B45F04',
                                                                 alpha=0.5, linewidth=0.1)  # edgecolor='black',
            # other colours that are compatible with colourblind readers: #8A084B Dark red, #B45F04 deep orange, reddish purple #4B088A
            # http://html-color-codes.info/
            # label the x-axis for each plot, based on the TMD
            axarr[row_nr, col_nr].set_xlabel('average conservation ratio (membranous over nonmembranous)',
                                             fontsize=fontsize)
            # move the x-axis label closer to the x-axis
            axarr[row_nr, col_nr].xaxis.set_label_coords(0.45, -0.085)
            xlim_min = 0.8
            xlim_max = 1.5
            axarr[row_nr, col_nr].set_xlim(xlim_min, xlim_max)
            # set x-axis ticks
            # use the slide selection to select every second item in the list as an xtick(axis label)
            axarr[row_nr, col_nr].set_xticks([float('%0.1f' % c) for c in centre_of_bar_in_x_axis[::3]])
            axarr[row_nr, col_nr].set_ylabel('freq', rotation='vertical', fontsize=fontsize)
            # change axis font size
            axarr[row_nr, col_nr].tick_params(labelsize=fontsize)
            # create legend?#http://stackoverflow.com/questions/9834452/how-do-i-make-a-single-legend-for-many-subplots-with-matplotlib
            legend_obj = axarr[row_nr, col_nr].legend(['AAIMON GPCR', 'AAIMON ALL'], loc='upper right',
                                                      fontsize=fontsize)
            # add figure number to top left of subplot
            axarr[row_nr, col_nr].annotate(s=str(Fig_Nr) + '.', xy=(0.04, 0.9), fontsize=5, xytext=None,
                                           xycoords='axes fraction', alpha=0.75)
            # add figure title to top left of subplot
            axarr[row_nr, col_nr].annotate(s=title, xy=(0.1, 0.9), fontsize=5, xytext=None, xycoords='axes fraction',
                                           alpha=0.75)
            # improve ggplot style for a canvas (fig) with 4 figures (plots)
            utils.improve_ggplot_for_4_plots(axarr, row_nr, col_nr, backgroundcolour, legend_obj)

            # save the figure if necessary (i.e., if the maximum number of plots per figure has been obtained)
            utils.savefig_if_necessary(savefig, fig, fig_nr, base_filepath=pathdict["base_filename_summaries"])

            '''
            Fig21: Boxplot only GPCRs
            '''
            Fig_Nr = 21
            title = 'Only GPCRs, boxplot for each TMD'
            sys.stdout.write(str(Fig_Nr) + ', ')
            newfig, savefig, fig_nr, plot_nr_in_fig, row_nr, col_nr = dict_organising_subplots[Fig_Nr]
            # if the plot is the last one, the figure should be saved
            # if a new figure should be created (either because the orig is full, or the last TMD is analysed)
            fig, axarr = utils.create_new_fig_if_necessary(newfig, fig, axarr, nrows_in_each_fig, ncols_in_each_fig,
                                                           dpi=300)

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
            for i in range(1, max_num_TMDs + 1):
                TM = 'TM%02d' % i
                hist_data_AAIMON_each_TM = df_GPCR['TM%02d_AAIMON_ratio_mean' % i].dropna()
                if len(hist_data_AAIMON_each_TM) > 0:
                    data_to_plot.append(hist_data_AAIMON_each_TM)
                    legend.append(TM)

            meanpointprops = dict(marker='o', markerfacecolor='black', markersize=3)  # markeredgecolor='0.75',

            # flierprops = dict(marker='o', color = 'black', markerfacecolor='black', markersize=1)

            # flierprops = dict(marker='o',color='0.1', alpha=0.1)
            flierprops = dict(marker='o', markerfacecolor='green', markersize=12,
                              linestyle='none')
            boxplotcontainer = axarr[row_nr, col_nr].boxplot(data_to_plot, sym='+', whis=1.5, showmeans=True,
                                                             meanprops=meanpointprops)
            axarr[row_nr, col_nr].tick_params(labelsize=fontsize)
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
            axarr[row_nr, col_nr].get_xaxis().tick_bottom()
            axarr[row_nr, col_nr].get_yaxis().tick_left()
            ## Custom x-axis labels
            axarr[row_nr, col_nr].set_xticklabels(legend, rotation=45)
            # add figure number to top left of subplot
            axarr[row_nr, col_nr].annotate(s=str(Fig_Nr) + '.', xy=(0.04, 0.9), fontsize=5, xytext=None,
                                           xycoords='axes fraction', alpha=0.75)
            # add figure title to top left of subplot
            axarr[row_nr, col_nr].annotate(s=title, xy=(0.1, 0.9), fontsize=5, xytext=None, xycoords='axes fraction',
                                           alpha=0.75)
            # improve ggplot style for a canvas (fig) with 4 figures (plots)
            utils.improve_ggplot_for_4_plots(axarr, row_nr, col_nr, backgroundcolour, legend_obj)

            '''
            SAVE GRAPH BEFORE KEYWORD ANALYSIS
            '''
            # save the figure as it is
            savefig = True
            # save the figure if necessary (i.e., if the maximum number of plots per figure has been obtained)
            utils.savefig_if_necessary(savefig, fig, fig_nr, base_filepath=pathdict["base_filename_summaries"])

            '''
            Fig25 AND ONWARDS: Histograms split by keywords
            '''
            starting_Fig_Nr = Fig_Nr + 1
            # add non-functional "if" line to aid document navigation in some IDEs (e.g. Spyder)
            if fig_title == "Fig25 AND ONWARDS: Histograms split by keywords":
                pass
            # convert the KW stringlist to a python list, if it hasn't been done already
            if isinstance(df['uniprot_KW'][0], str):
                df['uniprot_KW'] = df['uniprot_KW'].apply(lambda x: ast.literal_eval(x))
            # join all keywords together into a large list
            nested_list_all_KW = list(itertools.chain(*list(df['uniprot_KW'])))
            # convert list to pandas series
            all_KW_series = pd.Series(nested_list_all_KW)
            # obtain series of major keywords
            KW_counts = all_KW_series.value_counts()
            # exclude keywords with less than 50 applicable proteins
            KW_counts_major = KW_counts[KW_counts > 50]
            # create a list of keywords to be ignored
            list_ignored_KW = ['Transmembrane', 'Complete proteome', 'Reference proteome', 'Membrane',
                               'Transmembrane helix', 'Cell membrane', 'Repeat',
                               'Alternative splicing', 'Sodium', 'Potassium', 'Direct protein sequencing', 'Transducer',
                               'Polymorphism', 'Glycoprotein']
            # remove undesired keywords from the dataframe (Transmembrane, etc.)
            for KW in list_ignored_KW:
                if KW in KW_counts_major.index:
                    KW_counts_major = KW_counts_major.drop(KW)

            '''
            KW Cross-correlation analysis.
            Q: which uniprot keywords in this dataset are related to each other?
            '''
            KW_for_correlation_analysis = list(KW_counts.loc[KW_counts > 10].index)
            # remove undesired keywords from the dataframe (Transmembrane, etc.)
            for KW in list_ignored_KW:
                if KW in KW_for_correlation_analysis:
                    KW_for_correlation_analysis.remove(KW)

            # create an empty dataframe to hold the keyword correlation data
            len_df = len(KW_for_correlation_analysis)
            df_KW_corr = pd.DataFrame(index=KW_for_correlation_analysis, columns=KW_for_correlation_analysis)

            for KW in KW_for_correlation_analysis:
                # create a new columns describing if the KW is in the KW list of that protein
                df['contains_KW'] = df['uniprot_KW'].apply(lambda x: KW in x)
                # slice dataframe to view only entries with that keyword
                df_KW = df.loc[df['contains_KW'] == True]

                # now we want to count the number of proteins containing each other keyword (subKW)
                # iterate through the full keyword list agoun
                for subKW in KW_for_correlation_analysis:
                    # sys.stdout.write(subKW + ', ')
                    # create a new column describing whether the protein KW list also contains the subKW
                    df['contains_subKW'] = df_KW['uniprot_KW'].apply(lambda x: subKW in x)
                    # count how many of the proteins contain the subKW
                    val_counts = df['contains_subKW'].value_counts()
                    if True in val_counts.index:
                        num_prot_contain_subKW = val_counts[True]
                    else:
                        num_prot_contain_subKW = 0
                    # now add that number to the array of all KW against all KW
                    df_KW_corr.loc[KW, subKW] = num_prot_contain_subKW
            # normalise the array by dividing the number of common keywords by the total number with that keyword.
            # NOTE WELL: it it normalised by the COLUMN keyword, not the row
            # copy the dataframe
            df_KW_corr_norm_by_col = df_KW_corr.copy()
            for col_KW in df_KW_corr_norm_by_col.dropna(axis=1, how='all').columns:
                # obtain the total number with that keyword
                total_for_that_KW = df_KW_corr_norm_by_col.loc[col_KW, col_KW]
                # divide each value in the column by the total with that keyword. Present as a percentage.
                df_KW_corr_norm_by_col[col_KW] = df_KW_corr_norm_by_col[col_KW].dropna().apply(
                    lambda x: x / total_for_that_KW * 100)
                # round to whole numbers
                df_KW_corr_norm_by_col[col_KW] = df_KW_corr_norm_by_col[col_KW].dropna().apply(lambda x: np.round(x))

            # create a dictionary of the ten top correlated subkeywords for each keyword
            top_corr_KW_dict = {}
            for col_KW in df_KW_corr_norm_by_col:
                # a = df_KW_corr_norm_by_col[col_KW].order()
                array_top_corr_KW = df_KW_corr_norm_by_col[col_KW].order(ascending=False)[1:11]
                top_corr_KW_dict[col_KW] = array_top_corr_KW

            '''
            25 and onwards: Colour lists and Histograms
            '''
            # prepare colour lists
            colour_lists = utils.create_colour_lists()
            TUM_colours_list_with_greys = colour_lists['TUM_colours_list_with_greys']
            colourlist_greys = [(0.6, 0.7764705882352941, 0.9058823529411765), 'None']

            # iterate over the dataframe. Note that acc = uniprot accession here.
            linspace_binlist = np.linspace(s["mp_smallest_bin"],
                                           s["mp_largest_bin"],
                                           s["mp_number_of_bins"])

            # add 30 as the last bin, to make sure 100% of the data is added to the histogram, including major outliers
            binlist = np.append(linspace_binlist,
                                s["mp_final_highest_bin"])

            dict_ttest_pvalues = {}

            for m, KW in enumerate(KW_counts_major.index):
                Fig_Nr = starting_Fig_Nr + m
                sys.stdout.write(str(Fig_Nr) + ', ')
                newfig, savefig, fig_nr, plot_nr_in_fig, row_nr, col_nr = dict_organising_subplots[Fig_Nr]
                # if a new figure should be created (either because the orig is full, or the last TMD is analysed)
                fig, axarr = utils.create_new_fig_if_necessary(newfig, fig, axarr, nrows_in_each_fig, ncols_in_each_fig,
                                                               dpi=300)
                # KW = 'Receptor'
                # create a new column showing whether the protein contains that keyword
                df[KW] = df['uniprot_KW'].apply(lambda x: KW in x)
                # create dataframes with entries containing, or not containing the particular keyword
                df_KW = df.loc[df[KW] == True]
                df_nonKW = df.loc[df[KW] == False]
                # combine into a list for iteration
                df_list_KW = [df_KW, df_nonKW]

                data_dict = {}
                data_names_list = ['containing keyword', 'without keyword']

                for n, dfK in enumerate(df_list_KW):
                    data_column = 'AAIMON_ratio_mean_all_TMDs'
                    color = colourlist_greys[n]
                    alpha = 1.0
                    col_width_value = 0.95
                    hist2_shift_to_right = 0.01
                    backgroundcolour = '0.95'
                    hist_data_AAIMON_mean = np.array(dfK[data_column].dropna())
                    # add the data to the dictionary for later analysis, named according to the list above
                    data_dict[data_names_list[n]] = hist_data_AAIMON_mean
                    # move the second histogram bins slightly to the right, so that the bars do not overlap
                    binlist_for_this_data = binlist + (hist2_shift_to_right * n)
                    # use numpy to create a histogram
                    freq_counts, bin_array = np.histogram(hist_data_AAIMON_mean, bins=binlist_for_this_data)
                    # normalize the frequency counts
                    freq_counts_normalised = freq_counts / freq_counts.max()  # use numpy to create a histogram
                    # assuming all of the bins are exactly the same size, make the width of the column equal to a percentage of each bin
                    dist_between_bins = (bin_array_I[1] - bin_array_I[0])
                    col_width = float('%0.3f' % (col_width_value * dist_between_bins))
                    # when align='center', the central point of the bar in the x-axis is simply the middle of the bins ((bin_0-bin_1)/2, etc)
                    centre_of_bar_in_x_axis = (bin_array[:-2] + bin_array[1:-1]) / 2
                    # the last bin is open, and holds all larger datapoints than specified in the binlist. it needs to be added to the list of bars
                    centre_of_bar_in_x_axis = np.append(centre_of_bar_in_x_axis,
                                                        centre_of_bar_in_x_axis[-1] + dist_between_bins)
                    # print a warning if the last bar is not actually empty
                    if freq_counts[-1] != 0:
                        logging.warning(
                            'for %s, the final bin (%s) has %i entries. The settings "largest_bin" may need to be adjusted.' % (
                            KW, bin_array[-1], freq_counts[-1]))
                    # plot the data as a bar-chart
                    barcontainer_AAIMON_mean = axarr[row_nr, col_nr].bar(left=centre_of_bar_in_x_axis,
                                                                         height=freq_counts_normalised,
                                                                         align='center', width=col_width, color=color,
                                                                         alpha=alpha, linewidth=0.1)
                # obtain data from dictionary
                data1 = data_dict[data_names_list[0]]
                data2 = data_dict[data_names_list[1]]
                # run ttest, assuming unequal variance (http://stackoverflow.com/questions/22611446/perform-2-sample-t-test)
                t, p = ttest_ind(data1, data2, equal_var=False)
                # add the ttest results to a dictionary
                dict_ttest_pvalues[KW] = p
                # determine symbol describing stat significance
                signif_symbol = utils.get_signif_symbol(p)

                # label the x-axis for each plot, based on the TMD
                axarr[row_nr, col_nr].set_xlabel('average conservation ratio (membranous over nonmembranous)',
                                                 fontsize=fontsize)
                # move the x-axis label closer to the x-axis
                axarr[row_nr, col_nr].xaxis.set_label_coords(0.45, -0.085)
                # alter x-axis limits. Note that if data outside these limits the axis will be automatically adjusted)
                xlim_min = 0.5
                xlim_max = 1.6
                axarr[row_nr, col_nr].set_xlim(xlim_min, xlim_max)
                # set x-axis ticks
                # use the slide selection to select every second item in the list as an xtick(axis label)
                axarr[row_nr, col_nr].set_xticks([float('%0.1f' % c) for c in centre_of_bar_in_x_axis[::3]])
                axarr[row_nr, col_nr].set_ylabel('freq', rotation='vertical', fontsize=fontsize)
                # change axis font size
                axarr[row_nr, col_nr].tick_params(labelsize=fontsize)
                # remove background grid
                axarr[row_nr, col_nr].grid(False)
                # add figure number to top left of subplot
                axarr[row_nr, col_nr].annotate(s=str(Fig_Nr) + '.', xy=(0.04, 0.9), fontsize=5, xytext=None,
                                               xycoords='axes fraction', alpha=0.75)
                # add figure title to top left of subplot
                axarr[row_nr, col_nr].annotate(s=KW, xy=(0.1, 0.9), fontsize=5, xytext=None, xycoords='axes fraction',
                                               alpha=0.75)
                # add a legend describing which hist contains the KW, or doesn't contain the KW
                legend_obj = axarr[row_nr, col_nr].legend(['containing keyword', 'without keyword'], loc='upper right',
                                                          fontsize=fontsize, )
                # add ttest result to graph
                axarr[row_nr, col_nr].annotate(s='p = %0.2g (%s)' % (p, signif_symbol), xy=(0.70, 0.75), fontsize=5,
                                               xytext=None, xycoords='axes fraction', alpha=0.75)
                # improve ggplot style for a canvas (fig) with 4 figures (plots)
                utils.improve_ggplot_for_4_plots(axarr, row_nr, col_nr, backgroundcolour, legend_obj)
                # annotate the correlated keywords to the graph
                corr_KW = [xc[0:16] for xc in list(top_corr_KW_dict[KW].index)]
                corr_KW_value = [int(xv) for xv in list(top_corr_KW_dict[KW])]
                # add "correlated keywords" title
                axarr[row_nr, col_nr].annotate(s='correlated keywords', xy=(0.7, 0.65), fontsize=4,
                                               xycoords='axes fraction', alpha=0.75)
                for ann_num, cKW in enumerate(corr_KW):
                    # define y position for annotation
                    ypos_ann = 0.6 - (ann_num * 0.05)
                    # define KW name position
                    ann_pos_xy = (0.7, ypos_ann)
                    # determine position of KW correlation value
                    val_pos_xy = (0.95, ypos_ann)
                    # annotate the KW name
                    axarr[row_nr, col_nr].annotate(s=cKW, xy=ann_pos_xy, fontsize=4, xycoords='axes fraction', alpha=0.75)
                    # annotate correlation value
                    axarr[row_nr, col_nr].annotate(s=str(corr_KW_value[ann_num]), xy=val_pos_xy, fontsize=4,
                                                   xycoords='axes fraction', alpha=0.75)
                # calculate the odds ratio (TMD/rest conservation for proteins with KW) / (TMD/rest conservation for proteins without KW)
                odds_ratio_KW_over_nonKW = '%0.2f' % (data1.mean() / data2.mean())
                # add odds ratio to figure
                axarr[row_nr, col_nr].annotate(s='odds ratio = ' + odds_ratio_KW_over_nonKW, xy=(0.7, 0.70), fontsize=5,
                                               xytext=None, xycoords='axes fraction', alpha=0.75)

                # save the figure if necessary (i.e., if the maximum number of plots per figure has been obtained)
                utils.savefig_if_necessary(savefig, fig, fig_nr, base_filepath=pathdict["base_filename_summaries"])

            # only create the following for multi-pass proteins (at least 2 proteins in lost with 7 TMDs)
            if dataset_contains_multipass_prots:
                '''
                Histograms split by keywords, EXCLUDING GPCR
                '''
                starting_Fig_Nr = Fig_Nr + 1
                # add non-functional "if" line to aid document navigation in some IDEs (e.g. Spyder)
                if fig_title == "Histograms split by keywords, EXCLUDING GPCR":
                    pass

                # recreate the nonGPCR dataframe from the main dataframe
                df_nonGPCR = df.loc[df['G-protein coupled receptor'] == False]

                # join all keywords together into a large list
                nested_list_all_KW = list(itertools.chain(*list(df_nonGPCR['uniprot_KW'])))
                # convert list to pandas series
                all_KW_series = pd.Series(nested_list_all_KW)
                # obtain series of major keywords
                KW_counts = all_KW_series.value_counts()
                # exclude keywords with less than 50 applicable proteins
                KW_counts_major = KW_counts[KW_counts > 50]

                # remove undesired keywords from the dataframe (Transmembrane, etc.)
                for KW in list_ignored_KW:
                    if KW in KW_counts_major.index:
                        KW_counts_major = KW_counts_major.drop(KW)

                # remove 'G-protein coupled receptor' from the series of keywords to test
                if 'G-protein coupled receptor' in KW_counts_major.index:
                    KW_counts_major = KW_counts_major.drop('G-protein coupled receptor')

                # #iterate over the dataframe. Note that acc = uniprot accession here.
                # linspace_binlist = np.linspace(s["mp_smallest_bin"],
                #                                s["mp_largest_bin"],
                #                                s["mp_number_of_bins"])
                #
                # #add 30 as the last bin, to make sure 100% of the data is added to the histogram, including major outliers
                # binlist = np.append(linspace_binlist, s["mp_final_highest_bin"])

                dict_ttest_pvalues = {}

                for m, KW in enumerate(KW_counts_major.index):
                    Fig_Nr = starting_Fig_Nr + m
                    sys.stdout.write(str(Fig_Nr) + ', ')
                    newfig, savefig, fig_nr, plot_nr_in_fig, row_nr, col_nr = dict_organising_subplots[Fig_Nr]
                    # if the plot is the last one, the figure should be saved
                    # if a new figure should be created (either because the orig is full, or the last TMD is analysed)
                    fig, axarr = utils.create_new_fig_if_necessary(newfig, fig, axarr, nrows_in_each_fig, ncols_in_each_fig,
                                                                   dpi=300)
                    # The boolean for each KW should already be in the dataframe.
                    # Simply reselect the data from the smaller df excluding GPCRs
                    df_KW = df_nonGPCR.loc[df_nonGPCR[KW] == True]
                    df_nonKW = df_nonGPCR.loc[df_nonGPCR[KW] == False]
                    df_list_KW = [df_KW, df_nonKW]

                    data_dict = {}
                    data_names_list = ['containing keyword', 'without keyword']

                    for n, dfK in enumerate(df_list_KW):
                        data_column = 'AAIMON_ratio_mean_all_TMDs'
                        color = colourlist_greys[n]
                        alpha = 1.0
                        col_width_value = 0.95
                        hist2_shift_to_right = 0.01
                        backgroundcolour = '0.95'
                        hist_data_AAIMON_mean = np.array(dfK[data_column].dropna())
                        # add the data to the dictionary for later analysis, named according to the list above
                        data_dict[data_names_list[n]] = hist_data_AAIMON_mean
                        # move the second histogram bins slightly to the right, so that the bars do not overlap
                        binlist_for_this_data = binlist + (hist2_shift_to_right * n)
                        data_dict[data_names_list[n]] = hist_data_AAIMON_mean
                        # use numpy to create a histogram
                        freq_counts, bin_array = np.histogram(hist_data_AAIMON_mean, bins=binlist_for_this_data)
                        # check if there proteins in the final bin (should be 30, and for mult proteins have a freq of 0)
                        if freq_counts[-1] != 0:
                            logging.warning(
                                'for %s, the final bin (%s) has %i entries' % (KW, bin_array[-1], freq_counts[-1]))
                        # normalize the frequency counts
                        freq_counts_normalised = freq_counts / freq_counts.max()  # use numpy to create a histogram
                        # assuming all of the bins are exactly the same size, make the width of the column equal to 70% of each bin
                        col_width = float('%0.3f' % (col_width_value * (bin_array_I[1] - bin_array_I[0])))
                        # when align='center', the central point of the bar in the x-axis is simply the middle of the bins ((bin_0-bin_1)/2, etc)
                        centre_of_bar_in_x_axis = (bin_array[:-2] + bin_array[1:-1]) / 2
                        # add the final bin, which is physically located just after the last regular bin but represents all higher values
                        bar_width = centre_of_bar_in_x_axis[3] - centre_of_bar_in_x_axis[2]
                        centre_of_bar_in_x_axis = np.append(centre_of_bar_in_x_axis,
                                                            centre_of_bar_in_x_axis[-1] + bar_width)
                        barcontainer_AAIMON_mean = axarr[row_nr, col_nr].bar(left=centre_of_bar_in_x_axis,
                                                                             height=freq_counts_normalised,
                                                                             align='center', width=col_width, color=color,
                                                                             alpha=alpha, linewidth=0.1)
                    # obtain data from dictionary
                    data1 = data_dict[data_names_list[0]]
                    data2 = data_dict[data_names_list[1]]
                    # run ttest, assuming unequal variance (http://stackoverflow.com/questions/22611446/perform-2-sample-t-test)
                    t, p = ttest_ind(data1, data2, equal_var=False)
                    # add the ttest results to a dictionary
                    dict_ttest_pvalues[KW] = p
                    # determine symbol describing stat significance
                    signif_symbol = utils.get_signif_symbol(p)

                    # label the x-axis
                    axarr[row_nr, col_nr].set_xlabel('average conservation ratio (membranous over nonmembranous)',
                                                     fontsize=fontsize)
                    # move the x-axis label closer to the x-axis
                    axarr[row_nr, col_nr].xaxis.set_label_coords(0.45, -0.085)
                    xlim_min = 0.8
                    xlim_max = 1.5
                    axarr[row_nr, col_nr].set_xlim(xlim_min, xlim_max)
                    # set x-axis ticks
                    # use the slide selection to select every second item in the list as an xtick(axis label)
                    axarr[row_nr, col_nr].set_xticks([float('%0.1f' % c) for c in centre_of_bar_in_x_axis[::3]])
                    axarr[row_nr, col_nr].set_ylabel('freq', rotation='vertical', fontsize=fontsize)
                    # change axis font size
                    axarr[row_nr, col_nr].tick_params(labelsize=fontsize)
                    # add figure number to top left of subplot
                    axarr[row_nr, col_nr].annotate(s=str(Fig_Nr) + '.', xy=(0.04, 0.9), fontsize=5, xytext=None,
                                                   xycoords='axes fraction', alpha=0.75)
                    # add figure title to top left of subplot
                    axarr[row_nr, col_nr].annotate(s=KW + ' nonGPCR', xy=(0.1, 0.9), fontsize=5, xytext=None,
                                                   xycoords='axes fraction', alpha=0.75)
                    # add a legend describing which hist contains the KW, or doesn't contain the KW
                    legend_obj = axarr[row_nr, col_nr].legend(['containing keyword', 'without keyword'], loc='upper right',
                                                              fontsize=fontsize)
                    # add ttest result to graph
                    axarr[row_nr, col_nr].annotate(s='p = %0.2g (%s)' % (p, signif_symbol), xy=(0.70, 0.75), fontsize=5,
                                                   xytext=None, xycoords='axes fraction', alpha=0.75)
                    # improve ggplot style for a canvas (fig) with 4 figures (plots)
                    utils.improve_ggplot_for_4_plots(axarr, row_nr, col_nr, backgroundcolour, legend_obj)
                    # annotate the correlated keywords to the graph
                    corr_KW = [xc[0:16] for xc in list(top_corr_KW_dict[KW].index)]
                    corr_KW_value = [int(xv) for xv in list(top_corr_KW_dict[KW])]
                    # add "correlated keywords" title
                    axarr[row_nr, col_nr].annotate(s='correlated keywords', xy=(0.7, 0.65), fontsize=4,
                                                   xycoords='axes fraction', alpha=0.75)
                    for ann_num, cKW in enumerate(corr_KW):
                        # define y position for annotation
                        ypos_ann = 0.6 - (ann_num * 0.05)
                        # define KW name position
                        ann_pos_xy = (0.7, ypos_ann)
                        # determine position of KW correlation value
                        val_pos_xy = (0.95, ypos_ann)
                        # annotate the KW name
                        axarr[row_nr, col_nr].annotate(s=cKW, xy=ann_pos_xy, fontsize=4, xycoords='axes fraction',
                                                       alpha=0.75)
                        # annotate correlation value
                        axarr[row_nr, col_nr].annotate(s=str(corr_KW_value[ann_num]), xy=val_pos_xy, fontsize=4,
                                                       xycoords='axes fraction', alpha=0.75)
                    # calculate the odds ratio (TMD/rest conservation for proteins with KW) / (TMD/rest conservation for proteins without KW)
                    odds_ratio_KW_over_nonKW = '%0.2f' % (data1.mean() / data2.mean())
                    # add odds ratio to figure
                    axarr[row_nr, col_nr].annotate(s='odds ratio = ' + odds_ratio_KW_over_nonKW, xy=(0.7, 0.70), fontsize=5,
                                                   xytext=None, xycoords='axes fraction', alpha=0.75)

                    # save the figure if necessary (i.e., if the maximum number of plots per figure has been obtained)
                    utils.savefig_if_necessary(savefig, fig, fig_nr, base_filepath=pathdict["base_filename_summaries"])

            '''
            Enzymes vs NonEnzymes
            '''
            Fig_Nr = Fig_Nr + 1
            logging.info("Enzymes vs NonEnzymes Fig_Nr = %i" % Fig_Nr)
            # add non-functional "if" line to aid document navigation in some IDEs (e.g. Spyder)
            if fig_title == "Enzymes vs NonEnzymes":
                pass
            sys.stdout.write(str(Fig_Nr) + ', ')
            list_enzyme_KW = ['Transferase', 'Hydrolase', 'Glycosyltransferase', 'Protease', 'Kinase', 'Oxidoreductase',
                              'Metalloprotease',
                              'Serine protease', 'Protein phosphatase', 'Ligase', 'Acyltransferase',
                              'Serine/threonine-protein kinase', 'Glycosidase',
                              'Aminopeptidase', 'Isomerase', 'Methyltransferase', 'Carboxypeptidase', 'Hydroxylation',
                              'Aspartyl protease',
                              'Serine esterase', 'Lipid biosynthesis', 'GPI-anchor biosynthesis', 'Steroid biosynthesis',
                              'Melanin biosynthesis',
                              'Thyroid hormones biosynthesis', 'Phospholipid biosynthesis', 'Sterol biosynthesis',
                              'Glutathione biosynthesis',
                              'Cholesterol biosynthesis', 'Fatty acid biosynthesis', 'Prostaglandin biosynthesis',
                              'cGMP biosynthesis', 'Leukotriene biosynthesis',
                              'Catecholamine biosynthesis', 'Lipid metabolism', 'Carbohydrate metabolism',
                              'Steroid metabolism', 'Sterol metabolism',
                              'Sphingolipid metabolism', 'Cholesterol metabolism', 'Fatty acid metabolism',
                              'Phospholipid metabolism', 'Catecholamine metabolism', 'Prostaglandin metabolism',
                              'Glycogen metabolism', 'Fucose metabolism']

            df['enzyme'] = df['uniprot_KW'].apply(utils.KW_list_contains_any_desired_KW, args=(list_enzyme_KW,))

            newfig, savefig, fig_nr, plot_nr_in_fig, row_nr, col_nr = dict_organising_subplots[Fig_Nr]
            # if a new figure should be created (either because the orig is full, or the last TMD is analysed)
            fig, axarr = utils.create_new_fig_if_necessary(newfig, fig, axarr, nrows_in_each_fig, ncols_in_each_fig,
                                                           dpi=300)

            # create dataframes with entries containing, or not containing the particular keyword
            df_enzyme = df.loc[df['enzyme'] == True]
            df_nonenzyme = df.loc[df['enzyme'] == False]
            df_list_KW = [df_enzyme, df_nonenzyme]

            data_dict = {}
            data_names_list = ['enzyme', 'nonenzyme']

            for n, dfK in enumerate(df_list_KW):
                data_column = 'AAIMON_ratio_mean_all_TMDs'
                color = colourlist_greys[n]
                alpha = 1.0
                col_width_value = 0.95
                hist2_shift_to_right = 0.01
                backgroundcolour = '0.95'
                hist_data_AAIMON_mean = np.array(dfK[data_column].dropna())
                # add the data to the dictionary for later analysis, named according to the list above
                data_dict[data_names_list[n]] = hist_data_AAIMON_mean
                # move the second histogram bins slightly to the right, so that the bars do not overlap
                binlist_for_this_data = binlist + (hist2_shift_to_right * n)
                # use numpy to create a histogram
                freq_counts, bin_array = np.histogram(hist_data_AAIMON_mean, bins=binlist_for_this_data)
                # check if there proteins in the final bin (should be 30, and for mult proteins have a freq of 0)
                if freq_counts[-1] != 0:
                    logging.warning('for %s, the final bin (%s) has %i entries' % (KW, bin_array[-1], freq_counts[-1]))
                # normalize the frequency counts
                freq_counts_normalised = freq_counts / freq_counts.max()  # use numpy to create a histogram
                # assuming all of the bins are exactly the same size, make the width of the column equal to 70% of each bin
                col_width = float('%0.3f' % (col_width_value * (bin_array_I[1] - bin_array_I[0])))
                # when align='center', the central point of the bar in the x-axis is simply the middle of the bins ((bin_0-bin_1)/2, etc)
                centre_of_bar_in_x_axis = (bin_array[:-2] + bin_array[1:-1]) / 2
                # add the final bin, which is physically located just after the last regular bin but represents all higher values
                bar_width = centre_of_bar_in_x_axis[3] - centre_of_bar_in_x_axis[2]
                centre_of_bar_in_x_axis = np.append(centre_of_bar_in_x_axis, centre_of_bar_in_x_axis[-1] + bar_width)
                barcontainer_AAIMON_mean = axarr[row_nr, col_nr].bar(left=centre_of_bar_in_x_axis,
                                                                     height=freq_counts_normalised,
                                                                     align='center', width=col_width, color=color,
                                                                     alpha=alpha, linewidth=0.1)
            # obtain data from dictionary
            data1 = data_dict[data_names_list[0]]
            data2 = data_dict[data_names_list[1]]
            # run ttest, assuming unequal variance (http://stackoverflow.com/questions/22611446/perform-2-sample-t-test)
            t, p = ttest_ind(data1, data2, equal_var=False)
            # determine symbol describing stat significance
            signif_symbol = utils.get_signif_symbol(p)

            # label the x-axis for each plot, based on the TMD
            axarr[row_nr, col_nr].set_xlabel('average conservation ratio (membranous over nonmembranous)',
                                             fontsize=fontsize)
            # move the x-axis label closer to the x-axis
            axarr[row_nr, col_nr].xaxis.set_label_coords(0.45, -0.085)
            xlim_min = 0.8
            xlim_max = 1.5
            axarr[row_nr, col_nr].set_xlim(xlim_min, xlim_max)
            # set x-axis ticks
            # use the slide selection to select every second item in the list as an xtick(axis label)
            axarr[row_nr, col_nr].set_xticks([float('%0.1f' % c) for c in centre_of_bar_in_x_axis[::3]])
            axarr[row_nr, col_nr].set_ylabel('freq', rotation='vertical', fontsize=fontsize)
            # change axis font size
            axarr[row_nr, col_nr].tick_params(labelsize=fontsize)
            # add figure number to top left of subplot
            axarr[row_nr, col_nr].annotate(s=str(Fig_Nr) + '.', xy=(0.04, 0.9), fontsize=5, xytext=None,
                                           xycoords='axes fraction', alpha=0.75)
            # add figure title to top left of subplot
            axarr[row_nr, col_nr].annotate(s='. Enzymes', xy=(0.1, 0.9), fontsize=5, xytext=None, xycoords='axes fraction',
                                           alpha=0.75)
            # add a legend describing which hist contains the KW, or doesn't contain the KW
            legend_obj = axarr[row_nr, col_nr].legend(data_names_list, loc='upper right', fontsize=fontsize)
            # add ttest result to graph
            axarr[row_nr, col_nr].annotate(s='p = %0.2g (%s)' % (p, signif_symbol), xy=(0.70, 0.75), fontsize=5,
                                           xytext=None, xycoords='axes fraction', alpha=0.75)
            # improve ggplot style for a canvas (fig) with 4 figures (plots)
            utils.improve_ggplot_for_4_plots(axarr, row_nr, col_nr, backgroundcolour, legend_obj)
            # calculate the odds ratio (TMD/rest conservation for proteins with KW) / (TMD/rest conservation for proteins without KW)
            odds_ratio_KW_over_nonKW = '%0.2f' % (data1.mean() / data2.mean())
            # add odds ratio to figure
            axarr[row_nr, col_nr].annotate(s='odds ratio = ' + odds_ratio_KW_over_nonKW, xy=(0.7, 0.70), fontsize=5,
                                           xytext=None, xycoords='axes fraction', alpha=0.75)
            # calculate the odds ratio (TMD/rest conservation for proteins with KW) / (TMD/rest conservation for proteins without KW)
            odds_ratio_KW_over_nonKW = '%0.2f' % (data1.mean() / data2.mean())
            # add odds ratio to figure
            axarr[row_nr, col_nr].annotate(s='odds ratio = ' + odds_ratio_KW_over_nonKW, xy=(0.7, 0.70), fontsize=5,
                                           xytext=None, xycoords='axes fraction', alpha=0.75)
            # improve ggplot style for a canvas (fig) with 4 figures (plots)
            utils.improve_ggplot_for_4_plots(axarr, row_nr, col_nr, backgroundcolour, legend_obj)
            # save the figure if necessary (i.e., if the maximum number of plots per figure has been obtained)
            utils.savefig_if_necessary(savefig, fig, fig_nr, base_filepath=pathdict["base_filename_summaries"])

            if dataset_contains_multipass_prots:
                '''
                Enzymes vs NonEnzymes, Non-GPCR ONLY
                '''
                Fig_Nr = Fig_Nr + 1
                # add non-functional "if" line to aid document navigation in some IDEs (e.g. Spyder)
                if fig_title == "Enzymes vs NonEnzymes, Non-GPCR ONLY":
                    pass
                sys.stdout.write(str(Fig_Nr) + ', ')

                newfig, savefig, fig_nr, plot_nr_in_fig, row_nr, col_nr = dict_organising_subplots[Fig_Nr]
                logging.info("Enzymes vs NonEnzymes fig_nr = %i" % fig_nr)
                # if a new figure should be created (either because the orig is full, or the last TMD is analysed)
                fig, axarr = utils.create_new_fig_if_necessary(newfig, fig, axarr, nrows_in_each_fig, ncols_in_each_fig,
                                                               dpi=300)

                # create nonGPCR dataframe view
                df_nonGPCR = df.loc[df['G-protein coupled receptor'] == False]
                # create enzyme and nonenzyme dataframe views
                df_enzyme = df_nonGPCR.loc[df['enzyme'] == True]
                df_nonenzyme = df_nonGPCR.loc[df['enzyme'] == False]

                df_list_KW = [df_enzyme, df_nonenzyme]

                data_dict = {}
                data_names_list = ['enzyme', 'nonenzyme']

                for n, dfK in enumerate(df_list_KW):
                    data_column = 'AAIMON_ratio_mean_all_TMDs'
                    color = colourlist_greys[n]
                    alpha = 1.0
                    col_width_value = 0.95
                    hist2_shift_to_right = 0.01
                    backgroundcolour = '0.95'
                    hist_data_AAIMON_mean = np.array(dfK[data_column].dropna())
                    # add the data to the dictionary for later analysis, named according to the list above
                    data_dict[data_names_list[n]] = hist_data_AAIMON_mean
                    # move the second histogram bins slightly to the right, so that the bars do not overlap
                    binlist_for_this_data = binlist + (hist2_shift_to_right * n)
                    # use numpy to create a histogram
                    freq_counts, bin_array = np.histogram(hist_data_AAIMON_mean, bins=binlist_for_this_data)
                    # check if there proteins in the final bin (should be 30, and for mult proteins have a freq of 0)
                    if freq_counts[-1] != 0:
                        logging.warning('for %s, the final bin (%s) has %i entries' % (KW, bin_array[-1], freq_counts[-1]))
                    # normalize the frequency counts
                    freq_counts_normalised = freq_counts / freq_counts.max()  # use numpy to create a histogram
                    # assuming all of the bins are exactly the same size, make the width of the column equal to 70% of each bin
                    col_width = float('%0.3f' % (col_width_value * (bin_array_I[1] - bin_array_I[0])))
                    # when align='center', the central point of the bar in the x-axis is simply the middle of the bins ((bin_0-bin_1)/2, etc)
                    centre_of_bar_in_x_axis = (bin_array[:-2] + bin_array[1:-1]) / 2
                    # add the final bin, which is physically located just after the last regular bin but represents all higher values
                    bar_width = centre_of_bar_in_x_axis[3] - centre_of_bar_in_x_axis[2]
                    centre_of_bar_in_x_axis = np.append(centre_of_bar_in_x_axis, centre_of_bar_in_x_axis[-1] + bar_width)
                    barcontainer_AAIMON_mean = axarr[row_nr, col_nr].bar(left=centre_of_bar_in_x_axis,
                                                                         height=freq_counts_normalised,
                                                                         align='center', width=col_width, color=color,
                                                                         alpha=alpha, linewidth=0.1)
                # obtain data from dictionary
                data1 = data_dict[data_names_list[0]]
                data2 = data_dict[data_names_list[1]]
                # run ttest, assuming unequal variance (http://stackoverflow.com/questions/22611446/perform-2-sample-t-test)
                t, p = ttest_ind(data1, data2, equal_var=False)
                # determine symbol describing stat significance
                signif_symbol = utils.get_signif_symbol(p)

                # label the x-axis for each plot, based on the TMD
                axarr[row_nr, col_nr].set_xlabel('average conservation ratio (membranous over nonmembranous)',
                                                 fontsize=fontsize)
                # move the x-axis label closer to the x-axis
                axarr[row_nr, col_nr].xaxis.set_label_coords(0.45, -0.085)
                xlim_min = 0.8
                xlim_max = 1.5
                axarr[row_nr, col_nr].set_xlim(xlim_min, xlim_max)
                # set x-axis ticks
                # use the slide selection to select every second item in the list as an xtick(axis label)
                axarr[row_nr, col_nr].set_xticks([float('%0.1f' % c) for c in centre_of_bar_in_x_axis[::3]])
                axarr[row_nr, col_nr].set_ylabel('freq', rotation='vertical', fontsize=fontsize)
                # change axis font size
                axarr[row_nr, col_nr].tick_params(labelsize=fontsize)
                # add figure number to top left of subplot
                axarr[row_nr, col_nr].annotate(s=str(Fig_Nr) + '.', xy=(0.04, 0.9), fontsize=5, xytext=None,
                                               xycoords='axes fraction', alpha=0.75)
                # add figure title to top left of subplot
                axarr[row_nr, col_nr].annotate(s='. Enzymes (all NonGPCR)', xy=(0.1, 0.9), fontsize=5, xytext=None,
                                               xycoords='axes fraction', alpha=0.75)
                # add a legend describing which hist contains the KW, or doesn't contain the KW
                legend_obj = axarr[row_nr, col_nr].legend(data_names_list, loc='upper right', fontsize=fontsize)
                # add ttest result to graph
                axarr[row_nr, col_nr].annotate(s='p = %0.2g (%s)' % (p, signif_symbol), xy=(0.70, 0.75), fontsize=5,
                                               xytext=None, xycoords='axes fraction', alpha=0.75)
                # improve ggplot style for a canvas (fig) with 4 figures (plots)
                utils.improve_ggplot_for_4_plots(axarr, row_nr, col_nr, backgroundcolour, legend_obj)
                # calculate the odds ratio (TMD/rest conservation for proteins with KW) / (TMD/rest conservation for proteins without KW)
                odds_ratio_KW_over_nonKW = '%0.2f' % (data1.mean() / data2.mean())
                # add odds ratio to figure
                axarr[row_nr, col_nr].annotate(s='odds ratio = ' + odds_ratio_KW_over_nonKW, xy=(0.7, 0.70), fontsize=5,
                                               xytext=None, xycoords='axes fraction', alpha=0.75)
                # calculate the odds ratio (TMD/rest conservation for proteins with KW) / (TMD/rest conservation for proteins without KW)
                odds_ratio_KW_over_nonKW = '%0.2f' % (data1.mean() / data2.mean())
                # add odds ratio to figure
                axarr[row_nr, col_nr].annotate(s='odds ratio = ' + odds_ratio_KW_over_nonKW, xy=(0.7, 0.70), fontsize=5,
                                               xytext=None, xycoords='axes fraction', alpha=0.75)

                # save the figure if necessary (i.e., if the maximum number of plots per figure has been obtained)
                utils.savefig_if_necessary(savefig, fig, fig_nr, base_filepath=pathdict["base_filename_summaries"])
            '''
            Histograms split by keywords (NONENZYME ONLY)
            '''
            starting_Fig_Nr = Fig_Nr + 1
            # add non-functional "if" line to aid document navigation in some IDEs (e.g. Spyder)
            if fig_title == "Hist keywords (NONENZYME ONLY)":
                pass
            # convert the KW stringlist to a python list, if it hasn't been done already
            if isinstance(df['uniprot_KW'][0], str):
                df_nonenzyme['uniprot_KW'] = df_nonenzyme['uniprot_KW'].apply(lambda x: ast.literal_eval(x))
            # join all keywords together into a large list
            nested_list_all_KW = list(itertools.chain(*list(df_nonenzyme['uniprot_KW'])))
            # convert list to pandas series
            all_KW_series = pd.Series(nested_list_all_KW)
            # obtain series of major keywords
            KW_counts = all_KW_series.value_counts()
            # exclude keywords with less than 50 applicable proteins
            KW_counts_major = KW_counts[KW_counts > 50]
            # create a list of keywords to be ignored
            list_ignored_KW = ['Transmembrane', 'Complete proteome', 'Reference proteome', 'Membrane',
                               'Transmembrane helix', 'Cell membrane', 'Repeat',
                               'Alternative splicing', 'Sodium', 'Potassium', 'Direct protein sequencing']
            # remove undesired keywords from the dataframe (Transmembrane, etc.)
            for KW in list_ignored_KW:
                if KW in KW_counts_major.index:
                    KW_counts_major = KW_counts_major.drop(KW)

            # prepare colour lists
            colour_lists = utils.create_colour_lists()
            TUM_colours_list_with_greys = colour_lists['TUM_colours_list_with_greys']
            colourlist_greys = [(0.6, 0.7764705882352941, 0.9058823529411765), 'None']

            # #iterate over the dataframe. Note that acc = uniprot accession here.
            # linspace_binlist = np.linspace(s["mp_smallest_bin"],
            #                                s["mp_largest_bin"],
            #                                s["mp_number_of_bins"])
            #
            # #add 30 as the last bin, to make sure 100% of the data is added to the histogram, including major outliers
            # binlist = np.append(linspace_binlist, s["mp_final_highest_bin"]))

            dict_ttest_pvalues = {}

            # create a new column showing whether the protein contains that keyword
            df[KW] = df['uniprot_KW'].apply(lambda x: KW in x)
            # slice the dataframe to show only non-enzymes
            df_nonenzyme = df.loc[df['enzyme'] == False]

            for m, KW in enumerate(KW_counts_major.index):
                Fig_Nr = starting_Fig_Nr + m
                sys.stdout.write(str(Fig_Nr) + ', ')
                newfig, savefig, fig_nr, plot_nr_in_fig, row_nr, col_nr = dict_organising_subplots[Fig_Nr]
                # if a new figure should be created (either because the orig is full, or the last TMD is analysed)
                fig, axarr = utils.create_new_fig_if_necessary(newfig, fig, axarr, nrows_in_each_fig, ncols_in_each_fig,
                                                               dpi=300)
                # create dataframes with entries containing, or not containing the particular keyword
                df_KW = df_nonenzyme.loc[df_nonenzyme[KW] == True]
                df_nonKW = df_nonenzyme.loc[df_nonenzyme[KW] == False]
                # create list of dataframes for iteration
                df_list_KW = [df_KW, df_nonKW]

                data_dict = {}
                data_names_list = ['containing keyword', 'without keyword']

                for n, dfK in enumerate(df_list_KW):
                    data_column = 'AAIMON_ratio_mean_all_TMDs'
                    # color = TUM_colours_list_with_greys[n + 1]
                    color = colourlist_greys[n]
                    alpha = 1.0
                    col_width_value = 0.95
                    hist2_shift_to_right = 0.01
                    backgroundcolour = '0.95'
                    hist_data_AAIMON_mean = np.array(dfK[data_column].dropna())
                    # add the data to the dictionary for later analysis, named according to the list above
                    data_dict[data_names_list[n]] = hist_data_AAIMON_mean
                    # move the second histogram bins slightly to the right, so that the bars do not overlap
                    binlist_for_this_data = binlist + (hist2_shift_to_right * n)
                    # use numpy to create a histogram
                    freq_counts, bin_array = np.histogram(hist_data_AAIMON_mean, bins=binlist_for_this_data)
                    # normalize the frequency counts
                    freq_counts_normalised = freq_counts / freq_counts.max()  # use numpy to create a histogram
                    # assuming all of the bins are exactly the same size, make the width of the column equal to a percentage of each bin
                    dist_between_bins = (bin_array_I[1] - bin_array_I[0])
                    col_width = float('%0.3f' % (col_width_value * dist_between_bins))
                    # when align='center', the central point of the bar in the x-axis is simply the middle of the bins ((bin_0-bin_1)/2, etc)
                    centre_of_bar_in_x_axis = (bin_array[:-2] + bin_array[1:-1]) / 2
                    # the last bin is open, and holds all larger datapoints than specified in the binlist. it needs to be added to the list of bars
                    centre_of_bar_in_x_axis = np.append(centre_of_bar_in_x_axis,
                                                        centre_of_bar_in_x_axis[-1] + dist_between_bins)
                    # print a warning if the last bar is not actually empty
                    if freq_counts[-1] != 0:
                        logging.warning(
                            'for %s, the final bin (%s) has %i entries. The settings "largest_bin" may need to be adjusted.' % (
                            KW, bin_array[-1], freq_counts[-1]))
                    # plot the data as a bar-chart
                    barcontainer_AAIMON_mean = axarr[row_nr, col_nr].bar(left=centre_of_bar_in_x_axis,
                                                                         height=freq_counts_normalised,
                                                                         align='center', width=col_width, color=color,
                                                                         alpha=alpha, linewidth=0.1)
                # obtain data from dictionary
                data1 = data_dict[data_names_list[0]]
                data2 = data_dict[data_names_list[1]]
                # run ttest, assuming unequal variance (http://stackoverflow.com/questions/22611446/perform-2-sample-t-test)
                t, p = ttest_ind(data1, data2, equal_var=False)
                # add the ttest results to a dictionary
                dict_ttest_pvalues[KW] = p
                # determine symbol describing stat significance
                signif_symbol = utils.get_signif_symbol(p)

                # label the x-axis for each plot, based on the TMD
                axarr[row_nr, col_nr].set_xlabel('average conservation ratio (membranous over nonmembranous)',
                                                 fontsize=fontsize)
                # move the x-axis label closer to the x-axis
                axarr[row_nr, col_nr].xaxis.set_label_coords(0.45, -0.085)
                # alter x-axis limits. Note that if data outside these limits the axis will be automatically adjusted)
                xlim_min = 0.5
                xlim_max = 1.6
                axarr[row_nr, col_nr].set_xlim(xlim_min, xlim_max)
                # set x-axis ticks
                # use the slide selection to select every second item in the list as an xtick(axis label)
                axarr[row_nr, col_nr].set_xticks([float('%0.1f' % c) for c in centre_of_bar_in_x_axis[::3]])
                axarr[row_nr, col_nr].set_ylabel('freq', rotation='vertical', fontsize=fontsize)
                # change axis font size
                axarr[row_nr, col_nr].tick_params(labelsize=fontsize)
                # remove background grid
                axarr[row_nr, col_nr].grid(False)
                # add figure number to top left of subplot
                # axarr[row_nr, col_nr].annotate(s=str(Fig_Nr) + '.', xy = (0.04,0.9), fontsize=5, xytext=None, xycoords='axes fraction', alpha=0.75)
                # add figure title to top left of subplot
                axarr[row_nr, col_nr].annotate(s=str(Fig_Nr) + '.' + KW + ' nonenzyme', xy=(0.1, 0.9), fontsize=5,
                                               xytext=None, xycoords='axes fraction', alpha=0.75)
                # add a legend describing which hist contains the KW, or doesn't contain the KW
                legend_obj = axarr[row_nr, col_nr].legend(['containing keyword', 'without keyword'], loc='upper right',
                                                          fontsize=fontsize, )
                # add ttest result to graph
                axarr[row_nr, col_nr].annotate(s='p = %0.2g (%s)' % (p, signif_symbol), xy=(0.70, 0.75), fontsize=5,
                                               xytext=None, xycoords='axes fraction', alpha=0.75)
                # improve ggplot style for a canvas (fig) with 4 figures (plots)
                utils.improve_ggplot_for_4_plots(axarr, row_nr, col_nr, backgroundcolour, legend_obj)
                # annotate the correlated keywords to the graph
                corr_KW = [xc[0:16] for xc in list(top_corr_KW_dict[KW].index)]
                corr_KW_value = [int(xv) for xv in list(top_corr_KW_dict[KW])]
                # add "correlated keywords" title
                axarr[row_nr, col_nr].annotate(s='correlated keywords', xy=(0.7, 0.65), fontsize=4,
                                               xycoords='axes fraction', alpha=0.75)
                for ann_num, cKW in enumerate(corr_KW):
                    # define y position for annotation
                    ypos_ann = 0.6 - (ann_num * 0.05)
                    # define KW name position
                    ann_pos_xy = (0.7, ypos_ann)
                    # determine position of KW correlation value
                    val_pos_xy = (0.95, ypos_ann)
                    # annotate the KW name
                    axarr[row_nr, col_nr].annotate(s=cKW, xy=ann_pos_xy, fontsize=4, xycoords='axes fraction', alpha=0.75)
                    # annotate correlation value
                    axarr[row_nr, col_nr].annotate(s=str(corr_KW_value[ann_num]), xy=val_pos_xy, fontsize=4,
                                                   xycoords='axes fraction', alpha=0.75)
                # calculate the odds ratio (TMD/rest conservation for proteins with KW) / (TMD/rest conservation for proteins without KW)
                odds_ratio_KW_over_nonKW = '%0.2f' % (data1.mean() / data2.mean())
                # add odds ratio to figure
                axarr[row_nr, col_nr].annotate(s='odds ratio = ' + odds_ratio_KW_over_nonKW, xy=(0.7, 0.70), fontsize=5,
                                               xytext=None, xycoords='axes fraction', alpha=0.75)

                # save the figure if necessary (i.e., if the maximum number of plots per figure has been obtained)
                utils.savefig_if_necessary(savefig, fig, fig_nr, base_filepath=pathdict["base_filename_summaries"])

            # only create the following for multi-pass proteins (at least 2 proteins in lost with 7 TMDs)
            if dataset_contains_multipass_prots:
                '''
                Histograms split by keywords (NONENZYME AND nonGPCR)
                '''
                starting_Fig_Nr = Fig_Nr + 1
                # add non-functional "if" line to aid document navigation in some IDEs (e.g. Spyder)
                if fig_title == "Hist keywords (NONENZYME AND nonGPCR)":
                    pass

                df_nonenzyme_nonGPCR = df_nonenzyme.loc[df_nonenzyme['G-protein coupled receptor'] == False]

                # join all keywords together into a large list
                nested_list_all_KW = list(itertools.chain(*list(df_nonenzyme_nonGPCR['uniprot_KW'])))
                # convert list to pandas series
                all_KW_series = pd.Series(nested_list_all_KW)
                # obtain series of major keywords
                KW_counts = all_KW_series.value_counts()
                # exclude keywords with less than 50 applicable proteins
                KW_counts_major = KW_counts[KW_counts > 50]

                # remove undesired keywords from the dataframe (Transmembrane, etc.)
                for KW in list_ignored_KW:
                    if KW in KW_counts_major.index:
                        KW_counts_major = KW_counts_major.drop(KW)

                dict_ttest_pvalues = {}

                for m, KW in enumerate(KW_counts_major.index):
                    Fig_Nr = starting_Fig_Nr + m
                    sys.stdout.write(str(Fig_Nr) + ', ')
                    newfig, savefig, fig_nr, plot_nr_in_fig, row_nr, col_nr = dict_organising_subplots[Fig_Nr]
                    # if a new figure should be created (either because the orig is full, or the last TMD is analysed)
                    fig, axarr = utils.create_new_fig_if_necessary(newfig, fig, axarr, nrows_in_each_fig, ncols_in_each_fig,
                                                                   dpi=300)
                    # create dataframes with entries containing, or not containing the particular keyword
                    df_KW = df_nonenzyme_nonGPCR.loc[df_nonenzyme_nonGPCR[KW] == True]
                    df_nonKW = df_nonenzyme_nonGPCR.loc[df_nonenzyme_nonGPCR[KW] == False]
                    df_list_KW = [df_KW, df_nonKW]

                    data_dict = {}
                    data_names_list = ['containing keyword', 'without keyword']

                    for n, dfK in enumerate(df_list_KW):
                        data_column = 'AAIMON_ratio_mean_all_TMDs'
                        # color = TUM_colours_list_with_greys[n + 1]
                        color = colourlist_greys[n]
                        alpha = 1.0
                        col_width_value = 0.95
                        hist2_shift_to_right = 0.01
                        backgroundcolour = '0.95'
                        hist_data_AAIMON_mean = np.array(dfK[data_column].dropna())
                        # add the data to the dictionary for later analysis, named according to the list above
                        data_dict[data_names_list[n]] = hist_data_AAIMON_mean
                        # move the second histogram bins slightly to the right, so that the bars do not overlap
                        binlist_for_this_data = binlist + (hist2_shift_to_right * n)
                        # use numpy to create a histogram
                        freq_counts, bin_array = np.histogram(hist_data_AAIMON_mean, bins=binlist_for_this_data)
                        # normalize the frequency counts
                        freq_counts_normalised = freq_counts / freq_counts.max()  # use numpy to create a histogram
                        # assuming all of the bins are exactly the same size, make the width of the column equal to a percentage of each bin
                        dist_between_bins = (bin_array_I[1] - bin_array_I[0])
                        col_width = float('%0.3f' % (col_width_value * dist_between_bins))
                        # when align='center', the central point of the bar in the x-axis is simply the middle of the bins ((bin_0-bin_1)/2, etc)
                        centre_of_bar_in_x_axis = (bin_array[:-2] + bin_array[1:-1]) / 2
                        # the last bin is open, and holds all larger datapoints than specified in the binlist. it needs to be added to the list of bars
                        centre_of_bar_in_x_axis = np.append(centre_of_bar_in_x_axis,
                                                            centre_of_bar_in_x_axis[-1] + dist_between_bins)
                        # print a warning if the last bar is not actually empty
                        if freq_counts[-1] != 0:
                            logging.warning(
                                'for %s, the final bin (%s) has %i entries. The settings "largest_bin" may need to be adjusted.' % (
                                KW, bin_array[-1], freq_counts[-1]))
                        # plot the data as a bar-chart
                        barcontainer_AAIMON_mean = axarr[row_nr, col_nr].bar(left=centre_of_bar_in_x_axis,
                                                                             height=freq_counts_normalised,
                                                                             align='center', width=col_width, color=color,
                                                                             alpha=alpha, linewidth=0.1)
                    # obtain data from dictionary
                    data1 = data_dict[data_names_list[0]]
                    data2 = data_dict[data_names_list[1]]
                    # run ttest, assuming unequal variance (http://stackoverflow.com/questions/22611446/perform-2-sample-t-test)
                    t, p = ttest_ind(data1, data2, equal_var=False)
                    # add the ttest results to a dictionary
                    dict_ttest_pvalues[KW] = p
                    # determine symbol describing stat significance
                    signif_symbol = utils.get_signif_symbol(p)

                    # label the x-axis for each plot, based on the TMD
                    axarr[row_nr, col_nr].set_xlabel('average conservation ratio (membranous over nonmembranous)',
                                                     fontsize=fontsize)
                    # move the x-axis label closer to the x-axis
                    axarr[row_nr, col_nr].xaxis.set_label_coords(0.45, -0.085)
                    # alter x-axis limits. Note that if data outside these limits the axis will be automatically adjusted)
                    xlim_min = 0.5
                    xlim_max = 1.6
                    axarr[row_nr, col_nr].set_xlim(xlim_min, xlim_max)
                    # set x-axis ticks
                    # use the slide selection to select every second item in the list as an xtick(axis label)
                    axarr[row_nr, col_nr].set_xticks([float('%0.1f' % c) for c in centre_of_bar_in_x_axis[::3]])
                    axarr[row_nr, col_nr].set_ylabel('freq', rotation='vertical', fontsize=fontsize)
                    # change axis font size
                    axarr[row_nr, col_nr].tick_params(labelsize=fontsize)
                    # remove background grid
                    axarr[row_nr, col_nr].grid(False)
                    # add figure number to top left of subplot
                    # axarr[row_nr, col_nr].annotate(s=str(Fig_Nr) + '.', xy = (0.04,0.9), fontsize=5, xytext=None, xycoords='axes fraction', alpha=0.75)
                    # add figure title to top left of subplot
                    axarr[row_nr, col_nr].annotate(s=str(Fig_Nr) + '.' + KW, xy=(0.1, 0.9), fontsize=5, xytext=None,
                                                   xycoords='axes fraction', alpha=0.75)
                    # add figure title to top left of subplot
                    axarr[row_nr, col_nr].annotate(s='nonenzyme, nonGPCR', xy=(0.1, 0.8), fontsize=5, xytext=None,
                                                   xycoords='axes fraction', alpha=0.75)

                    # add a legend describing which hist contains the KW, or doesn't contain the KW
                    legend_obj = axarr[row_nr, col_nr].legend(['containing keyword', 'without keyword'], loc='upper right',
                                                              fontsize=fontsize, )
                    # add ttest result to graph
                    axarr[row_nr, col_nr].annotate(s='p = %0.2g (%s)' % (p, signif_symbol), xy=(0.70, 0.75), fontsize=5,
                                                   xytext=None, xycoords='axes fraction', alpha=0.75)
                    # Remove top axes and right axes ticks
                    axarr[row_nr, col_nr].get_xaxis().tick_bottom()
                    axarr[row_nr, col_nr].get_yaxis().tick_left()
                    # chane the position of the axis ticklabels so they are closer to the axis
                    axarr[row_nr, col_nr].tick_params(direction='out', pad=0.4)
                    # change background colour of graph
                    axarr[row_nr, col_nr].set_axis_bgcolor(backgroundcolour)
                    # change background colour of legend box
                    legend_obj.get_frame().set_facecolor(backgroundcolour)
                    # annotate the correlated keywords to the graph
                    corr_KW = [xc[0:16] for xc in list(top_corr_KW_dict[KW].index)]
                    corr_KW_value = [int(xv) for xv in list(top_corr_KW_dict[KW])]
                    # add "correlated keywords" title
                    axarr[row_nr, col_nr].annotate(s='correlated keywords', xy=(0.7, 0.65), fontsize=4,
                                                   xycoords='axes fraction', alpha=0.75)
                    for ann_num, cKW in enumerate(corr_KW):
                        # define y position for annotation
                        ypos_ann = 0.6 - (ann_num * 0.05)
                        # define KW name position
                        ann_pos_xy = (0.7, ypos_ann)
                        # determine position of KW correlation value
                        val_pos_xy = (0.95, ypos_ann)
                        # annotate the KW name
                        axarr[row_nr, col_nr].annotate(s=cKW, xy=ann_pos_xy, fontsize=4, xycoords='axes fraction',
                                                       alpha=0.75)
                        # annotate correlation value
                        axarr[row_nr, col_nr].annotate(s=str(corr_KW_value[ann_num]), xy=val_pos_xy, fontsize=4,
                                                       xycoords='axes fraction', alpha=0.75)
                    # calculate the odds ratio (TMD/rest conservation for proteins with KW) / (TMD/rest conservation for proteins without KW)
                    odds_ratio_KW_over_nonKW = '%0.2f' % (data1.mean() / data2.mean())
                    # add odds ratio to figure
                    axarr[row_nr, col_nr].annotate(s='odds ratio = ' + odds_ratio_KW_over_nonKW, xy=(0.7, 0.70), fontsize=5,
                                                   xytext=None, xycoords='axes fraction', alpha=0.75)

                    # save the figure if necessary (i.e., if the maximum number of plots per figure has been obtained)
                    utils.savefig_if_necessary(savefig, fig, fig_nr, base_filepath=pathdict["base_filename_summaries"])

    '''
    FINAL - MAKE SURE LAST GRAPHS ARE SAVED
    '''
    # save the figure as it is
    savefig = True
    # save the figure if necessary (i.e., if the maximum number of plots per figure has been obtained)

    # add a column indicating that save_figures_describing_proteins_in_list has been run
    # df["save_figures_describing_proteins_in_list"] = True
    '''
    save the updated dataframe, containing the various extra columns used for the figure
    '''
    # df.to_csv(pathdict["list_summary_csv"], sep=",", quoting=csv.QUOTE_NONNUMERIC)
    logging.info('run_save_figures_describing_proteins_in_list is finished')
