import matplotlib.pyplot as plt
import numpy as np
import os


def save_hist_AAIMON_ratio_single_protein (fig_nr, fig, axarr, df_cr, s, TMD, binarray, zipout, row_nr, col_nr, fontsize, savefig, AAIMON_hist_path_prefix):
    """Save histogram showing the AAIMON ratio for homologues of a single TMD, for a single protein.

    Parameters
    ----------
    fig_nr : int
        Figure number (canvas number). Used to save Fig01, Fig02 etc, each with 4 plots (e.g. for 4 TMDs)
    fig : matplotlib.figure
        Figure (canvas) object, containing plots for AAIMON of up to 4 TMDs.
    axarr : array
        Array used to select the plots in the figure object. E.g. axarr[0,0] refers to the plot on the top left.
    df_cr : pd.DataFrame
        Dataframe with conservation ratios for a particular TMD (or region).
    s : dict
        Settings dictionary extracted from excel settings file.
    TMD : str
        String denoting transmembrane domain number (e.g. "TM01")
    binarray : np.ndarray
        List of bins used for the histogram.
    zipout : zipfile.Zipfile handle
        Handle for zipfile, open for writing.
    row_nr : int
        Row number in Figure canvas. E.g. when row_nr & col_nr are both 0, the axarr[row_nr, col_nr] is axarr[0, 0], which refers to the plot on the top left.
    col_nr : int
        Column number in Figure canvas. E.g. when row_nr & col_nr are both 0, the axarr[row_nr, col_nr] is axarr[0, 0], which refers to the plot on the top left.
    fontsize : int
        Fontsize in plots.
    savefig : bool
        Whether the figure/canvas needs to be saved. Since there are 4 plots in each figure, this is True for every 4th plot (4th fig_nr).
    AAIMON_hist_path_prefix : str
        Path and beginning of filename for histogram figure. (df.homol_base + '_AAIMON_hist')
    Returns
    -------

    """
    # use linspace to get a fixid number of points between tha min and the max for the histogram
    # set up evenly distributed bins between the chosen min and max
    # if possible, 1.0 should be in the centre of a bin, to catch cases where a lot of homologues have a ratio that approximates 1

    #with tarfile.open(df.loc[acc, 'output_tarfile_path'], mode='w:gz') as tar_out:
    #with zipfile.ZipFile(homol_cr_ratios_zip, mode="w", compression=zipfile.ZIP_DEFLATED) as zipout:
    # calculate ratio of Amino acid Identity Membranous Over Nonmembranous  (AAIMON ratio)
    #for TMD in list_of_TMDs:

    # following the general filters, filter to only analyse sequences with TMD identity above cutoff, and a nonTMD_perc_ident above zero ,to avoid a divide by zero error
    #RESET SO THAT IT TAKES IT FROM df_cr AGAIN
    #df_cr_filt_AAIMON = df_cr_filt.loc[df_cr['%s_perc_ident'%TMD] >= min_identity_of_TMD_initial_filter]

    # create numpy array of membranous over nonmembranous conservation ratios (identity)
    hist_data_I = np.array(df_cr['%s_AAIMON_ratio'%TMD].dropna())
    # use numpy to create a histogram
    freq_counts, bin_array = np.histogram(hist_data_I, bins=binarray)
    # assuming all of the bins are exactly the same size, make the width of the column equal to 70% of each bin
    col_width = float('%0.3f' % (0.7 * (bin_array[1] - bin_array[0])))
    # when align='center', the central point of the bar in the x-axis is simply the middle of the bins ((bin_0-bin_1)/2, etc)
    centre_of_bar_in_x_axis = (bin_array[:-2] + bin_array[1:-1]) / 2
    # add the final bin, which is physically located just after the last regular bin but represents all higher values
    centre_of_bar_in_x_axis = np.append(centre_of_bar_in_x_axis,
                                        centre_of_bar_in_x_axis[-1] +
                                        centre_of_bar_in_x_axis[0])
    barcontainer_I = axarr[row_nr, col_nr].bar(left=centre_of_bar_in_x_axis,
                                               height=freq_counts, align='center',
                                               width=col_width, color="#0489B1",
                                               alpha=0.5)  # edgecolor='black',

    # SHENGER TO INSERT TM01_AAIMON_ratio_n normalised stuff here
    hist_data_N = np.array(df_cr['%s_AAIMON_ratio_n' % TMD].dropna())
    # use numpy to create a histogram
    freq_counts, bin_array = np.histogram(hist_data_N, bins=binarray)
    # assuming all of the bins are exactly the same size, make the width of the column equal to 70% of each bin
    col_width = float('%0.3f' % (0.7 * (bin_array[1] - bin_array[0])))
    # when align='center', the central point of the bar in the x-axis is simply the middle of the bins ((bin_0-bin_1)/2, etc)
    centre_of_bar_in_x_axis = (bin_array[:-2] + bin_array[1:-1]) / 2
    # add the final bin, which is physically located just after the last regular bin but represents all higher values
    centre_of_bar_in_x_axis = np.append(centre_of_bar_in_x_axis,
                                        centre_of_bar_in_x_axis[-1] +
                                        centre_of_bar_in_x_axis[0])
    barcontainer_N = axarr[row_nr, col_nr].bar(left=centre_of_bar_in_x_axis,
                                               height=freq_counts, align='center',
                                               width=col_width, color="green",
                                               alpha=0.2)  # edgecolor='black',

    # create numpy array of membranous over nonmembranous conservation ratios (identity + similarity)
    hist_data_S = np.array(df_cr['%s_AASMON_ratio'%TMD].dropna())
    # use numpy to create a histogram
    freq_counts, bin_array = np.histogram(hist_data_S, bins=binarray)
    # create a line graph rather than a bar graph for the AASMON (ident + similarity)
    linecontainer_S = axarr[row_nr, col_nr].plot(centre_of_bar_in_x_axis, freq_counts,color="#0101DF", alpha=0.5)
    # other colours that are compatible with colourblind readers: #8A084B Dark red, #B45F04 deep orange, reddish purple #4B088A
    # http://html-color-codes.info/
    # label the x-axis for each plot, based on the TMD
    axarr[row_nr, col_nr].set_xlabel('%s conservation ratio (membranous over nonmembranous)'%TMD,fontsize=fontsize)

    if savefig:
        # take x-axis min from settings
        xlim_min = s["1p_smallest_bin"]
        # take x-axis max from settings
        xlim_max = s["1p_largest_bin"]
        # apply the following formatting changes to all plots in the figure
        for ax in axarr.flat:
            # set x-axis min
            ax.set_xlim(xlim_min, xlim_max)
            # set x-axis ticks
            # use the slide selection to select every second item in the list as an xtick(axis label)
            ax.set_xticks([float('%0.1f' % c) for c in centre_of_bar_in_x_axis[::3]])
            ax.set_ylabel('freq', rotation='vertical', fontsize=fontsize)
            # change axis font size
            ax.tick_params(labelsize=fontsize)
            # create legend?#http://stackoverflow.com/questions/9834452/how-do-i-make-a-single-legend-for-many-subplots-with-matplotlib
            ax.legend(['AASMON (identity + similarity)', 'AAIMON (identity)','AAIMON_n (correction)'],loc='upper right', fontsize=fontsize)
            # add background grid
            ax.grid(True, color='0.75', alpha=0.5)
        # automatically tighten the layout of plots in the figure
        fig.tight_layout()
        # save files
        fig.savefig(AAIMON_hist_path_prefix + '_%01d.png' % fig_nr,format='png', dpi=200)
        #fig.savefig(AAIMON_hist_path_prefix + '_%01d.pdf' % fig_nr,format='pdf')
        # close figure
        plt.close('all')
        # add to zipfile
        zipout.write(AAIMON_hist_path_prefix + '_%01d.png' % fig_nr,arcname=os.path.basename(AAIMON_hist_path_prefix) + '_%01d.png' % fig_nr)
        #zipout.write(AAIMON_hist_path_prefix + '_%01d.pdf' % fig_nr,arcname=os.path.basename(AAIMON_hist_path_prefix) + '_%01d.pdf' % fig_nr)
        # delete temporory files
        os.remove(AAIMON_hist_path_prefix + '_%01d.png' % fig_nr)
        #os.remove(AAIMON_hist_path_prefix + '_%01d.pdf' % fig_nr)
