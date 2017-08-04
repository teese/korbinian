import matplotlib.pyplot as plt
import numpy as np
import os
import math
# import debugging tools
from korbinian.utils import pr, pc, pn, aaa

def save_hist_AAIMON_single_protein (fig_nr, fig, axarr, df_cr, s, TMD, binarray, zipout, row_nr, col_nr, fontsize, savefig, AAIMON_hist_path_prefix):
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
    hist_data_I = np.array(df_cr['%s_AAIMON'%TMD].dropna())
    # use numpy to create a histogram
    freq_counts, bin_array = np.histogram(hist_data_I, bins=binarray)
    # assuming all of the bins are exactly the same size, make the width of the column equal to 70% of each bin
#    col_width = float('%0.3f' % (0.7 * (bin_array[1] - bin_array[0])))
    # when align='center', the central point of the bar in the x-axis is simply the middle of the bins ((bin_0-bin_1)/2, etc)
    centre_of_bar_in_x_axis = (bin_array[:-2] + bin_array[1:-1]) / 2
    # add the final bin, which is physically located just after the last regular bin but represents all higher values
    centre_of_bar_in_x_axis = np.append(centre_of_bar_in_x_axis,
                                        centre_of_bar_in_x_axis[-1] +
                                        centre_of_bar_in_x_axis[0])
    linecontainer_I = axarr[row_nr, col_nr].plot(centre_of_bar_in_x_axis, freq_counts,color="blue", alpha=0.5)
#    barcontainer_I = axarr[row_nr, col_nr].bar(left=centre_of_bar_in_x_axis,
#                                               height=freq_counts, align='center',
#                                               width=col_width, color="#0489B1",
#                                               alpha=0.5)  # edgecolor='black',

    # create linegraph for AAIMON after correction
    hist_data_N = np.array(df_cr['%s_AAIMON_n' % TMD].dropna())
    # use numpy to create a histogram
    freq_counts, bin_array = np.histogram(hist_data_N, bins=binarray)
    # assuming all of the bins are exactly the same size, make the width of the column equal to 70% of each bin
#    col_width = float('%0.3f' % (0.7 * (bin_array[1] - bin_array[0])))
    # when align='center', the central point of the bar in the x-axis is simply the middle of the bins ((bin_0-bin_1)/2, etc)
#    centre_of_bar_in_x_axis = (bin_array[:-2] + bin_array[1:-1]) / 2
    # add the final bin, which is physically located just after the last regular bin but represents all higher values
#    centre_of_bar_in_x_axis = np.append(centre_of_bar_in_x_axis,
#                                        centre_of_bar_in_x_axis[-1] +
#                                        centre_of_bar_in_x_axis[0])
    linecontainer_N = axarr[row_nr, col_nr].plot(centre_of_bar_in_x_axis, freq_counts,color="green", alpha=0.5)
#    linecontainer_N = axarr[row_nr, col_nr].bar(left=centre_of_bar_in_x_axis,
#                                               height=freq_counts, align='center',
#                                               width=col_width, color="green",
#                                               alpha=0.2)  # edgecolor='black',

    # create numpy array of membranous over nonmembranous conservation ratios (identity + similarity)
    hist_data_S = np.array(df_cr['%s_AASMON'%TMD].dropna())
    # use numpy to create a histogram
    freq_counts, bin_array = np.histogram(hist_data_S, bins=binarray)
    # create a line graph rather than a bar graph for the AASMON (ident + similarity)
    linecontainer_S = axarr[row_nr, col_nr].plot(centre_of_bar_in_x_axis, freq_counts,color="red", alpha=0.3)
    # other colours that are compatible with colourblind readers: #8A084B Dark red, #B45F04 deep orange, reddish purple #4B088A
    # http://html-color-codes.info/
    # label the x-axis for each plot, based on the TMD
    axarr[row_nr, col_nr].set_xlabel('%s conservation ratio (membranous over nonmembranous)'%TMD,fontsize=fontsize)

    if savefig:
        # take x-axis min from settings
        #xlim_min = s["1p_smallest_bin"]
        # take x-axis max from settings
        #xlim_max = s["1p_largest_bin"]
        # apply the following formatting changes to all plots in the figure
        for ax in axarr.flat:
            # set x-axis min
            #ax.set_xlim(xlim_min, xlim_max)
            ax.set_xlim(0)
            # set x-axis ticks
            # use the slide selection to select every second item in the list as an xtick(axis label)
            ax.set_xticks([float('%0.1f' % c) for c in centre_of_bar_in_x_axis[::3]])
            ax.set_ylabel('freq', rotation='vertical', fontsize=fontsize)
            # change axis font size
            ax.tick_params(labelsize=fontsize)
            # create legend?#http://stackoverflow.com/questions/9834452/how-do-i-make-a-single-legend-for-many-subplots-with-matplotlib
            ax.legend(['AAIMON (identity)','AAIMON_n (correction)','AASMON (identity + similarity)'],loc='upper right', fontsize=fontsize)
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

# code from http://stackoverflow.com/questions/13226038/calculating-angle-between-two-lines-in-python
def angle_between_slopes(pt1, pt2):
    x1, y1 = pt1
    x2, y2 = pt2
    inner_product = x1*x2 + y1*y2
    len1 = math.hypot(x1, y1)
    len2 = math.hypot(x2, y2)
    return math.acos(inner_product/(len1*len2))

def save_scatter_AAIMON_norm_and_AAIMON_slope_single_protein (fig_nr, fig, axarr, df_cr, x_data, y_data, y_data_n, AAIMON_slope, AAIMON_n_slope,
                                                                    TMD, zipout, row_nr, col_nr, fontsize, savefig, norm_scatter_path_prefix):
    # define data to plot
    datapointsize = 0.5
    # avoid the matplotlib nan/rgb bug by dropping the nan values, if present
    # https://github.com/matplotlib/matplotlib/issues/8658
    # ValueError: Invalid RGBA argument: 0.098000000000000004 (this is a perfectly good python rgb colour value!)
    df_scatter = df_cr.loc[:, ['obs_changes', '%s_AAIMON'%TMD, '%s_AAIMON_n'%TMD]]
    df_scatter.dropna(inplace=True)

    x_data_obs_changes = df_scatter['obs_changes']
    scatter_data_AAIMON = df_scatter['%s_AAIMON'%TMD]
    scatter_data_AAIMON_n = df_scatter['%s_AAIMON_n'%TMD]

    xlim_min = 0
    xlim_max = 60

    axarr[row_nr, col_nr].scatter(x_data_obs_changes, scatter_data_AAIMON, color="k", alpha=0.3, s=datapointsize)
    axarr[row_nr, col_nr].scatter(x_data_obs_changes, scatter_data_AAIMON_n, color=(0.843, 0.098, 0.1098), marker='^', alpha=0.3, s=datapointsize)
    axarr[row_nr, col_nr].plot(x_data, y_data, color="k", alpha=0.3)
    axarr[row_nr, col_nr].plot(x_data, y_data_n, color=(0.843, 0.098, 0.1098), alpha=0.3)
    # calculate angle between AAIMON_slope and AAIMON_n_slope
    angle = angle_between_slopes((x_data[0], y_data[0]), (x_data[1], y_data[1]))

    axarr[row_nr, col_nr].set_ylabel('%s AAIMON' % TMD, rotation='vertical', fontsize=fontsize)
    #axarr[row_nr, col_nr].set_xlabel('% identity')
    axarr[row_nr, col_nr].set_ylim(0.0, 3)
    axarr[row_nr, col_nr].annotate(s='AAIMON slope: {a:0.3f}, AAIMON_n slope: {b:0.3f}, angle {c:.2f}°'.format(a=AAIMON_slope, b=AAIMON_n_slope, c=angle),
                                   xy=(0.01, 1.01), xytext=None, xycoords='axes fraction', alpha=0.75, fontsize=fontsize)
    #axarr[row_nr, col_nr].set_xticks(range(xlim_min,xlim_max+1,10))

    if savefig:

        # apply the following formatting changes to all plots in the figure
        for ax in axarr.flat:
            # set x-axis min
            ax.set_xlim(xlim_min, xlim_max)
            # set x-axis ticks
            # use the slide selection to select every second item in the list as an xtick(axis label)
            ax.set_xticks(range(xlim_min,xlim_max+1,10))
            ax.set_xlabel('% observed aa substitutions', fontsize=fontsize)
            # change axis font size
            ax.tick_params(labelsize=fontsize)
            # create legend?#http://stackoverflow.com/questions/9834452/how-do-i-make-a-single-legend-for-many-subplots-with-matplotlib
            ax.legend(['AAIMON slope', 'AAIMON_n slope', 'AAIMON','AAIMON_n'],loc='upper right', fontsize=fontsize)
            # add background grid
            ax.grid(True, color='0.75', alpha=0.5)
        # automatically tighten the layout of plots in the figure
        fig.tight_layout()
        # save files
        fig.savefig(norm_scatter_path_prefix + '_%01d.png' % fig_nr,format='png', dpi=200)
        # close figure
        plt.close('all')
        # add to zipfile
        zipout.write(norm_scatter_path_prefix + '_%01d.png' % fig_nr,arcname=os.path.basename(norm_scatter_path_prefix) + '_%01d.png' % fig_nr)
        # delete temporory files
        os.remove(norm_scatter_path_prefix + '_%01d.png' % fig_nr)

    return angle



