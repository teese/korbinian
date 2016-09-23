import korbinian.mtutils as utils
import matplotlib.pyplot as plt
import numpy as np
import os
import zipfile

def save_hist_AAIMON_ratio_single_protein (fig_nr, fig, axarr, df_cr, set_, TMD, binlist, zipout, row_nr, col_nr, fontsize, savefig, AAIMON_hist_path_prefix):
    """

    Parameters
    ----------
    fig_nr
    fig
    axarr
    df_cr
    set_
    TMD
    binlist
    zipout
    row_nr
    col_nr
    fontsize
    savefig
    AAIMON_hist_path_prefix

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
    min_identity_of_TMD_initial_filter = set_['cr_min_identity_of_TMD_initial_filter']
    #RESET SO THAT IT TAKES IT FROM df_cr AGAIN
    #df_cr_filt_AAIMON = df_cr_filt.loc[df_cr['%s_perc_ident'%TMD] >= min_identity_of_TMD_initial_filter]

    # create numpy array of membranous over nonmembranous conservation ratios (identity)
    hist_data_I = np.array(df_cr['%s_AAIMON_ratio'%TMD].dropna())
    # use numpy to create a histogram
    freq_counts, bin_array = np.histogram(hist_data_I, bins=binlist)
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
    # create numpy array of membranous over nonmembranous conservation ratios (identity + similarity)
    hist_data_S = np.array(df_cr['%s_AASMON_ratio'%TMD].dropna())
    # use numpy to create a histogram
    freq_counts, bin_array = np.histogram(hist_data_S, bins=binlist)
    # create a line graph rather than a bar graph for the AASMON (ident + similarity)
    linecontainer_S = axarr[row_nr, col_nr].plot(centre_of_bar_in_x_axis, freq_counts,
                                                 color="#0101DF", alpha=0.5)
    # other colours that are compatible with colourblind readers: #8A084B Dark red, #B45F04 deep orange, reddish purple #4B088A
    # http://html-color-codes.info/
    # label the x-axis for each plot, based on the TMD
    axarr[row_nr, col_nr].set_xlabel('%s conservation ratio (membranous over nonmembranous)'%TMD,
                                     fontsize=fontsize)
    if savefig:
        # take x-axis min from settings
        xlim_min = set_["1p_smallest_bin"]
        # take x-axis max from settings
        xlim_max = set_["1p_largest_bin"]
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
            ax.legend(['AASMON (identity + similarity)', 'AAIMON (identity)'],
                      loc='upper right', fontsize=fontsize)
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
