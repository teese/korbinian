from scipy.stats import ttest_ind
import ast
import csv
import itertools
import korbinian
import korbinian.utils as utils
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import pandas as pd
import sys
import zipfile
import os
import pickle
import scipy
# import debugging tools
from korbinian.utils import pr, pc, pn, aaa

def save_figures_describing_proteins_in_list(pathdict, s, logging):
    logging.info("~~~~~~~~~~~~         starting run_save_figures_describing_proteins_in_list          ~~~~~~~~~~~~")
    plt.style.use('seaborn-whitegrid')
    plt.rcParams["legend.framealpha"] = 1
    plt.rcParams['errorbar.capsize'] = 3
    # set resolution for plots in png format
    dpi = 300
    #plt.rcParams.update(plt.rcParamsDefault)

    # define save settings
    save_png = s["save_png"]
    save_pdf = s["save_pdf"]

    # get a color list (HTML works best). Make it a long list, to accept list numbers from 1-1000
    color_list = utils.create_colour_lists()['HTML_list01'] * 1000
    # set default font size for plot
    fontsize = 8
    datapointsize = 8

    cdict = utils.create_colour_lists()
    #TUMblues = ['#0F3750', '#0076B8', '#9ECEEC']
    TUMblues = [cdict['TUM_colours']["TUM2"],cdict['TUM_colours']["TUMBlue"],cdict['TUM_colours']["TUM1"]]

    # letters for saving variations o0f a figure
    letters = list("abcdefghijk")

    try:
        base_filepath = pathdict["single_list_fig_path"]
    except:
        logging.info("Rare bug related to 'not subscriptable' has occurred. Logging full pathdict and pathdict type:")
        logging.info(pathdict)
        logging.info(type(pathdict))

    list_number = s["list_number"]

    # for xlim, use the min and max evolutionary distance settings for the full dataset
    # this is used for subsequent figures
    min_evol_distance = int((1 - s["max_ident"]) * 100)
    max_evol_distance = int((1 - s["min_ident"]) * 100)

    '''Prepare data for the following plots'''
    # load cr_summary file
    df_cr_summary = pd.read_csv(pathdict["list_cr_summary_csv"], sep=",", quoting=csv.QUOTE_NONNUMERIC, index_col=0, low_memory=False)
    if df_cr_summary.shape[0] < s["min_n_proteins_in_list"]:
        return "~~~~~~~~~~~~            run_save_figures skipped, only {} proteins in list           ~~~~~~~~~~~~".format(df_cr_summary.shape[0])

    sys.stdout.write('Preparing data for plotting'), sys.stdout.flush()
    # load summary file
    df_list = pd.read_csv(pathdict["list_csv"], sep=",", quoting=csv.QUOTE_NONNUMERIC, index_col=0, low_memory=False)

    # if 'uniprot_KW' in df_list.columns:
    #     if not "uniprot_KW_for_analysis" in df_list.columns:
    #         raise ValueError("Please run keyword analysis.")

    # merge cr_summary and summary file, if columns are equal in both files, suffix _dfc will be added in cr_summary column names for backwards compatibility
    df = pd.merge(df_cr_summary, df_list, left_index=True, right_index=True, suffixes=('_dfc', ''))

    # consider the dataset to be multipass if it contains at least 3 proteins with TM02
    if "TM02_AAIMON_slope" in df.columns and df["TM02_AAIMON_slope"].dropna().shape[0] >= 3:
        dataset_is_multipass = True
    else:
        dataset_is_multipass = False

    # create number of datapoint dependent alpha_dpd
    alpha_dpd = utils.calc_alpha_from_datapoints(df['AAIMON_mean_all_TM_res'])
    #sys.stdout.write('\nopacity of datapoints: {a:.2f}\n'.format(a=alpha_dpd))

    # filter to remove proteins that have less than ~5 homologues
    # this is only important for the beta-barrel dataset, which has a lot of these proteins!
    min_n_homol = s["min_homol"]
    n_prot_before_n_homol_cutoff = df.shape[0]
    df = df.loc[df['AAIMON_n_homol'] >= min_n_homol]
    n_prot_after_n_homol_cutoff = df.shape[0]
    n_removed = n_prot_before_n_homol_cutoff - n_prot_after_n_homol_cutoff
    # if any proteins have been removed, then print the exact number.
    if n_removed >= 1:
        sys.stdout.write("-- {}/{} -- proteins were removed, as they contained less than {} valid homologues. "
              "\nFinal number of proteins = {}\n".format(n_removed, n_prot_before_n_homol_cutoff, min_n_homol, n_prot_after_n_homol_cutoff))
        sys.stdout.flush()

    # open list_csv file
    #df_uniprot = pd.read_csv(pathdict["list_csv"], sep=",", quoting=csv.QUOTE_NONNUMERIC, index_col=0)

    prot_family_df_dict = {}
    # only split into families if there are uniprot annotations are are multipass
    if 'uniprot_KW' in df.columns and df.number_of_TMDs_excl_SP.mean() > 2:
        if df.number_of_TMDs_excl_SP.mean() > 2:
            # create a new column showing whether the protein is a GPCR
            if "GPCR" not in df.columns:
                # convert the keywords from a stringlist to a python list
                if isinstance(df['uniprot_KW'][0], str):
                    df['uniprot_KW'] = df['uniprot_KW'].apply(lambda x: ast.literal_eval(x))
                df['GPCR'] = df['uniprot_KW'].apply(lambda x: 'G-protein coupled receptor' in x)
                df['olfactory_receptor'] = df['prot_descr'].apply(korbinian.cons_ratio.keywords.KW_list_contains_any_desired_KW, args=(['Olfactory receptor'],))

            df_GPCR = df.loc[df['GPCR'] == True]
            # add the dataframe segments to a dictionary for easy access?
            prot_family_df_dict["df_GPCR"] = df_GPCR
            prot_family_df_dict["df_nonGPCR"] = df.loc[df['GPCR'] == False]
            prot_family_df_dict["df_olfactory_receptorGPCR"] = df_GPCR.loc[df_GPCR['olfactory_receptor'] == True]
            prot_family_df_dict["df_non_olfactory_receptorGPCR"] = df_GPCR.loc[df_GPCR['olfactory_receptor'] == False]
            list_prot_families = ["df_GPCR", "df_nonGPCR", "df_olfactory_receptorGPCR", "df_non_olfactory_receptorGPCR"]
        else:
            sys.stdout.write('Not a multipass dataset. GPCR figs will be skipped.\n')
    else:
        sys.stdout.write('No multipass proteins found with UniProt keywords. Cannot create figures 19-21 \n')
        list_prot_families = []

    # print mean aaimon slope for the full dataset
    mean_AAIMON_slope = df['AAIMON_slope_all_TM_res'].mean()
    logging.info("\n{p} mean AAIMON slope = {m:0.02f} * 10^-3".format(p="full dataset", m=mean_AAIMON_slope))

    # print mean values for protein families
    for prot_family in list_prot_families:
        mean_AAIMON_slope = prot_family_df_dict[prot_family]['AAIMON_slope_all_TM_res'].mean()
        logging.info("{p} mean AAIMON slope = {m:0.02f} * 10^-3".format(p=prot_family[3:], m=mean_AAIMON_slope))

    # # save dataframe
    # df.to_csv(pathdict["base_filename_summaries"] + '_df_figs.csv', sep=",", quoting=csv.QUOTE_NONNUMERIC)

    # create binlist
    linspace_binlist = np.linspace(s["mp_smallest_bin"],
                                   s["mp_largest_bin"],
                                   s["mp_number_of_bins"])

    # add 30 as the last bin, to make sure 100% of the data is added to the histogram, including major outliers
    binlist = np.append(linspace_binlist, s["mp_final_highest_bin"])

    # create list of colours to use in figures
    colour_lists = utils.create_colour_lists()
    tableau20 = colour_lists['tableau20']

    # DEPRECATED?
    # # create dataframe mean_AAIMON_each_TM
    # df_mean_AAIMON_each_TM = pd.DataFrame()
    # # add AAIMON each TMD to dataframe
    # for acc in df.index:
    #     for TMD in ast.literal_eval(df.loc[acc, 'list_of_TMDs']):
    #         df_mean_AAIMON_each_TM.loc[acc, '{a}_AAIMON_mean'.format(a=TMD)] = df.loc[acc, '{b}_AAIMON_mean'.format(b=TMD)]

    # count the maximum number of TMDs (excluding signal peptides) in the dataset
    max_num_TMDs = df.number_of_TMDs_excl_SP.max()

    # make list_of_TMDs a python list
    df['list_of_TMDs'] = df['list_of_TMDs'].apply(ast.literal_eval)
    df['list_of_TMDs_excl_SP'] = df['list_of_TMDs_excl_SP'].apply(ast.literal_eval)
    # logging saved data types
    sys.stdout.write('\nSaving figures as: ')
    if s['save_pdf']:
        sys.stdout.write(' .pdf ')
    if s['save_png']:
        sys.stdout.write(' .png ')
    sys.stdout.write('\n')

    # for Figs 97 and 98 set data to 'None' not to load data twice
    #data = False
    #binned_data = False

    if s['Fig01_Hist_AAIMON_and_AASMON']:
        Fig_Nr = 1
        title = 'Mean ratios'
        Fig_name = 'List{:02d}_Fig01_Hist_AAIMON_and_AASMON'.format(list_number)
        # create a new figure
        fig, ax = plt.subplots()
        # create numpy array of membranous over nonmembranous conservation ratios (identity)
        hist_data_AAIMON_mean = np.array(df['AAIMON_mean_all_TM_res'].dropna())
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
        hist_data_AAIMON_n_mean = np.array(df['AAIMON_n_mean_all_TM_res'].dropna())
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
        hist_data_AASMON_mean = np.array(df['AASMON_mean_all_TMDs_mean'].dropna())
        # use numpy to create a histogram
        freq_counts_S, bin_array_S = np.histogram(hist_data_AASMON_mean, bins=binlist)
        # barcontainer_S = axarr[row_nr,col_nr].bar(left=centre_of_bar_in_x_axis, height=freq_counts_S, align='center', width=col_width, color="#0101DF", edgecolor="#0101DF", alpha = 0.5)
        # create a line graph rather than a bar graph for the AASMON (ident + similarity)
        linecontainer_AASMON_mean = ax.plot(centre_of_bar_in_x_axis, freq_counts_S, color="#0101DF",
                                            alpha=0.5)
        # other colours that are compatible with colourblind readers: #8A084B Dark red, #B45F04 deep orange, reddish purple #4B088A
        # http://html-color-codes.info/
        # label the x-axis for each plot, based on the TMD
        ax.set_xlabel('TM/EM conservation ratio', fontsize=fontsize)
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

    if s['Fig02_Density_AAIMON_vs_evol_dist']:
        Fig_Nr = 2
        Fig_name = 'List{:02d}_Fig02_Density_AAIMON_vs_evol_dist'.format(list_number)

        # read data from disk
        in_zipfile = pathdict["save_df_characterising_each_homol_TMD"]
        if os.path.isfile(in_zipfile):
            with zipfile.ZipFile(in_zipfile, "r", zipfile.ZIP_DEFLATED) as openzip:
                data = pickle.load(openzip.open("data_characterising_each_homol_TMD.pickle", "r"))
                binned_data = pickle.load(openzip.open("binned_data_characterising_each_homol_TMD.pickle", "r"))
        else:
            raise FileNotFoundError("{} not found".format(in_zipfile))

        vmax = s['vmax']

        fig, ax = plt.subplots(figsize=(5, 5))

        x = data[:, 0]  # FASTA_gapped_identity
        y = data[:, 1]  # AAIMON for each TMD

        # histogram definition
        # data range
        xyrange = [[0, max_evol_distance], [0, 3]]
        # number of bins
        bins = [max_evol_distance*2, 120]
        # density threshold
        thresh = 1

        # histogram the data
        hh, locx, locy = scipy.histogram2d(x, y, range=xyrange, bins=bins)

        # fill the areas with low density by NaNs
        hh[hh < thresh] = np.nan

        # ax.scatter(x=x, y=y, color="#EE762C", alpha=0.2, s=0.008, marker='x', linewidths=0.003)

        im = ax.imshow(np.flipud(hh.T), cmap='Oranges', extent=np.array(xyrange).flatten(),
                       interpolation='none', origin='upper', aspect='auto', vmin=0, vmax=vmax)

        """Plot the binned data (line graph)
        The binned_data typically has the shape 50, 8
        binned_data[:, 0] = % substitutions bins, eg array([ 49.5,  48.5,  47.5,  46.5,  45.5]
        binned_data[:, 1] = AAIMON data, e.g. array([ 0.97796493,  0.98020355,  0.98166374,  0.98894771,  0.99462542,
        binned_data[:, 2-6] as above. 2 is the normalised data. I am not sure what 3-6 are.
        binned_data[:, 7] = homologue counts, e.g. array([ 214263.,  238795.,  215564.,  191406.,  154984.,

        """
        ax.plot(binned_data[:, 0], binned_data[:, 1], color='#0F3750', label='non-normalised')
        ax.plot(binned_data[:, 0], binned_data[:, 2], color='#9ECEEC', label='normalised')
        # ax.grid(False, which='both')
        ax.tick_params(axis='both', which='major', length=3, width=1, color='#CCCCCC')
        ax.set_xlabel('evolutionary distance (% substitutions)', fontsize=fontsize)
        #ax.set_ylabel(r'$\mathrm{\mathsf{\frac {TM}{EM}  conservation}}}}$', rotation=90, fontsize=fontsize)
        ax.set_ylabel('TM/EM  conservation ratio', rotation=90, fontsize=fontsize)

        # get colorbar from latest imshow element (color scale should be the same for all subplots)
        # fig.subplots_adjust(right=0.8)
        cbar_ax = fig.add_axes([0.12, 0.89, 0.78, 0.02])
        fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
        cbar_ax.xaxis.set_ticks_position('top')
        labels = cbar_ax.get_xmajorticklabels()
        labels[-1] = '>{}    '.format(vmax)
        cbar_ax.set_xticklabels(labels)
        cbar_ax.tick_params(pad=0, labelsize=fontsize)
        ax.tick_params(pad=1, labelsize=fontsize)
        ax.legend(frameon=True, loc='upper left', fontsize=fontsize)

        # set the xlim based on the chosen settings for evolutionary distance of those homologues
        ax.set_xlim(min_evol_distance, max_evol_distance)

        utils.save_figure(fig, Fig_name, base_filepath, save_png, save_pdf, dpi)

    if s['Fig03_Density_lipo_vs_TM_conservation']:
        Fig_Nr = 3

        # get the maximum mnumber of TMDs in the full dataset (e.g. 32)
        max_n_TMDs = int(df.number_of_TMDs_excl_SP.max())
        # create a large list of columns, e.g. ['TM01_AAIMON_slope',  'TM02_AAIMON_slope',  'TM03_AAIMON_slope', ...
        col_list_AAIMON_slope = ['TM{:02d}_AAIMON_slope'.format(TM_nr) for TM_nr in range(1, max_n_TMDs + 1)]
        col_list_lipo = ['TM{:02d}_lipo'.format(TM_nr) for TM_nr in range(1, max_n_TMDs + 1)]

        #plot the original dataset (excluding signal peptides)
        Fig03_Density_lipo_vs_TM_conservation(list_number, df, "", "", col_list_AAIMON_slope, col_list_lipo, max_evol_distance, base_filepath, save_png, save_pdf, dpi, fontsize)

        # create a separate graph of signal peptids if available
        if "SP01_start" in df.columns:
            df_contains_SP = df.loc[df.SP01_AAIMON_slope.notnull()]
            col_list_AAIMON_slope_SP = ["SP01_AAIMON_slope"]
            col_list_lipo_SP = ["SP01_lipo"]
            suffix_SP = "_signal_peptide"
            Fig03_Density_lipo_vs_TM_conservation(list_number, df_contains_SP, "", suffix_SP, col_list_AAIMON_slope_SP, col_list_lipo_SP, max_evol_distance, base_filepath, save_png, save_pdf, dpi, fontsize)

        for i, prot_family in enumerate(list_prot_families):
            # a, b, c, etc
            letter = letters[i]
            # _GPCR, _nonGPCR, etc
            suffix = "_{}".format(prot_family[3:])
            # get appropriate dataframe subset for analysis (above)
            df_Fig03 = prot_family_df_dict[prot_family]
            if not df_Fig03.empty:
                # plot
                Fig03_Density_lipo_vs_TM_conservation(list_number, df_Fig03, letter, suffix, col_list_AAIMON_slope, col_list_lipo, max_evol_distance, base_filepath, save_png, save_pdf, dpi, fontsize)

        # for human multipass, test GPCR TM01 and TM07 only
        if list_number in [2,5]:
            # redefine as only the first and 7th TM (first and lost GPCR TM)
            col_list_AAIMON_slope = ['TM01_AAIMON_slope', 'TM07_AAIMON_slope']
            col_list_lipo = ['TM01_lipo', 'TM07_lipo']
            for i, prot_family in enumerate(list_prot_families):
                letter = letters[i]
                suffix = "_{}".format(prot_family[3:]) + "_TM01_and_TM07_only"
                df_Fig03 = prot_family_df_dict[prot_family]
                if not df_Fig03.empty:
                    Fig03_Density_lipo_vs_TM_conservation(list_number, df_Fig03, letter, suffix, col_list_AAIMON_slope, col_list_lipo, max_evol_distance, base_filepath, save_png, save_pdf, dpi, fontsize)

    if s['Fig04_Boxplot_AAIMON_slope_each_TMD']:
        Fig_Nr = 4
        title = 'Boxplot of all TMDs'
        Fig_name = 'List{:02d}_Fig04_Boxplot_AAIMON_slope_each_TMD'.format(list_number)
        fig, ax = plt.subplots()
        # "#0489B1"
        alpha = 0.25
        col_width_value = 0.95
        ylabel = 'freq'
        xlabel = 'average conservation ratio (membranous over nonmembranous)'
        # legend =

        max_number_of_TMDs_excl_SP = int(df.number_of_TMDs_excl_SP.max())
        legend = []
        data_to_plot = []
        for i in range(1, max_number_of_TMDs_excl_SP + 1):
            TM = 'TM%02d' % i
            hist_data_AAIMON_each_TM = df['TM%02d_AAIMON_slope' % i].dropna()*1000
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

        ax.set_ylabel(r'm$_{\rm TM/EM} *10^{\rm -3}$', rotation='vertical', fontsize=fontsize)

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

    if s['Fig05_Boxplot_lipo_each_TMDs']:
        Fig_Nr = 5
        title = 'Boxplot lipophilicity of all TMDs'
        Fig_name = 'List{:02d}_Fig05_Boxplot_lipo_each_TMDs'.format(list_number)

        ### boxplot of all TMDs

        ### this section specifies the last bin to avoid bins containing only one TMD
        # join all numbers of TMDs together into a large list
        nested_list_all_TMDs = list(df['number_of_TMDs_excl_SP'])
        # convert list to pandas series
        all_TMDs_series = pd.Series(nested_list_all_TMDs)
        # obtain series of TMD_counts
        TMD_counts = all_TMDs_series.value_counts()
        # exclude TMD numbers with less than x applic0able proteins from boxplot max detection
        boxplot_cutoff_number_of_TMDs = 20
        #TMD_counts_major = TMD_counts[TMD_counts >= boxplot_cutoff_number_of_TMDs]
        max_num_TMDs = int(TMD_counts.index.max())

        if max_num_TMDs > boxplot_cutoff_number_of_TMDs:
            # title = str(keyword) + '_Boxplot'
            # Fig_name = str(str(Fig_Nr) + '._' + 'Keyword_' + title)
            fig, ax = plt.subplots()
            ax2 = plt.twinx()

            legend = []
            data_to_plot = []
            for i in range(1, max_num_TMDs + 1):
                TM = 'TM%02d' % i
                hist_data_AAIMON_each_TM = df['TM%02d_lipo' % i].dropna()
                if len(hist_data_AAIMON_each_TM) > 0:
                    data_to_plot.append(hist_data_AAIMON_each_TM)
                    legend.append(TM)

            # add values of every TMD number that is larger than the boxplot_cutoff_number_of_TMDs to final bin
            data_for_final_bin = []
            for i in range(max_num_TMDs + 1, df.number_of_TMDs_excl_SP.max().astype('int') + 1):
                # TM_final = 'TM%02d' % i
                hist_data_AAIMON_each_TM_final_bin = df['TM%02d_lipo' % i].dropna()
                # if len(hist_data_AAIMON_each_TM) > 0:
                data_for_final_bin.append(hist_data_AAIMON_each_TM_final_bin)
            final_bin = list(itertools.chain.from_iterable(data_for_final_bin))
            data_to_plot.append(final_bin)
            legend.append('>{}'.format(TM))

            n_elements_in_bin = []
            for element in data_to_plot:
                n_elements_in_bin.append(len(element))

            x = range(1, len(legend) + 1)
            ax2.plot(x, n_elements_in_bin, color='#0076B8', alpha=0.5)
            ax2.grid(b=False)
            ax2.set_ylabel('number of TMDs in bin', rotation='vertical', fontsize=fontsize)
            ax2.tick_params(labelsize=fontsize)

            meanpointprops = dict(marker='o', markerfacecolor='black', markersize=2, markeredgecolor='black')  # markeredgecolor='0.75',

            flierprops = dict(marker='o', markerfacecolor='green', markersize=12,
                              linestyle='none')
            # plot boxplot
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

            ax.set_ylabel('lipophilicity (Hessa scale)', rotation='vertical', fontsize=fontsize)
            # ax.set_ylim(-20, 30)

            ## Remove top axes and right axes ticks
            ax.get_xaxis().tick_bottom()
            ax.get_yaxis().tick_left()
            ## Custom x-axis labels
            ax.set_xticklabels(legend, rotation=45)
            ax.set_ylim(-0.5, 1)

            # add figure number to top left of subplot
            ax.annotate(s=str(Fig_Nr) + '.', xy=(0.04, 0.9), fontsize=fontsize, xytext=None,
                        xycoords='axes fraction', alpha=0.75)
            # add figure title to top left of subplot
            ax.annotate(s=title, xy=(0.1, 0.9), fontsize=fontsize, xytext=None, xycoords='axes fraction',
                        alpha=0.75)

            utils.save_figure(fig, Fig_name, base_filepath, save_png, save_pdf, dpi)


    if s['Fig06_Boxplot_AAIMON_slope_by_number_of_TMDs']:
        Fig_Nr = 6
        title = 'num_TMDs vs AAIMON'
        Fig_name = 'List{:02d}_Fig06_Boxplot_AAIMON_slope_by_number_of_TMDs'.format(list_number)
        fig, ax = plt.subplots()

        alpha = 0.25
        col_width_value = 0.95
        legend = []
        data_to_plot = []
        # iterate through df and get all AAIMONs with specified number of TMD
        for i in range(1, max_num_TMDs + 1):
            hist_data = []
            for acc in df.loc[df['list_of_TMDs'].notnull()].loc[df['list_of_TMDs'] != 'nan'].index:
                if df.loc[acc, 'number_of_TMDs'] == i:
                    hist_data.append(df.loc[acc, 'AAIMON_slope_all_TM_res']*1000)
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
        ax.set_ylabel(r'm$_{\rm TM/EM} *10^{\rm -3}$', fontsize=fontsize)
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


    if s['Fig07_Density_AAIMON_or_slope_vs_evol_distance']:
        Fig_Nr = 7
        title = 'compare AAIMON with AAIMON_slope'
        Fig_name = 'List{:02d}_Fig07_Density_AAIMON_or_slope_vs_evol_distance'.format(list_number)

        vmax = 3
        fig, (ax1, ax2) = plt.subplots(2, sharex=True, sharey=False, figsize=(5, 10))

        # histogram definition
        # data range
        xyrange = [[0, max_evol_distance], [0.2, 1.8]]
        # number of bins
        bins = [max_evol_distance*2, 120]
        # density threshold
        thresh = 1

        # plot AAIMON data
        x = df.obs_changes_mean
        y = df.AAIMON_mean_all_TM_res
        # histogram the data
        hh, locx, locy = scipy.histogram2d(x, y, range=xyrange, bins=bins)
        # fill the areas with low density by NaNs
        hh[hh < thresh] = np.nan
        im = ax1.imshow(np.flipud(hh.T), cmap='Oranges', extent=np.array(xyrange).flatten(),
                        interpolation='none', origin='upper', aspect='auto', vmin=0, vmax=vmax)

        # plot AAIMON_slope data with changed xyrange
        xyrange = [[0, max_evol_distance], [-20, 20]]
        # plot AAIMON_slope data
        x = df.obs_changes_mean
        y = df.AAIMON_slope_all_TMDs_mean * 1000
        # histogram the data
        hh, locx, locy = scipy.histogram2d(x, y, range=xyrange, bins=bins)
        # fill the areas with low density by NaNs
        hh[hh < thresh] = np.nan
        im = ax2.imshow(np.flipud(hh.T), cmap='Oranges', extent=np.array(xyrange).flatten(),
                        interpolation='none', origin='upper', aspect='auto', vmin=0, vmax=vmax)

        # define axis limits and parameters, set labels
        ax1.set_ylim(0.2, 1.8)
        ax2.set_ylim(-20, 20)
        ax1.tick_params(labelsize=fontsize, pad=2)
        ax2.tick_params(labelsize=fontsize, pad=2)
        ax1.set_ylabel('TM/EM conservation ratio', fontsize=fontsize)
        ax2.set_ylabel(r'm$_{\rm TM/EM} *10^{\rm -3}$', fontsize=fontsize)
        ax1.set(adjustable='box-forced')
        ax2.set(adjustable='box-forced')
        plt.xlabel('average evolutionary distance of homologues (% substitutions)', fontsize=fontsize)

        # add colorbar
        cbar_ax = fig.add_axes([0.12, 0.89, 0.78, 0.01])
        fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
        cbar_ax.xaxis.set_ticks_position('top')
        labels = cbar_ax.get_xmajorticklabels()
        labels[-1] = '>{}'.format(vmax)
        cbar_ax.set_xticklabels(labels)
        cbar_ax.tick_params(pad=0, labelsize=fontsize)
        # remove white space between subplots
        fig.subplots_adjust(wspace=0.2, hspace=0.075)

        utils.save_figure(fig, Fig_name, base_filepath, save_png, save_pdf, dpi)

    if s['Fig08_Hist_AAIMON_slope_TM01_vs_lastTM']:
        ###### for backwards compatibility ##### can be removed if all data is re-processed after march 5 2017
        if not 'AAIMON_slope_last_TMD' in df.columns:
            sys.stdout.write('AAIMON_slope_last_TMD not in dataframe -> older version of data, re-run "gather_AAIMON_ratios"; adding data for figure')
            for n, acc in enumerate(df.index):
                if n % 200 == 0:
                    sys.stdout.write('. '), sys.stdout.flush()
                last_TMD = df.loc[acc, 'last_TMD']
                df.loc[acc, 'AAIMON_slope_last_TMD'] = df.loc[acc, '%s_AAIMON_slope' % last_TMD]

        Fig_Nr = 8
        title = 'AAIMON_slope TM01 vs lastTM'
        Fig_name = 'List{:02d}_Fig08_Hist_AAIMON_slope_TM01_vs_lastTM'.format(list_number)
        binlist = np.linspace(-40, 40, 61)
        linewidth = 1
        fig, ax = plt.subplots()

        ###   TM01 AAIMON_slope   ###
        # create numpy array of membranous over nonmembranous conservation ratios (identity)
        hist_data = (df['TM01_AAIMON_slope'] * 1000).dropna()
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
        linecontainer_AAIMON_mean = ax.plot(centre_of_bar_in_x_axis, freq_counts_normalised, color='k',
                                            alpha=0.9, linewidth=linewidth)

        ###   last TMD AAIMON_slope   ###
        # create numpy array of membranous over nonmembranous conservation ratios (identity)
        hist_data = (df['AAIMON_slope_last_TMD'] * 1000).dropna()
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
        linecontainer_AAIMON_mean = ax.plot(centre_of_bar_in_x_axis, freq_counts_normalised, ':', color='k',
                                            alpha=0.9,
                                            linewidth=linewidth)

        ax.set_xlabel(r'm$_{\rm TM/EM} *10^{\rm -3}$', fontsize=fontsize)
        # move the x-axis label closer to the x-axis
        ax.xaxis.set_label_coords(0.5, -0.085)
        # x and y axes min and max
        xlim_min = -30
        xlim_max = 30
        ax.set_xlim(xlim_min, xlim_max)
        ax.set_ylabel('freq', rotation='vertical', fontsize=fontsize)
        # change axis font size
        ax.tick_params(labelsize=fontsize)
        # create legend
        ax.xaxis.set_label_coords(0.5, -0.07)
        # ax.yaxis.set_label_coords(-0.005, 0.5)

        # add legend
        ax.legend(['TM01', 'last TM'], fontsize=fontsize, frameon=True)

        # add annotations
        ax.annotate(s="TM less conserved", xy=(0, -0.09), fontsize=fontsize, xytext=None, xycoords='axes fraction')
        ax.annotate(s="TM more conserved", xy=(1.0, -0.09), fontsize=fontsize, xytext=None, horizontalalignment='right', xycoords='axes fraction')

        utils.save_figure(fig, Fig_name, base_filepath, save_png, save_pdf, dpi)

    if s['Fig09_Hist_lipo_TM01_vs_lastTM']:
        if dataset_is_multipass:
            Fig_Nr = 9
            title = 'Lipo TM01 vs lastTM'
            Fig_name = 'List{:02d}_Fig09_Hist_lipo_TM01_vs_lastTM'.format(list_number)
            min_ = -0.5
            max_ = 0.8
            binlist = np.linspace(min_, max_, 21)
            fig, ax = plt.subplots()
            # offset = len(protein_lists) - 1

            # create numpy array of membranous over nonmembranous conservation ratios (identity)
            hist_data = np.array(df['lipo_mean_excl_TM01'].dropna())
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
            linecontainer_AAIMON_mean = ax.plot(centre_of_bar_in_x_axis, freq_counts_normalised, color=TUMblues[0],
                                                alpha=1,
                                                linewidth=1)

            # create numpy array of membranous over nonmembranous conservation ratios (identity)
            hist_data = np.array(df['TM01_lipo'].dropna())
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
            linecontainer_AAIMON_mean = ax.plot(centre_of_bar_in_x_axis, freq_counts_normalised, color=TUMblues[1],
                                                alpha=1, linewidth=1)

            # # create numpy array of membranous over nonmembranous conservation ratios (identity)
            # hist_data = np.array(df['lipo_last_TMD'].dropna())
            # # use numpy to create a histogram
            # freq_counts, bin_array = np.histogram(hist_data, bins=binlist)
            # freq_counts_normalised = freq_counts / freq_counts.max()
            # # assuming all of the bins are exactly the same size, make the width of the column equal to XX% (e.g. 95%) of each bin
            # col_width = float('%0.3f' % (0.95 * (bin_array[1] - bin_array[0])))
            # # when align='center', the central point of the bar in the x-axis is simply the middle of the bins ((bin_0-bin_1)/2, etc)
            # centre_of_bar_in_x_axis = (bin_array[:-2] + bin_array[1:-1]) / 2
            # # add the final bin, which is physically located just after the last regular bin but represents all higher values
            # bar_width = centre_of_bar_in_x_axis[3] - centre_of_bar_in_x_axis[2]
            # centre_of_bar_in_x_axis = np.append(centre_of_bar_in_x_axis, centre_of_bar_in_x_axis[-1] + bar_width)
            # linecontainer_AAIMON_mean = ax.plot(centre_of_bar_in_x_axis, freq_counts_normalised, color=TUMblues[2],
            #                                     alpha=1,
            #                                     linewidth=1)

            ###############################################################
            #                                                             #
            #                       set up plot style                     #
            #                                                             #
            ###############################################################

            ax.set_xlabel('lipophilicity (Hessa scale)', fontsize=fontsize)
            # move the x-axis label closer to the x-axis
            ax.xaxis.set_label_coords(0.5, -0.085)
            # x axes min and max
            xlim_min = min_
            xlim_max = max_
            ax.set_xlim(xlim_min, xlim_max)
            ax.set_ylabel('freq', rotation='vertical', fontsize=fontsize)
            # change axis font size
            ax.tick_params(labelsize=fontsize)

            # add annotations
            ax.annotate(s="more lipophilic", xy=(0, -0.1), fontsize=fontsize, xytext=None, xycoords='axes fraction')
            ax.annotate(s="less lipophilic", xy=(1.0, -0.1), fontsize=fontsize, xytext=None, horizontalalignment='right', xycoords='axes fraction')

            ax.legend(['mean excl. TM01', 'TM01'], fontsize=fontsize, frameon=True)

            utils.save_figure(fig, Fig_name, base_filepath, save_png, save_pdf, dpi)


    if s['Fig10_Boxplot_AAIMON_by_number_of_TMDs'] and s["max_TMDs"] >= 2:
        Fig_Nr = 10
        title = 'num_TMDs vs AAIMON'
        Fig_name = 'List{:02d}_Fig10_Boxplot_AAIMON_by_number_of_TMDs'.format(list_number)
        fig, ax = plt.subplots()

        # data that is binned
        column_for_bins = 'number_of_TMDs'
        # data that is plotted in bin
        column_for_data = 'AAIMON_mean_all_TM_res'
        hist_data = []
        legend = []
        TMD_number = list(range(1, 16, 1))
        for element in TMD_number:
            select = df[column_for_bins] == element
            data = df.loc[select, column_for_data].values
            hist_data.append(data)
            legend.append(element)
        select = df[column_for_bins] > TMD_number[-1]
        data = df.loc[select, column_for_data].values
        hist_data.append(data)
        legend.append('>15')

        fig, ax = plt.subplots()

        meanpointprops = dict(marker='o', markerfacecolor='black', markersize=3, markeredgecolor='black')

        flierprops = dict(marker='o', markerfacecolor='black', markersize=12, linestyle='none')
        boxplotcontainer = ax.boxplot(hist_data, sym='+', whis=1.5, showmeans=True,
                                      meanprops=meanpointprops)

        list_n_datapoints = [len(x) for x in hist_data]
        x_for_list_n_datapoints = list(range(1, len(list_n_datapoints) + 1, 1))
        ax2 = ax.twinx()
        line_graph_container = ax2.plot(x_for_list_n_datapoints, list_n_datapoints, color="#53A7D5", alpha=0.8)

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

        ax.set_xlabel('number of TMDs in protein', fontsize=fontsize)
        ax.tick_params(labelsize=fontsize)
        ax2.tick_params(labelsize=fontsize)
        ax.set_ylim(ymin=0, ymax=2)
        ax2.set_ylim(ymin=0, ymax=800)
        ax.set_xlim(xmin=0, xmax=17)
        # move the x-axis label closer to the x-axis
        # ax.xaxis.set_label_coords(0.45, -0.085)
        ax.set_ylabel('Average TM/EM conservation ratio for all TMDs', fontsize=fontsize)
        ax2.set_ylabel('Number of proteins in bin', fontsize=fontsize)
        ## Remove top axes and right axes ticks
        ax.get_xaxis().tick_bottom()
        ax.get_yaxis().tick_left()
        ## Custom x-axis labels
        ax.set_xticklabels(legend)

        utils.save_figure(fig, Fig_name, base_filepath, save_png, save_pdf, dpi)

    if s['Fig11_Boxplot_AAIMON_slope_by_seqlen']:
        Fig_Nr = 11
        title = 'seqlen vs AAIMON'
        Fig_name = 'List{:02d}_Fig11_Boxplot_AAIMON_by_seqlen'.format(list_number)

        fig, ax = plt.subplots()

        # data that is binned
        column_for_bins = 'seqlen'
        # data that is plotted in bin
        column_for_data = 'AAIMON_slope_all_TM_res'
        # specify variable for binning function
        x = df[column_for_bins]
        # specify number of bins
        nbin = 10
        npt = len(x)
        # get bin-borders
        borders = np.interp(np.linspace(0, npt, nbin + 1), np.arange(npt), np.sort(x)).astype(int).tolist()
        # extend the bin borders to catch up every value
        borders[0] = 0
        borders[-1] = borders[-1] + 1
        # initialise lists for legend and data in bin
        legend = []
        hist_data = []
        # generate data in bin via selecting rows from pandas dataframe, create legend
        for n in range(1, len(borders), 1):
            legend.append('-'.join([str(borders[n - 1]), str(borders[n])]))
            select = (df[column_for_bins] > borders[n - 1]) & (df[column_for_bins] <= borders[n])
            data = df.loc[select, column_for_data].values
            hist_data.append(data*1000)

        meanpointprops = dict(marker='o', markerfacecolor='black', markersize=3)  # markeredgecolor='0.75',
        flierprops = dict(marker='o', markerfacecolor='green', markersize=12, linestyle='none')
        boxplotcontainer = ax.boxplot(hist_data, sym='+', whis=1.5, showmeans=True,
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

        ax.set_xlabel('Length of protein in bins', fontsize=fontsize)
        # move the x-axis label closer to the x-axis
        # ax.xaxis.set_label_coords(0.45, -0.085)
        ax.set_ylabel(r'm$_{\rm TM/EM} *10^{\rm -3}$', fontsize=fontsize)
        ## Remove top axes and right axes ticks
        ax.get_xaxis().tick_bottom()
        ax.get_yaxis().tick_left()
        ## Custom x-axis labels
        ax.set_xticklabels(legend, rotation=25)
        fig.tight_layout()
        utils.save_figure(fig, Fig_name, base_filepath, save_png, save_pdf, dpi)

    if s['Fig12_Boxplot_AAIMON_slope_by_nonTMD_len']:
        Fig_Nr = 12
        title = 'seqlen vs AAIMON'
        Fig_name = 'List{:02d}_Fig12_Fig12_Boxplot_AAIMON_slope_by_nonTMD_len'.format(list_number)

        fig, ax = plt.subplots()

        # data that is binned
        column_for_bins = 'nonTMD_SW_align_len_mean'
        # data that is plotted in bin
        column_for_data = 'AAIMON_slope_all_TM_res'
        # specify variable for binning function
        x = df[column_for_bins]
        # specify number of bins
        nbin = 10
        npt = len(x)
        # get bin-borders
        borders = np.interp(np.linspace(0, npt, nbin + 1), np.arange(npt), np.sort(x)).astype(int).tolist()
        # extend the bin borders to catch up every value
        borders[0] = 0
        borders[-1] = borders[-1] + 1
        # initialise lists for legend and data in bin
        legend = []
        hist_data = []
        # generate data in bin via selecting rows from pandas dataframe, create legend
        for n in range(1, len(borders), 1):
            legend.append('-'.join([str(borders[n - 1]), str(borders[n])]))
            select = (df[column_for_bins] > borders[n - 1]) & (df[column_for_bins] <= borders[n])
            data = df.loc[select, column_for_data].values
            hist_data.append(data*1000)

        meanpointprops = dict(marker='o', markerfacecolor='black', markersize=3)  # markeredgecolor='0.75',

        flierprops = dict(marker='o', markerfacecolor='green', markersize=12, linestyle='none')
        boxplotcontainer = ax.boxplot(hist_data, sym='+', whis=1.5, showmeans=True,
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

        ax.set_xlabel('Average length of nonTMD region in homologues', fontsize=fontsize)
        #ax.set_ylim(ymin=0, ymax=2)
        # move the x-axis label closer to the x-axis
        # ax.xaxis.set_label_coords(0.45, -0.085)
        ax.set_ylabel(r'm$_{\rm TM/EM} *10^{\rm -3}$', fontsize=fontsize)
        ## Remove top axes and right axes ticks
        ax.get_xaxis().tick_bottom()
        ax.get_yaxis().tick_left()
        ## Custom x-axis labels
        ax.set_xticklabels(legend, rotation=25)
        fig.tight_layout()
        utils.save_figure(fig, Fig_name, base_filepath, save_png, save_pdf, dpi)

    if s['Fig13_Boxplot_AAIMON_slope_by_num_of_simap_hits']:
        Fig_Nr = 13
        title = 'number SIMAP hits'
        Fig_name = 'List{:02d}_Fig13_Boxplot_AAIMON_slope_by_num_of_simap_hits'.format(list_number)

        fig, ax = plt.subplots()

        # data that is binned
        column_for_bins = 'AAIMON_n_homol'
        # data that is plotted in bin
        column_for_data = 'AAIMON_slope_all_TM_res'
        # specify variable for binning function
        x = df[column_for_bins]
        # specify number of bins
        nbin = 10
        npt = len(x)
        # get bin-borders
        borders = np.interp(np.linspace(0, npt, nbin + 1), np.arange(npt), np.sort(x)).astype(int).tolist()
        # extend the bin borders to catch up every value
        borders[0] = 0
        borders[-1] = borders[-1] + 1
        # initialise lists for legend and data in bin
        legend = []
        hist_data = []
        # generate data in bin via selecting rows from pandas dataframe, create legend
        for n in range(1, len(borders), 1):
            legend.append('-'.join([str(borders[n - 1]), str(borders[n])]))
            select = (df[column_for_bins] > borders[n - 1]) & (df[column_for_bins] <= borders[n])
            data = df.loc[select, column_for_data].values
            hist_data.append(data*1000)

        meanpointprops = dict(marker='o', markerfacecolor='black', markersize=3)  # markeredgecolor='0.75',

        flierprops = dict(marker='o', markerfacecolor='green', markersize=12, linestyle='none')
        boxplotcontainer = ax.boxplot(hist_data, sym='+', whis=1.5, showmeans=True,
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

        ax.set_xlabel('total number of homologues', fontsize=fontsize)
        #ax.set_ylim(ymin=0, ymax=2)
        # move the x-axis label closer to the x-axis
        # ax.xaxis.set_label_coords(0.45, -0.085)
        ax.set_ylabel(r'm$_{\rm TM/EM} *10^{\rm -3}$', fontsize=fontsize)
        ## Remove top axes and right axes ticks
        ax.get_xaxis().tick_bottom()
        ax.get_yaxis().tick_left()
        ## Custom x-axis labels
        ax.set_xticklabels(legend, rotation=25)
        fig.tight_layout()
        utils.save_figure(fig, Fig_name, base_filepath, save_png, save_pdf, dpi)

    if 'uniprot_KW' in df.columns:

        if s['Fig14_Hist_AAIMON_GPCRs_vs_nonGPCRs']:
            if True in df.GPCR.tolist():
                Fig_Nr = 14
                title = 'only GPCR in uniprot KW, NORM'
                Fig_name = 'List{:02d}_Fig14_Hist_AAIMON_GPCRs_vs_nonGPCRs'.format(list_number)
                fig, ax = plt.subplots()
                # create numpy array of membranous over nonmembranous conservation ratios (identity)
                hist_data_AAIMON_mean = np.array(df_GPCR['AAIMON_mean_all_TM_res'].dropna())
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
                df_nonGPCR = df.loc[df['GPCR'] == False]

                # create numpy array of membranous over nonmembranous conservation ratios (identity)
                hist_data_AAIMON_mean = np.array(df_nonGPCR['AAIMON_mean_all_TM_res'].dropna())
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

        if s['Fig15_Boxplot_AAIMON_by_number_of_TMDs_GPCRs_only']:
            if True in df.GPCR.tolist():
                Fig_Nr = 15
                title = 'Only GPCRs, boxplot for each TMD'
                Fig_name = 'List{:02d}_Fig15_Boxplot_AAIMON_by_number_of_TMDs_GPCRs_only'.format(list_number)
                fig, ax = plt.subplots()
                # "#0489B1"
                alpha = 0.25
                col_width_value = 0.95
                ylabel = 'freq'
                xlabel = 'average conservation ratio (membranous over nonmembranous)'
                # legend =

                number_of_TMDs_excl_SP = df_GPCR.number_of_TMDs_excl_SP.max()
                legend = []
                data_to_plot = []
                for i in range(1, max_num_TMDs + 1):
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
                ax.set_ylabel('TM/EM conservation ratio', fontsize=fontsize)
                # add figure number to top left of subplot
                ax.annotate(s=str(Fig_Nr) + '.', xy=(0.04, 0.9), fontsize=fontsize, xytext=None,
                            xycoords='axes fraction', alpha=0.75)
                # add figure title to top left of subplot
                ax.annotate(s=title, xy=(0.1, 0.9), fontsize=fontsize, xytext=None, xycoords='axes fraction',
                            alpha=0.75)

                utils.save_figure(fig, Fig_name, base_filepath, save_png, save_pdf, dpi)
            else:
                sys.stdout.write('Dataset does not contain GPCRs; cannot create figure 21 \n')

    if s['Fig16_Scatterplot_AAIMON_n_vs_slope']:
        Fig_Nr = "16a"
        title = 'AAIMON vs. AAIMON_slope'
        Fig_name = 'List{:02d}_Fig16a_Scatterplot_AAIMON_vs_AAIMON_slope'.format(list_number)
        fig, ax = plt.subplots()
        colour = color_list[2]

        x = df['AAIMON_mean_all_TM_res']
        y = df['AAIMON_slope_all_TMDs_mean']*1000

        if len(x) > 5:
            # calculate linear regression for fitted line
            linear_regression = np.polyfit(x, y, 1)
            fit_fn = np.poly1d(linear_regression)
            x_first_last_dp = [x.min(), x.max()]
            y_fitted = fit_fn(x_first_last_dp)
            ax.plot(x_first_last_dp, y_fitted, "--", alpha=0.75, color=colour, label="fitted")
            ax.annotate(s='y = {a:.5f}x + {b:.5f}'.format(a=linear_regression[0], b=linear_regression[1]), xy=(0.85, 0.95),
                        fontsize=fontsize-2, xytext=None, xycoords='axes fraction', alpha=0.75)
        else:
            logging.info("The dataset has less than 5 proteins. Lines of best fit will not be calculated.")

        ax.scatter(x, y, color = colour, alpha=alpha_dpd, s=datapointsize, label="data")
        ax.set_ylabel(r'm$_{\rm TM/EM} *10^{\rm -3}$', rotation='vertical', fontsize=fontsize)
        ax.set_xlabel('TM/EM conservation ratio', fontsize=fontsize)
        ax.annotate(s=str(Fig_Nr) + '.', xy=(0.04, 0.9), fontsize=fontsize, xytext=None,
                    xycoords='axes fraction', alpha=0.75)
        # add figure title to top left of subplot
        ax.annotate(s=title, xy=(0.1, 0.9), fontsize=fontsize, xytext=None, xycoords='axes fraction',
                    alpha=0.75)

        # change axis font size
        ax.tick_params(labelsize=fontsize)
        ax.legend(loc="lower right")
        fig.tight_layout()
        utils.save_figure(fig, Fig_name, base_filepath, save_png, save_pdf, dpi)

        ##################################################################
        Fig_Nr = "16b"
        title = 'AAIMON_slope_all_TM_res vs AAIMON_slope_all_TMDs_mean'
        Fig_name = 'List{:02d}_Fig16b_Scatterplot_AAIMON_slope_all_TM_res_vs_AAIMON_slope_all_TMDs_mean'.format(list_number)
        fig, ax = plt.subplots()

        x = df['AAIMON_slope_all_TM_res']*1000
        y = df['AAIMON_slope_all_TMDs_mean']*1000

        if len(x) > 5:
            # calculate linear regression for fitted line
            linear_regression = np.polyfit(x, y, 1)
            fit_fn = np.poly1d(linear_regression)
            x_first_last_dp = [x.min(), x.max()]
            y_fitted = fit_fn(x_first_last_dp)
            ax.plot(x_first_last_dp, y_fitted, "--", alpha=0.75, color=colour, label="fitted")
            ax.plot(x_first_last_dp, x_first_last_dp, "--", alpha=0.75, color="k", label="symmetrical")
            ax.annotate(s='y = {a:.5f}x + {b:.5f}'.format(a=linear_regression[0], b=linear_regression[1]), xy=(0.85, 0.95),
                        fontsize=fontsize-2, xytext=None, xycoords='axes fraction', alpha=0.75)
        else:
            logging.info("The dataset has less than 5 proteins. Lines of best fit will not be calculated.")

        ax.scatter(x, y, color = colour, alpha=alpha_dpd, s=datapointsize, label="data")
        ax.set_xlabel(r'm$_{\rm TM/EM} *10^{\rm -3}$ (all TM residues in protein)', fontsize=fontsize)
        ax.set_ylabel(r'm$_{\rm TM/EM} *10^{\rm -3}$ (mean of all TMDs calculated separately)', rotation='vertical', fontsize=fontsize)

        ax.annotate(s=str(Fig_Nr) + '.', xy=(0.04, 0.9), fontsize=fontsize, xytext=None,
                    xycoords='axes fraction', alpha=0.75)
        # add figure title to top left of subplot
        ax.annotate(s=title, xy=(0.1, 0.9), fontsize=fontsize, xytext=None, xycoords='axes fraction',
                    alpha=0.75)

        # change axis font size
        ax.tick_params(labelsize=fontsize)
        ax.legend(loc = "lower right")
        fig.tight_layout()
        utils.save_figure(fig, Fig_name, base_filepath, save_png, save_pdf, dpi)


    if s['Fig17_Scatterplot_perc_identity_nonTMD_vs_TMD']:
        Fig_Nr = 17
        title = 'perc identity all TMDs  vs. perc identity nonTMD'
        Fig_name = 'List{:02d}_Fig17_Scatterplot_perc_identity_nonTMD_vs_TMD'.format(list_number)
        fig, ax = plt.subplots()

        x = df['TMD_perc_ident_mean'] * 100
        y = df['nonTMD_perc_ident_mean'] * 100

        if len(x) > 5:
            linear_regression = np.polyfit(x, y, 1)
            fit_fn = np.poly1d(linear_regression)
            x_first_last_dp = [x.min(), x.max()]
            y_fitted = fit_fn(x_first_last_dp)
            ax.plot(x_first_last_dp, y_fitted, ":", color="k", alpha=0.75, label="fitted")
            ax.annotate(s='y = {a:.5f}x + {b:.5f}'.format(a=linear_regression[0], b=linear_regression[1]),
                        xy=(0.65, 0.95), fontsize=fontsize - 2, xytext=None, xycoords='axes fraction',alpha=0.75)

        ax.scatter(x, y, s=datapointsize, alpha=alpha_dpd, color=TUMblues, label="data")
        #symmetrical = [s["min_ident"]*100, s["max_ident"]*100]
        symmetrical = x_first_last_dp
        ax.plot(symmetrical, symmetrical, "--", color=TUMblues[0], alpha=0.5, label="symmetrical")
        ax.set_xlabel('Average % identity in TM region, all homologues', fontsize=fontsize)
        ax.set_ylabel('Average % identity in EM region, all homologues', rotation='vertical', fontsize=fontsize)
        ax.tick_params(labelsize=fontsize)
        #ax.set_xlim(100-max_evol_distance, 100)
        #ax.set_ylim(s["min_ident"]*100, 100)

        ax.annotate(s=str(Fig_Nr) + '.', xy=(0.04, 0.9), fontsize=fontsize, xytext=None,
                    xycoords='axes fraction', alpha=0.75)
        # add figure title to top left of subplot
        ax.annotate(s=title, xy=(0.1, 0.9), fontsize=fontsize, xytext=None, xycoords='axes fraction',
                    alpha=0.75)

        # change axis font size
        ax.tick_params(labelsize=fontsize)
        ax.legend(loc="lower right")
        fig.tight_layout()
        utils.save_figure(fig, Fig_name, base_filepath, save_png, save_pdf, dpi)

    if 'uniprot_KW' in df.columns and "uniprot_KW_for_analysis" in df.columns:

        # convert stringlists to python lists
        #if isinstance(df.uniprot_KW_for_analysis.dropna().iloc[0], str):
        df["uniprot_KW_for_analysis"] = df["uniprot_KW_for_analysis"].fillna("[]")
        df["uniprot_KW_for_analysis"] = df.uniprot_KW_for_analysis.dropna().apply(lambda x: ast.literal_eval(x))


        if s['Fig18_KW_assoc_with_large_number_of_homol']:
            Fig_Nr = 18
            title = 'keywords associated with many valid homologues'
            Fig_name = 'List{:02d}_Fig18_KW_assoc_with_large_number_of_homol'.format(list_number)

            # create cutoff, mean + 1 std
            cutoff = df.AAIMON_n_homol.mean() + df.AAIMON_n_homol.std()
            cutoff_int = int(np.round(cutoff))

            # get index of proteins with many or few homologues
            many_homol_index = df.AAIMON_n_homol.loc[df.AAIMON_n_homol > cutoff].index
            few_homol_index = df.AAIMON_n_homol.loc[df.AAIMON_n_homol <= cutoff].index
            # select subset of dataframe with many or few homologues
            df_few = df.loc[few_homol_index, :]
            df_many = df.loc[many_homol_index, :]
            # get a large list of keywords
            many_KW_list = df_many.uniprot_KW_for_analysis.dropna().tolist()
            few_KW_list = df_few.uniprot_KW_for_analysis.dropna().tolist()

            # convert list of keywords into pandas series, and use value_counts to count
            many_ser = pd.Series(utils.flatten(many_KW_list))
            few_ser = pd.Series(utils.flatten(few_KW_list))
            df_KW = pd.DataFrame()
            df_KW["many"] = many_ser.dropna().value_counts()
            df_KW["few"] = few_ser.dropna().value_counts()
            # total number of proteins with keyword (used as a cutoff)
            df_KW["total"] = df_KW["many"] + df_KW["few"]
            # fraction of proteins containing the keyword, in the many or few dataset
            df_KW["many_frac_containing_KW"] = df_KW.many / len(many_ser)
            df_KW["few_frac_containing_KW"] = df_KW.few / len(few_ser)
            # relative abundance of the keyword in the fraction with many homologues
            df_KW["relative_abundance_KW"] = df_KW["many_frac_containing_KW"] / df_KW["few_frac_containing_KW"]
            df_KW.sort_values("relative_abundance_KW", inplace=True, ascending=False)

            # only examine keywords with a significant number of proteins
            cutoff_min_num_prot_with_KW = 100
            df_KW_signif = df_KW.loc[df_KW.total > cutoff_min_num_prot_with_KW]
            if not df_KW_signif.empty:
                fig, ax = plt.subplots()
                x = df_KW_signif.index
                y = df_KW_signif.many_frac_containing_KW
                y2 = df_KW_signif.few_frac_containing_KW
                width = 0.4
                x_ind = np.array(range(len(x))) + width
                ax.bar(left=x_ind, height=y, width=width, color="#0489B1", label=">{} homologues".format(cutoff_int))
                ax.bar(left=x_ind + width, height=y2, width=width, color="0.9", label="<{} homologues".format(cutoff_int))
                ax.set_xlim(0, x_ind[-1] + 1 + width / 2)
                ax.set_xticks(x_ind + width)
                ax.set_xticklabels(x, rotation=90)
                ax.set_ylabel("Fraction of proteins containing keyword")
                ax.legend(frameon=True)
                fig.tight_layout()
                utils.save_figure(fig, Fig_name, base_filepath, save_png, save_pdf, dpi)
    else:
        logging.info("keyword figures were not created. 'uniprot_KW_for_analysis' was not in columns")

    if s['Fig19_barchart_keywords_signif_TM_cons']:
        Fig_Nr = 19
        title = 'keywords with significant difference in AAIMON slope'
        Fig_name = 'List{:02d}_Fig19_barchart_keywords_signif_TM_cons'.format(list_number)
        fig, ax = plt.subplots()

        list_enzyme_KW, list_ignored_KW, PFAM_dict = korbinian.cons_ratio.keywords.get_list_enzyme_KW_and_list_ignored_KW()


        # open up the keywords csv
        kw_csv = os.path.join(pathdict["keywords"], 'List%02d_keywords.csv' % s["list_number"])
        if os.path.isfile(kw_csv):
            dfKW = pd.read_csv(kw_csv, index_col=0)
            # drop any annotations where the difference was not significant from the rest of the dataset
            dfKW = dfKW.loc[dfKW["p-value_AAIMON_slope"] < 0.05]
            # drop any ignored keywords
            index_excl_ignored_KW = set(dfKW.index) - set(list_ignored_KW)
            # replace the specific betabarrel PFAM ID

            dropped_KW = set(dfKW.index).intersection(set(list_ignored_KW))
            if len(dropped_KW) != 0:
                logging.info("dropped KW from list_ignored_KW : {}".format(dropped_KW))
            dfKW = dfKW.reindex(index_excl_ignored_KW)

            # DEPRECATED
            #dfKW.index = pd.Series(dfKW.index).replace(PFAM_dict)

            # re-sort by p-value
            dfKW.sort_values("p-value_AAIMON_slope", inplace=True, ascending=True)
            # ONLY CONTINUE IF THERE ARE ACTUALLY SOME KW WITH DATAPOINTS
            if not dfKW.AAIMON_slope_whole_dataset_mean.empty:
                # get the mean value for the whole dataset
                AAIMON_slope_whole_dataset_mean = dfKW.AAIMON_slope_whole_dataset_mean.iloc[0]

                # if it is a betabarrel dataset, replace the PFAM ids with the description
                if True in df.betabarrel.tolist():
                    dfKW["orig_index"] = dfKW.index
                    dfKW.index = dfKW.index.astype(str) + "\n(" + dfKW["PFAM_desc"] + ")"

                # create a new dataframe with the mean values, including the "All {}".format(s["list_description"]) of the full dataset
                orig_index = dfKW.index.tolist()

                new_index = ["All {}".format(s["list_description"])] + orig_index
                df_barchart = pd.DataFrame(index=new_index)
                # add the means
                df_barchart["AAIMON_slope_mean"] = dfKW.AAIMON_slope_keyword_mean
                df_barchart.loc["All {}".format(s["list_description"]), "AAIMON_slope_mean"] = AAIMON_slope_whole_dataset_mean
                # calculate the SEM from the std
                df_barchart["AAIMON_slope_keyword_SEM"] = dfKW.AAIMON_slope_keyword_std / np.sqrt(dfKW.number_of_proteins_keyword)
                # get the SEM using the original dataframe derived from df_cr
                sqrt_n = np.sqrt(df.AAIMON_slope_all_TMDs_mean.dropna().shape[0])
                df_barchart.loc["All {}".format(s["list_description"]), "AAIMON_slope_keyword_SEM"] = scipy.stats.sem(df.AAIMON_slope_all_TMDs_mean)

                color_first_bar_all_proteins = "0.35"

                ind = np.arange(df_barchart.shape[0]) + 0.5
                width = 0.4
                # add first color to be different ["first"] + ["chosen"]*10
                colours = [color_first_bar_all_proteins] + [color_list[list_number - 1]] * df_barchart.shape[0]
                ax.bar(ind, df_barchart["AAIMON_slope_mean"] * 1000, color=colours)
                ax.errorbar(ind, df_barchart["AAIMON_slope_mean"] * 1000, yerr=df_barchart["AAIMON_slope_keyword_SEM"] * 1000, fmt="none",  ecolor="k", ls="none", capthick=1, elinewidth=1, capsize=4)
                #ax.set_xticks(ind + width)
                ax.set_xticks(ind)
                # take only first 20 characters for the x-axis label
                ax.set_xticklabels(pd.Series(df_barchart.index).str[0:20], rotation=90)
                ax.set_ylabel(r'm$_{\rm TM/EM} *10^{\rm -3}$', rotation='vertical', fontsize=fontsize + 3)
                if df_barchart.shape[0] < 10:
                    ax.set_xlim(0, 10)
                # add annotations
                ax.annotate(s="TM less\nconserved", xy=(-0.08, 0.1), fontsize=fontsize, xytext=None, xycoords='axes fraction', rotation=90)
                ax.annotate(s="TM more\nconserved", xy=(-0.08, 0.9), fontsize=fontsize, xytext=None, xycoords='axes fraction', rotation=90)
                fig.tight_layout()
                utils.save_figure(fig, Fig_name, base_filepath, save_png, save_pdf, dpi)
            else:
                sys.stdout.write("Fig19_barchart_keywords_signif_TM_cons skipped, dfKW.AAIMON_slope_whole_dataset_mean is empty, no KW were significant?")

    if s['Fig20_BB_linechart_lipo_patterns_orig_dataset']:

        if True in df.betabarrel.tolist():

            Fig_Nr = "20a"
            title = 'keywords with significant difference in AAIMON slope'
            Fig_name = 'List{:02d}_Fig20a_BB_linechart_lipo_last_TM'.format(list_number)
            fig, ax = plt.subplots()

            final_TM_aligned_fasta = os.path.join(base_filepath, "BB_terminal_TM_aligned.fas")
            final_TM_aligned_fasta_finalF = os.path.join(base_filepath, "BB_terminal_TM_aligned_subset_final_AA_is_F.fas")

            def get_pos_dep_lipo(seq_ser):
                """Get the position-dependent lipophilicity

                Parameters
                ----------
                seq_ser : pd.Series
                    Series of amino acid sequences

                Returns
                -------
                n_residue_list : list
                    number of residues at that position
                mean_lipo_list : list
                    mean lipophilicity at that position
                """
                # drop empty
                seq_ser = seq_ser.dropna()
                # get number of proteins
                n_prot = seq_ser.shape[0]
                # convert the series of sequences to a 2-D numpy array
                arr = np.array(seq_ser.apply(lambda x: list(x)).tolist())
                n_residue_list = []
                mean_lipo_list = []
                # iterate through the rows of the array, corresponding to positions
                for i in range(arr.shape[1]):
                    # get amino acids as a list
                    rowlist = arr[:, i]
                    # join to make a string
                    joined = "".join(rowlist)
                    # count gaps, and calculate number of residues without gaps
                    n_gaps = joined.count("-")
                    n_residues = n_prot - n_gaps
                    n_residue_list.append(n_residues)
                    # get lipophilicity
                    mean_lipo = korbinian.utils.calc_lipophilicity(joined)
                    mean_lipo_list.append(mean_lipo)
                return n_residue_list, mean_lipo_list

            def get_altern_pattern_lipo_and_save_fasta(dfp, TM, base_filepath):
                TM_seqs = dfp["{}_seq".format(TM)].dropna()
                len_longest_TM_strand = dfp["{}_seq".format(TM)].str.len().max()
                for acc in dfp["{}_seq".format(TM)].dropna().index:
                    seq = dfp.loc[acc, "{}_seq".format(TM)]
                    f0_seq = seq[::2]
                    f1_seq = seq[1::2]
                    if korbinian.utils.calc_lipophilicity(f0_seq) < korbinian.utils.calc_lipophilicity(f1_seq):
                        frame = 0
                        dfp.loc[acc, "{}_lipo_frame".format(TM)] = 0
                    else:
                        dfp.loc[acc, "{}_lipo_frame".format(TM)] = 1
                        frame = 1
                    endpadding = len_longest_TM_strand - frame + 1
                    frontgap = "-" * frame
                    dfp.loc[acc, "{}_padded".format(TM)] = ("{}{:-<%d}" % endpadding).format(frontgap, seq)

                n_residue_list, mean_lipo_list = get_pos_dep_lipo(dfp["{}_padded".format(TM)])

                aligned_fasta = os.path.join(base_filepath, "BB_{}_aligned.fas".format(TM))

                with open(aligned_fasta, "w") as f:
                    for acc in dfp["{}_padded".format(TM)].dropna().index:
                        f.write(">{}\n{}\n".format(acc,  dfp["{}_padded".format(TM)][acc]))

                return n_residue_list, mean_lipo_list

            # get length of longest last TMD
            longest_TM_strand = df["last_TMD_seq"].str.len().max()
            # pad with - to create alignment
            df["last_TMD_seq_padded"] = df.last_TMD_seq.dropna().apply(lambda x: ("{:->%d}" % longest_TM_strand).format(x))
            df["last_TMD_seq_padded"].head()

            last_TMD_seq_padded_ser = df["last_TMD_seq_padded"].dropna()
            with open(final_TM_aligned_fasta, "w") as f:
                for acc in last_TMD_seq_padded_ser.index:
                    f.write(">{}\n{}\n".format(acc, last_TMD_seq_padded_ser[acc]))

            final_F_ser = last_TMD_seq_padded_ser.loc[last_TMD_seq_padded_ser.str[-1] == "F"]

            # final_TM_aligned_fasta_finalF
            with open(final_TM_aligned_fasta_finalF, "w") as f:
                for acc in final_F_ser.index:
                    f.write(">{}\n{}\n".format(acc, final_F_ser[acc]))

            n_residue_list, mean_lipo_list = get_pos_dep_lipo(last_TMD_seq_padded_ser)

            color_list = utils.create_colour_lists()['HTML_list01']
            cdict = utils.create_colour_lists()
            c1 = cdict["TUM_colours"]["TUM1"]
            c2 = cdict["TUM_colours"]["TUM2"]

            fig, ax = plt.subplots()
            ax2 = ax.twinx()

            range_0_13 = np.arange(1, len(n_residue_list) + 1)
            x = range_0_13 - len(n_residue_list)
            ax.plot(x, mean_lipo_list, color=c2)
            ax.grid(b=False)
            ax.set_ylabel("mean hydrophobicity, Hessa scale", color=c2)

            ax2.plot(x, n_residue_list, "--", color=c1)
            ax2.grid(b=False)
            ax2.set_ylabel("number of residues at position", color=c1)

            ax2.spines["left"].set_color(c2)
            ax2.spines["right"].set_color(c1)
            ax.yaxis.label.set_color(c2)
            ax2.yaxis.label.set_color(c1)
            ax.tick_params(axis="y", colors=c2)
            ax2.tick_params(axis="y", colors=c1)

            ax.set_xlabel("residue position (aligned to terminal position, 0)")
            ax.set_title("C-terminal strand hydrophobicity pattern")

            fig.tight_layout()
            utils.save_figure(fig, Fig_name, base_filepath, save_png, save_pdf, dpi)



            ###############################################################
            #                                                             #
            #            linechart hydrophobicity central TM              #
            #                                                             #
            ###############################################################
            Fig_Nr = "20b"
            title = 'linechart hydrophobicity central TM'

            for TM in ["TM03", "TM04"]:
                n_residue_list, mean_lipo_list = get_altern_pattern_lipo_and_save_fasta(df, TM, base_filepath)

                Fig_name = 'List{:02d}_Fig20b_BB_linechart_lipo_{}'.format(list_number, TM)
                fig, ax = plt.subplots()
                ax2 = ax.twinx()

                range_0_13 = np.arange(1, len(n_residue_list) + 1)
                x = range_0_13
                ax.plot(x, mean_lipo_list, color=c2)
                ax.grid(b=False)
                ax.set_ylabel("mean hydrophobicity, Hessa scale", color=c2)

                ax2.plot(x, n_residue_list, color=c1, linestyle="--")

                ax2.grid(b=False)
                ax2.set_ylabel("number of residues at position", color=c1)

                ax2.spines["left"].set_color(c2)
                ax2.spines["right"].set_color(c1)
                ax.yaxis.label.set_color(c2)
                ax2.yaxis.label.set_color(c1)
                ax.tick_params(axis="y", colors=c2)
                ax2.tick_params(axis="y", colors=c1)

                ax.set_xlabel("residue position (aligned to initial position, 1)")
                ax.set_title("hydrophobicity pattern")
                #fig_res_numbers = r"D:\Schweris\Projects\TMD_Conservation\20170612 beta signal analysis\weblogo\residue_numbers_and_hydrophobicity_TM03.png"
                fig.tight_layout()
                utils.save_figure(fig, Fig_name, base_filepath, save_png, save_pdf, dpi)

            ###############################################################
            #                                                             #
            #            AA propensity at terminal position               #
            #                                                             #
            ###############################################################
            Fig_Nr = "20c"
            title = 'AA propensity at terminal position'
            Fig_name = 'List{:02d}_Fig20c_BB_barchart_AA_prop_terminal_pos'.format(list_number)
            # number of proteins with a last TMD seq
            n_prot = df["last_TMD_seq"].dropna().shape[0]
            # value counts of each amino acid as a percentage
            vc_last = df["last_TMD_seq"].str[-1].value_counts() / n_prot * 100
            vc_minus_2 = df["last_TMD_seq"].str[-3].value_counts() / n_prot * 100
            vc_minus_4 = df["last_TMD_seq"].str[-5].value_counts() / n_prot * 100
            vc_minus_6 = df["last_TMD_seq"].str[-7].value_counts() / n_prot * 100
            # combine into a single dataframe
            dfa = pd.DataFrame()
            dfa["terminal"] = vc_last
            dfa["pos -2"] = vc_minus_2
            dfa["pos -4"] = vc_minus_4
            dfa["pos -6"] = vc_minus_6
            # plot directly from dataframe
            dfa.plot(kind="bar", color=color_list[1:])
            plt.title("amino acid propensity at terminal position")
            plt.ylabel("frequency (%)")
            plt.tight_layout()
            if save_png:
                plt.savefig(os.path.join(base_filepath, Fig_name + ".png"),dpi=dpi)
            if save_pdf:
                plt.savefig(os.path.join(base_filepath, "pdf", Fig_name + ".pdf"))


    # dictionary showing where the bitopic data is stored (e.g. list02 (human polytopic) matches list01 (human bitopic)
    bitopic_list_dict = {"02":"01", "08":"07", "12":"01", "31":"30", "48":"47"}


    if s["Fig21_linechart_lipo_f_c_l_vs_number_of_TMDs"]:
        Fig_Nr = 21
        title = 'linechart_lipo_f_c_l_vs_number_of_TMDs'
        Fig_name = 'List{:02d}_Fig21a_linechart_lipo_f_c_l_vs_number_of_TMDs'.format(list_number)

        # minimum number of proteins in bin based on number of TMDs
        min_n_prot = 8

        if max_num_TMDs >= 2:

            fig, ax = plt.subplots()

            if True in df.betabarrel.tolist():
                min_n_TMDs = 8
                max_num_TMDs_fig21 = 31
            else:
                min_n_TMDs = 1
                max_num_TMDs_fig21 = int(max_num_TMDs)

            """DEPRECATED: ORIGINAL LISTS NOW GENERALLY EXCLUDE GPCRs.
            # exclude GPCRs from multipass datasets
            # NOTE: source of data is the original list csv,
            # including proteins with insufficient homologues for conservation analyses
            if "GPCR" in df.columns:
                df_filt = df_list.loc[df_list.GPCR == False]
            else:
                df_filt = df
            """

            list_number_str = "{:02d}".format(list_number)
            if list_number_str in bitopic_list_dict:
                bitopic_list_number_str = bitopic_list_dict[list_number_str]
                df_filt = add_bitopic_proteins_to_df_from_another_list(df, list_number_str, bitopic_list_number_str, min_n_homol, pathdict, logging)
            else:
                df_filt = df

            # DEPRECATED, use add_bitopic_proteins_to_df_from_another_list function instead
            # # for list 2, add the singlepass data to the analysis for this plot
            # if list_number in [2, 12]:
            #     list1_csv = pathdict["list_csv"].replace("{:02d}".format(list_number), "01")
            #     list1_cr_csv = pathdict["list_cr_summary_csv"].replace("02", "01")
            #     df_list1 = pd.read_csv(list1_csv, index_col=0)
            #     df_cr1 = pd.read_csv(list1_cr_csv, index_col=0)
            #     df_list1_merged = pd.merge(df_list1, df_cr1, left_index=True, right_index=True, suffixes=('_dfc', ''))
            #     df_list1_merged = df_list1_merged.loc[df_list1_merged['AAIMON_n_homol'] >= min_n_homol]
            #     df_filt = pd.merge(df, df_list1_merged, how="outer")
            #     logging.info("{} proteins added from List01 for Figure 21, linechart_lipo_f_c_l_vs_number_of_TMDs".format(df_list1_merged.shape[0]))
            #
            # # for list 31, add the singlepass data to the analysis for this plot
            # if list_number == 31:
            #     list30_csv = pathdict["list_csv"].replace("31", "30")
            #     list30_cr_csv = pathdict["list_cr_summary_csv"].replace("31", "30")
            #     df_list30 = pd.read_csv(list30_csv, index_col=0)
            #     df_cr30 = pd.read_csv(list30_cr_csv, index_col=0)
            #     df_list30_merged = pd.merge(df_list30, df_cr30, left_index=True, right_index=True, suffixes=('_dfc', ''))
            #     df_list30_merged = df_list30_merged.loc[df_list30_merged['AAIMON_n_homol'] >= min_n_homol]
            #     #df_filt = pd.concat([df, df_list30_merged])
            #     df_filt = pd.merge(df, df_list30_merged, how="outer")
            #     logging.info("{} proteins added from List01 for Figure 21, linechart_lipo_f_c_l_vs_number_of_TMDs".format(df_list30_merged.shape[0]))
            #
            # # for list 31, add the singlepass data to the analysis for this plot
            # if list_number == 31:
            #     list30_csv = pathdict["list_csv"].replace("31", "30")
            #     list30_cr_csv = pathdict["list_cr_summary_csv"].replace("31", "30")
            #     df_list30 = pd.read_csv(list30_csv, index_col=0)
            #     df_cr30 = pd.read_csv(list30_cr_csv, index_col=0)
            #     df_list30_merged = pd.merge(df_list30, df_cr30, left_index=True, right_index=True, suffixes=('_dfc', ''))
            #     df_list30_merged = df_list30_merged.loc[df_list30_merged['AAIMON_n_homol'] >= min_n_homol]
            #     #df_filt = pd.concat([df, df_list30_merged])
            #     df_filt = pd.merge(df, df_list30_merged, how="outer")
            #     logging.info("{} proteins added from List01 for Figure 21, linechart_lipo_f_c_l_vs_number_of_TMDs".format(df_list30_merged.shape[0]))


            sys.stdout.write("\ncode gives RuntimeWarning\n")
            # NOTE: this code gives a RuntimeWarning: Degrees of freedom <= 0 for slice
            # warnings.warn("Degrees of freedom <= 0 for slice", RuntimeWarning)
            # It's likely to be related to SEM of an empty series or nan values,
            # however adding a check for dff["lipo_mean_central_TMDs"].dropna().empty did not get rid of the warning
            dfn = pd.DataFrame()
            for i in range(min_n_TMDs, max_num_TMDs_fig21):
                # create a filtered selection with just that number of TMDs
                dff = df_filt.loc[df_filt.number_of_TMDs == i]
                dfn.loc[i, "TM01_lipo_mean"] = dff["TM01_lipo"].mean()
                dfn.loc[i, "TM01_lipo_sem"] = scipy.stats.sem(dff["TM01_lipo"])
                # lipo_mean_central_TMDs
                dfn.loc[i, "central_lipo_mean"] = dff["lipo_mean_central_TMDs"].mean()
                dfn.loc[i, "central_lipo_sem"] = scipy.stats.sem(dff["lipo_mean_central_TMDs"])
                # last TM
                dfn.loc[i, "last_lipo_mean"] = dff["lipo_last_TMD"].mean()
                dfn.loc[i, "last_lipo_sem"] = scipy.stats.sem(dff["lipo_last_TMD"])
                # number of proteins
                dfn.loc[i, "n_prot"] = dff["lipo_last_TMD"].dropna().shape[0]

            fig, ax = plt.subplots()
            ax2 = ax.twinx()
            fontsize = 8
            c0 = "0.5"
            c1 = cdict["TUM_colours"]["TUM1"]
            c2 = cdict["TUM_colours"]["TUM2"]

            dfn = dfn.loc[dfn.n_prot >= min_n_prot]

            if not dfn.empty:

                # plot central
                ax.plot(dfn.index, dfn.central_lipo_mean, color=c0, label="central TMs")
                ax.errorbar(x=dfn.index, y=dfn.central_lipo_mean, yerr=np.array(dfn.central_lipo_sem),
                            color=c0, capthick=1, elinewidth=1.5, capsize=3, label=None)
                # plot first
                ax.plot(dfn.index, dfn.TM01_lipo_mean, color=c1, label="TM01", linestyle="--")
                ax.errorbar(x=dfn.index, y=dfn.TM01_lipo_mean, yerr=np.array(dfn.TM01_lipo_sem),
                            color=c1, capthick=1, elinewidth=1.5, capsize=3, linestyle="--", label=None)

                # plot last
                ax.plot(dfn.index, dfn.last_lipo_mean, color=c2, label="last TM", linestyle=":", )
                ax.errorbar(x=dfn.index, y=dfn.last_lipo_mean, yerr=np.array(dfn.last_lipo_sem),
                            color=c2, capthick=1, elinewidth=1.5, capsize=3, linestyle=":", label=None)

                # plot number of proteins
                ax2.plot(dfn.index, dfn.n_prot, color=cdict['TUM_oranges']["TUM0"], label="number of proteins")
                ax2.grid(b=False)
                ax2.set_ylabel("number of proteins", color=cdict['TUM_oranges']["TUM0"])
                #ax2.legend(loc="upper right")
                ylim_max = dfn.n_prot.max()
                if not np.isnan(ylim_max):
                    ax2.set_ylim(0, ylim_max * 5)
                ax2.spines["right"].set_color(cdict['TUM_oranges']["TUM0"])
                ax2.yaxis.label.set_color(cdict['TUM_oranges']["TUM0"])
                ax2.tick_params(axis="y", colors=cdict['TUM_oranges']["TUM0"])

                ax.set_xlabel("number of TM regions")
                ax.set_ylabel("mean lipophilicity\n(Hessa scale)")
                #ax.set_ylim(0.45, 0.75)
                if not np.isnan(ylim_max):
                    ax.set_xlim(dfn.index.min() - 1, dfn.index.max() + 1 + 1)

                # add annotations
                ax.annotate(s="TM more\nlipophilic", xy=(-0.15, 0.1), fontsize=fontsize, xytext=None, xycoords='axes fraction', rotation=90)
                ax.annotate(s="TM less\nlipophilic", xy=(-0.15, 0.9), fontsize=fontsize, xytext=None, xycoords='axes fraction', rotation=90)

                ax.legend(frameon=True, loc="upper left")
                #fig.tight_layout()
                utils.save_figure(fig, Fig_name, base_filepath, save_png, save_pdf, dpi)

                ##################################################################
                #                                                                #
                #        linechart_AAIMON_slope_f_c_l_vs_number_of_TMDs          #
                #                                                                #
                ##################################################################
                Fig_name = 'List{:02d}_Fig21b_linechart_AAIMON_slope_f_c_l_vs_number_of_TMDs'.format(list_number)
                fig, ax = plt.subplots()

                list_number_str = "{:02d}".format(list_number)
                if list_number_str in bitopic_list_dict:
                    bitopic_list_number_str = bitopic_list_dict[list_number_str]
                    df_filt = add_bitopic_proteins_to_df_from_another_list(df, list_number_str, bitopic_list_number_str, min_n_homol, pathdict, logging)
                else:
                    df_filt = df
                # 1
                # # exclude GPCRs from multipass datasets
                # # NOTE: source of data is the merged cr_summary and list csv, filtered by sufficient homologues
                # if "GPCR" in df.columns:
                #     df_filt = df.loc[df.GPCR == False]
                # else:
                #     df_filt = df
                #
                # # HUMAN MULTIPASS: for list 2 or 12 , add the singlepass data to the analysis for this plot
                # if list_number in [2, 12]:
                #     # limit columns to avoid errors from python lists during merging
                #     df_list1_merged = df_list1_merged.loc[:, ["number_of_TMDs", "TM01_AAIMON_slope", "AAIMON_slope_central_TMDs", "AAIMON_slope_last_TMD"]]
                #     df_list1_merged_shape = df_list1_merged.shape
                #     df_filt = pd.merge(df_filt, df_list1_merged, how="outer")
                #     #logging.info("{} proteins added from List01 for Figure 21, linechart_lipo_f_c_l_vs_number_of_TMDs".format(df_list1_merged.shape[0]))
                #
                # # CRYSTAL STRUCTURES: for list 31, add singlepass data from list 30
                # if list_number == 31:
                #     # limit columns to avoid errors from python lists during merging
                #     df_list30_merged = df_list30_merged.loc[:, ["number_of_TMDs", "TM01_AAIMON_slope", "AAIMON_slope_central_TMDs", "AAIMON_slope_last_TMD"]]
                #     df_list30_merged_shape = df_list30_merged.shape
                #     df_filt = pd.merge(df_filt, df_list30_merged, how="outer")
                #
                # # YEAST: for list 48, add singlepass data from list 47
                # if list_number == 48:
                #     # limit columns to avoid errors from python lists during merging
                #     df_list30_merged = df_list47_merged.loc[:, ["number_of_TMDs", "TM01_AAIMON_slope", "AAIMON_slope_central_TMDs", "AAIMON_slope_last_TMD"]]
                #     df_list30_merged_shape = df_list30_merged.shape
                #     df_filt = pd.merge(df_filt, df_list30_merged, how="outer")


                dfn = pd.DataFrame()
                for i in range(min_n_TMDs, max_num_TMDs_fig21):
                    # create a filtered selection with just that number of TMDs
                    dff = df_filt.loc[df_filt.number_of_TMDs == i]
                    dfn.loc[i, "TM01_AAIMON_slope_mean"] = dff["TM01_AAIMON_slope"].mean()*1000
                    dfn.loc[i, "TM01_AAIMON_slope_sem"] = scipy.stats.sem(dff["TM01_AAIMON_slope"])*1000
                    # lipo_mean_central_TMDs
                    dfn.loc[i, "central_AAIMON_slope_mean"] = dff["AAIMON_slope_central_TMDs"].mean()*1000
                    dfn.loc[i, "central_AAIMON_slope_sem"] = scipy.stats.sem(dff["AAIMON_slope_central_TMDs"])*1000
                    # last TM
                    dfn.loc[i, "last_AAIMON_slope_mean"] = dff["AAIMON_slope_last_TMD"].mean()*1000
                    dfn.loc[i, "last_AAIMON_slope_sem"] = scipy.stats.sem(dff["AAIMON_slope_last_TMD"])*1000
                    # number of proteins
                    dfn.loc[i, "n_prot"] = dff["AAIMON_slope_last_TMD"].dropna().shape[0]

                dfn = dfn.loc[dfn.n_prot >= min_n_prot]

                if not dfn.empty:

                    fig, ax = plt.subplots()
                    ax2 = ax.twinx()
                    fontsize = 8
                    c0 = "0.5"
                    c1 = cdict["TUM_colours"]["TUM1"]
                    c2 = cdict["TUM_colours"]["TUM2"]
                    # plot central
                    ax.plot(dfn.index, dfn.central_AAIMON_slope_mean, color=c0, label="central TMs")
                    ax.errorbar(x=dfn.index, y=dfn.central_AAIMON_slope_mean, yerr=np.array(dfn.central_AAIMON_slope_sem),
                                color=c0, capthick=1, elinewidth=1.5, capsize=3, label=None)
                    # plot first
                    ax.plot(dfn.index, dfn.TM01_AAIMON_slope_mean, color=c1, label="TM01", linestyle="--")
                    ax.errorbar(x=dfn.index, y=dfn.TM01_AAIMON_slope_mean, yerr=np.array(dfn.TM01_AAIMON_slope_sem),
                                color=c1, capthick=1, elinewidth=1.5, capsize=3, linestyle="--", label=None)

                    # plot last
                    ax.plot(dfn.index, dfn.last_AAIMON_slope_mean, color=c2, label="last TM", linestyle=":", )
                    ax.errorbar(x=dfn.index, y=dfn.last_AAIMON_slope_mean, yerr=np.array(dfn.last_AAIMON_slope_sem),
                                color=c2, capthick=1, elinewidth=1.5, capsize=3, linestyle=":", label=None)

                    # plot number of proteins
                    ax2.plot(dfn.index, dfn.n_prot, color=cdict['TUM_oranges']["TUM0"], label="number of proteins")
                    ax2.grid(b=False)
                    ax2.set_ylabel("number of proteins", color=cdict['TUM_oranges']["TUM0"])
                    #ax2.legend(loc="upper right")
                    ax2.set_ylim(0, dfn.n_prot.max()*5)
                    ax2.spines["right"].set_color(cdict['TUM_oranges']["TUM0"])
                    ax2.yaxis.label.set_color(cdict['TUM_oranges']["TUM0"])
                    ax2.tick_params(axis="y", colors=cdict['TUM_oranges']["TUM0"])

                    ax.set_xlabel("number of TM regions")
                    #ax.set_ylabel("m$_{\rm TM/EM} *10^{\rm -3}$")
                    ax.set_ylabel(r'm$_{\rm TM/EM} *10^{\rm -3}$', rotation='vertical', fontsize=fontsize + 3)
                    #ax.set_ylim(0.45, 0.75)
                    ax.set_xlim(dfn.index.min() - 1, dfn.index.max() + 1 + 1)

                    # add annotations
                    ax.annotate(s="TM less\nconserved", xy=(-0.1, 0.1), fontsize=fontsize, xytext=None, xycoords='axes fraction', rotation=90)
                    ax.annotate(s="TM more\nconserved", xy=(-0.1, 0.9), fontsize=fontsize, xytext=None, xycoords='axes fraction', rotation=90)

                    ax.legend(frameon=True, loc="upper left")

                    #fig.tight_layout()
                    utils.save_figure(fig, Fig_name, base_filepath, save_png, save_pdf, dpi)


    # plot only for human multipass?
    if s["Fig22_multipass_linechart_f_c_l_cons_lipo_protein_subgroups"] and list_number in [2, 12]:
        if 'uniprot_KW' in df.columns and "uniprot_KW_for_analysis" in df.columns:

            # select only those proteins with at least 3 TMDs, to give first, central and last TM conservation
            df_min_3 = df.loc[df.number_of_TMDs_excl_SP >= 3].copy()

            # select the solute carriers
            df_solute = df_min_3.loc[df_min_3.prot_descr.str.contains("Solute carrier")].copy()
            # select only those with 11, 12 and 13 TM helices, which show a highly conserved TM01
            n_TMDs_to_include = [11, 12, 13]
            df_solute["selected_n_TMDs"] = df_solute.number_of_TMDs.apply(lambda x: x in n_TMDs_to_include)
            df_solute_sel = df_solute.loc[df_solute["selected_n_TMDs"]]
            # create df with subset of ion channels
            df_min_3["ion channel"] = df_min_3.uniprot_KW_for_analysis.apply(lambda x: "Ion channel" in x)
            df_ion_channel = df_min_3.loc[df_min_3["ion channel"]]
            # create df with subset of enzymes
            df_min_3["enzyme"] = df_min_3.uniprot_KW_for_analysis.apply(lambda x: "Enzyme" in x)
            df_enzyme = df_min_3.loc[df_min_3["enzyme"]]
            min_ = df_enzyme.number_of_TMDs.min()
            # create df with three transporters, who show a conserved TM01
            sel_transporter_kw = ["Amino-acid transport", "Neurotransmitter transport", "Lipid transport"]
            def iter_in_list(input_list):
                item_found = False
                for item in input_list:
                    if item in sel_transporter_kw:
                        item_found = True
                        break
                return item_found
            df_min_3["sel_transporter"] = df_min_3.uniprot_KW_for_analysis.apply(iter_in_list)
            df_sel_transporter = df_min_3.loc[df_min_3["sel_transporter"]]
            # create a selection that excludes all above
            others_index = set(df_min_3.index) - set(df_GPCR.index) - set(df_solute_sel.index) - set(df_ion_channel.index) - set(df_enzyme.index) - set(df_sel_transporter.index)
            df_others = df_min_3.loc[others_index, :]

            min_ = df_others.number_of_TMDs.min()
            df_others = df_others.loc[df_others.number_of_TMDs >= 3]

            # define colours
            c1 = cdict["TUM_oranges"]["TUM3"]
            c2 = cdict["TUM_oranges"]["TUM2"]
            c3 = cdict["TUM_colours"]["TUM2"]
            c4 = cdict["TUM_colours"]["TUM1"]
            c5 = "k"
            c6 = "0.5"

            ##################################################################
            #                                                                #
            #                 linechart for CONSERVATION                     #
            #             (first, central, last for protein subgroups)       #
            #                                                                #
            ##################################################################
            Fig_Nr = 22
            title = 'linechart_cons_f_c_l_multipass_prot_subgroups'
            Fig_name = 'List{:02d}_Fig22a_multipass_linechart_f_c_l_cons_protein_subgroups'.format(list_number)

            fig = plt.figure()
            ax = fig.add_axes([0.1,0.1,0.6,0.6])
            x = np.arange(3)

            y_solute_sel = [df_solute_sel.TM01_AAIMON_slope.mean() * 1000, df_solute_sel.AAIMON_slope_central_TMDs.mean() * 1000, df_solute_sel.AAIMON_slope_last_TMD.mean() * 1000]
            ax.plot(x, y_solute_sel, color=c1, label="solute carriers with 11-13 TM helices, n={}".format(df_solute_sel.shape[0]))

            y_sel_transporter = [df_sel_transporter.TM01_AAIMON_slope.mean() * 1000, df_sel_transporter.AAIMON_slope_central_TMDs.mean() * 1000, df_sel_transporter.AAIMON_slope_last_TMD.mean() * 1000]
            ax.plot(x, y_sel_transporter, color=c2, label="amino acid, lipid, or neurotransmitter transporters, n={}".format(df_sel_transporter.shape[0]), linestyle="-.")

            y_GPCR = [df_GPCR.TM01_AAIMON_slope.mean() * 1000, df_GPCR.AAIMON_slope_central_TMDs.mean() * 1000, df_GPCR.AAIMON_slope_last_TMD.mean() * 1000]
            ax.plot(x, y_GPCR, color=c3, label="GPCRs, n={}".format(df_GPCR.shape[0]))

            # y_transport = [df_transport.TM01_AAIMON_slope.mean()*1000, df_transport.AAIMON_slope_central_TMDs.mean()*1000, df_transport.AAIMON_slope_last_TMD.mean()*1000]
            # ax.plot(x, y_transport, color = "b", label="other transporters, n={}".format(df_transport.shape[0]), linestyle="--")

            y_ionchannel = [df_ion_channel.TM01_AAIMON_slope.mean() * 1000, df_ion_channel.AAIMON_slope_central_TMDs.mean() * 1000, df_ion_channel.AAIMON_slope_last_TMD.mean() * 1000]
            ax.plot(x, y_ionchannel, color=c4, label="ion channels, n={}".format(df_ion_channel.shape[0]), linestyle="--")

            y_enzyme = [df_enzyme.TM01_AAIMON_slope.mean() * 1000, df_enzyme.AAIMON_slope_central_TMDs.mean() * 1000, df_enzyme.AAIMON_slope_last_TMD.mean() * 1000]
            ax.plot(x, y_enzyme, color=c5, label="enzymes, n={}".format(df_enzyme.shape[0]), linestyle="--")

            y_others = [df_others.TM01_AAIMON_slope.mean() * 1000, df_others.AAIMON_slope_central_TMDs.mean() * 1000, df_others.AAIMON_slope_last_TMD.mean() * 1000]
            ax.plot(x, y_others, color=c6, label="other human multipass proteins, n={}".format(df_others.shape[0]), linestyle=":")

            x_labels = ["first", "central", "last"]
            ax.set_xticks(x)
            ax.set_xticklabels(x_labels)
            ax.set_ylabel(r'm$_{\rm TM/EM} *10^{\rm -3}$')
            ax.set_xlabel("transmembrane domain(s)")
            lgd = ax.legend(frameon=True, bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                      ncol=1, mode="expand", borderaxespad=0.)
            ax.set_xlim(-0.3, 2.3)
            # add annotations
            fontsize = 10
            ax.annotate(s="TM less\nconserved", xy=(-0.12, 0.1), fontsize=fontsize, xytext=None, xycoords='axes fraction', rotation=90)
            ax.annotate(s="TM more\nconserved", xy=(-0.12, 0.9), fontsize=fontsize, xytext=None, xycoords='axes fraction', rotation=90)

            utils.save_figure(fig, Fig_name, base_filepath, save_png, save_pdf, dpi)

            ##################################################################
            #                                                                #
            #                    linechart for LIPOPHILICITY                 #
            #             (first, central, last for protein subgroups)       #
            #                                                                #
            ##################################################################
            # plot the lipophilicity line graphs
            Fig_Nr = 22
            title = 'linechart_cons_f_c_l_multipass_prot_subgroups'
            Fig_name = 'List{:02d}_Fig22b_multipass_linechart_f_c_l_lipo_protein_subgroups'.format(list_number)

            fig = plt.figure()
            ax = fig.add_axes([0.1, 0.1, 0.6, 0.6])
            x = np.arange(3)

            y_solute_sel = [df_solute_sel.TM01_lipo.mean(), df_solute_sel.lipo_mean_central_TMDs.mean(), df_solute_sel.lipo_last_TMD.mean()]
            ax.plot(x, y_solute_sel, color=c1, label="solute carriers with 11-13 TM helices, n={}".format(df_solute_sel.shape[0]))

            y_sel_transporter = [df_sel_transporter.TM01_lipo.mean(), df_sel_transporter.lipo_mean_central_TMDs.mean(), df_sel_transporter.lipo_last_TMD.mean()]
            ax.plot(x, y_sel_transporter, color=c2, label="amino acid, lipid, or neurotransmitter transporters, n={}".format(df_sel_transporter.shape[0]), linestyle="-.")

            y_GPCR = [df_GPCR.TM01_lipo.mean(), df_GPCR.lipo_mean_central_TMDs.mean(), df_GPCR.lipo_last_TMD.mean()]
            ax.plot(x, y_GPCR, color=c3, label="GPCRs, n={}".format(df_GPCR.shape[0]))

            # y_transport = [df_transport.TM01_lipo.mean()*1000, df_transport.lipo_mean_central_TMDs.mean()*1000, df_transport.lipo_last_TMD.mean()*1000]
            # ax.plot(x, y_transport, color = "b", label="other transporters, n={}".format(df_transport.shape[0]), linestyle="--")

            y_ionchannel = [df_ion_channel.TM01_lipo.mean(), df_ion_channel.lipo_mean_central_TMDs.mean(), df_ion_channel.lipo_last_TMD.mean()]
            ax.plot(x, y_ionchannel, color=c4, label="ion channels, n={}".format(df_ion_channel.shape[0]), linestyle="--")

            y_enzyme = [df_enzyme.TM01_lipo.mean(), df_enzyme.lipo_mean_central_TMDs.mean(), df_enzyme.lipo_last_TMD.mean()]
            ax.plot(x, y_enzyme, color=c5, label="enzymes, n={}".format(df_enzyme.shape[0]), linestyle="--")

            y_others = [df_others.TM01_lipo.mean(), df_others.lipo_mean_central_TMDs.mean(), df_others.lipo_last_TMD.mean()]
            ax.plot(x, y_others, color=c6, label="other human multipass proteins, n={}".format(df_others.shape[0]), linestyle=":")

            x_labels = ["first", "central", "last"]
            ax.set_xticks(x)
            ax.set_xticklabels(x_labels)
            ax.set_ylabel('mean lipophilicity\n(Hessa scale)')
            ax.set_xlabel("transmembrane domain(s)")
            ax.legend(frameon=True, bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                      ncol=1, mode="expand", borderaxespad=0.)
            ax.set_xlim(-0.3, 2.3)
            ax.set_ylim(0.05, 0.35)
            # add annotations
            fontsize = 10
            ax.annotate(s="TM more\nlipophilic", xy=(-0.14, 0.1), fontsize=fontsize, xytext=None, xycoords='axes fraction', rotation=90)
            ax.annotate(s="TM less\nlipophilic", xy=(-0.14, 0.9), fontsize=fontsize, xytext=None, xycoords='axes fraction', rotation=90)

            utils.save_figure(fig, Fig_name, base_filepath, save_png, save_pdf, dpi)
    return "~~~~~~~~~~~~        run_save_figures_describing_proteins_in_list is finished        ~~~~~~~~~~~~"

def Fig03_Density_lipo_vs_TM_conservation(list_number, df, letter, suffix, col_list_AAIMON_slope, col_list_lipo, max_evol_distance, base_filepath, save_png, save_pdf, dpi, fontsize):
    Fig_name = 'List{:02d}_Fig03{}_Density_lipo_vs_TM_conservation{}'.format(list_number, letter, suffix)
    title = suffix[1:]
    '''
    data array columns:
    |   0   |   1  |
    | Slope | lipo |
    '''
    data = np.empty([0, 2])

    # select all AAIMON slopes or lipo data
    df_slopes = df.loc[:, col_list_AAIMON_slope]
    df_lipos = df.loc[:, col_list_lipo]
    # check that .stack drops nans, and that there were exactly equal number of nans in the lipo and slope datasets
    if df_slopes.stack().shape != df_lipos.stack().shape:
        raise ValueError("There must be a nan in the lipo or AAIMON slopes. Check code, revert to orig if necessary.")

    # convert slopes and lipos to a 1D numpy array
    y_slopes = df_slopes.stack().values*1000
    x_lipos = df_lipos.stack().values
    x_first_last_dp = [x_lipos.min(), x_lipos.max()]

    # # join to a single numpy array
    # data = np.array([y_slopes, x_lipos]).T

    fig, (cbar_ax, ax) = plt.subplots(2, 1, figsize=(5, 5.5), gridspec_kw={'height_ratios': [0.2, 12]})
    #fontsize = 16
    # number of bins
    n_bins_x = int(max_evol_distance*2)
    n_bins_y = 120
    bins = [n_bins_x, n_bins_y]
    # density threshold
    thresh = 1

    x_border = 1.5
    y_border = 30
    xyrange = [[-x_border, x_border], [-y_border, y_border]]

    # histogram the data
    hh, locx, locy = scipy.histogram2d(x_lipos, y_slopes, range=xyrange, bins=bins)
    hh1 = hh.reshape(1, n_bins_x * n_bins_y)
    hh1 = hh1[hh1 > 0]
    vmax = np.percentile(hh1, 99)
    if vmax % 2 == True:
        vmax = vmax - 1
    # fill the areas with low density by NaNs
    hh[hh < thresh] = np.nan
    im = ax.imshow(np.flipud(hh.T), cmap='Oranges', extent=np.array(xyrange).flatten(),
                   interpolation='none', origin='upper', aspect='auto', vmin=0, vmax=vmax)

    cbar = matplotlib.colorbar.ColorbarBase(cbar_ax, cmap='Oranges', orientation='horizontal')
    if vmax < 10:
        cbar.set_ticks(np.linspace(0, 1, vmax + 1))
        labels = list(range(0, int(vmax), 1))
    else:
        cbar.set_ticks(np.linspace(0, 1, vmax / 2 + 1))
        labels = list(range(0, int(vmax), 2))

    labels.append('>{}'.format(int(vmax)))
    cbar.set_ticklabels(labels)
    cbar_ax.xaxis.set_ticks_position('top')

    # linear regression for data
    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x_lipos, y_slopes)
    fit_fn = np.poly1d(np.array([slope, intercept]))
    y_fitted = fit_fn(x_first_last_dp)
    # plot regression line
    ax.plot(x_first_last_dp, y_fitted, color='k', linewidth=1)
    # annotate regression line formula and R2 value
    ax.annotate(s='$y = {s:.03f}x + {i:.03f}$\n$R^2 = {r_sq:.05f}$'.format(s=slope, i=intercept, r_sq=r_value ** 2),
                xy=(-1.45, 29), fontsize=fontsize - 2, verticalalignment='top')

    ax.set_title(title, fontsize=fontsize)
    ax.set_xlabel('lipophilicity (Hessa scale)', fontsize=fontsize)
    ax.set_ylabel(r'm$_{\rm TM/EM} *10^{\rm -3}$', fontsize=fontsize)
    # set fontsize for axis labels and specify their separation from axis ticks
    ax.tick_params(labelsize=fontsize, pad=3)
    cbar_ax.tick_params(labelsize=fontsize, pad=0)
    # set x and y axis limits to avoid weird limits caused by linear regression
    ax.set_ylim(-y_border, y_border)
    ax.set_xlim(-x_border, x_border)
    # specify space between plot and colorbar
    plt.subplots_adjust(hspace=0.03)
    plt.tight_layout()

    utils.save_figure(fig, Fig_name, base_filepath, save_png, save_pdf, dpi)


def add_bitopic_proteins_to_df_from_another_list(df, list_number_str, bitopic_list_number_str, min_n_homol, pathdict, logging):
    """Add bitopic proteins to the dataframe from another list.

    For the creation of graphs that go from 1 to many TMDs.

    Graph of both bitopic and polytopic proteins is saved under the "polytopic" list of proteins in the summaries folder.

    Previously filtered out GPCRs, but now this is usually handled by altering the protein list.

    Parameters
    ----------
    list_number_str : str
        Protein list number, e.g. "01"
    bitopic_list_number_str : str
        Protein list containing the bitopic data. (e.g. "01" which contains bitopics to be added to list02)
    min_n_homol : int
        Minimum number of homologues
    pathdict : dict
        Dictionary with filepaths
    logging : logging.Logger
        Logger for printing to console and logfile.

    Returns
    -------
    df_filt : pd.DataFrame
        Dataframe of both bitopic and polytopic proteins.
    """
    bitopic_list_csv = pathdict["list_csv"].replace(list_number_str, bitopic_list_number_str)
    bitopic_cr_csv = pathdict["list_cr_summary_csv"].replace(list_number_str, bitopic_list_number_str)
    df_bitopic = pd.read_csv(bitopic_list_csv, index_col=0)
    df_bitopic_cr = pd.read_csv(bitopic_cr_csv, index_col=0)
    df_bitopic_merged = pd.merge(df_bitopic, df_bitopic_cr, left_index=True, right_index=True, suffixes=('_dfc', ''))
    # limit columns to avoid errors from python lists during merging
    df_bitopic_merged = df_bitopic_merged.loc[:, ["number_of_TMDs", "TM01_AAIMON_slope", "AAIMON_slope_central_TMDs", "AAIMON_slope_last_TMD",
                                                  "TM01_lipo", "lipo_mean_central_TMDs", "lipo_last_TMD", 'AAIMON_n_homol']]
    df_bitopic_merged = df_bitopic_merged.loc[df_bitopic_merged['AAIMON_n_homol'] >= min_n_homol]

    df_filt = pd.merge(df, df_bitopic_merged, how="outer")
    logging.info("{} proteins added from List{} for Figure 21, linechart_lipo_f_c_l_vs_number_of_TMDs".format(df_bitopic_merged.shape[0], bitopic_list_number_str))
    return df_filt

