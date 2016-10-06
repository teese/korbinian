import ast
import csv
import korbinian.utils as utils
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os

def compare_rel_con_lists(pathdict, s, logging):
    protein_lists = s["protein_lists"]
    protein_list_names = s["protein_list_names"]

    protein_lists_joined = '_'.join(['%02d' % n for n in protein_lists])
    # output will be saved under data_dir/summaries/compare_lists
    base_path_summ_two_lists = os.path.join(pathdict["data_dir"], "summaries", "compare_lists")
    if os.path.isdir(base_path_summ_two_lists) == False:
        os.mkdir(base_path_summ_two_lists)
    base_filename_summ_two_lists = os.path.join(base_path_summ_two_lists, 'Lists_%s' % protein_lists_joined)

    df_list = []
    for index, list_num in enumerate(protein_lists):
        base_filename_summ = os.path.join(pathdict["data_dir"], "summaries", 'List%02d' % list_num)
        pathdict["dfout09_simap_AAIMON_02"] = '%s_simap_AAIMON_02.csv' % base_filename_summ
        if index == 0:
            df1 = pd.read_csv(pathdict["dfout09_simap_AAIMON_02"], sep=",", quoting=csv.QUOTE_NONNUMERIC, index_col=0)
            df_list.append(df1)
        if index == 1:
            df2 = pd.read_csv(pathdict["dfout09_simap_AAIMON_02"], sep=",", quoting=csv.QUOTE_NONNUMERIC, index_col=0)
            df_list.append(df2)
        if index == 2:
            df3 = pd.read_csv(pathdict["dfout09_simap_AAIMON_02"], sep=",", quoting=csv.QUOTE_NONNUMERIC, index_col=0)
            df_list.append(df3)
        if index == 3:
            df4 = pd.read_csv(pathdict["dfout09_simap_AAIMON_02"], sep=",", quoting=csv.QUOTE_NONNUMERIC, index_col=0)
            df_list.append(df4)
        if index == 4:
            df5 = pd.read_csv(pathdict["dfout09_simap_AAIMON_02"], sep=",", quoting=csv.QUOTE_NONNUMERIC, index_col=0)
            df_list.append(df5)
        if index == 5:
            df6 = pd.read_csv(pathdict["dfout09_simap_AAIMON_02"], sep=",", quoting=csv.QUOTE_NONNUMERIC, index_col=0)
            df_list.append(df6)
        if index == 6:
            logging.warning('ERROR! Too many lists included for analysis in compare_rel_con_lists')

    '''
    The beta-barrel dataset contained a lot of proteins with an average AAIMON of 1.000000. This can only mean that there are not enough homologues.
    The singlepass dataset contained only 2 proteins, the alpha-helicas multipass only 8 with 1.000000.
    All these are excluded from the dataset (and all following graphs). Note that it would be better to exactly count the valid homologues, rather than exclude them afterwards like this.
    '''
    # create new list of dataframes
    df_list_excluding_AAIMON_ones = []
    for dfl in df_list:
        # count numbers of each AAIMON ratio
        vc_AAIMON = dfl.AAIMON_ratio_mean_all_TMDs.value_counts()
        # replace dfl with a filtered dataframe, with all rows excluded where AAIMON_ratio_mean_all_TMDs is 1.000000
        if 1.000000 in vc_AAIMON:
            num_proteins_with_AAIMON_of_ONE = vc_AAIMON[1.000000]
            total_num_prot_with_data = len(dfl['AAIMON_ratio_mean_all_TMDs'].dropna())
            logging.info('num_proteins_with_AAIMON_of_ONE in orig dataframe : %i from %i total' % (
            num_proteins_with_AAIMON_of_ONE, total_num_prot_with_data))
            dfl = dfl.loc[dfl['AAIMON_ratio_mean_all_TMDs'] != 1.000000]
        df_list_excluding_AAIMON_ones.append(dfl)
    df_list = df_list_excluding_AAIMON_ones

    '''
    Prepare fonts, colours etc for following figures
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

    colour_lists = utils.create_colour_lists()
    TUM_colours_list_with_greys = colour_lists['TUM_colours_list_with_greys']

    '''
    Fig01: Histogram of mean AAIMON ratios
    '''
    Fig_Nr = 1
    title = 'Histogram of mean AAIMON ratios'
    newfig, savefig, fig_nr, plot_nr_in_fig, row_nr, col_nr = dict_organising_subplots[Fig_Nr]
    # if the plot is the last one, the figure should be saved
    # if a new figure should be created (either because the orig is full, or the last TMD is analysed)
    if newfig:
        # clear any previous plots
        plt.close('all')
        # create a new figure
        fig, axarr = plt.subplots(nrows=nrows_in_each_fig, ncols=ncols_in_each_fig, dpi=300)  # sharex=True

    for n, dfl in enumerate(df_list):
        data_column = 'AAIMON_ratio_mean_all_TMDs'
        color = TUM_colours_list_with_greys[n + 2]  # "#0489B1"
        alpha = 0.7
        col_width_value = 0.95
        ylabel = 'freq'
        xlabel = 'average conservation ratio (membranous over nonmembranous)'
        utils.create_hist_from_df_col(df=dfl,
                                      title=title, row_nr=row_nr, col_nr=col_nr, axarr=axarr,
                                      settings=s,
                                      data_column=data_column,
                                      color=color,
                                      alpha=alpha,
                                      col_width_value=col_width_value,
                                      fontsize=fontsize,
                                      xlabel=xlabel,
                                      ylabel=ylabel, legend=protein_list_names
                                      )

    Fig_Nr = 2
    title = 'Histogram of mean AASMON ratios'
    newfig, savefig, fig_nr, plot_nr_in_fig, row_nr, col_nr = dict_organising_subplots[Fig_Nr]
    # if the plot is the last one, the figure should be saved
    # if a new figure should be created (either because the orig is full, or the last TMD is analysed)
    if newfig:
        # clear any previous plots
        plt.close('all')
        # create a new figure
        fig, axarr = plt.subplots(nrows=nrows_in_each_fig, ncols=ncols_in_each_fig, dpi=300)  # sharex=True

    for n, dfl in enumerate(df_list):
        data_column = 'AASMON_ratio_mean_all_TMDs'
        color = TUM_colours_list_with_greys[n + 2]  # "#0489B1"
        alpha = 0.7
        col_width_value = 0.95
        ylabel = 'freq'
        xlabel = 'average conservation ratio (membranous over nonmembranous)'
        utils.create_hist_from_df_col(df=dfl,
                                      title=title, row_nr=row_nr, col_nr=col_nr, axarr=axarr,
                                      settings=s,
                                      data_column=data_column,
                                      color=color,
                                      alpha=alpha,
                                      col_width_value=col_width_value,
                                      fontsize=fontsize,
                                      xlabel=xlabel,
                                      ylabel=ylabel, legend=protein_list_names
                                      )
    Fig_Nr = 3
    title = 'len_nonTMD_seq_match_mean'
    newfig, savefig, fig_nr, plot_nr_in_fig, row_nr, col_nr = dict_organising_subplots[Fig_Nr]
    # if the plot is the last one, the figure should be saved
    # if a new figure should be created (either because the orig is full, or the last TMD is analysed)
    if newfig:
        # clear any previous plots
        plt.close('all')
        # create a new figure
        fig, axarr = plt.subplots(nrows=nrows_in_each_fig, ncols=ncols_in_each_fig, dpi=300)  # sharex=True

    for n, dfl in enumerate(df_list):
        data_column = 'len_nonTMD_seq_match_mean'
        color = TUM_colours_list_with_greys[n + 2]  # "#0489B1"
        alpha = 0.7
        col_width_value = 0.95
        ylabel = 'freq'
        xlabel = data_column

        utils.create_hist_from_df_col_with_auto_binlist(df=dfl,
                                                        title=title, row_nr=row_nr, col_nr=col_nr, axarr=axarr, num_bins=60,
                                                        settings=s,
                                                        data_column=data_column,
                                                        color=color,
                                                        alpha=alpha,
                                                        col_width_value=col_width_value,
                                                        fontsize=fontsize,
                                                        xlabel=xlabel,
                                                        ylabel=ylabel, legend=protein_list_names
                                                        )
    Fig_Nr = 4
    title = 'len_nonTMD_align_mean'
    newfig, savefig, fig_nr, plot_nr_in_fig, row_nr, col_nr = dict_organising_subplots[Fig_Nr]
    # if the plot is the last one, the figure should be saved
    # if a new figure should be created (either because the orig is full, or the last TMD is analysed)
    if newfig:
        # clear any previous plots
        plt.close('all')
        # create a new figure
        fig, axarr = plt.subplots(nrows=nrows_in_each_fig, ncols=ncols_in_each_fig, dpi=300)  # sharex=True

    for n, dfl in enumerate(df_list):
        data_column = 'len_nonTMD_align_mean'
        color = TUM_colours_list_with_greys[n + 2]  # "#0489B1"
        alpha = 0.7
        col_width_value = 0.95
        ylabel = 'freq'
        xlabel = data_column

        utils.create_hist_from_df_col_with_auto_binlist(df=dfl,
                                                        title=title, row_nr=row_nr, col_nr=col_nr, axarr=axarr,
                                                        num_bins=100,
                                                        settings=s,
                                                        data_column=data_column,
                                                        color=color,
                                                        alpha=alpha,
                                                        col_width_value=col_width_value,
                                                        fontsize=fontsize,
                                                        xlabel=xlabel,
                                                        ylabel=ylabel, legend=protein_list_names
                                                        )
    # save the figure as it is
    if savefig:
        utils.save_fig_with_subplots(fig=fig, axarr=axarr, base_filename=base_filename_summ_two_lists,
                                     fig_nr=fig_nr, fontsize=fontsize)

    Fig_Nr = 5
    title = 'nonTMD_perc_ident_mean'
    newfig, savefig, fig_nr, plot_nr_in_fig, row_nr, col_nr = dict_organising_subplots[Fig_Nr]
    # if the plot is the last one, the figure should be saved
    # if a new figure should be created (either because the orig is full, or the last TMD is analysed)
    if newfig:
        # close any old figures
        plt.close('all')
        # create a new figure
        fig, axarr = plt.subplots(nrows=nrows_in_each_fig, ncols=ncols_in_each_fig, dpi=300)  # sharex=True

    for n, dfl in enumerate(df_list):
        data_column = 'nonTMD_perc_ident_mean'
        color = TUM_colours_list_with_greys[n + 2]  # "#0489B1"
        alpha = 0.7
        col_width_value = 0.95
        ylabel = 'freq'
        xlabel = data_column

        utils.create_hist_from_df_col_with_auto_binlist(df=dfl,
                                                        title=title, row_nr=row_nr, col_nr=col_nr, axarr=axarr, num_bins=60,
                                                        settings=s,
                                                        data_column=data_column,
                                                        color=color,
                                                        alpha=alpha,
                                                        col_width_value=col_width_value,
                                                        fontsize=fontsize,
                                                        xlabel=xlabel,
                                                        ylabel=ylabel, legend=protein_list_names
                                                        )

    Fig_Nr = Fig_Nr + 1
    title = 'TMD_perc_ident_mean'
    newfig, savefig, fig_nr, plot_nr_in_fig, row_nr, col_nr = dict_organising_subplots[Fig_Nr]
    # if the plot is the last one, the figure should be saved
    # if a new figure should be created (either because the orig is full, or the last TMD is analysed)
    if newfig:
        # close any old figures
        plt.close('all')
        # create a new figure
        fig, axarr = plt.subplots(nrows=nrows_in_each_fig, ncols=ncols_in_each_fig, dpi=300)  # sharex=True

    for n, dfl in enumerate(df_list):
        data_column = 'nonTMD_perc_ident_mean'
        color = TUM_colours_list_with_greys[n + 2]  # "#0489B1"
        alpha = 0.7
        col_width_value = 0.95
        ylabel = 'freq'
        xlabel = data_column

        utils.create_hist_from_df_col_with_auto_binlist(df=dfl,
                                                        title=title, row_nr=row_nr, col_nr=col_nr, axarr=axarr, num_bins=60,
                                                        settings=s,
                                                        data_column=data_column,
                                                        color=color,
                                                        alpha=alpha,
                                                        col_width_value=col_width_value,
                                                        fontsize=fontsize,
                                                        xlabel=xlabel,
                                                        ylabel=ylabel, legend=protein_list_names
                                                        )
    Fig_Nr = Fig_Nr + 1
    title = 'nonTMD_perc_sim_mean'
    newfig, savefig, fig_nr, plot_nr_in_fig, row_nr, col_nr = dict_organising_subplots[Fig_Nr]
    # if the plot is the last one, the figure should be saved
    # if a new figure should be created (either because the orig is full, or the last TMD is analysed)
    if newfig:
        # clear any previous plots
        plt.close('all')
        # create a new figure
        fig, axarr = plt.subplots(nrows=nrows_in_each_fig, ncols=ncols_in_each_fig, dpi=300)  # sharex=True

    for n, dfl in enumerate(df_list):
        data_column = 'nonTMD_perc_sim_mean'
        color = TUM_colours_list_with_greys[n + 2]  # "#0489B1"
        alpha = 0.7
        col_width_value = 0.95
        ylabel = 'freq'
        xlabel = data_column

        utils.create_hist_from_df_col_with_auto_binlist(df=dfl,
                                                        title=title, row_nr=row_nr, col_nr=col_nr, axarr=axarr, num_bins=60,
                                                        settings=s,
                                                        data_column=data_column,
                                                        color=color,
                                                        alpha=alpha,
                                                        col_width_value=col_width_value,
                                                        fontsize=fontsize,
                                                        xlabel=xlabel,
                                                        ylabel=ylabel, legend=protein_list_names
                                                        )

    Fig_Nr = Fig_Nr + 1
    title = 'nonTMD_perc_sim_plus_ident_mean'
    newfig, savefig, fig_nr, plot_nr_in_fig, row_nr, col_nr = dict_organising_subplots[Fig_Nr]
    # if the plot is the last one, the figure should be saved
    # if a new figure should be created (either because the orig is full, or the last TMD is analysed)
    if newfig:
        # clear any previous plots
        plt.close('all')
        # create a new figure
        fig, axarr = plt.subplots(nrows=nrows_in_each_fig, ncols=ncols_in_each_fig, dpi=300)  # sharex=True

    for n, dfl in enumerate(df_list):
        data_column = 'nonTMD_perc_sim_plus_ident_mean'
        color = TUM_colours_list_with_greys[n + 2]  # "#0489B1"
        alpha = 0.7
        col_width_value = 0.95
        ylabel = 'freq'
        xlabel = data_column

        utils.create_hist_from_df_col_with_auto_binlist(df=dfl,
                                                        title=title, row_nr=row_nr, col_nr=col_nr, axarr=axarr, num_bins=60,
                                                        settings=s,
                                                        data_column=data_column,
                                                        color=color,
                                                        alpha=alpha,
                                                        col_width_value=col_width_value,
                                                        fontsize=fontsize,
                                                        xlabel=xlabel,
                                                        ylabel=ylabel, legend=protein_list_names
                                                        )

    Fig_Nr = Fig_Nr + 1
    title = 'nonTMD_qm_gaps_per_q_residue_mean'
    newfig, savefig, fig_nr, plot_nr_in_fig, row_nr, col_nr = dict_organising_subplots[Fig_Nr]
    # if the plot is the last one, the figure should be saved
    # if a new figure should be created (either because the orig is full, or the last TMD is analysed)
    if newfig:
        # clear any previous plots
        plt.close('all')
        # create a new figure
        fig, axarr = plt.subplots(nrows=nrows_in_each_fig, ncols=ncols_in_each_fig, dpi=300)  # sharex=True

    for n, dfl in enumerate(df_list):
        data_column = 'nonTMD_qm_gaps_per_q_residue_mean'
        color = TUM_colours_list_with_greys[n + 2]  # "#0489B1"
        alpha = 0.7
        col_width_value = 0.95
        ylabel = 'freq'
        xlabel = data_column

        utils.create_hist_from_df_col_with_auto_binlist(df=dfl,
                                                        title=title, row_nr=row_nr, col_nr=col_nr, axarr=axarr, num_bins=60,
                                                        settings=s,
                                                        data_column=data_column,
                                                        color=color,
                                                        alpha=alpha,
                                                        col_width_value=col_width_value,
                                                        fontsize=fontsize,
                                                        xlabel=xlabel,
                                                        ylabel=ylabel, legend=protein_list_names
                                                        )

    # save the figure as it is
    if savefig:
        utils.save_fig_with_subplots(fig=fig, axarr=axarr, base_filename=base_filename_summ_two_lists,
                                     fig_nr=fig_nr, fontsize=fontsize)

    Fig_Nr = Fig_Nr + 1
    title = 'TM01_AAIMON_ratio_mean'
    newfig, savefig, fig_nr, plot_nr_in_fig, row_nr, col_nr = dict_organising_subplots[Fig_Nr]
    # if the plot is the last one, the figure should be saved
    # if a new figure should be created (either because the orig is full, or the last TMD is analysed)
    if newfig:
        # clear any previous plots
        plt.close('all')
        # create a new figure
        fig, axarr = plt.subplots(nrows=nrows_in_each_fig, ncols=ncols_in_each_fig, dpi=300)  # sharex=True

    for n, dfl in enumerate(df_list):
        data_column = 'TM01_AAIMON_ratio_mean'
        color = TUM_colours_list_with_greys[n + 2]  # "#0489B1"
        alpha = 0.7
        col_width_value = 0.95
        ylabel = 'freq'
        xlabel = data_column

        utils.create_hist_from_df_col_with_auto_binlist(df=dfl,
                                                        title=title, row_nr=row_nr, col_nr=col_nr, axarr=axarr, num_bins=60,
                                                        settings=s,
                                                        data_column=data_column,
                                                        color=color,
                                                        alpha=alpha,
                                                        col_width_value=col_width_value,
                                                        fontsize=fontsize,
                                                        xlabel=xlabel,
                                                        ylabel=ylabel, legend=protein_list_names
                                                        )

    Fig_Nr = Fig_Nr + 1
    title = 'TM01_perc_ident_mean'
    newfig, savefig, fig_nr, plot_nr_in_fig, row_nr, col_nr = dict_organising_subplots[Fig_Nr]
    # if the plot is the last one, the figure should be saved
    # if a new figure should be created (either because the orig is full, or the last TMD is analysed)
    if newfig:
        # clear any previous plots
        plt.close('all')
        # create a new figure
        fig, axarr = plt.subplots(nrows=nrows_in_each_fig, ncols=ncols_in_each_fig, dpi=300)  # sharex=True

    for n, dfl in enumerate(df_list):
        data_column = 'TM01_perc_ident_mean'
        if data_column in df.columns:
            color = TUM_colours_list_with_greys[n + 2]  # "#0489B1"
            alpha = 0.7
            col_width_value = 0.95
            ylabel = 'freq'
            xlabel = data_column

            utils.create_hist_from_df_col_with_auto_binlist(df=dfl,
                                                            title=title, row_nr=row_nr, col_nr=col_nr, axarr=axarr,
                                                            num_bins=60,
                                                            settings=s,
                                                            data_column=data_column,
                                                            color=color,
                                                            alpha=alpha,
                                                            col_width_value=col_width_value,
                                                            fontsize=fontsize,
                                                            xlabel=xlabel,
                                                            ylabel=ylabel, legend=protein_list_names
                                                            )
    Fig_Nr = Fig_Nr + 1
    title = 'TM01_start'
    newfig, savefig, fig_nr, plot_nr_in_fig, row_nr, col_nr = dict_organising_subplots[Fig_Nr]
    # if the plot is the last one, the figure should be saved
    # if a new figure should be created (either because the orig is full, or the last TMD is analysed)
    if newfig:
        # clear any previous plots
        plt.close('all')
        # create a new figure
        fig, axarr = plt.subplots(nrows=nrows_in_each_fig, ncols=ncols_in_each_fig, dpi=300)  # sharex=True

    for n, df in enumerate(df_list):
        if 'TM01_start' in df.columns:
            data_column = 'TM01_start'
        elif 'TM01_start_in_SW_alignment' in df.columns:
            data_column = 'TM01_start_in_SW_alignment'
        else:
            ValueError('Neither TM01_start nor TM01_start_in_SW_alignment are in df.columns')

        color = TUM_colours_list_with_greys[n + 2]  # "#0489B1"
        alpha = 0.7
        col_width_value = 0.95
        ylabel = 'freq'
        xlabel = data_column

        utils.create_hist_from_df_col_with_auto_binlist(df=dfl,
                                                        title=title, row_nr=row_nr, col_nr=col_nr, axarr=axarr, num_bins=60,
                                                        settings=s,
                                                        data_column=data_column,
                                                        color=color,
                                                        alpha=alpha,
                                                        col_width_value=col_width_value,
                                                        fontsize=fontsize,
                                                        xlabel=xlabel,
                                                        ylabel=ylabel, legend=protein_list_names
                                                        )

    Fig_Nr = Fig_Nr + 1
    title = 'number_of_TMDs'
    newfig, savefig, fig_nr, plot_nr_in_fig, row_nr, col_nr = dict_organising_subplots[Fig_Nr]
    # if the plot is the last one, the figure should be saved
    # if a new figure should be created (either because the orig is full, or the last TMD is analysed)
    if newfig:
        # clear any previous plots
        plt.close('all')
        # create a new figure
        fig, axarr = plt.subplots(nrows=nrows_in_each_fig, ncols=ncols_in_each_fig, dpi=300)  # sharex=True

    for n, df in enumerate(df_list):
        data_column = 'number_of_TMDs'
        if data_column not in df.columns:
            data_column = 'number_of_TMDs'
        color = TUM_colours_list_with_greys[n + 2]  # "#0489B1"
        alpha = 0.7
        col_width_value = 0.95
        ylabel = 'freq'
        xlabel = data_column

        utils.create_hist_from_df_col_with_auto_binlist(df=dfl,
                                                        title=title, row_nr=row_nr, col_nr=col_nr, axarr=axarr, num_bins=60,
                                                        settings=s,
                                                        data_column=data_column,
                                                        color=color,
                                                        alpha=alpha,
                                                        col_width_value=col_width_value,
                                                        fontsize=fontsize,
                                                        xlabel=xlabel,
                                                        ylabel=ylabel, legend=protein_list_names
                                                        )

    Fig_Nr = Fig_Nr + 1
    title = 'TM01_ratio_length_of_query_TMD_to_rest_of_match_protein_mean'
    newfig, savefig, fig_nr, plot_nr_in_fig, row_nr, col_nr = dict_organising_subplots[Fig_Nr]
    # if the plot is the last one, the figure should be saved
    # if a new figure should be created (either because the orig is full, or the last TMD is analysed)
    if newfig:
        # clear any previous plots
        plt.close('all')
        # create a new figure
        fig, axarr = plt.subplots(nrows=nrows_in_each_fig, ncols=ncols_in_each_fig, dpi=300)  # sharex=True

    for n, df in enumerate(df_list):
        data_column = 'TM01_ratio_length_of_query_TMD_to_rest_of_match_protein_mean'
        color = TUM_colours_list_with_greys[n + 2]  # "#0489B1"
        alpha = 0.7
        col_width_value = 0.95
        ylabel = 'freq'
        xlabel = data_column

        utils.create_hist_from_df_col_with_auto_binlist(df=dfl,
                                                        title=title, row_nr=row_nr, col_nr=col_nr, axarr=axarr, num_bins=60,
                                                        settings=s,
                                                        data_column=data_column,
                                                        color=color,
                                                        alpha=alpha,
                                                        col_width_value=col_width_value,
                                                        fontsize=fontsize,
                                                        xlabel=xlabel,
                                                        ylabel=ylabel, legend=protein_list_names
                                                        )
    # save the figure as it is
    if savefig:
        utils.save_fig_with_subplots(fig=fig, axarr=axarr, base_filename=base_filename_summ_two_lists,
                                     fig_nr=fig_nr, fontsize=fontsize)
    Fig_Nr = Fig_Nr + 1
    title = 'SIMAP_total_hits'
    newfig, savefig, fig_nr, plot_nr_in_fig, row_nr, col_nr = dict_organising_subplots[Fig_Nr]
    # if the plot is the last one, the figure should be saved
    # if a new figure should be created (either because the orig is full, or the last TMD is analysed)
    if newfig:
        # clear any previous plots
        plt.close('all')
        # create a new figure
        fig, axarr = plt.subplots(nrows=nrows_in_each_fig, ncols=ncols_in_each_fig, dpi=300)  # sharex=True

    for n, df in enumerate(df_list):
        data_column = 'SIMAP_total_hits'
        color = TUM_colours_list_with_greys[n + 2]  # "#0489B1"
        alpha = 0.7
        col_width_value = 0.95
        ylabel = 'freq'
        xlabel = data_column

        utils.create_hist_from_df_col_with_auto_binlist(df=dfl,
                                                        title=title, row_nr=row_nr, col_nr=col_nr, axarr=axarr, num_bins=60,
                                                        settings=s,
                                                        data_column=data_column,
                                                        color=color,
                                                        alpha=alpha,
                                                        col_width_value=col_width_value,
                                                        fontsize=fontsize,
                                                        xlabel=xlabel,
                                                        ylabel=ylabel, legend=protein_list_names
                                                        )

    #
    #    Fig_Nr = 14
    #    title = 'number_of_valid_hits'
    #    newfig, savefig, fig_nr, plot_nr_in_fig, row_nr, col_nr = dict_organising_subplots[Fig_Nr]
    #    #if the plot is the last one, the figure should be saved
    #    #if a new figure should be created (either because the orig is full, or the last TMD is analysed)
    #    if newfig:
    #        #clear any previous plots
    #        plt.close('all')
    #        #create a new figure
    #        fig, axarr = plt.subplots(nrows=nrows_in_each_fig, ncols=ncols_in_each_fig, dpi=300)  # sharex=True
    #
    #    for n, df in enumerate(df_list):
    #        data_column = 'number_of_valid_hits'
    #        color = TUM_colours_list_with_greys[n+2]#"#0489B1"
    #        alpha = 0.7
    #        col_width_value = 0.95
    #        ylabel = 'freq'
    #        xlabel = data_column
    #
    #        utils.create_hist_from_df_col_with_auto_binlist(df=dfl,
    #                                                title=title,row_nr=row_nr,col_nr=col_nr,axarr=axarr,num_bins=60,
    #                                                settings=s,
    #                                                data_column = data_column,
    #                                                color=color,
    #                                                alpha=alpha,
    #                                                col_width_value=col_width_value,
    #                                                fontsize=fontsize,
    #                                                xlabel=xlabel,
    #                                                ylabel=ylabel,legend=protein_list_names
    #                                                )

    #    #save the figure as it is
    #    if savefig:
    #        utils.save_fig_with_subplots(fig=fig, axarr = axarr, base_filename = base_filename_summ_two_lists,
    #                                     fig_nr = fig_nr, fontsize=fontsize)
    #


    # save the figure as it is
    savefig = True
    if savefig:
        utils.save_fig_with_subplots(fig=fig, axarr=axarr, base_filename=base_filename_summ_two_lists,
                                     fig_nr=fig_nr, fontsize=fontsize)
    logging.info('A10_compare_lists is finished')
    logging.info('df1 AAIMON_ratio_mean_all_TMDs : %0.5f' % df_list[0]['AAIMON_ratio_mean_all_TMDs'].mean())
    logging.info('df1 AASMON_ratio_mean_all_TMDs : %0.5f' % df_list[0]['AASMON_ratio_mean_all_TMDs'].mean())
    logging.info('df2 AAIMON_ratio_mean_all_TMDs : %0.5f' % df_list[1]['AAIMON_ratio_mean_all_TMDs'].mean())
    logging.info('df2 AASMON_ratio_mean_all_TMDs : %0.5f' % df_list[1]['AASMON_ratio_mean_all_TMDs'].mean())
    if len(df_list) > 2:
        logging.info('df3 AAIMON_ratio_mean_all_TMDs : %0.5f' % df_list[2]['AAIMON_ratio_mean_all_TMDs'].mean())
        logging.info('df3 AASMON_ratio_mean_all_TMDs : %0.5f' % df_list[2]['AASMON_ratio_mean_all_TMDs'].mean())
    if len(df_list) > 3:
        logging.info('df4 AAIMON_ratio_mean_all_TMDs : %0.5f' % df_list[3]['AAIMON_ratio_mean_all_TMDs'].mean())
        logging.info('df4 AASMON_ratio_mean_all_TMDs : %0.5f' % df_list[3]['AASMON_ratio_mean_all_TMDs'].mean())

    '''
    SUMMARY of proportions of TM, JM regions, etc.
    SHOULD BE SHIFTED TO AN EARLIER SCRIPT, WHICH DOES NOT REQUIRE A LIST OF DATAFRAMES
    GIVES A SettingWithCopyWarning: WHICH MAY, OR MAY NOT BE JUSTIFIED.
    '''
    logging.info("starting mean calculations")
    for n, dfl in enumerate(df_list):
        for acc in dfl.index:
            dict_len_TMD = {}
            for TMD in ast.literal_eval(dfl.loc[acc, 'list_of_TMDs']):
                if '%s_len'%TMD in dfl.columns:
                    dict_len_TMD[TMD] = dfl.loc[acc, '%s_len'%TMD]
                elif '%s_start'%TMD in dfl.columns:
                    dict_len_TMD[TMD] = dfl.loc[acc, '%s_end'%TMD] - dfl.loc[acc, '%s_start'%TMD]
                else:
                    raise IndexError('neither TM01_len nor TM01_start are in dfl.columns')
            dfl.loc[acc, 'len_all_TMD_region'] = np.sum(list(dict_len_TMD.values()))
            # series_len_all_TMD_region = np.sum(list(dict_len_TMD.values()))
            dfl.loc[acc, 'average_len_all_TMDs'] = np.mean(list(dict_len_TMD.values()))
            # series_average_len_all_TMDs = np.mean(list(dict_len_TMD.values()))

        dfl.loc[:, 'ratio_len_all_TMD_to_seqlen'] = dfl['len_all_TMD_region'] / dfl.len_query_align_seq
        dfl.loc[:,
        'proportion_JM_region_3aa'] = dfl.number_of_TMDs * 3 * 2 / dfl.len_query_align_seq
        dfl.loc[:,
        'proportion_JM_region_6aa'] = dfl.number_of_TMDs * 6 * 2 / dfl.len_query_align_seq
        print(
            "dfl%i\nAverage AAIMON ratio all TMDs = %0.3f\nAverage len TMD = %0.2f\nAverage len full sequence = %0.2f\nAverage TMD proportion = %0.4f (%0.2f%%)\nAverage proportion JM (3aa each side) = %0.4f (%0.2f%%)\nAverage proportion JM (6aa each side) = %0.4f (%0.2f%%)" % (
                n,
                dfl['AAIMON_ratio_mean_all_TMDs'].mean(),
                dfl.len_all_TMD_region.mean(),
                dfl.len_query_align_seq.mean(),
                dfl['ratio_len_all_TMD_to_seqlen'].mean(),
                dfl['ratio_len_all_TMD_to_seqlen'].mean() * 100,
                dfl['proportion_JM_region_3aa'].mean(),
                dfl['proportion_JM_region_3aa'].mean() * 100,
                dfl['proportion_JM_region_6aa'].mean(),
                dfl['proportion_JM_region_6aa'].mean() * 100
            ))

    # for ax in axarr.flat:
    #            #change axis font size
    #            ax.tick_params(labelsize = fontsize)
    #            #hide spines
    #            ax.spines["top"].set_visible(False)
    #            ax.spines["right"].set_visible(False)
    #            ax.tick_params(
    #                axis='x',          # changes apply to the x-axis
    #                which='both',      # both major and minor ticks are affected
    #                bottom='off',      # ticks along the bottom edge are off
    #                top='off',         # ticks along the top edge are off
    #                labelbottom='on') # labels along the bottom edge are off
    #        #automatically tighten the layout of plots in the figure
    #        fig.subplots_adjust(bottom = 0)
    #        fig.subplots_adjust(top = 1)
    #        fig.subplots_adjust(right = 1)
    #        fig.subplots_adjust(left = 0)
    #        #save files
    #        fig.savefig(base_filename_summ_two_lists + '_%01d.png' % fig_nr, format='png', dpi=400)
    #        fig.savefig(base_filename_summ_two_lists + '_%01d.pdf' % fig_nr, format='pdf')

    '''+++++++++++++++TMD CONSERVATION (OLD)++++++++++++++++++'''