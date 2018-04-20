import pandas as pd
import numpy as np
import ast
import csv
import sys
import itertools
from scipy.stats import ttest_ind
import korbinian.utils as utils
import matplotlib.pyplot as plt
import os
# import debugging tools
from korbinian.utils import pr, pc, pn, aaa

def keyword_analysis(pathdict, s, logging):
    """


    Martin to write docstring and drop proteins with not enough homologues.

    Parameters
    ----------
    pathdict
    s
    logging

    Returns
    -------

    """
    logging.info("~~~~~~~~~~~~                      starting keyword_analysis                         ~~~~~~~~~~~~")
    # load summary file
    df_list = pd.read_csv(pathdict["list_csv"], sep=",", quoting=csv.QUOTE_NONNUMERIC, index_col=0, low_memory=False)
    if df_list.shape[0] < s["min_n_proteins_in_list"]:
        return "~~~~~~~~~~~~           keyword_analysis skipped, only {} proteins in list            ~~~~~~~~~~~~".format(df_list.shape[0])

    if utils.file_is_old(pathdict["list_cr_summary_csv"], s["oldest_acceptable_file_date"]):
        raise ValueError("{} is too old for analysis".format(pathdict["list_cr_summary_csv"]))

    # load cr_summary file
    df_cr_summary = pd.read_csv(pathdict["list_cr_summary_csv"], sep=",", quoting=csv.QUOTE_NONNUMERIC, index_col=0, low_memory=False)
    # merge cr_summary and summary file, if columns are equal in both files, suffix _dfc will be added in cr_summary column names for backwards compatibility
    df = pd.merge(df_cr_summary, df_list, left_index=True, right_index=True, suffixes=('_dfc', ''))

    # skip the keyword analysis if there are no keywords
    # if it's a betabarrel dataset, take the Pfam IDs as keywords
    is_betabarrel_dataset = False
    if 'uniprot_KW' not in df.columns:
        if "betabarrel" in df.columns and "Pfam_ID" in df.columns:
            if True in df.betabarrel.tolist():
                is_betabarrel_dataset = True
                df["uniprot_KW"] = df["Pfam_ID"]
        else:
            return "Keyword analysis not conducted. No keywords found in protein summary file."

    # drop proteins with less then x homologues
    n_before = df.shape[0]
    df = df[df.AAIMON_n_homol >= s['min_homol']]
    n_after = df.shape[0]
    n_dropped = n_before - n_after
    logging.info('{} proteins dropped due to insufficient homologues\nOrig = {}, Final = {}, cutoff = {} homologues\n'.format(n_dropped, n_before, n_after, s['min_homol']))

    # create folder in list summary directory to hold keyword data
    if not os.path.exists(pathdict["keywords"]):
        os.makedirs(pathdict["keywords"])
    # base filepath for figures
    base_filepath = os.path.join(pathdict["keywords"], 'histograms')
    #save figures to .pdf or .png
    save_png = s["save_png"]
    save_pdf = s["save_pdf"]
    # max number of correlated keywords in keywords summary dataframe
    max_n_kw_in_df = 3 #s['cutoff_number_correlated_keywords']
    # create binlist for histograms
    smallest_bin = -0.04
    largest_bin = 0.04
    number_of_bins = 51
    binlist = np.linspace(smallest_bin, largest_bin, number_of_bins)
    # initialise basic settings for figures
    plt.style.use('seaborn-whitegrid')
    fontsize = 12
    alpha = 0.8


    # define list of ignored and enzyme keywords
    list_enzyme_KW, list_ignored_KW, PFAM_dict = get_list_enzyme_KW_and_list_ignored_KW()

    # remove proteins containing nan in AAIMON_slope_all_TMDs_mean
    list_of_acc_without_nan = []
    for acc in df.index:
        if not pd.isnull(df.loc[acc, 'AAIMON_slope_all_TMDs_mean']):
            list_of_acc_without_nan.append(acc)
    df = df.loc[list_of_acc_without_nan, :]

    # to avoid reprocessing of all keywords, a new column with removed ignored and replaced enzyme keywords is created and saved
    # if this column is already present in any dataframe, re-creation is skipped
    if not 'uniprot_KW_for_analysis' in df.columns:
        df['uniprot_KW_for_analysis'] = df['uniprot_KW']

        # apply ast.literal_eval to every item in df['uniprot_KW_for_analysis']
        if isinstance(df['uniprot_KW_for_analysis'][0], str):
            df['uniprot_KW_for_analysis'] = df['uniprot_KW_for_analysis'].dropna().apply(lambda x: ast.literal_eval(x))

        if not True in df.betabarrel.tolist():
            # check if protein is an Enzyme or GPCR
            df['enzyme'] = df['uniprot_KW_for_analysis'].dropna().apply(KW_list_contains_any_desired_KW, args=(list_enzyme_KW,))

            # DEPRECATED. SHOULD BE IN PROT_LIST
            # # check if protein is a GPCR
            # list_GPCR_KW = ['G-protein coupled receptor']
            # df['GPCR'] = df['uniprot_KW_for_analysis'].apply(KW_list_contains_any_desired_KW, args=(list_GPCR_KW,))

            # remove ignored keywords; replace Enzyme keywords with single keyword 'Enzyme
            sys.stdout.write('removing ignored keywords; replacing enzyme associated keywords with "Enzyme"\n')
            #n = 0

            for n, acc in enumerate(df['uniprot_KW_for_analysis'].dropna().index):
                n += 1
                if n % 20 == 0:
                    sys.stdout.write('.'), sys.stdout.flush()
                    if n % 600 == 0:
                        sys.stdout.write('\n'), sys.stdout.flush()
                # remove ignored keywords from dataframe 'uniprot_KW_for_analysis'
                for element in list_ignored_KW:
                    if element in df.loc[acc, 'uniprot_KW_for_analysis']:
                        df.loc[acc, 'uniprot_KW_for_analysis'].remove(element)
                # replace keywords associated with enzymes with keyword 'Enzyme'
                for element in list_enzyme_KW:
                    if element in df.loc[acc, 'uniprot_KW_for_analysis']:
                        df.loc[acc, 'uniprot_KW_for_analysis'].remove(element)
                if df.loc[acc, 'enzyme']:
                    df.loc[acc, 'uniprot_KW_for_analysis'].append('Enzyme')
            else:
                df['enzyme'] = False

        df_list['uniprot_KW_for_analysis'] = df['uniprot_KW_for_analysis']
        df_list.to_csv(pathdict["list_csv"], sep=",", quoting=csv.QUOTE_NONNUMERIC)

    else:
        # apply ast.literal_eval to every item in df['uniprot_KW_for_analysis']
        if isinstance(df['uniprot_KW_for_analysis'][0], str):
            df['uniprot_KW_for_analysis'] = df['uniprot_KW_for_analysis'].dropna().apply(lambda x: ast.literal_eval(x))

    n_before = df.shape[0]
    df.dropna(subset=['uniprot_KW_for_analysis'], inplace=True)
    n_after = df.shape[0]
    n_dropped = n_before - n_after
    if n_dropped > 0:
        logging.info("\n{} proteins dropped due to empty uniprot_KW_for_analysis".format(n_dropped))

    # check if dataset is SP or MP
    dataset = 'none'
    if df.number_of_TMDs_excl_SP.max() == 1:
        dataset = 'SP'
    elif df.number_of_TMDs_excl_SP.max() > 1:
        dataset = 'MP'
    else:
        dataset = 'Unknown'

    # join all keywords together into a large list
    nested_list_all_KW = list(itertools.chain(*list(df['uniprot_KW_for_analysis'].dropna())))
    # convert list to pandas series
    all_KW_series = pd.Series(nested_list_all_KW)
    # obtain series of major keywords
    KW_counts = all_KW_series.value_counts()
    # exclude keywords with less than x applicable proteins
    cutoff_major_keywords = s['cutoff_major_keywords']
    KW_counts_major = KW_counts[KW_counts > cutoff_major_keywords]
    # extract series indices and make them a python list
    list_KW_counts_major = sorted(list(KW_counts_major.index))
    logging.info ('valid keywords for analysis (n = {a}):\n{b}\n'.format(a = len(list_KW_counts_major), b = list_KW_counts_major))

    # define list of keywords that get excluded from analysis related to other keywords, if keyword is analysed, it gets included
    # move to settings file !?!?
    keywords_for_exclusion = ['Enzyme'] # Ion channel taken out on Feb 17 17 - number 7 not so important any more
    # symbols for annotation
    symbols = ['+', '#']
    # create bool in column for keyword to remove
    for element in keywords_for_exclusion:
        excl_list = ['{}'.format(element)]
        df[element] = df['uniprot_KW_for_analysis'].dropna().apply(KW_list_contains_any_desired_KW, args=(excl_list,))

    if list_KW_counts_major:
        # initialise pandas dataframe with keywords as index
        dfk = pd.DataFrame(index=list_KW_counts_major)
        # initialise pandas dataframe holding data from correlation analysis
        df_correlation = pd.DataFrame(index=list_KW_counts_major, columns=list_KW_counts_major)
        # initialise pandas dataframe that holds significant raw data for histogram re-creation
        dfr = pd.DataFrame(index=df.index)
        # initialise pretty dataframe
        dfp = pd.DataFrame(index=list_KW_counts_major + ['annotations'])
        # initialise figure number
        Fig_Nr = 0
        Fig_Nr_box = 0
        # add mean and std of whole dataset (AAIMON and AAIMON_slope)
        dfk['AAIMON_whole_dataset_mean'] = np.mean(df['AAIMON_mean_all_TM_res'])
        dfk['AAIMON_whole_dataset_std'] = np.std(df['AAIMON_mean_all_TM_res'])
        dfk['AAIMON_slope_whole_dataset_mean'] = np.mean(df['AAIMON_slope_all_TMDs_mean'])
        dfk['AAIMON_slope_whole_dataset_std'] = np.std(df['AAIMON_slope_all_TMDs_mean'])

        # exclude KW from analysis based on analysed KW, i.e. if Enzymes are excluded and the KW Enzyme is analysed, Enzymes get included again
        for keyword in list_KW_counts_major:
            # copy initial dataframe to drop enzymes and GPCRs dependent on keyword
            dfq = df.copy()
            acc_to_keep = []
            list_to_exclude = keywords_for_exclusion.copy()
            if keyword in list_to_exclude:
                list_to_exclude.remove(keyword)
            for element in list_to_exclude:
                dfq = dfq[dfq[element] == False]

            # create a new columns describing if the KW is in the KW list of that protein
            dfq['contains_KW'] = dfq['uniprot_KW_for_analysis'].apply(lambda x: keyword in x)
            # slice dataframe to view only entries with that keyword
            df_keyword = dfq.loc[dfq['contains_KW'] == True]
            df_no_keyword = dfq.loc[dfq['contains_KW'] == False]
            # check if keyword-dataframe still matches cutoff_major_keywords requirements
            if len(df_keyword) < cutoff_major_keywords:
                dfk = dfk.drop(keyword)
                dfp = dfp.drop(keyword)
                logging.info('\n\nkeyword "{}" does not match the requirements after excluding {}\n'.format(keyword, list_to_exclude))
                continue

            ###############################################################
            #                                                             #
            #                 calculate data for keyword                  #
            #                                                             #
            ###############################################################

            # calculate mean and std of AAIMON_mean_all_TM_res
            dfk.loc[keyword, 'AAIMON_keyword_mean'] = np.mean(df_keyword['AAIMON_mean_all_TM_res'])
            dfk.loc[keyword, 'AAIMON_keyword_std'] = np.std(df_keyword['AAIMON_mean_all_TM_res'])
            AAIMON_slope_keyword_mean = np.mean(df_keyword['AAIMON_slope_all_TMDs_mean'])
            dfk.loc[keyword, 'AAIMON_slope_keyword_mean'] = AAIMON_slope_keyword_mean
            dfk.loc[keyword, 'AAIMON_slope_keyword_std'] = np.std(df_keyword['AAIMON_slope_all_TMDs_mean'])
            number_of_proteins_keyword = len(df_keyword.index)
            dfk.loc[keyword, 'number_of_proteins_keyword'] = number_of_proteins_keyword
            dfk.loc[keyword, 'obs_changes_mean_keyword'] = np.mean(df_keyword['obs_changes_mean'])

            dfk.loc[keyword, 'AAIMON_no_keyword_mean'] = np.mean(df_no_keyword['AAIMON_mean_all_TM_res'])
            dfk.loc[keyword, 'AAIMON_no_keyword_std'] = np.std(df_no_keyword['AAIMON_mean_all_TM_res'])
            AAIMON_slope_no_keyword_mean = np.mean(df_no_keyword['AAIMON_slope_all_TMDs_mean'])
            dfk.loc[keyword, 'AAIMON_slope_no_keyword_mean'] = AAIMON_slope_no_keyword_mean
            dfk.loc[keyword, 'AAIMON_slope_no_keyword_std'] = np.std(df_no_keyword['AAIMON_slope_all_TMDs_mean'])
            number_of_proteins_no_keyword = len(df_no_keyword.index)
            dfk.loc[keyword, 'number_of_proteins_no_keyword'] = number_of_proteins_no_keyword
            dfk.loc[keyword, 'obs_changes_mean_no_keyword'] = np.mean(df_no_keyword['obs_changes_mean'])

            # calculate odds ratio, p- and t-values for AAIMONs
            dfk.loc[keyword, 'odds_ratio_AAIMON'] = dfk.loc[keyword, 'AAIMON_keyword_mean'] / dfk.loc[keyword, 'AAIMON_no_keyword_mean']
            data1 = df_keyword['AAIMON_mean_all_TM_res'].dropna()
            data2 = df_no_keyword['AAIMON_mean_all_TM_res'].dropna()
            t_AAIMON, p_AAIMON = ttest_ind(data1, data2, equal_var=True)             # equal_var True or False ?!?!
            dfk.loc[keyword, 't-value_AAIMON'] = t_AAIMON
            dfk.loc[keyword, 'p-value_AAIMON'] = p_AAIMON

            # calculate difference, p- and t-values for AAIMON_slopes
            difference_AAIMON_slope = dfk.loc[keyword, 'AAIMON_slope_keyword_mean'] - dfk.loc[keyword, 'AAIMON_slope_no_keyword_mean']
            dfk.loc[keyword, 'difference_AAIMON_slope'] = difference_AAIMON_slope
            KW_slope = df_keyword['AAIMON_slope_all_TMDs_mean'].dropna()
            no_KW_slope = df_no_keyword['AAIMON_slope_all_TMDs_mean'].dropna()
            t_AAIMON_slope, p_AAIMON_slope = ttest_ind(KW_slope, no_KW_slope, equal_var=True)             # equal_var True or False ?!?!
            dfk.loc[keyword, 't-value_AAIMON_slope'] = t_AAIMON_slope
            dfk.loc[keyword, 'p-value_AAIMON_slope'] = p_AAIMON_slope

            logging.info('\nmean AAIMON_slope containing keyword  "{a}" : {k:.3f} ± {l:.3f}, n = {d:.0f}\n'
                             'mean AAIMON_slope   without  keyword  "{a}" : {m:.3f} ± {n:.3f}, n = {g:.0f}\n'
                             'difference_AAIMON_slope = {o:.3f}, t-value_AAIMON_slope = {p:.3f}, p-value_AAIMON_slope = {q:.3f}'
                             .format(a=keyword,
                                     d=dfk.loc[keyword, 'number_of_proteins_keyword'],
                                     g=dfk.loc[keyword, 'number_of_proteins_no_keyword'],
                                     k=dfk.loc[keyword, 'AAIMON_slope_keyword_mean'],
                                     l=dfk.loc[keyword, 'AAIMON_slope_keyword_std'],
                                     m=dfk.loc[keyword, 'AAIMON_slope_no_keyword_mean'],
                                     n=dfk.loc[keyword, 'AAIMON_slope_no_keyword_std'],
                                     o=dfk.loc[keyword, 'difference_AAIMON_slope'],
                                     p=dfk.loc[keyword, 't-value_AAIMON_slope'],
                                     q=dfk.loc[keyword, 'p-value_AAIMON_slope']))


            ###############################################################
            #                                                             #
            #                  find correlated keywords                   #
            #                                                             #
            ###############################################################
            if len(df_keyword) != 0:
                for subKW in list_KW_counts_major:
                    # create a new column describing whether the protein KW list also contains the subKW
                    dfq['contains_subKW'] = df_keyword['uniprot_KW_for_analysis'].apply(lambda x: subKW in x)
                    # count how many of the proteins contain the subKW
                    val_counts = dfq['contains_subKW'].value_counts()
                    if True in val_counts.index:
                        num_prot_contain_subKW = val_counts[True]
                    else:
                        num_prot_contain_subKW = 0
                    # now add that number to the array of all KW against all KW
                    df_correlation.loc[subKW, keyword] = num_prot_contain_subKW
                # add correlated keywords to dfk
                correlated_keywords = df_correlation[keyword]
                # remove rows with 0
                correlated_keywords = correlated_keywords[(correlated_keywords != 0)]
                correlated_keywords_pretty = correlated_keywords.drop(keyword).sort_values(ascending=False).head(max_n_kw_in_df).index.tolist()
                correlated_keywords = correlated_keywords.drop(keyword).sort_values(ascending=False).head(10).to_string()
                dfk.loc[keyword, 'correlated_KW'] = correlated_keywords
            else:
                dfk.loc[keyword, 'correlated_KW'] = 'no correlated keywords found'
                correlated_keywords_pretty = 'no correlated keywords found'

            ###############################################################
            #                                                             #
            #        create normalised line histogram per keyword         #
            #                 save significant raw data                   #
            #                                                             #
            ###############################################################
            if dfk.loc[keyword, 'p-value_AAIMON_slope'] <= s['p_value_cutoff_for_histograms']:
                Fig_Nr += 1
                title = str(keyword)
                Fig_name = str(str(Fig_Nr) + '._' + 'Keyword_' + title)
                fig, ax = plt.subplots()

                ###############################################################
                #                                                             #
                #         plot histogram data containing keyword              #
                #                                                             #
                ###############################################################

                data_keyword = df_keyword['AAIMON_slope_all_TMDs_mean'].values
                # create numpy array of membranous over nonmembranous conservation ratios (identity)
                hist_data_keyword = np.array(data_keyword)
                # use numpy to create a histogram
                freq_counts, bin_array = np.histogram(hist_data_keyword, bins=binlist)
                freq_counts_normalised = freq_counts / freq_counts.max()
                # assuming all of the bins are exactly the same size, make the width of the column equal to XX% (e.g. 95%) of each bin
                col_width = float('%0.3f' % (0.95 * (bin_array[1] - bin_array[0])))
                # when align='center', the central point of the bar in the x-axis is simply the middle of the bins ((bin_0-bin_1)/2, etc)
                centre_of_bar_in_x_axis = (bin_array[:-2] + bin_array[1:-1]) / 2
                # add the final bin, which is physically located just after the last regular bin but represents all higher values
                bar_width = centre_of_bar_in_x_axis[3] - centre_of_bar_in_x_axis[2]
                centre_of_bar_in_x_axis = np.append(centre_of_bar_in_x_axis, centre_of_bar_in_x_axis[-1] + bar_width)
                linecontainer_AAIMON_mean = ax.plot(centre_of_bar_in_x_axis, freq_counts_normalised, color="#0489B1", alpha=alpha)

                ###############################################################
                #                                                             #
                #            plot histogram data without keyword              #
                #                                                             #
                ###############################################################

                data_no_keyword = df_no_keyword['AAIMON_slope_all_TMDs_mean'].values
                # create numpy array of membranous over nonmembranous conservation ratios (identity)
                hist_data_no_keyword = np.array(data_no_keyword)
                # use numpy to create a histogram
                freq_counts, bin_array = np.histogram(hist_data_no_keyword, bins=binlist)
                freq_counts_normalised = freq_counts / freq_counts.max()
                # assuming all of the bins are exactly the same size, make the width of the column equal to XX% (e.g. 95%) of each bin
                col_width = float('%0.3f' % (0.95 * (bin_array[1] - bin_array[0])))
                # when align='center', the central point of the bar in the x-axis is simply the middle of the bins ((bin_0-bin_1)/2, etc)
                centre_of_bar_in_x_axis = (bin_array[:-2] + bin_array[1:-1]) / 2
                # add the final bin, which is physically located just after the last regular bin but represents all higher values
                bar_width = centre_of_bar_in_x_axis[3] - centre_of_bar_in_x_axis[2]
                centre_of_bar_in_x_axis = np.append(centre_of_bar_in_x_axis, centre_of_bar_in_x_axis[-1] + bar_width)
                linecontainer_AAIMON_mean = ax.plot(centre_of_bar_in_x_axis, freq_counts_normalised, color="#EE762C", alpha=alpha)

                ###############################################################
                #                                                             #
                #                       set up plot style                     #
                #                                                             #
                ###############################################################

                ax.set_xlabel('AAIMON slope', fontsize=fontsize)
                # move the x-axis label closer to the x-axis
                #ax.xaxis.set_label_coords(0.45, -0.085)
                # x and y axes min and max
                xlim_min = -0.04
                xlim_max = 0.04
                ax.set_xlim(xlim_min, xlim_max)
                ylim_min = 0
                ylim_max = 1.2
                ax.set_ylim(ylim_min, ylim_max)
                # set x-axis ticks
                # use the slide selection to select every second item in the list as an xtick(axis label)
                #ax.set_xticks([float('%0.1f' % c) for c in centre_of_bar_in_x_axis[::3]])
                ax.set_ylabel('freq', rotation='vertical', fontsize=fontsize)
                # change axis font size
                ax.tick_params(labelsize=fontsize)
                # create legend
                legend_obj = ax.legend(['containing KW n={}'.format(number_of_proteins_keyword), 'without KW n={}'.format(number_of_proteins_no_keyword)], loc='upper right', fontsize=fontsize, frameon=True)
                # add figure number to top left of subplot
                ax.annotate(s=str(Fig_Nr) + '.', xy=(0.04, 0.9), fontsize=fontsize, xytext=None, xycoords='axes fraction', alpha=0.75)
                # add figure title to top left of subplot
                ax.annotate(s=title, xy=(0.1, 0.9), fontsize=fontsize, xytext=None, xycoords='axes fraction', alpha=0.75)
                # add p-value to top left of subplot
                ax.annotate(s='p-value = {:.5f}'.format(p_AAIMON_slope), xy=(0.1, 0.85), fontsize=fontsize, xytext=None, xycoords='axes fraction', alpha=0.75)
                # save every individual figure
                utils.save_figure(fig, Fig_name, base_filepath, save_png, save_pdf)

                if dataset == 'MP':
                    ### boxplot of all TMDs

                    ### this section specifies the last bin to avoid bins containing only one TMD
                    # join all numbers of TMDs together into a large list
                    nested_list_all_TMDs = list(df_keyword['number_of_TMDs_excl_SP'])
                    # convert list to pandas series
                    all_TMDs_series = pd.Series(nested_list_all_TMDs)
                    # obtain series of TMD_counts
                    TMD_counts = all_TMDs_series.value_counts()
                    # exclude TMD numbers with less than x applicable proteins from boxplot max detection
                    boxplot_cutoff_number_of_TMDs = 5
                    TMD_counts_major = TMD_counts[TMD_counts >= boxplot_cutoff_number_of_TMDs]
                    max_num_TMDs = TMD_counts_major.index.max()


                    if pd.notnull(max_num_TMDs):
                        title = str(keyword) + '_AAIMON_Boxplot'
                        Fig_Nr_box += 1
                        Fig_name = str(str(Fig_Nr_box) + '._' + 'Keyword_' + title)
                        fig, ax = plt.subplots()
                        ax2 = plt.twinx()

                        legend = []
                        data_to_plot = []
                        for i in range(1, max_num_TMDs.astype('int') + 1):
                            TM = 'TM%02d' % i
                            hist_data_AAIMON_each_TM = df_keyword['TM%02d_AAIMON_slope' % i].dropna() * 1000
                            if len(hist_data_AAIMON_each_TM) > 0:
                                data_to_plot.append(hist_data_AAIMON_each_TM)
                                legend.append(TM)

                        # add values of every TMD number that is larger than the boxplot_cutoff_number_of_TMDs to final bin
                        data_for_final_bin = []
                        for i in range(max_num_TMDs.astype('int') + 1, df_keyword.number_of_TMDs_excl_SP.max().astype('int') + 1):
                            # TM_final = 'TM%02d' % i
                            hist_data_AAIMON_each_TM_final_bin = df_keyword['TM%02d_AAIMON_slope' % i].dropna() * 1000
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
                        ax2.set_ylabel('number of TMDs in bin', rotation='vertical', fontsize=fontsize, color='#0076B8')
                        ax2.tick_params(labelsize=fontsize, labelcolor='#0076B8')

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

                        ax.set_ylabel(r'm$_{\rm TM/nonTM} *10^{\rm -3}$', rotation='vertical', fontsize=fontsize)
                        ax.set_ylim(-20, 30)

                        ## Remove top axes and right axes ticks
                        ax.get_xaxis().tick_bottom()
                        ax.get_yaxis().tick_left()
                        ## Custom x-axis labels
                        ax.set_xticklabels(legend, rotation=45)
                        # add figure number to top left of subplot
                        ax.annotate(s=str(Fig_Nr_box) + '. ' + title + ' ; p-value = {:.5f}'.format(p_AAIMON_slope), xy=(0.02, 0.95), fontsize=fontsize, xytext=None, xycoords='axes fraction', alpha=0.75)
                        plt.tight_layout()

                        utils.save_figure(fig, Fig_name, base_filepath, save_png, save_pdf)

                        ##############################
                        #       plot lipo data       #
                        ##############################

                        title = str(keyword) + '_lipo_Boxplot'
                        Fig_name = str(str(Fig_Nr_box) + '._' + 'Keyword_' + title)
                        fig, ax = plt.subplots()
                        ax2 = plt.twinx()

                        legend = []
                        data_to_plot = []
                        for i in range(1, max_num_TMDs.astype('int') + 1):
                            TM = 'TM%02d' % i
                            hist_data_AAIMON_each_TM = df_keyword['TM%02d_lipo' % i].dropna()
                            if len(hist_data_AAIMON_each_TM) > 0:
                                data_to_plot.append(hist_data_AAIMON_each_TM)
                                legend.append(TM)

                        # add values of every TMD number that is larger than the boxplot_cutoff_number_of_TMDs to final bin
                        data_for_final_bin = []
                        for i in range(max_num_TMDs.astype('int') + 1, df_keyword.number_of_TMDs_excl_SP.max().astype('int') + 1):
                            # TM_final = 'TM%02d' % i
                            hist_data_AAIMON_each_TM_final_bin = df_keyword['TM%02d_lipo' % i].dropna()
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
                        ax2.set_ylabel('number of TMDs in bin', rotation='vertical', fontsize=fontsize, color='#0076B8')
                        ax2.tick_params(labelsize=fontsize, labelcolor='#0076B8')

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
                        ax.set_ylim(-0.4, 1)

                        ## Remove top axes and right axes ticks
                        ax.get_xaxis().tick_bottom()
                        ax.get_yaxis().tick_left()
                        ## Custom x-axis labels
                        ax.set_xticklabels(legend, rotation=45)
                        # add figure number to top left of subplot
                        ax.annotate(s=str(Fig_Nr_box) + '. ' + title + ' ; p-value = {:.5f}'.format(p_AAIMON_slope), xy=(0.02, 0.95), fontsize=fontsize, xytext=None, xycoords='axes fraction', alpha=0.75)
                        plt.tight_layout()

                        utils.save_figure(fig, Fig_name, base_filepath, save_png, save_pdf)

            # add selected stuff to pretty dataframe dfp if significant, if not, drop keyword from dfp, create annotations
            if p_AAIMON_slope <= s['p_value_cutoff_for_histograms']:
                dfr['RAW_KW_{}'.format(keyword)] = KW_slope
                dfr['RAW_no_KW_{}'.format(keyword)] = no_KW_slope

                if not keyword in keywords_for_exclusion:
                    replace = keyword
                    keyword = '{}*'.format(keyword)
                    dfp = dfp.rename(index={replace : keyword})
                # else:
                #     if not keyword in keywords_for_exclusion[0]:
                #         replace = keyword
                #         keyword = '{}{}'.format(keyword, symbols[0])
                #         dfp = dfp.rename(index={replace: keyword})
                    # elif not keyword in keywords_for_exclusion[1]:
                    #     replace = keyword
                    #     keyword = '{}{}'.format(keyword, symbols[1])
                    #     dfp = dfp.rename(index={replace: keyword})
                #dfp.loc[keyword, 'Dataset'] = dataset
                dfp.loc[keyword, 'p-value'] = p_AAIMON_slope
                dfp.loc[keyword, 'mean AAIMON slope *10^-3'] = AAIMON_slope_keyword_mean * 1000
                dfp.loc[keyword, 'd AAIMON slope *10^-3'] = difference_AAIMON_slope * 1000
                dfp.loc[keyword, '# proteins'] = '{}/{}'.format(number_of_proteins_keyword, number_of_proteins_keyword+number_of_proteins_no_keyword)
                dfp.loc[keyword, 'Top correlated keywords'] = ', '.join(correlated_keywords_pretty)
            else:
                dfp = dfp.drop(keyword)

        # sort dataframes by p-value of AAIMON slopes
        dfk = dfk.sort_values('p-value_AAIMON_slope')
        if not dfp.empty:
            dfp.sort_values('p-value', inplace=True)

        # for betabarrel dataset, the keywords are PFAM acc. Find relevant accession names
        if is_betabarrel_dataset:
            # open pfam mapping csv. See PFAM "FTP" site. Downlead and save at database_folder\PFAM\pdb_pfam_mapping.txt
            pfam_mapping_csv = os.path.join(s["data_dir"], "PFAM\pdb_pfam_mapping.txt")
            df_pfam = pd.read_csv(pfam_mapping_csv, sep="\t")
            # delete the .15, .17, etc at end of PFAM accessions
            df_pfam["pfam_acc_base"] = df_pfam.PFAM_ACC.apply(lambda x: x.split(".")[0])
            # set as index
            df_pfam.set_index("pfam_acc_base", inplace=True)
            for acc in dfk.index:
                if acc in df_pfam.index:
                    # get either a series of redundant records, or a single PFAM description
                    PFAM_desc_ser_or_str = df_pfam.loc[acc, "PFAM_desc"]
                    # if it's a series, grab the first one
                    if isinstance(PFAM_desc_ser_or_str, pd.Series):
                        PFAM_desc = PFAM_desc_ser_or_str.iloc[0]
                    else:
                        # assume it's a string, grab the single description
                        PFAM_desc = PFAM_desc_ser_or_str
                else:
                    if acc in PFAM_dict:
                        PFAM_desc = PFAM_dict[acc]
                    else:
                        PFAM_desc = "PFAM not found"
                # add as a new column to the dataframe
                dfk.loc[acc, "PFAM_desc"] = PFAM_desc
                # for dfp, assume the star is already in the index
                dfp.loc[acc + "*", "PFAM_desc"] = PFAM_desc

        annotate = ['*excluding {}'.format(', '.join(list_to_exclude))]
        # for n, element in enumerate(list_to_exclude):
        #     annotate.append('{}excluding {}'.format(symbols[n], element))

        dfp.loc['annotations', 'p-value'] = '; '.join(annotate)

        # save pandas dataframes with values
        dfk.to_csv(os.path.join(pathdict["keywords"], 'List%02d_keywords.csv' % s["list_number"]), sep=",", quoting=csv.QUOTE_NONNUMERIC)
        df_correlation.to_csv(os.path.join(pathdict["keywords"], 'List%02d_KW_cross_correlation.csv' % s["list_number"]), sep=",", quoting=csv.QUOTE_NONNUMERIC)
        dfr.to_csv(os.path.join(pathdict["keywords"], 'List%02d_keywords_significant_RAW_data.csv' % s["list_number"]), sep=",", quoting=csv.QUOTE_NONNUMERIC)
        dfp.to_csv(os.path.join(pathdict["keywords"], 'List%02d_keywords_pretty.csv' % s["list_number"]), sep=",", quoting=csv.QUOTE_NONNUMERIC)
    else:
        return 'no valid keywords found! change "cutoff_major_keywords" setting! \ncurrent value: {}'.format(s['cutoff_major_keywords'])

    return "\n~~~~~~~~~~~~                      finished keyword_analysis                         ~~~~~~~~~~~~"


def get_list_enzyme_KW_and_list_ignored_KW():
    ''' defines keywords that are Enzyme associated and keywords that are ignored

        Parameters
    ----------
    none

        Return
    ----------
    list_enzyme_KW: list
        list of enzyme associated keywords
    list_ignored_KW: list
        list of keywords that can be ignored

    '''
    list_enzyme_KW = ['Transferase', 'Hydrolase', 'Glycosyltransferase', 'Protease', 'Kinase', 'Oxidoreductase', 'Metalloprotease', 'Serine protease',
                      'Protein phosphatase', 'Ligase', 'Acyltransferase', 'Serine/threonine-protein kinase', 'Glycosidase', 'Aminopeptidase',
                      'Isomerase', 'Methyltransferase', 'Carboxypeptidase', 'Hydroxylation', 'Aspartyl protease', 'Serine esterase',
                      'Lipid biosynthesis', 'GPI-anchor biosynthesis', 'Steroid biosynthesis', 'Melanin biosynthesis', 'Thyroid hormones biosynthesis',
                      'Phospholipid biosynthesis', 'Sterol biosynthesis', 'Glutathione biosynthesis', 'Cholesterol biosynthesis',
                      'Fatty acid biosynthesis', 'Prostaglandin biosynthesis', 'cGMP biosynthesis', 'Leukotriene biosynthesis', 'Catecholamine biosynthesis',
                      'Lipid metabolism', 'Carbohydrate metabolism', 'Steroid metabolism', 'Sterol metabolism', 'Sphingolipid metabolism',
                      'Cholesterol metabolism', 'Fatty acid metabolism', 'Phospholipid metabolism', 'Catecholamine metabolism', 'Prostaglandin metabolism',
                      'Glycogen metabolism', 'Fucose metabolism']

    list_ignored_KW = ['Transmembrane', 'Complete proteome', 'Reference proteome', 'Membrane',
                       'Transmembrane helix', 'Repeat', 'Alternative splicing', 'Sodium', 'Potassium', 'Direct protein sequencing',
                       'Transducer', 'Polymorphism', 'Glycoprotein', 'Calcium transport', 'Ion transport', 'Transport', 'Protein transport',
                       'Voltage-gated channel', 'ATP-binding', 'Calcium', 'Zinc', 'Synapse', 'Signal', 'Disulfide bond', '3D-structure',
                       'Host-virus interaction', 'Palmitate', 'Potassium transport', 'Endoplasmic reticulum', 'G-protein coupled receptor',
                       'Plasma membrane', 'Golgi apparatus','Sensory transduction', 'Calcium channel', 'Sugar transport', 'Metal-binding']

    list_ignored_PFAM_acc = ["PF13505", "PF13568", "DUF3078", "PF04338"] # Outer membrane protein beta-barrel domain,  Outer membrane protein beta-barrel domain : This domain is found in a wide range of outer membrane proteins. This domain assumes a membrane bound beta-barrel fold.

    list_ignored_KW += list_ignored_PFAM_acc

    PFAM_dict = {"PF06629": "MltA-interacting protein MipA"}


    return list_enzyme_KW, list_ignored_KW, PFAM_dict


def KW_list_contains_any_desired_KW(KW_list_to_search,list_desired_KW):
    ''' Determine if two lists contain any common values.
    Used to determine, for example, if a list of keywords contains any
    enzyme related words, from another list
    input:
    KW_list_to_search
    list_desired_KW
    note: in theory, this function could be updated to use set(), which should be slightly quicker
    '''
    is_enzyme = False
    for KW in list_desired_KW:
        if KW in KW_list_to_search:
            is_enzyme = True
            break
    return is_enzyme