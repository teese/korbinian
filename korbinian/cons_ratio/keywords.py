import pandas as pd
import numpy as np
import ast
import csv
import sys
import itertools
from scipy.stats import ttest_ind
import korbinian.utils as utils
import matplotlib.pyplot as plt

def keyword_analysis(pathdict, s, logging):
    logging.info("~~~~~~~~~~~~         starting keyword_analysis           ~~~~~~~~~~~~")

    # base filepath for figures
    base_filepath = pathdict["figs_keywords"]
    #save figures to .pdf or .png
    save_png = s["save_png"]
    save_pdf = s["save_pdf"]
    # create binlist for histograms
    linspace_binlist = np.linspace(s["mp_smallest_bin"], s["mp_largest_bin"], s["mp_number_of_bins"])
    # add 30 as the last bin, to make sure 100% of the data is added to the histogram, including major outliers
    binlist = np.append(linspace_binlist, s["mp_final_highest_bin"])
    # initialise basic settings for figures
    plt.style.use('bmh')
    fontsize = 12
    alpha = 0.8

    # load cr_summary file
    dfc = pd.read_csv(pathdict["list_cr_summary_csv"], sep=",", quoting=csv.QUOTE_NONNUMERIC, index_col=0)
    # load summary file
    dfu = pd.read_csv(pathdict["list_summary_csv"], sep=",", quoting=csv.QUOTE_NONNUMERIC, index_col=0)
    # merge cr_summary and summary file, if columns are equal in both files, suffix _dfc will be added in cr_summary column names for backwards compatibility
    df = pd.merge(dfc, dfu, left_index=True, right_index=True, suffixes=('_dfc', ''))

    # remove proteins containing nan in AAIMON_ratio_mean_all_TMDs
    list_of_acc_without_nan = []
    for acc in df.index:
        if not pd.isnull(df.loc[acc, 'AAIMON_ratio_mean_all_TMDs']):
            list_of_acc_without_nan.append(acc)
    df = df.loc[list_of_acc_without_nan, :]
    # apply ast.literal_eval to every item in df['uniprot_KW']
    if isinstance(df['uniprot_KW'][0], str):
        df['uniprot_KW'] = df['uniprot_KW'].apply(lambda x: ast.literal_eval(x))
    # join all keywords together into a large list
    nested_list_all_KW = list(itertools.chain(*list(df['uniprot_KW'])))
    # convert list to pandas series
    all_KW_series = pd.Series(nested_list_all_KW)
    # obtain series of major keywords
    KW_counts = all_KW_series.value_counts()
    # exclude keywords with less than x applicable proteins
    cutoff_major_keywords = s['cutoff_major_keywords']
    KW_counts_major = KW_counts[KW_counts > cutoff_major_keywords]
    # create a list of keywords to be ignored
    list_ignored_KW = ['Transmembrane', 'Complete proteome', 'Reference proteome', 'Membrane',
                       'Transmembrane helix', 'Cell membrane', 'Repeat',
                       'Alternative splicing', 'Sodium', 'Potassium', 'Direct protein sequencing', 'Transducer',
                       'Polymorphism', 'Glycoprotein']
    # check if protein is an enzyme
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
    # check if protein is a GPCR
    list_GPCR_KW = ['G-protein coupled receptor']
    df['GPCR'] = df['uniprot_KW'].apply(utils.KW_list_contains_any_desired_KW, args=(list_GPCR_KW,))

    for KW in list_ignored_KW:
        if KW in KW_counts_major.index:
            KW_counts_major = KW_counts_major.drop(KW)
    # extract series indices and make them a python list
    list_KW_counts_major = sorted(list(KW_counts_major.index))
    sys.stdout.write ('valid keywords for analysis (n = {a}):\n{b}\n\n'.format(a = len(list_KW_counts_major), b = list_KW_counts_major))
    # check if keywords are present
    if list_KW_counts_major:
        # initialise pandas dataframe with keywords as index
        dfk = pd.DataFrame(index=list_KW_counts_major)
        # initialise figure number
        Fig_Nr = 0
        for keyword in list_KW_counts_major:
            # check if keyword is GPCR or enzyme - important for exclusion analysis
            if keyword in list_enzyme_KW:
                KW_is_enzyme = True
            else:
                KW_is_enzyme = False

            if keyword == 'G-protein coupled receptor':
                KW_is_GPCR = True
            else:
                KW_is_GPCR = False

            # initialise lists of acc that do or do not contain the keyword
            list_of_acc_containing_kw = []
            list_of_acc_without_kw = []
            for acc in df.index:
                if keyword in df.loc[acc, 'uniprot_KW']:  # here, no ast.literal_eval() is necessary as it was applied as lambda function above
                    list_of_acc_containing_kw.append(acc)
                else:
                    list_of_acc_without_kw.append(acc)
            # dropping all proteins that do or do not contain the keyword - create two dataframes with and without keyword
            df_keyword = df.loc[list_of_acc_containing_kw, :]
            df_no_keyword = df.loc[list_of_acc_without_kw, :]

            # removing enzymes from dataframes
            df_keyword_nonenzyme = df_keyword.query('enzyme == False')
            df_no_keyword_nonenzyme = df_no_keyword.query('enzyme == False')
            # removing GPCRs from dataframes
            df_keyword_nonGPCR = df_keyword.query('GPCR == False')
            df_no_keyword_nonGPCR = df_no_keyword.query('GPCR == False')
            # removing Enzymes and GPCRs from dataframes
            df_keyword_nonenzyme_nonGPCR = df_keyword_nonenzyme.query('GPCR == False')
            df_no_keyword_nonenzyme_nonGPCR = df_no_keyword_nonenzyme.query('GPCR == False')

            ###############################################################
            #                                                             #
            #                     whole dataset data                      #
            #                                                             #
            ###############################################################

            # calculate mean and std of AAIMON_ratio_mean_all_TMDs
            dfk.loc[keyword, 'AAIMON_keyword_mean'] = np.mean(df_keyword['AAIMON_ratio_mean_all_TMDs'])
            dfk.loc[keyword, 'AAIMON_keyword_std'] = np.std(df_keyword['AAIMON_ratio_mean_all_TMDs'])
            number_of_proteins_keyword = len(df_keyword.index)
            dfk.loc[keyword, 'number_of_proteins_keyword'] = number_of_proteins_keyword
            dfk.loc[keyword, 'AAIMON_no_keyword_mean'] = np.mean(df_no_keyword['AAIMON_ratio_mean_all_TMDs'])
            dfk.loc[keyword, 'AAIMON_no_keyword_std'] = np.std(df_no_keyword['AAIMON_ratio_mean_all_TMDs'])
            number_of_proteins_no_keyword = len(df_no_keyword.index)
            dfk.loc[keyword, 'number_of_proteins_no_keyword'] = number_of_proteins_no_keyword
            # calculate odds ratio
            dfk.loc[keyword, 'odds_ratio_whole_dataset'] = dfk.loc[keyword, 'AAIMON_keyword_mean'] / dfk.loc[keyword, 'AAIMON_no_keyword_mean']
            # ttest p- and t- values calculation
            data1 = df_keyword['AAIMON_ratio_mean_all_TMDs']
            data2 = df_no_keyword['AAIMON_ratio_mean_all_TMDs']
            t, p = ttest_ind(data1, data2, equal_var=False)             # equal_var True or False ?!?!
            dfk.loc[keyword, 't-value_whole_dataset'] = t
            dfk.loc[keyword, 'p-value_whole_dataset'] = p

            if not KW_is_enzyme:
                ###############################################################
                #                                                             #
                #                       without enzymes                       #
                #                                                             #
                ###############################################################

                # calculate mean and std of AAIMON_ratio_mean_all_TMDs
                dfk.loc[keyword, 'AAIMON_keyword_nonenzyme_mean'] = np.mean(df_keyword_nonenzyme['AAIMON_ratio_mean_all_TMDs'])
                dfk.loc[keyword, 'AAIMON_keyword_nonenzyme_std'] = np.std(df_keyword_nonenzyme['AAIMON_ratio_mean_all_TMDs'])
                number_of_proteins_keyword_nonenzyme = len(df_keyword_nonenzyme.index)
                dfk.loc[keyword, 'number_of_proteins_keyword_nonenzyme'] = number_of_proteins_keyword_nonenzyme

                dfk.loc[keyword, 'AAIMON_no_keyword_nonenzyme_mean'] = np.mean(df_no_keyword_nonenzyme['AAIMON_ratio_mean_all_TMDs'])
                dfk.loc[keyword, 'AAIMON_no_keyword_nonenzyme_std'] = np.std(df_no_keyword_nonenzyme['AAIMON_ratio_mean_all_TMDs'])
                number_of_proteins_no_keyword_nonenzyme = len(df_no_keyword_nonenzyme.index)
                dfk.loc[keyword, 'number_of_proteins_no_keyword_nonenzyme'] = number_of_proteins_no_keyword_nonenzyme

                # calculate odds ratio
                dfk.loc[keyword, 'odds_ratio_excluding_enzymes'] = dfk.loc[keyword, 'AAIMON_keyword_nonenzyme_mean'] / dfk.loc[keyword, 'AAIMON_no_keyword_nonenzyme_mean']
                # ttest p- and t- values calculation
                data1 = df_keyword_nonenzyme['AAIMON_ratio_mean_all_TMDs']
                data2 = df_no_keyword_nonenzyme['AAIMON_ratio_mean_all_TMDs']
                t, p = ttest_ind(data1, data2, equal_var=False)
                dfk.loc[keyword, 't-value_excluding_enzymes'] = t
                dfk.loc[keyword, 'p-value_excluding_enzymes'] = p

            if not KW_is_GPCR:
                ###############################################################
                #                                                             #
                #                        without GPCRs                        #
                #                                                             #
                ###############################################################

                # calculate mean and std of AAIMON_ratio_mean_all_TMDs
                dfk.loc[keyword, 'AAIMON_keyword_nonGPCR_mean'] = np.mean(df_keyword_nonGPCR['AAIMON_ratio_mean_all_TMDs'])
                dfk.loc[keyword, 'AAIMON_keyword_nonGPCR_std'] = np.std(df_keyword_nonGPCR['AAIMON_ratio_mean_all_TMDs'])
                number_of_proteins_keyword_nonGPCR = len(df_keyword_nonGPCR.index)
                dfk.loc[keyword, 'number_of_proteins_keyword_nonGPCR'] = number_of_proteins_keyword_nonGPCR

                dfk.loc[keyword, 'AAIMON_no_keyword_nonGPCR_mean'] = np.mean(df_no_keyword_nonGPCR['AAIMON_ratio_mean_all_TMDs'])
                dfk.loc[keyword, 'AAIMON_no_keyword_nonGPCR_std'] = np.std(df_no_keyword_nonGPCR['AAIMON_ratio_mean_all_TMDs'])
                number_of_proteins_no_keyword_nonGPCR = len(df_no_keyword_nonGPCR.index)
                dfk.loc[keyword, 'number_of_proteins_no_keyword_nonGPCR'] = number_of_proteins_no_keyword_nonGPCR

                # calculate odds ratio
                dfk.loc[keyword, 'odds_ratio_excluding_GPCRs'] = dfk.loc[keyword, 'AAIMON_keyword_nonGPCR_mean'] / dfk.loc[keyword, 'AAIMON_no_keyword_nonGPCR_mean']
                # ttest p- and t- values calculation
                data1 = df_keyword_nonGPCR['AAIMON_ratio_mean_all_TMDs']
                data2 = df_no_keyword_nonGPCR['AAIMON_ratio_mean_all_TMDs']
                t, p = ttest_ind(data1, data2, equal_var=False)
                dfk.loc[keyword, 't-value_excluding_GPCRs'] = t
                dfk.loc[keyword, 'p-value_excluding_GPCRs'] = p

            if not KW_is_enzyme and not KW_is_GPCR:
                ###############################################################
                #                                                             #
                #                  without enzymes and GPCRs                  #
                #                                                             #
                ###############################################################

                # calculate mean and std of AAIMON_ratio_mean_all_TMDs
                dfk.loc[keyword, 'AAIMON_keyword_nonenzyme_nonGPCR_mean'] = np.mean(df_keyword_nonenzyme_nonGPCR['AAIMON_ratio_mean_all_TMDs'])
                dfk.loc[keyword, 'AAIMON_keyword_nonenzyme_nonGPCR_std'] = np.std(df_keyword_nonenzyme_nonGPCR['AAIMON_ratio_mean_all_TMDs'])
                number_of_proteins_keyword_nonenzyme_nonGPCR = len(df_keyword_nonenzyme_nonGPCR.index)
                dfk.loc[keyword, 'number_of_proteins_keyword_nonenzyme_nonGPCR'] = number_of_proteins_keyword_nonenzyme_nonGPCR

                dfk.loc[keyword, 'AAIMON_no_keyword_nonenzyme_nonGPCR_mean'] = np.mean(df_no_keyword_nonenzyme_nonGPCR['AAIMON_ratio_mean_all_TMDs'])
                dfk.loc[keyword, 'AAIMON_no_keyword_nonenzyme_nonGPCR_std'] = np.std(df_no_keyword_nonenzyme_nonGPCR['AAIMON_ratio_mean_all_TMDs'])
                number_of_proteins_no_keyword_nonenzyme_nonGPCR = len(df_no_keyword_nonenzyme_nonGPCR.index)
                dfk.loc[keyword, 'number_of_proteins_no_keyword_nonenzyme_nonGPCR'] = number_of_proteins_no_keyword_nonenzyme_nonGPCR

                # calculate odds ratio
                dfk.loc[keyword, 'odds_ratio_excluding_Enzymes+GPCRs'] = dfk.loc[keyword, 'AAIMON_keyword_nonenzyme_nonGPCR_mean'] / dfk.loc[keyword, 'AAIMON_no_keyword_nonenzyme_nonGPCR_mean']
                # ttest p- and t- values calculation
                data1 = df_keyword_nonenzyme_nonGPCR['AAIMON_ratio_mean_all_TMDs']
                data2 = df_no_keyword_nonenzyme_nonGPCR['AAIMON_ratio_mean_all_TMDs']
                t, p = ttest_ind(data1, data2, equal_var=False)
                dfk.loc[keyword, 't-value_excluding_Enzymes+GPCRs'] = t
                dfk.loc[keyword, 'p-value_excluding_Enzymes+GPCRs'] = p



            sys.stdout.write('mean AAIMON containing keyword  "{a}" : {b:.3f} ± {c:.3f}, n = {d:.0f}\n'
                             'mean AAIMON   without  keyword  "{a}" : {e:.3f} ± {f:.3f}, n = {g:.0f}\n'
                             'odds_ratio = {h:.3f}, t-value = {i:.3f}, p-value = {j:.3f}\n'
                             .format(a=keyword, b=dfk.loc[keyword, 'AAIMON_keyword_mean'], c=dfk.loc[keyword, 'AAIMON_keyword_std'], d=dfk.loc[keyword, 'number_of_proteins_keyword'],
                                     e=dfk.loc[keyword, 'AAIMON_no_keyword_mean'], f=dfk.loc[keyword, 'AAIMON_no_keyword_std'], g=dfk.loc[keyword, 'number_of_proteins_no_keyword'],
                                     h=dfk.loc[keyword, 'odds_ratio_whole_dataset'], i=dfk.loc[keyword, 't-value_whole_dataset'], j=dfk.loc[keyword, 'p-value_whole_dataset']))

            ###############################################################
            #                                                             #
            #        create normalised line histogram per keyword         #
            #                                                             #
            ###############################################################

            Fig_Nr += 1
            title = str(keyword)
            Fig_name = str(str(Fig_Nr) + '._' + 'Keyword_' + title)

            fig, ax = plt.subplots()

            ###############################################################
            #                                                             #
            #         plot histogram data containing keyword              #
            #                                                             #
            ###############################################################

            data_keyword = df_keyword['AAIMON_ratio_mean_all_TMDs'].values
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

            data_no_keyword = df_no_keyword['AAIMON_ratio_mean_all_TMDs'].values
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

            ax.set_xlabel('average conservation ratio (membranous over nonmembranous)', fontsize=fontsize)
            # move the x-axis label closer to the x-axis
            ax.xaxis.set_label_coords(0.45, -0.085)
            # x and y axes min and max
            xlim_min = 0.5
            xlim_max = 1.75
            ax.set_xlim(xlim_min, xlim_max)
            ylim_min = 0
            ylim_max = 1.2
            ax.set_ylim(ylim_min, ylim_max)
            # set x-axis ticks
            # use the slide selection to select every second item in the list as an xtick(axis label)
            ax.set_xticks([float('%0.1f' % c) for c in centre_of_bar_in_x_axis[::3]])
            ax.set_ylabel('freq', rotation='vertical', fontsize=fontsize)
            # change axis font size
            ax.tick_params(labelsize=fontsize)
            # create legend
            legend_obj = ax.legend(['containing KW n={}'.format(number_of_proteins_keyword), 'without KW n={}'.format(number_of_proteins_no_keyword)], loc='upper right', fontsize=fontsize)
            # add figure number to top left of subplot
            ax.annotate(s=str(Fig_Nr) + '.', xy=(0.04, 0.9), fontsize=fontsize, xytext=None, xycoords='axes fraction', alpha=0.75)
            # add figure title to top left of subplot
            ax.annotate(s=title, xy=(0.1, 0.9), fontsize=fontsize, xytext=None, xycoords='axes fraction', alpha=0.75)
            # save every individual figure
            utils.save_figure(fig, Fig_name, base_filepath, save_png, save_pdf)

        # # add mean and std of whole dataset
        # dfk['AAIMON_whole_dataset_mean'] = np.mean(df['AAIMON_ratio_mean_all_TMDs'])
        # dfk['AAIMON_whole_dataset_std'] = np.std(df['AAIMON_ratio_mean_all_TMDs'])
        # save pandas dataframe with values
        dfk.to_csv(pathdict["list_keywords_csv"], sep=",", quoting=csv.QUOTE_NONNUMERIC)   # - transpose dataframe here ?!
    else:
        sys.stdout.write ('no valid keywords found! change "cutoff_major_keywords" setting! current value: {}'.format(s['cutoff_major_keywords']))

    logging.info("\n~~~~~~~~~~~~        keyword_analysis is finished         ~~~~~~~~~~~~")