import pandas as pd
import numpy as np
import ast
import csv
import sys
import itertools
from scipy.stats import ttest_ind
import matplotlib.pyplot as plt

def keyword_analysis(pathdict, s, logging):
    logging.info("~~~~~~~~~~~~         starting keyword_analysis           ~~~~~~~~~~~~")

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
    for KW in list_ignored_KW:
        if KW in KW_counts_major.index:
            KW_counts_major = KW_counts_major.drop(KW)
    # extract series indices and make them a python list
    list_KW_counts_major = list(KW_counts_major.index)
    sys.stdout.write ('valid keywords for analysis (n = {a}):\n{b}\n\n'.format(a = len(list_KW_counts_major), b = list_KW_counts_major))
    # check if keywords are present
    if list_KW_counts_major:
        # initialise pandas dataframe with keywords as index
        dfk = pd.DataFrame(index=list_KW_counts_major)
        for keyword in list_KW_counts_major:
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
            # calculate mean and std of AAIMON_ratio_mean_all_TMDs
            dfk.loc[keyword, 'AAIMON_keyword_mean'] = np.mean(df_keyword['AAIMON_ratio_mean_all_TMDs'])
            dfk.loc[keyword, 'AAIMON_keyword_std'] = np.std(df_keyword['AAIMON_ratio_mean_all_TMDs'])
            dfk.loc[keyword, 'AAIMON_no_keyword_mean'] = np.mean(df_no_keyword['AAIMON_ratio_mean_all_TMDs'])
            dfk.loc[keyword, 'AAIMON_no_keyword_std'] = np.std(df_no_keyword['AAIMON_ratio_mean_all_TMDs'])
            dfk.loc[keyword, 'number_of_proteins_keyword'] = len(df_keyword.index)
            dfk.loc[keyword, 'number_of_proteins_no_keyword'] = len(df_no_keyword.index)
            # calculate odds ratio
            dfk.loc[keyword, 'odds_ratio'] = dfk.loc[keyword, 'AAIMON_keyword_mean'] / dfk.loc[keyword, 'AAIMON_no_keyword_mean']
            # ttest p- and t- values calculation
            data1 = df_keyword['AAIMON_ratio_mean_all_TMDs']
            data2 = df_no_keyword['AAIMON_ratio_mean_all_TMDs']
            t, p = ttest_ind(data1, data2, equal_var=False)
            dfk.loc[keyword, 't-value'] = t
            dfk.loc[keyword, 'p-value'] = p
            sys.stdout.write('mean AAIMON containing keyword  "{a}" : {b:.3f} ± {c:.3f}, number of proteins: {d:.0f}\n'
                             'mean AAIMON   without  keyword  "{a}" : {e:.3f} ± {f:.3f}, number of proteins: {g:.0f}\n'
                             'odds_ratio = {h:.3f}, t-value = {i:.3f}, p-value = {j:.3f}\n'
                             .format(a=keyword, b=dfk.loc[keyword, 'AAIMON_keyword_mean'], c=dfk.loc[keyword, 'AAIMON_keyword_std'], d=dfk.loc[keyword, 'number_of_proteins_keyword'],
                                     e=dfk.loc[keyword, 'AAIMON_no_keyword_mean'], f=dfk.loc[keyword, 'AAIMON_no_keyword_std'], g=dfk.loc[keyword, 'number_of_proteins_no_keyword'],
                                     h=dfk.loc[keyword, 'odds_ratio'], i=dfk.loc[keyword, 't-value'], j=dfk.loc[keyword, 'p-value']))
        # add mean and std of whole dataset
        dfk['AAIMON_whole_dataset_mean'] = np.mean(df['AAIMON_ratio_mean_all_TMDs'])
        dfk['AAIMON_whole_dataset_std'] = np.std(df['AAIMON_ratio_mean_all_TMDs'])
        # save pandas dataframe with values
        dfk.to_csv(pathdict["list_keywords_csv"], sep=",", quoting=csv.QUOTE_NONNUMERIC)
    else:
        sys.stdout.write ('no valid keywords found! change "cutoff_major_keywords" setting! current value: {}'.format(s['cutoff_major_keywords']))

    logging.info("\n~~~~~~~~~~~~        keyword_analysis is finished         ~~~~~~~~~~~~")