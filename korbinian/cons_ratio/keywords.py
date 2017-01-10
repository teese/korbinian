import pandas as pd
import numpy as np
import ast
import csv
import sys
import itertools
import matplotlib.pyplot as plt

def keyword_analysis(pathdict, s, logging):
    logging.info("~~~~~~~~~~~~         starting keyword_analysis           ~~~~~~~~~~~~")
    # load cr_summary file
    dfc = pd.read_csv(pathdict["list_cr_summary_csv"], sep=",", quoting=csv.QUOTE_NONNUMERIC, index_col=0)
    # load summary file
    dfu = pd.read_csv(pathdict["list_summary_csv"], sep=",", quoting=csv.QUOTE_NONNUMERIC, index_col=0)

    # merge cr_summary and summary file
    df = pd.merge(dfc, dfu, left_index=True, right_index=True, suffixes=('', '_dfu'))


    # from old keyword analysis
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
    for KW in list_ignored_KW:
        if KW in KW_counts_major.index:
            KW_counts_major = KW_counts_major.drop(KW)

    dict_KW_counts_major = pd.Series.to_dict(KW_counts_major)
    sys.stdout.write('\navailable keywords - number of proteins containing keyword\n\n')
    for key in dict_KW_counts_major:
        sys.stdout.write('{}    -   {}\n'.format(key, dict_KW_counts_major[key]))


    # keyword to analyze
    keyword = s['keyword_to_analyse']
    # initialise lists of acc that do or do not contain the keyword
    list_of_acc_containing_kw = []
    list_of_acc_without_kw = []
    for acc in df.index:
        if keyword in df.loc[acc, 'uniprot_KW']:
            list_of_acc_containing_kw.append(acc)
        else:
            list_of_acc_without_kw.append(acc)

    # initialise empty pandas dataframe to hold values for each keyword
    dfk = pd.DataFrame

    # dropping all proteins that do or do not contain the keyword - create two dataframes with and without keyword
    df_keyword = df.loc[list_of_acc_containing_kw, :]
    df_no_keyword = df.loc[list_of_acc_without_kw, :]

    # calculate mean and std of AAIMON_ratio_mean_all_TMDs
    AAIMON_keyword_mean = np.mean(df_keyword['AAIMON_ratio_mean_all_TMDs'])
    AAIMON_keyword_std = np.std(df_keyword['AAIMON_ratio_mean_all_TMDs'])
    AAIMON_no_keyword_mean = np.mean(df_no_keyword['AAIMON_ratio_mean_all_TMDs'])
    AAIMON_no_keyword_std = np.std(df_no_keyword['AAIMON_ratio_mean_all_TMDs'])
    number_of_proteins_keyword = len(df_keyword.index)
    number_of_proteins_no_keyword = len(df_no_keyword.index)
    # calculate odds ratio
    odds_ratio = AAIMON_keyword_mean / AAIMON_no_keyword_mean
    sys.stdout.write ('mean AAIMON containing keyword  "{a}" : {b:.5f} ± {c:.5f}, number of proteins: {d}\n'
                      'mean AAIMON   without  keyword  "{a}" : {e:.5f} ± {f:.5f}, number of proteins: {g}\n'
           'odds_ratio = {h:.5f}'.format(a=keyword, b=AAIMON_keyword_mean, c=AAIMON_keyword_std, d=number_of_proteins_keyword,
                                         e=AAIMON_no_keyword_mean, f=AAIMON_no_keyword_std, g=number_of_proteins_no_keyword, h=odds_ratio))





    logging.info("\n~~~~~~~~~~~~        keyword_analysis is finished         ~~~~~~~~~~~~")