import pandas as pd
import numpy as np
import ast
import csv
import matplotlib.pyplot as plt

def keyword_analysis(pathdict, s, logging):
    logging.info("~~~~~~~~~~~~         starting keyword_analysis           ~~~~~~~~~~~~")
    # load cr_summary file
    dfc = pd.read_csv(pathdict["list_cr_summary_csv"], sep=",", quoting=csv.QUOTE_NONNUMERIC, index_col=0)
    # load summary file
    dfu = pd.read_csv(pathdict["list_summary_csv"], sep=",", quoting=csv.QUOTE_NONNUMERIC, index_col=0)

    # merge cr_summary and summary file
    df = pd.merge(dfc, dfu, left_index=True, right_index=True, suffixes=('', '_dfu'))

    # keyword to analyze
    keyword = s['keyword_to_analyse']
    # initialise lists of acc that do or do not contain the keyword
    list_of_acc_containing_kw = []
    list_of_acc_without_kw = []
    for acc in df.index:
        if keyword in ast.literal_eval(df.loc[acc, 'uniprot_KW']):
            list_of_acc_containing_kw.append(acc)
        else:
            list_of_acc_without_kw.append(acc)
    # log if keyword is not available
    if not list_of_acc_containing_kw:
        print ('your specified keyword  > {} <  is not available - cannot perform keyword analysis' .format(keyword))
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
    print ('mean AAIMON containing keyword "{a}": {b:.5f} ± {c:.5f}, number of proteins: {d}\nmean AAIMON without keyword "{a}":  {e:.5f} ± {f:.5f}, number of proteins: {g}\nodds_ratio = {h:.5f}'.format(
        a=keyword, b=AAIMON_keyword_mean, c=AAIMON_keyword_std, d=number_of_proteins_keyword, e=AAIMON_no_keyword_mean, f=AAIMON_no_keyword_std, g=number_of_proteins_no_keyword, h=odds_ratio))

    logging.info("\n~~~~~~~~~~~~        keyword_analysis is finished         ~~~~~~~~~~~~")