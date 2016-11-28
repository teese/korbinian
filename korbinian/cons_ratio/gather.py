import ast
import csv
import numpy as np
import os
import korbinian.utils as utils
import pandas as pd

def gather_AAIMON_ratios(pathdict, logging):
    logging.info("~~~~~~~~~~~~         starting gather_AAIMON_ratios           ~~~~~~~~~~~~")
    df = pd.read_csv(pathdict["list_summary_csv"], sep=",", quoting=csv.QUOTE_NONNUMERIC, index_col=0)

    dfg = pd.DataFrame()

    # iterate over the dataframe for proteins with an existing list_of_TMDs. acc = uniprot accession.
    for acc in df.loc[df['list_of_TMDs'].notnull()].loc[df['list_of_TMDs'] != 'nan'].index:
        protein_name = df.loc[acc, 'protein_name']
        logging.info(protein_name)
        if not os.path.exists(df.loc[acc, 'homol_cr_ratios_zip']):
            logging.info("{} Protein skipped. File does not exist".format(df.loc[acc, 'homol_cr_ratios_zip']))
            continue
        # open csv as pandas dataframe (note, it was originally a series, and contains only one column and an index)
        mean_ser_filename = "{}_cr_mean.csv".format(protein_name)
        mean_ser = utils.open_df_from_csv_zip(df.loc[acc, 'homol_cr_ratios_zip'], filename=mean_ser_filename)
        dfg = pd.concat([dfg,mean_ser], axis=1)

    # transpose dataframe dfg
    dfg = dfg.T

    # calculate mean AAIMON for all TMDs
    for acc in dfg.loc[dfg['list_of_TMDs'].notnull()].loc[dfg['list_of_TMDs'] != 'nan'].index:
        dict_AAIMON_ratio_mean = {}
        for TMD in ast.literal_eval(dfg.loc[acc, 'list_of_TMDs']):
            dict_AAIMON_ratio_mean[TMD] = dfg.loc[acc, '%s_AAIMON_ratio_mean' % TMD]
        dfg.loc[acc, 'AAIMON_ratio_mean_all_TMDs'] = np.mean(pd.to_numeric(pd.Series(list(dict_AAIMON_ratio_mean.values()))))

    # calculate mean normalised AAIMON_n for all TMDs
    for acc in dfg.loc[dfg['list_of_TMDs'].notnull()].loc[dfg['list_of_TMDs'] != 'nan'].index:
        dict_AAIMON_ratio_mean_n = {}
        for TMD in ast.literal_eval(dfg.loc[acc, 'list_of_TMDs']):
            dict_AAIMON_ratio_mean_n[TMD] = dfg.loc[acc, '%s_AAIMON_ratio_mean_n' % TMD]
        dfg.loc[acc, 'AAIMON_ratio_mean_all_TMDs_n'] = np.mean(
            pd.to_numeric(pd.Series(list(dict_AAIMON_ratio_mean_n.values()))))

    # count the number of TMDs for each protein
    for acc in dfg.loc[dfg['list_of_TMDs'].notnull()].loc[dfg['list_of_TMDs'] != 'nan'].index:
        dfg.loc[acc, 'number_of_TMDs'] = len(dfg.loc[acc, 'list_of_TMDs'].split(','))

    # add sequence length to dfg
    for acc in dfg.loc[dfg['list_of_TMDs'].notnull()].loc[dfg['list_of_TMDs'] != 'nan'].index:
        dfg.loc[acc, 'seqlen'] = df.loc[acc, 'seqlen']

    # add total_number_of_simap_hits
    for acc in dfg.loc[dfg['list_of_TMDs'].notnull()].loc[dfg['list_of_TMDs'] != 'nan'].index:
        dfg.loc[acc, 'total_number_of_simap_hits'] = dfg.loc[acc, 'TM01_AAIMON_n_homol']

    # add 'uniprot_entry_name'
    for acc in dfg.loc[dfg['list_of_TMDs'].notnull()].loc[dfg['list_of_TMDs'] != 'nan'].index:
        dfg.loc[acc, 'uniprot_entry_name'] = df.loc[acc, 'uniprot_entry_name']


    dfg.to_csv(pathdict["list_cr_summary_csv"], sep=",", quoting=csv.QUOTE_NONNUMERIC)
    logging.info("~~~~~~~~~~~~        gather_AAIMON_ratios is finished         ~~~~~~~~~~~~")

