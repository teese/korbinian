import ast
import numpy as np
import csv
import os
import pickle
import korbinian
import korbinian.mtutils as utils
import pandas as pd
import zipfile

def gather_AAIMON_ratios(pathdict, set_, logging):
    logging.info("~~~~~~~~~~~~         starting gather_AAIMON_ratios           ~~~~~~~~~~~~")
    df = pd.read_csv(pathdict["list_summary_csv"], sep=",", quoting=csv.QUOTE_NONNUMERIC, index_col=0)

    dfg = pd.DataFrame()

    # iterate over the dataframe for proteins with an existing list_of_TMDs. acc = uniprot accession.
    for acc in df.loc[df['list_of_TMDs'].notnull()].loc[df['list_of_TMDs'] != 'nan'].index:
        logging.info(df.loc[acc, 'protein_name'])
        if not os.path.exists(df.loc[acc, 'homol_cr_ratios_zip']):
            logging.info("{} Protein skipped. File does not exist".format(df.loc[acc, 'homol_cr_ratios_zip']))
            continue
        # open csv as pandas dataframe (note, it was originally a series, and contains only one column and an index)
        mean_ser_filename = "{}_cr_mean.csv".format(acc)
        mean_ser = utils.open_df_from_csv_zip(df.loc[acc, 'homol_cr_ratios_zip'], filename=mean_ser_filename)
        dfg = pd.concat([dfg,mean_ser], axis=1)

    dfg.T.to_csv(pathdict["list_cr_summary_csv"], sep=",", quoting=csv.QUOTE_NONNUMERIC)
    logging.info("~~~~~~~~~~~~        gather_AAIMON_ratios is finished         ~~~~~~~~~~~~")
