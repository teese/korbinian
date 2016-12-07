import ast
import csv
import numpy as np
import os
import korbinian.utils as utils
import pandas as pd
import statsmodels.stats.api as sms
import pickle
import sys
import zipfile

def gather_AAIMON_ratios(pathdict, logging, s):
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

    # for the OMPdb dataset, there is no uniprot_entry_name
    uniprot_entry_name_in_df = "uniprot_entry_name" in df.columns
    if not uniprot_entry_name_in_df:
        dfg['uniprot_entry_name'] = "OMPdb_dataset"

    # calculate mean AAIMON for all TMDs
    for acc in dfg.loc[dfg['list_of_TMDs'].notnull()].loc[dfg['list_of_TMDs'] != 'nan'].index:

        dict_AAIMON_ratio_mean = {}
        for TMD in ast.literal_eval(dfg.loc[acc, 'list_of_TMDs']):
            dict_AAIMON_ratio_mean[TMD] = dfg.loc[acc, '%s_AAIMON_ratio_mean' % TMD]
        dfg.loc[acc, 'AAIMON_ratio_mean_all_TMDs'] = np.mean(pd.to_numeric(pd.Series(list(dict_AAIMON_ratio_mean.values()))))

        # calculate mean normalised AAIMON_n for all TMDs
        dict_AAIMON_ratio_mean_n = {}
        for TMD in ast.literal_eval(dfg.loc[acc, 'list_of_TMDs']):
            dict_AAIMON_ratio_mean_n[TMD] = dfg.loc[acc, '%s_AAIMON_ratio_mean_n' % TMD]
        dfg.loc[acc, 'AAIMON_ratio_mean_all_TMDs_n'] = np.mean(
            pd.to_numeric(pd.Series(list(dict_AAIMON_ratio_mean_n.values()))))

        # count the number of TMDs for each protein
        dfg.loc[acc, 'number_of_TMDs'] = len(dfg.loc[acc, 'list_of_TMDs'].split(','))

        # add sequence length to dfg
        dfg.loc[acc, 'seqlen'] = df.loc[acc, 'seqlen']

        # add total_number_of_simap_hits
        dfg.loc[acc, 'total_number_of_simap_hits'] = dfg.loc[acc, 'TM01_AAIMON_n_homol']

        # add 'uniprot_entry_name'
        if uniprot_entry_name_in_df:
            dfg.loc[acc, 'uniprot_entry_name'] = df.loc[acc, 'uniprot_entry_name']

        # add 'uniprot_KW'
        if 'uniprot_KW' in df.columns:
            dfg.loc[acc, 'uniprot_KW'] = df.loc[acc, 'uniprot_KW']


    dfg.copy().to_csv(pathdict["list_cr_summary_csv"], sep=",", quoting=csv.QUOTE_NONNUMERIC)

    ########################################################################################
    #                                                                                      #
    #               Save a huge dataframe with the AAIMONs for                             #
    #                all homologues of all TMDs of all proteins                            #
    #                                                                                      #
    ########################################################################################

    if s['save_df_characterising_each_homol_TMD']:

        # defining cutoff for max and min number of homologues for each protein
        max_num_homologues = s['cutoff_max_characterising_each_homol_TMD']
        min_num_homologues = s['cutoff_min_characterising_each_homol_TMD']

        # filter summary file for min and max number of homologues based on TM01 number of homologues
        sys.stdout.write('Dropped homologues after filtering: \n')
        list_of_acc_to_keep = []
        for acc in dfg.loc[dfg['list_of_TMDs'].notnull()].loc[dfg['list_of_TMDs'] != 'nan'].index:
            TM01_AAIMON_n_homol = pd.to_numeric(dfg.loc[acc, 'TM01_AAIMON_n_homol'])
            if max_num_homologues > TM01_AAIMON_n_homol > min_num_homologues:
                list_of_acc_to_keep.append(acc)

        # keep only proteins that have the desired number of homologues
        dfg = dfg.loc[list_of_acc_to_keep, :]
        df = df.loc[list_of_acc_to_keep, :]

        # generate column names necessary for current file
        columns = ['FASTA_gapped_identity', '{}_AAIMON_ratio'.format(TMD), '{}_AAIMON_ratio_n'.format(TMD)]

        sys.stdout.write("\nLoading data\n")
        # initiate empty numpy array
        data = np.empty([0, 3])
        # navigate through filesystem and open pickles from .zip
        for acc in dfg.index:
            sys.stdout.write('.'), sys.stdout.flush()
            protein_name = df.loc[acc, "protein_name"]
            homol_cr_ratios_zip = df.loc[acc, "homol_cr_ratios_zip"]
            if not os.path.isfile(homol_cr_ratios_zip):
                # skip to next protein
                continue
            for TMD in ast.literal_eval(df.loc[acc, "list_of_TMDs"]):
                TM_cr_pickle = "{}_{}_cr_df.pickle".format(protein_name, TMD)
                # open dataframe  with function from korbinian, extract required columns, convert to np array
                df_TMD = utils.open_df_from_pickle_zip(homol_cr_ratios_zip, TM_cr_pickle)
                if columns[2] not in df_TMD.columns:
                    # file is old, and should be deleted
                    #os.remove(homol_cr_ratios_zip)
                    logging.info("{} file is presumed out of date, and WILL IN THE FUTURE been deleted".format(homol_cr_ratios_zip))
                    # skip to next protein
                    break
                df_TMD = df_TMD[columns].as_matrix()
                # join output data file with currently opened dataframe
                data = np.concatenate((data, df_TMD))
        # drop every row with nan
        data = data[~np.isnan(data).any(axis=1)]
        # create real percentage values, multiply column 1 with 100
        data[:, 0] = data[:, 0] * 100

    # create bins, calculate mean and 95% confidence interval
    sys.stdout.write('\nBinning data - calculating 95% confidence interval\n')
    number_of_bins = s['specify_number_of_bins_characterising_TMDs']
    linspace_binlist = np.linspace(1, 100, number_of_bins)
    binwidth = 100/number_of_bins
    binned_data = np.empty([0, 7])
    conf_95 = np.array([1, 2])
    conf95_norm = np.array([1, 2])
    for percentage in linspace_binlist:
        sys.stdout.write("."), sys.stdout.flush()
        bin_for_mean = np.empty([0, 3])
        for row in data:
            if row[0] < percentage and row[0] > percentage - binwidth:
                bin_for_mean = np.concatenate((bin_for_mean, row.reshape(1, 3)))
                # HERE DROP THE NANs????
        #bin_for_mean = bin_for_mean[~np.isnan(bin_for_mean).any(axis=1)]
        # calculate 95% conf. interv. in bin
        conf_95 = sms.DescrStatsW(bin_for_mean[:, 1]).tconfint_mean()
        # calculate 95% conf. interv. in bin _n
        conf95_norm = sms.DescrStatsW(bin_for_mean[:, 2]).tconfint_mean()
        mean_data_in_bin = np.array([percentage,
                                     # calculate mean in bin
                                     bin_for_mean[:, 1].mean(),
                                     # calculate mean in bin _n
                                     bin_for_mean[:, 2].mean(),
                                     # add 95% conf. interv. results to np array
                                     conf_95[0], conf_95[1], conf95_norm[0], conf95_norm[1]])
        # merge data from bin to the others
        binned_data = np.concatenate((mean_data_in_bin.reshape(1, 7), binned_data))
    # drop every row containing nan in array
    binned_data = binned_data[~np.isnan(binned_data).any(axis=1)]

    # save data and binned_data as zipped pickle

    with zipfile.ZipFile(pathdict['save_df_characterising_each_homol_TMD'], mode="w", compression=zipfile.ZIP_DEFLATED) as zipout:

        # save dataframe "data" as pickle
        with open('data_characterising_each_homol_TMD.pickle', "wb") as f:
            pickle.dump(data, f, protocol=pickle.HIGHEST_PROTOCOL)
        zipout.write('data_characterising_each_homol_TMD.pickle', arcname='data_characterising_each_homol_TMD.pickle')
        os.remove('data_characterising_each_homol_TMD.pickle')

        # save dataframe "binned_data" as pickle
        with open('binned_data_characterising_each_homol_TMD.pickle', "wb") as f:
            pickle.dump(binned_data, f, protocol=pickle.HIGHEST_PROTOCOL)
        zipout.write('binned_data_characterising_each_homol_TMD.pickle', arcname='binned_data_characterising_each_homol_TMD.pickle')
        os.remove('binned_data_characterising_each_homol_TMD.pickle')

    logging.info("\n~~~~~~~~~~~~        gather_AAIMON_ratios is finished         ~~~~~~~~~~~~")

