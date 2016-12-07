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

    if s['save_df_characterising_each_homol_TMD']:

        # dfg contains all informations for initialisation of navigation through file system
        df_summary = dfg

        # defining cutoff for max and min number of homologues for each protein
        max_num_homologues = s['cutoff_max_characterising_each_homol_TMD']
        min_num_homologues = s['cutoff_min_characterising_each_homol_TMD']

        # filter summary file for min and max number of homologues based on TM01 number of homologues
        sys.stdout.write('Dropped homologues after filtering: \n')
        list_of_acc_to_keep = []
        for acc in df_summary.loc[df_summary['list_of_TMDs'].notnull()].loc[df_summary['list_of_TMDs'] != 'nan'].index:
            a = pd.to_numeric(df_summary.loc[acc, 'TM01_AAIMON_n_homol'])
            if max_num_homologues > a > min_num_homologues:
                list_of_acc_to_keep.append(acc)
            if not acc in list_of_acc_to_keep:
                df_summary = df_summary.drop([acc])
                print (acc)

        # save relevant parts for navigation through file system (database) in dictionaries
        dict_TMDs = {}
        for acc in df_summary.loc[df_summary['list_of_TMDs'].notnull()].loc[df_summary['list_of_TMDs'] != 'nan'].index:
            dict_TMDs[acc] = df_summary.loc[acc, 'list_of_TMDs']

        dict_folder = {}
        for key in dict_TMDs.keys():
            dict_folder[key] = key[:2]

        dict_uniprot_entry = {}
        for acc in df_summary.loc[df_summary['list_of_TMDs'].notnull()].loc[df_summary['list_of_TMDs'] != 'nan'].index:
            dict_uniprot_entry[acc] = df_summary.loc[acc, 'uniprot_entry_name']

        data_dir = s['data_dir']
        print ("\nLoading data")
        # initiate empty numpy array
        data = np.empty([0, 3])
        # navigate through filesystem and open pickles from .zip
        for key in dict_folder.keys():
            sys.stdout.write('.')
            sys.stdout.flush()
            in_zipfile = os.path.join(data_dir, 'homol', '{a}'.format(a=dict_folder[key]), '{b}_{c}_cr_ratios.zip'.format(b=key,
                                                                                           c=dict_uniprot_entry[key]))
            #print ('current directory: {a}' .format(a=in_zipfile))
            for TMD in ast.literal_eval(dict_TMDs[key]):
                filename = "{a}_{c}_{b}_cr_df.pickle".format(a=key, b=TMD, c=dict_uniprot_entry[key])
                # filename_csv = "{a}_{c}_AAIMON_normalisation_data.csv".format(a=key, c=dict_uniprot_entry[key])
                #print ('current file: {a}' .format(a=filename))
                # generate column names necessary for current file
                columns = ['FASTA_gapped_identity', '{a}_AAIMON_ratio'.format(a=TMD),'{a}_AAIMON_ratio_n'.format(a=TMD)]
                # open dataframe  with function from korbinian, extract required columns, convert to np array
                if not os.path.isfile(in_zipfile):
                    # skip to next TMD or protein
                    continue # SHOULD THIS BE A BREAK?
                df = utils.open_df_from_pickle_zip(in_zipfile, filename)
                if columns[2] not in df.columns:
                    # file is old, and should be deleted
                    os.remove(in_zipfile)
                    logging.info("{} file is presumed out of date, and has been deleted".format(in_zipfile))
                    # skip to next TMD or protein
                    continue # SHOULD THIS BE A BREAK?
                df = df[columns].as_matrix()
                # join output data file with currently opened dataframe
                data = np.concatenate((data, df))
        # drop every row with nan
        data = data[~np.isnan(data).any(axis=1)]
        # create real percentage values, multiply column 1 with 100
        data[:, 0] = data[:, 0] * 100

    # create bins, calculate mean and 95% confidence interval
    sys.stdout.write('\nBinning data - calculate 95% confidence interval\n')
    number_of_bins = s['specify_number_of_bins_characterising_TMDs']
    linspace_binlist = np.linspace(1, 100, number_of_bins)
    binned_data = np.empty([0, 7])
    conf_95 = np.array([1, 2])
    conf95_norm = np.array([1, 2])
    for percentage in linspace_binlist:
        sys.stdout.write(".")
        sys.stdout.flush()
        bin_for_mean = np.empty([0, 3])
        for line in data:
            if line[0] < percentage and line[0] > percentage - 1:
                bin_for_mean = np.concatenate((bin_for_mean, line.reshape(1, 3)))
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

