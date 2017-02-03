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

def gather_AAIMONs(pathdict, logging, s):
    logging.info("~~~~~~~~~~~~         starting gather_AAIMONs           ~~~~~~~~~~~~")
    df = pd.read_csv(pathdict["list_csv"], sep=",", quoting=csv.QUOTE_NONNUMERIC, index_col=0)

    dfg = pd.DataFrame()

    # iterate over the dataframe for proteins with an existing list_of_TMDs. acc = uniprot accession.
    for acc in df.loc[df['list_of_TMDs'].notnull()].loc[df['list_of_TMDs'] != 'nan'].index:
        protein_name = df.loc[acc, 'protein_name']
        logging.info(protein_name)
        if not os.path.exists(df.loc[acc, 'homol_cr_ratios_zip']):
            logging.info("{} Protein skipped. File does not exist".format(df.loc[acc, 'homol_cr_ratios_zip']))
            continue
        # open csv as pandas dataframe (note, it was originally a series, and contains only one column and an index)
        # set delete_corrupt=True so that if the expected csv is not in the zip, the wholezipfile will be deleted
        mean_ser_filename = "{}_cr_mean.csv".format(protein_name)
        mean_ser = utils.open_df_from_csv_zip(df.loc[acc, 'homol_cr_ratios_zip'], filename=mean_ser_filename, delete_corrupt=True)
        dfg = pd.concat([dfg,mean_ser], axis=1)

    # transpose dataframe dfg
    dfg = dfg.T

    # for the OMPdb dataset, there is no uniprot_entry_name
    uniprot_entry_name_in_df = "uniprot_entry_name" in df.columns
    if not uniprot_entry_name_in_df:
        dfg['uniprot_entry_name'] = "OMPdb_dataset"

    # calculate mean AAIMON for all TMDs
    for acc in dfg.loc[dfg['list_of_TMDs'].notnull()].loc[dfg['list_of_TMDs'] != 'nan'].index:

        # dict_AAIMON_mean = {}
        # for TMD in ast.literal_eval(dfg.loc[acc, 'list_of_TMDs']):
        #     dict_AAIMON_mean[TMD] = dfg.loc[acc, '%s_AAIMON_mean' % TMD]
        # dfg.loc[acc, 'AAIMON_mean_all_TMDs'] = np.mean(pd.to_numeric(pd.Series(list(dict_AAIMON_mean.values()))))
        #
        # # calculate mean normalised AAIMON_n for all TMDs
        # dict_AAIMON_mean_n = {}
        # for TMD in ast.literal_eval(dfg.loc[acc, 'list_of_TMDs']):
        #     dict_AAIMON_mean_n[TMD] = dfg.loc[acc, '%s_AAIMON_mean_n' % TMD]
        # dfg.loc[acc, 'AAIMON_mean_all_TMDs_n'] = np.mean(pd.to_numeric(pd.Series(list(dict_AAIMON_mean_n.values()))))

        # CODE COPIED FROM FIGS.PY
        dict_AAIMON_mean = {}
        dict_AAIMON_std = {}
        dict_AAIMON_mean_n = {}
        dict_AAIMON_std_n = {}
        dict_AASMON_ratio_mean = {}
        dict_AASMON_ratio_std = {}
        dict_AAIMON_slope_mean = {}
        dict_AAIMON_n_slope_mean = {}
        dict_TMD_perc_identity_mean_all_TMDs = {}
        for TMD in ast.literal_eval(dfg.loc[acc, 'list_of_TMDs']):
            dict_AAIMON_mean[TMD] = dfg.loc[acc, '%s_AAIMON_mean' % TMD]
            dict_AAIMON_std[TMD] = dfg.loc[acc, '%s_AAIMON_std' % TMD]
            dict_AAIMON_mean_n[TMD] = dfg.loc[acc, '%s_AAIMON_mean_n' % TMD]
            #dict_AAIMON_std_n[TMD] = dfg.loc[acc, '%s_AAIMON_std_n' % TMD]
            dict_AASMON_ratio_mean[TMD] = dfg.loc[acc, '%s_AASMON_ratio_mean' % TMD]
            dict_AASMON_ratio_std[TMD] = dfg.loc[acc, '%s_AASMON_ratio_std' % TMD]
            dict_AAIMON_slope_mean[TMD] = dfg.loc[acc, '%s_AAIMON_slope' %TMD]
            dict_AAIMON_n_slope_mean[TMD] = dfg.loc[acc, '%s_AAIMON_n_slope' % TMD]
            dict_TMD_perc_identity_mean_all_TMDs[TMD] = dfg.loc[acc, '%s_perc_ident_mean' % TMD]

        dfg.loc[acc, 'AAIMON_mean_all_TMDs'] = np.mean(pd.to_numeric(pd.Series(list(dict_AAIMON_mean.values()))))
        dfg.loc[acc, 'AAIMON_mean_all_TMDs_n'] = np.mean(pd.to_numeric(pd.Series(list(dict_AAIMON_mean_n.values()))))
        #dfg.loc[acc, 'AAIMON_std_all_TMDs_n'] = np.mean(pd.to_numeric(pd.Series(list(dict_AAIMON_std_n.values()))))
        dfg.loc[acc, 'AAIMON_std_all_TMDs'] = np.mean(pd.to_numeric(pd.Series(list(dict_AAIMON_std.values()))))
        dfg.loc[acc, 'AASMON_ratio_mean_all_TMDs'] = np.mean(pd.to_numeric(pd.Series(list(dict_AASMON_ratio_mean.values()))))
        dfg.loc[acc, 'AASMON_ratio_std_all_TMDs'] = np.mean(pd.to_numeric(pd.Series(list(dict_AASMON_ratio_std.values()))))
        dfg.loc[acc, 'AAIMON_slope_mean_all_TMDs'] = np.mean(pd.to_numeric(pd.Series(list(dict_AAIMON_slope_mean.values()))))
        dfg.loc[acc, 'AAIMON_n_slope_mean_all_TMDs'] = np.mean(pd.to_numeric(pd.Series(list(dict_AAIMON_n_slope_mean.values()))))
        dfg.loc[acc, 'TMD_perc_identity_mean_all_TMDs'] = np.mean(pd.to_numeric(pd.Series(list(dict_TMD_perc_identity_mean_all_TMDs.values()))))

        # count the number of TMDs for each protein
        dfg.loc[acc, 'number_of_TMDs'] = len(dfg.loc[acc, 'list_of_TMDs'].split(','))

        # add sequence length to dfg
        dfg.loc[acc, 'seqlen'] = df.loc[acc, 'seqlen']

        # # add total_number_of_simap_hits
        # dfg.loc[acc, 'total_number_of_simap_hits'] = dfg.loc[acc, 'TM01_AAIMON_n_homol']
        #
        # # add 'uniprot_entry_name'
        # if uniprot_entry_name_in_df:
        #     dfg.loc[acc, 'uniprot_entry_name'] = df.loc[acc, 'uniprot_entry_name']
        #
        # # add 'uniprot_KW'
        # if 'uniprot_KW' in df.columns:
        #     dfg.loc[acc, 'uniprot_KW'] = df.loc[acc, 'uniprot_KW']


    dfg.copy().to_csv(pathdict["list_cr_summary_csv"], sep=",", quoting=csv.QUOTE_NONNUMERIC)

    ########################################################################################
    #                                                                                      #
    #                    Save a huge dataframe with the AAIMONs for                        #
    #                    all homologues of all TMDs of all proteins                        #
    #                                                                                      #
    ########################################################################################

    if s['save_df_characterising_each_homol_TMD']:

        # defining cutoff for max and min number of homologues for each protein
        max_num_homologues = s['cutoff_max_characterising_each_homol_TMD']
        min_num_homologues = s['cutoff_min_characterising_each_homol_TMD']
        min_match_to_query_ratio = 0.9 # allowed length of truncated alignment is 90% coverage

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

        sys.stdout.write("\nLoading data\n")
        # initiate empty numpy array
        data = np.empty([0, 3])
        # navigate through filesystem and open pickles from .zip
        n=0
        for acc in dfg.index:
            n += 1
            if n % 20 == 0:
                sys.stdout.write('.'), sys.stdout.flush()
                if n % 600 == 0:
                    sys.stdout.write('\n'), sys.stdout.flush()
            protein_name = df.loc[acc, "protein_name"]
            homol_cr_ratios_zip = df.loc[acc, "homol_cr_ratios_zip"]
            if not os.path.isfile(homol_cr_ratios_zip):
                # skip to next protein
                continue
            for TMD in ast.literal_eval(df.loc[acc, "list_of_TMDs"]):
                # generate column names necessary for current file
                columns = ['obs_changes', '{}_AAIMON'.format(TMD), '{}_AAIMON_n'.format(TMD)]
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
        #data[:, 0] = data[:, 0] * 100

        # create bins, calculate mean and 95% confidence interval
        sys.stdout.write('\nBinning data - calculating 95% confidence interval\n')
        number_of_bins = s['specify_number_of_bins_characterising_TMDs']
        linspace_binlist = np.linspace(1, 100, number_of_bins)
        binwidth = 100/number_of_bins
        binned_data = np.empty([0, 8])
        conf_95 = np.array([1, 2])
        conf95_norm = np.array([1, 2])
        for percentage in linspace_binlist:
            if percentage % 5 == 0:
                sys.stdout.write('{}%, '.format(int(percentage))), sys.stdout.flush()
            bin_for_mean = np.empty([0, 3])
            for row in data:
                if row[0] < percentage and row[0] > percentage - binwidth:
                    bin_for_mean = np.concatenate((bin_for_mean, row.reshape(1, 3)))
            # calculate 95% conf. interv. in bin
            if bin_for_mean.size != 0:
                conf_95 = sms.DescrStatsW(bin_for_mean[:, 1]).tconfint_mean()
                # calculate 95% conf. interv. in bin _n
                conf95_norm = sms.DescrStatsW(bin_for_mean[:, 2]).tconfint_mean()
                mean_data_in_bin = np.array([percentage - binwidth/2,
                                             # calculate mean in bin
                                             bin_for_mean[:, 1].mean(),
                                             # calculate mean in bin _n
                                             bin_for_mean[:, 2].mean(),
                                             # add 95% conf. interv. results to np array
                                             conf_95[0], conf_95[1], conf95_norm[0], conf95_norm[1],
                                             # add the number of TMDs in bin to bin
                                             len(bin_for_mean[:, 0])])
                # merge data from bin to the others
                binned_data = np.concatenate((mean_data_in_bin.reshape(1, 8), binned_data))
        # drop every row containing nan in array
        binned_data = binned_data[~np.isnan(binned_data).any(axis=1)]
        '''
        description of columns in numpy arrays:

        numpy array data:
        |       0       |   1    |    2     |
        | % obs_changes | AAIMON | AAIMON_n |

        numpy array binned_Data:
        |       0       |      1      |       2       |     3    |    4    |      5     |     6     |         7          |
        | % obs_changes | mean AAIMON | mean AAIMON_n | CI95_low | CI95_hi | CI95_low_n | CI95_hi_n | number of Proteins |
        '''
        # save data and binned_data as zipped pickle
        with zipfile.ZipFile(pathdict['save_df_characterising_each_homol_TMD'], mode="w", compression=zipfile.ZIP_DEFLATED) as zipout:

            # save dataframe "data_filt" as pickle
            with open('data_characterising_each_homol_TMD.pickle', "wb") as f:
                pickle.dump(data, f, protocol=pickle.HIGHEST_PROTOCOL)
            zipout.write('data_characterising_each_homol_TMD.pickle', arcname='data_characterising_each_homol_TMD.pickle')
            os.remove('data_characterising_each_homol_TMD.pickle')

            # save dataframe "binned_data" as pickle
            with open('binned_data_characterising_each_homol_TMD.pickle', "wb") as f:
                pickle.dump(binned_data, f, protocol=pickle.HIGHEST_PROTOCOL)
            zipout.write('binned_data_characterising_each_homol_TMD.pickle', arcname='binned_data_characterising_each_homol_TMD.pickle')
            os.remove('binned_data_characterising_each_homol_TMD.pickle')

    logging.info("\n~~~~~~~~~~~~        gather_AAIMONs is finished         ~~~~~~~~~~~~")

def gather_pretty_alignments(pathdict, logging, s):
    """

    Parameters
    ----------
    pathdict
    logging
    s

    Returns
    -------

    """
    logging.info("\n~~~~~~~~~~~~    starting gather_pretty_alignments      ~~~~~~~~~~~~")
    #create an empty output file
    open(pathdict["pretty_alignments_csv"], 'w').close()
    #reopen to add match details iteratively from dictionary
    csvfile = open(pathdict["pretty_alignments_csv"], 'a')

    df = pd.read_csv(pathdict["list_csv"], sep=",", quoting=csv.QUOTE_NONNUMERIC, index_col=0)
    # drop proteins that have no list of TMDs
    df = df.loc[df['list_of_TMDs'].notnull()].loc[df['list_of_TMDs'] != 'nan']
    # iterate over the dataframe for proteins with an existing list_of_TMDs. acc = uniprot accession.
    num_TMDs_in_all_proteins_processed = 0
    for num, acc in enumerate(df.index):
        homol_cr_ratios_zip = df.loc[acc, "homol_cr_ratios_zip"]
        sys.stdout.write("{}, ".format(acc))
        if num % 20 == 0:
            sys.stdout.write("\n")
        sys.stdout.flush()
        # create an output dictionary, d, to hold the data for that TMD of that protein
        d = {}
        protein_name = df.loc[acc, "protein_name"]
        if os.path.isfile(homol_cr_ratios_zip):
            homol_df_orig_zip = df.loc[acc, "homol_df_orig_zip"]
            if os.path.isfile(homol_df_orig_zip):
                SIMAP_align_pretty_csv_filename = os.path.basename(df.loc[acc, "SIMAP_align_pretty_csv"])
                dfp = utils.open_df_from_csv_zip(homol_df_orig_zip, SIMAP_align_pretty_csv_filename)
                for TMD in ast.literal_eval(df.loc[acc, "list_of_TMDs"]):
                    TM_cr_pickle = "{}_{}_cr_df.pickle".format(protein_name, TMD)
                    # open dataframe  with function from korbinian, extract required columns, convert to np array
                    df_TMD = utils.open_df_from_pickle_zip(homol_cr_ratios_zip, TM_cr_pickle)
                    # if the pickle doesn't seem to hold the correct file
                    if df_TMD.empty:
                        # skip to next TMD
                        continue

                    max_gaps = s["cr_max_n_gaps_in_TMD"]
                    max_lipo_homol = s["max_lipo_homol"]
                    min_ident = s["cr_min_identity_of_TMD"]

                    """This is used as a filter in filter_and_save_fasta, therefore is conducted earlier in the slicing function. """
                    ## count the number of gaps in the query and match sequences
                    cr_TMD_query_str = '{TMD}_perc_ident >= {min_ident} & ' \
                                       '{TMD}_SW_query_num_gaps <= {max_gaps} & ' \
                                       '{TMD}_SW_match_num_gaps <= {max_gaps} & ' \
                                       '{TMD}_SW_match_lipo <= {max_lipo_homol}'.format(TMD=TMD, max_gaps=max_gaps,
                                                                                        max_lipo_homol=max_lipo_homol,
                                                                                        min_ident=min_ident)
                    # n_homol_before_filter = df_cr.shape[0]
                    # filter by the above query
                    df_TMD.query(cr_TMD_query_str, inplace=True)

                    # add protein and TMD-specific values
                    d["protein_name"] = protein_name
                    d["TMD"] = TMD
                    ########################################################################################
                    #                                                                                      #
                    #                Getting the hit number with the median AAIMON is                      #
                    #                    rather painful in pandas, but it works                            #
                    #                                                                                      #
                    ########################################################################################

                    AAIMON_ser = df_TMD["{}_AAIMON".format(TMD)].dropna()
                    if len(AAIMON_ser) == 0:
                        # skip to next TMD
                        continue

                    min_ = AAIMON_ser.idxmin()
                    max_ = AAIMON_ser.idxmax()
                    median_value = AAIMON_ser.median()
                    AAMON_minus_median = AAIMON_ser - median_value


                    # convert all to positive numbers
                    AAMON_minus_median = AAMON_minus_median.abs()
                    # sort so that the 0.00 is the first hit
                    AAMON_minus_median = AAMON_minus_median.sort_values()
                    # get index of first value (when AAIMON - median == 0, must be a hit with an AAIMON very close to the median)
                    med_ = AAMON_minus_median.index[0]
                    list_outliers = [min_, max_, med_]
                    list_outlier_names = ["min", "max", "med"]
                    # add/or overwrite values specific for that particular outlier
                    for m, outlier_index in enumerate(list_outliers):
                        outlier_name = list_outlier_names[m]
                        d["outlier"] = outlier_name
                        d["hit"] = outlier_index
                        columns = ['FASTA_gapped_identity', 'obs_changes', "{}_AAIMON", '{}_perc_ident', 'nonTMD_perc_ident', '{}_start_in_SW_alignment', '{}_SW_query_seq', '{}_SW_markup_seq',
                                   '{}_SW_match_seq', '{}_ratio_len_TMD_to_len_nonTMD', '{}_SW_align_len', "{}_SW_match_lipo"]
                        col_names = ['FASTA_gapped_identity', 'obs_changes', "AAIMON", 'TM_perc_ident', 'nonTMD_perc_ident', 'TM_start_in_SW_alignment', 'SW_query_seq', 'SW_markup_seq', 'SW_match_seq',
                                     'ratio_len_TMD_to_len_nonTMD', 'SW_align_len', "SW_match_lipo"]
                        for n, col in enumerate(columns):
                            col_name = col_names[n]
                            value = df_TMD.loc[outlier_index, col.format(TMD)]
                            # add each one to the dictionary
                            d[col_name] = value
                        d["TM_align"] = "{}\r\r\n\r\r{}\r\r\n\r\r{}".format(d['SW_query_seq'], d['SW_markup_seq'], d['SW_match_seq'])
                        # # remove individual seqs from dictionary
                        # del d['SW_query_seq']
                        # del d['SW_markup_seq']
                        # del d['SW_match_seq']
                        # add the pretty alignment
                        d["align_pretty"] = dfp.loc[min_, "align_pretty"]


                        if num_TMDs_in_all_proteins_processed == 0:
                            # sort
                            csv_header = ["protein_name", "TMD", "outlier", "TM_align","SW_match_lipo", "align_pretty", 'FASTA_gapped_identity', 'obs_changes', "AAIMON", "hit", 'TM_perc_ident', 'nonTMD_perc_ident', 'TM_start_in_SW_alignment', 'SW_query_seq', 'SW_markup_seq', 'SW_match_seq',
                                     'ratio_len_TMD_to_len_nonTMD', 'SW_align_len']
                            # make sure that the csv header is up-to-date, and isn't missing items from dict
                            assert len(csv_header) is len(d)
                            # save the csv header to the csv file
                            writer = csv.writer(csvfile, delimiter=',', quotechar='"', lineterminator='\n', quoting=csv.QUOTE_NONNUMERIC, doublequote=True)
                            writer.writerow(csv_header)
                        # save the dict values as a line in the csv file
                        writer = csv.DictWriter(csvfile, fieldnames=csv_header,
                                                extrasaction='ignore', delimiter=',', quotechar='"',
                                                lineterminator='\n', quoting=csv.QUOTE_NONNUMERIC,
                                                doublequote=True)
                        writer.writerow(d)
                        num_TMDs_in_all_proteins_processed += 1

    csvfile.close()
    logging.info("\n~~~~~~~~~~~~    finished gather_pretty_alignments      ~~~~~~~~~~~~")