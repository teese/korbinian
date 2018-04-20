import ast
import csv
import numpy as np
import os
import korbinian.utils as utils
import korbinian
import pandas as pd
#from statsmodels.stats.api import DescrStatsW
from statsmodels.stats import weightstats as sms
import pickle
import re
import sys
import zipfile
# import debugging tools
from korbinian.utils import pr, pc, pn, aaa

def gather_AAIMONs(pathdict, logging, s):
    """ Gathers the AAIMON ratios and slopes for each protein, created by the run_calculate_AAIMONs scripts.

    To be compatible with multiprocessing, the run_calculate_AAIMONs script creates a separate output summary file
    for each protein. The gather_AAIMONs simply concatenates all of these files together

    Note that the gather_AAIMONs script does not do ANY filtering. This is all done earlier by run_calculate_AAIMONs.
    It is assumed that for homologues that did not pass the filter (e.g., because they had X in the sequence),
    that no AAIMON or AAIMON slope was calculated.

    Parameters
    ----------
    pathdict : dict
        Dictionary of the key paths and files associated with that List number.
    logging : logging.Logger
        Logger for printing to console and logfile.
    s : dict
        Settings dictionary extracted from excel settings file.

    Saved Files and Figures
    -----------------------
    list_cr_summary_csv : csv
        comma separated csv file with the AAIMON ratios etc
        contains all data within the {}_cr_mean.csv summary file for each protein
    pretty_alignments_csv : csv
        comma separated csv file with the pretty alignments of all the outliers
    data_characterising_each_homol_TMD.pickle : pickle
        Raw AAIMON and % identity (or % aa sub rate) datapoints for all TM of all homologues of all proteins
        Used to create large scatterplot of all datapoints.

    Returns
    -------
    pathdict : dict
        Dictionary of the key paths and files associated with that List number.
        In special cases, the pathdict is modified.
    """
    logging.info("~~~~~~~~~~~~                           starting gather_AAIMONs                      ~~~~~~~~~~~~")
    df = pd.read_csv(pathdict["list_csv"], sep=",", quoting=csv.QUOTE_NONNUMERIC, index_col=0)
    # drop any proteins without a list of TMDs
    df = df.loc[df['list_of_TMDs'].notnull()].loc[df['list_of_TMDs'] != 'nan']
    # convert list_of_TMDs from string to python list
    df['list_of_TMDs'] = df.list_of_TMDs.apply(lambda x: ast.literal_eval(x))

    ###############################################################
    #                                                             #
    #                        Filter keywords                      #
    #                                                             #
    ###############################################################
    if s['filter_keywords_in_gather']:
        # filter list file by keywords for exclusion analysis, e.g. enzyme only
        list_number = s['list_number']
        # specify allowed and disallowed keywords
        allowed_KW = ast.literal_eval(s['gather_filter_allowed_keywords'])
        disallowed_KW = ast.literal_eval(s['gather_filter_forbidden_keywords'])
        # generate new pathdict
        base_filename_summaries = os.path.join(s["data_dir"], "summaries", '%02d' % list_number, 'List%02d_filtered' % list_number, ' - '.join(allowed_KW),'List%02d' % list_number)
        pathdict = korbinian.common.create_pathdict(base_filename_summaries, s)
        # create new folder with new pathdict
        if not os.path.exists(base_filename_summaries[:-7]):
            os.makedirs(base_filename_summaries[:-7])

        # copy keyword column, apply ast.literal_eval to the copied column
        df['KW'] = df['uniprot_KW']
        # apply ast.literal_eval to every item in df['uniprot_KW']
        if isinstance(df['KW'][0], str):
            df['KW'] = df['KW'].apply(lambda x: ast.literal_eval(x))
        # get list of enzyme keywords
        list_enzyme_KW, list_ignored_KW, PFAM_dict = korbinian.cons_ratio.keywords.get_list_enzyme_KW_and_list_ignored_KW()
        # create column per allowed keyword that holds a bool if keyword is present in that protein
        for KW in allowed_KW:
            if KW == 'Enzyme':
                df['Enzyme'] = df['KW'].apply(korbinian.cons_ratio.keywords.KW_list_contains_any_desired_KW, args=(list_enzyme_KW,))
            else:
                df[KW] = df['KW'].apply(korbinian.cons_ratio.keywords.KW_list_contains_any_desired_KW, args=([KW],))
        # create column for every protein holding bool if protein contains at least one of the allowed keywords
        for acc in df.index:
            df.loc[acc, 'keep'] = df.loc[acc, allowed_KW].any()
        # drop all proteins whose keywords do not match the requirements
        df = df.loc[df['keep'] == True]

        # drop all proteins that contain one of the disallowed keywords
        for KW in disallowed_KW:
            if KW == 'Enzyme':
                df['Enzyme'] = df['KW'].apply(korbinian.cons_ratio.keywords.KW_list_contains_any_desired_KW, args=(list_enzyme_KW,))
            else:
                df[KW] = df['KW'].apply(korbinian.cons_ratio.keywords.KW_list_contains_any_desired_KW, args=([KW],))
            df = df.loc[df[KW] == False]
        # remove copied and edited keyword list
        df = df.drop('KW', 1)

        df.to_csv(pathdict["list_csv"], sep=",", quoting=csv.QUOTE_NONNUMERIC)

    #############################################################################
    #                                                                           #
    #       Collate all the "_cr_mean.csv" files into a single dataframe        #
    #                                                                           #
    #############################################################################

    dfg = pd.DataFrame()
    # iterate over the dataframe for proteins with an existing list_of_TMDs
    for acc in df.index:
        protein_name = df.loc[acc, 'protein_name']
        #logging.info(protein_name)
        sys.stdout.write("{}, ".format(acc)), sys.stdout.flush()
        if not os.path.exists(df.loc[acc, 'homol_cr_ratios_zip']):
            logging.info("{} skipped. homol_cr_ratios_zip does not exist".format(acc))
            continue

        if utils.file_is_old(df.loc[acc, 'homol_cr_ratios_zip'], s["oldest_acceptable_file_date"]):
            os.remove(df.loc[acc, 'homol_cr_ratios_zip']),
            logging.info("{} skipped, file is old and has been deleted".format(acc))
            continue

        # open csv as pandas dataframe (note, it was originally a series, and contains only one column and an index)
        # set delete_corrupt=True so that if the expected csv is not in the zip, the wholezipfile will be deleted
        mean_ser_filename = "{}_cr_mean.csv".format(acc)
        mean_ser = utils.open_df_from_csv_zip(df.loc[acc, 'homol_cr_ratios_zip'], filename=mean_ser_filename, delete_corrupt=True)
        dfg = pd.concat([dfg,mean_ser], axis=1)

    if dfg.empty:
        raise ValueError("\n\ndfg is an empty dataframe.\nThis means that none of the proteins had any correctly processed conservation ratios.\nSuggest checking the output of all previous steps.")

    # transpose dataframe (flip index and columns)
    dfg = dfg.T.copy()

    # for the OMPdb dataset, there is no uniprot_entry_name
    uniprot_entry_name_in_df = "uniprot_entry_name" in df.columns
    if not uniprot_entry_name_in_df:
        dfg['uniprot_entry_name'] = "OMPdb_dataset"

    # drop any proteins in dfg without a list of TMDs
    dfg = dfg.loc[df['list_of_TMDs'].notnull()].loc[dfg['list_of_TMDs'] != 'nan']

    # if the list_of_TMDs is a stringlist, convert to a python list
    dfg['list_of_TMDs'] = dfg['list_of_TMDs'].dropna().apply(lambda x : ast.literal_eval(x))

    # # for singlepass datasets, leave row blank by default
    # dfg['AAIMON_slope_central_TMDs'] = np.nan

    # CONVERT STRINGS TO FLOATS FOR SELECTED COLUMNS
    # note that after saving dfg to CSV, pandas then gets the dtype correct upon reopening for figs.py etc
    cols_to_convert = ["AAIMON_mean_all_TM_res", "AAIMON_n_mean_all_TM_res", "AAIMON_slope_all_TM_res", "AAIMON_n_slope_all_TM_res", 'AAIMON_n_homol']
    for col in cols_to_convert:
        dfg[col] = pd.to_numeric(dfg[col])
    # print out mean AAIMON values in dataset
    mean_AAIMON_in_dataset = dfg['AAIMON_mean_all_TM_res'].mean()
    mean_AAIMON_n_in_dataset = dfg['AAIMON_n_mean_all_TM_res'].mean()
    mean_AAIMON_slope_in_dataset = dfg['AAIMON_slope_all_TM_res'].mean()
    mean_AAIMON_n_slope_in_dataset = dfg['AAIMON_n_slope_all_TM_res'].mean()
    sys.stdout.write('\n\nmean AAIMON         in dataset: {a:.05f}'
                     '\nmean AAIMON_n       in dataset: {b:.05f}'
                     '\nmean AAIMON_slope   in dataset: {c:.07f}'
                     '\nmean AAIMON_n_slope in dataset: {d:.07f}\n'
                     .format(a=mean_AAIMON_in_dataset,
                             b=mean_AAIMON_n_in_dataset,
                             c=mean_AAIMON_slope_in_dataset,
                             d=mean_AAIMON_n_slope_in_dataset))

    dfg.to_csv(pathdict["list_cr_summary_csv"], sep=",", quoting=csv.QUOTE_NONNUMERIC)

    ########################################################################################
    #                                                                                      #
    #                    Save a huge dataframe with the AAIMONs for                        #
    #                    all homologues of all TMDs of all proteins                        #
    #                                                                                      #
    ########################################################################################

    if s['save_df_characterising_each_homol_TMD']:

        # defining cutoff for max and min number of homologues for each protein
        min_num_homologues = s['min_homol']

        # filter summary file for min and max number of homologues based on TM01 number of homologues
        #sys.stdout.write('Dropped homologues after filtering: \n')
        list_of_acc_to_keep = []
        for acc in dfg.index:
            AAIMON_n_homol = dfg.loc[acc, 'AAIMON_n_homol']
            if AAIMON_n_homol > min_num_homologues:
                list_of_acc_to_keep.append(acc)

        # keep only proteins that have the desired number of homologues
        dfg = dfg.loc[list_of_acc_to_keep, :]
        df = df.loc[list_of_acc_to_keep, :]

        # # convert from string to python list
        # if isinstance(dfg['list_of_TMDs'][0], str):
        #     dfg['list_of_TMDs'] = dfg['list_of_TMDs'].dropna().apply(lambda x: ast.literal_eval(x))

        #sys.stdout.write("\nLoading data\n")
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

            # Here we filter to take only datapoints where all TMDs were in the alignment
            AAIMON_all_TMD = protein_name + '_AAIMON_all_TMD.csv'
            df_AAIMON_all_TMD = utils.open_df_from_csv_zip(homol_cr_ratios_zip, filename=AAIMON_all_TMD, delete_corrupt=False)

            ########################################################################################
            #                                                                                      #
            #       CODE COPIED FROM cons_ratio.py. Delete the following two lines after           #
            #         re-running all calculated cons ratios                                         #
            #                                                                                      #
            ########################################################################################
            # first get a list of all the homologues that have AAIMON ratios for all TMDs
            df_AAIMON_all_TMD["AAIMON_avail_all_TMDs"] = df_AAIMON_all_TMD.n_TMDs_with_measurable_AAIMON == df.loc[acc, "number_of_TMDs"]
            filt_index = df_AAIMON_all_TMD.loc[df_AAIMON_all_TMD["AAIMON_avail_all_TMDs"] == True].index.tolist()
            #filt_index = [int(x) for x in filt_index]

            if not os.path.isfile(homol_cr_ratios_zip):
                # skip to next protein
                continue
            for TMD in df.loc[acc, "list_of_TMDs"]:
                # generate column names necessary for current file
                columns = ['obs_changes', '{}_AAIMON'.format(TMD), '{}_AAIMON_n'.format(TMD)]
                # Open pickle file with conservation-ratios.
                # NOTE that these have already been filtered according to cons_ratio.py.
                # if the homologues were not acceptable, AAIMON ratios WERE NOT CALCULATED
                TM_cr_pickle = "{}_{}_cr_df.pickle".format(protein_name, TMD)
                # open dataframe  with function from korbinian, extract required columns, convert to np array
                df_TMD = utils.open_df_from_pickle_zip(homol_cr_ratios_zip, TM_cr_pickle)
                if columns[2] not in df_TMD.columns:
                    # file is old, and should be deleted
                    #os.remove(homol_cr_ratios_zip)
                    logging.info("{} file is presumed out of date, and has been deleted".format(homol_cr_ratios_zip))
                    os.remove(homol_cr_ratios_zip)
                    # skip to next protein
                    break

                if set(filt_index).intersection(set(df_TMD.index)) == set():
                    # there is a mismatch between the filt_index for df_AAIMON_all_TMD, and the columns in df_TMD
                    # replace filt_index with empty list
                    logging.warning("Indexing Error in gather script. set(filt_index).intersection(set(df_TMD.index)) == set(). Try re-running calculate_AAIMON_ratios")
                    filt_index = []
                # use the filt_index above that shows homologues with AAIMON available for all TMDs
                df_TMD = df_TMD.loc[filt_index, :]
                # convert to numpy array
                df_TMD = df_TMD[columns].as_matrix()
                # join output data file with currently opened dataframe
                data = np.concatenate((data, df_TMD))
        # drop every row with nan
        data = data[~np.isnan(data).any(axis=1)]

        # create bins, calculate mean and confidence interval in bin - use multiprocessing if possible
        sys.stdout.write('\nBinning data - calculating confidence interval\n')
        number_of_bins = s['specify_number_of_bins_characterising_TMDs']
        # process confidence interval value to appropriate input format for function
        confidence_interval = (100 - s['CI']) / 100
        linspace_binlist = np.linspace(1, 100, number_of_bins)
        binned_data = np.empty([0, 8])
        binwidth = 100 / number_of_bins
        for percentage in linspace_binlist:
            bin_for_mean = data[(percentage >= data[:, 0]) & (data[:, 0] > percentage - binwidth)]
            if bin_for_mean.size != 0:
                # calculate conf. interv. in bin, alpha describes the significance level in the style 1-alpha
                conf = sms.DescrStatsW(bin_for_mean[:, 1]).tconfint_mean(alpha=confidence_interval)
                # calculate conf. interv. in bin _n, alpha describes the significance level in the style 1-alpha
                conf_norm = sms.DescrStatsW(bin_for_mean[:, 2]).tconfint_mean(alpha=confidence_interval)
                mean_data_in_bin = np.array([percentage - binwidth / 2,
                                             # calculate mean in bin
                                             bin_for_mean[:, 1].mean(),
                                             # calculate mean in bin _n
                                             bin_for_mean[:, 2].mean(),
                                             # add conf. interv. results to np array
                                             conf[0], conf[1], conf_norm[0], conf_norm[1],
                                             # add the number of TMDs in bin to bin
                                             len(bin_for_mean[:, 0])])
                mean_data_in_bin = mean_data_in_bin.reshape(1, 8)
                sys.stdout.write('.'), sys.stdout.flush()
                binned_data = np.concatenate((mean_data_in_bin.reshape(1, 8), binned_data))


        # # create bins, calculate mean and confidence interval in bin - use multiprocessing if possible
        # sys.stdout.write('\nBinning data - calculating confidence interval\n')
        #
        # use_multiprocessing = s['use_multiprocessing']
        # n_processes = s['multiprocessing_cores']
        # remove_from_binlist = int((1 - s['fa_min_identity_of_full_protein']) * 100)
        # number_of_bins = s['specify_number_of_bins_characterising_TMDs']
        # confidence_interval = (100 - s['CI']) / 100
        # linspace_binlist = np.linspace(1, 100, number_of_bins)[:remove_from_binlist]
        # binned_data = np.empty([0, 8])
        # binwidth = 100 / number_of_bins
        # list_p = []
        # for percentage in linspace_binlist:
        #     data_as_dict = {'data': data, 'percentage': percentage, 'binwidth': binwidth, 'confidence_interval': confidence_interval}
        #     list_p.append(data_as_dict)
        #
        # if use_multiprocessing:
        #     with Pool(processes=n_processes) as pool:
        #             mean_data_in_bin = pool.map(binning_data_multiprocessing, list_p)
        # else:
        #     mean_data_in_bin = []
        #     for p in list_p:
        #         output = binning_data_multiprocessing(p)
        #         if type(output) is np.ndarray:
        #             mean_data_in_bin.append(output)
        #
        # for n, element in enumerate(mean_data_in_bin):
        #     if type(mean_data_in_bin[n]) is np.ndarray:
        #         binned_data = np.concatenate((mean_data_in_bin[n].reshape(1, 8), binned_data))


        # create bins, calculate mean and 95% confidence interval
        # sys.stdout.write('\nBinning data - calculating confidence interval\n')
        # confidence_interval = (100 - s['CI'])/100
        # number_of_bins = s['specify_number_of_bins_characterising_TMDs']
        # linspace_binlist = np.linspace(1, 100, number_of_bins)
        # binwidth = 100/number_of_bins
        # binned_data = np.empty([0, 8])
        # # conf_95 = np.array([1, 2])
        # # conf95_norm = np.array([1, 2])
        # for percentage in linspace_binlist:
        #     if percentage % 5 == 0:
        #         sys.stdout.write('{}%, '.format(int(percentage))), sys.stdout.flush()
        #     bin_for_mean = np.empty([0, 3])
        #     for row in data:
        #         if row[0] < percentage and row[0] > percentage - binwidth:
        #             bin_for_mean = np.concatenate((bin_for_mean, row.reshape(1, 3)))
        #     if bin_for_mean.size != 0:
        #         # calculate conf. interv. in bin, alpha describes the significance level in the style 1-alpha
        #         conf = sms.DescrStatsW(bin_for_mean[:, 1]).tconfint_mean(alpha=confidence_interval)
        #         # calculate conf. interv. in bin _n, alpha describes the significance level in the style 1-alpha
        #         conf_norm = sms.DescrStatsW(bin_for_mean[:, 2]).tconfint_mean(alpha=confidence_interval)
        #         mean_data_in_bin = np.array([percentage - binwidth/2,
        #                                      # calculate mean in bin
        #                                      bin_for_mean[:, 1].mean(),
        #                                      # calculate mean in bin _n
        #                                      bin_for_mean[:, 2].mean(),
        #                                      # add conf. interv. results to np array
        #                                      conf[0], conf[1], conf_norm[0], conf_norm[1],
        #                                      # add the number of TMDs in bin to bin
        #                                      len(bin_for_mean[:, 0])])
        #         # merge data from bin to the others
        #         binned_data = np.concatenate((mean_data_in_bin.reshape(1, 8), binned_data))
        # # drop every row containing nan in array
        # binned_data = binned_data[~np.isnan(binned_data).any(axis=1)]
        '''
        description of columns in numpy arrays:

        numpy array data:
        |       0       |   1    |    2     |
        | % obs_changes | AAIMON | AAIMON_n |

        numpy array binned_Data:
        |       0       |      1      |       2       |     3    |    4    |      5     |     6     |         7          |
        | % obs_changes | mean AAIMON | mean AAIMON_n |  CI_low  |  CI_hi  |  CI_low_n  |  CI_hi_n  | number of Proteins |
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

    logging.info("\n~~~~~~~~~~~~                           finished gather_AAIMONs                      ~~~~~~~~~~~~")
    return pathdict

# ########################################################################################
# #                                                                                      #
# #                function for binning data using multiprocessing                       #
# #                                                                                      #
# ########################################################################################
#
# def binning_data_multiprocessing(list_p):
#     data = list_p['data']
#     percentage = list_p['percentage']
#     binwidth = list_p['binwidth']
#     confidence_interval = list_p['confidence_interval']
#     bin_for_mean = np.empty([0, 3])
#     mean_data_in_bin = np.empty([0, 8])
#     sys.stdout.write('.'), sys.stdout.flush()
#     for row in data:
#         if row[0] < percentage and row[0] > percentage - binwidth:
#             bin_for_mean = np.concatenate((bin_for_mean, row.reshape(1, 3)))
#     if bin_for_mean.size != 0:
#         # calculate conf. interv. in bin, alpha describes the significance level in the style 1-alpha
#         conf = sms.DescrStatsW(bin_for_mean[:, 1]).tconfint_mean(alpha=confidence_interval)
#         # calculate conf. interv. in bin _n, alpha describes the significance level in the style 1-alpha
#         conf_norm = sms.DescrStatsW(bin_for_mean[:, 2]).tconfint_mean(alpha=confidence_interval)
#         mean_data_in_bin = np.array([percentage - binwidth/2,
#                                      # calculate mean in bin
#                                      bin_for_mean[:, 1].mean(),
#                                      # calculate mean in bin _n
#                                      bin_for_mean[:, 2].mean(),
#                                      # add conf. interv. results to np array
#                                      conf[0], conf[1], conf_norm[0], conf_norm[1],
#                                      # add the number of TMDs in bin to bin
#                                      len(bin_for_mean[:, 0])])
#         mean_data_in_bin = mean_data_in_bin.reshape(1, 8)
#     if mean_data_in_bin.any():
#         return mean_data_in_bin



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
    # convert from string to python list
    df['list_of_TMDs'] = df['list_of_TMDs'].apply(lambda x: ast.literal_eval(x))

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

                # open parsed data
                #homol_parsed_zip = df.loc[acc, "homol_parsed"]
                parsed_pickle = "{}_df_orig.pickle".format(protein_name)
                dfh = utils.open_df_from_pickle_zip(homol_df_orig_zip, parsed_pickle)

                # open collected pretty alignments
                SIMAP_align_pretty_csv_filename = os.path.basename(df.loc[acc, "SIMAP_align_pretty_csv"])
                dfp = utils.open_df_from_csv_zip(homol_df_orig_zip, SIMAP_align_pretty_csv_filename)

                # open dataframe with the nonTMD calculations
                nonTMD_cr_pickle = "{}_nonTMD_cr_df.pickle".format(protein_name)
                df_nonTMD = utils.open_df_from_pickle_zip(homol_cr_ratios_zip, nonTMD_cr_pickle)


                for TMD in df.loc[acc, "list_of_TMDs"]:
                    TM_cr_pickle = "{}_{}_cr_df.pickle".format(protein_name, TMD)
                    # open dataframe  with function from korbinian, extract required columns, convert to np array
                    df_TMD = utils.open_df_from_pickle_zip(homol_cr_ratios_zip, TM_cr_pickle)
                    # if the pickle doesn't seem to hold the correct file
                    if df_TMD.empty:
                        # skip to next TMD
                        continue

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
                        TMD_cols = ['obs_changes', "norm_factor", "{}_AAIMON", '{}_perc_ident', '{}_start_in_SW_alignment', '{}_SW_query_seq', '{}_SW_markup_seq',
                                   '{}_SW_match_seq', '{}_ratio_len_TMD_to_len_nonTMD', '{}_SW_align_len', "{}_SW_match_lipo",
                                    '{}_SW_query_num_gaps', '{}_SW_match_num_gaps','{}_SW_align_len_excl_gaps'] # 'FASTA_gapped_identity',
                        TMD_col_names = ['obs_changes', "norm_factor", "AAIMON", 'TM_perc_ident', 'TM_start_in_SW_alignment', 'SW_query_seq', 'SW_markup_seq', 'SW_match_seq',
                                     'ratio_len_TMD_to_len_nonTMD', 'SW_align_len', "SW_match_lipo",
                                     'SW_query_num_gaps', 'SW_match_num_gaps', 'SW_align_len_excl_gaps'] # 'FASTA_gapped_identity',

                        ########################################################################################
                        #                                                                                      #
                        #         nonTMD and dfh stuff is very inefficient! cycles through for each TM         #
                        #         also cycles through each outlier!!
                        #                                                                                      #
                        ########################################################################################

                        # list of columns from which to obtain data
                        nonTMD_cols = ['perc_nonTMD_coverage', 'nonTMD_perc_ident', 'nonTMD_SW_align_len_excl_gaps']#, "nonTMD_seq_query"
                        nonTMD_col_names = nonTMD_cols

                        dfh_cols = ['SW_identity', 'SW_coverage_ratio', 'FASTA_identity', 'match_align_seq', 'query_align_seq', 'align_markup_seq']
                        dfh_col_names = dfh_cols

                        TMD_tuple = (TMD_cols, TMD_col_names, df_TMD)
                        nonTMD_tuple = (nonTMD_cols, nonTMD_col_names, df_nonTMD)
                        dfh_tuple = (dfh_cols, dfh_col_names, dfh)
                        tuples_to_process = (TMD_tuple, nonTMD_tuple, dfh_tuple)

                        # iterate through the TMD and nonTMD dataframes, extracting the relevant information
                        for tup in tuples_to_process:
                            columns = tup[0]
                            col_names = tup[1]
                            dfx = tup[2]
                            for n, col in enumerate(columns):
                                col_name = col_names[n]
                                if col.format(TMD) not in dfx.columns:
                                    # skip
                                    logging.info("{} skipped, {} not in columns".format(acc, col))
                                    continue
                                value = dfx.loc[outlier_index, col.format(TMD)]
                                # add each one to the dictionary
                                d[col_name] = value

                        # add the pretty alignment, which is extracted from the csv in the homol_df_orig_zip, and not from the pickle file with the other calculated variables
                        d["align_pretty"] = dfp.loc[outlier_index, "align_pretty"]

                        # create a "TMD alignment" by joining the query, markup and match sequences together.
                        # leave the individual columns in the dataframe, as they could be useful for other analyses (e.g. aa abundance tests)
                        d["TM_align"] = "{}\r\r\n\r\r{}\r\r\n\r\r{}".format(d['SW_query_seq'], d['SW_markup_seq'], d['SW_match_seq'])

                        # calculate an alternative observed_changes
                        d["fl_aln_len"] = len(d['match_align_seq'])
                        d["fl_gaps_qm"] = d['match_align_seq'].count("-") + d['query_align_seq'].count("-")
                        d["fl_ident"] = len(re.findall("\|", d['align_markup_seq']))
                        d["fl_obs_changes"] = 1 - d["fl_ident"] / (d["fl_aln_len"] - d["fl_gaps_qm"])

                        # remove the large sequences from the dictionary
                        for key in ['match_align_seq', 'query_align_seq', 'align_markup_seq']:
                            del d[key]

                        if num_TMDs_in_all_proteins_processed == 0:
                            # sort
                            csv_header = ["protein_name", "TMD", "outlier", "TM_align","SW_match_lipo", "align_pretty", 'obs_changes', "AAIMON", "norm_factor", 'perc_nonTMD_coverage', "hit",
                                          'TM_perc_ident', 'nonTMD_perc_ident', 'TM_start_in_SW_alignment', 'SW_query_seq', 'SW_markup_seq', 'SW_match_seq',
                                          'ratio_len_TMD_to_len_nonTMD', 'SW_align_len', 'SW_query_num_gaps', 'SW_match_num_gaps', 'SW_align_len_excl_gaps','nonTMD_SW_align_len_excl_gaps',
                                          'SW_identity', 'SW_coverage_ratio', 'FASTA_identity',
                                          "fl_aln_len", "fl_gaps_qm", "fl_ident", "fl_obs_changes"] # 'FASTA_gapped_identity','match_align_seq', 'query_align_seq', 'align_markup_seq', "nonTMD_seq_query"

                            # make sure that the csv header is up-to-date, and isn't missing items from dict
                            if len(csv_header) != len(d):
                                raise ValueError("Columns in CSV header and dictionary don't match.\nSuggest double-checking added columns.")
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