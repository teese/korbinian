import ast
import os
import pickle
import zipfile
import korbinian
import korbinian.cons_ratio.calc
import korbinian.cons_ratio.norm
import korbinian.utils as utils
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from multiprocessing import Pool
import sys
from scipy.optimize import leastsq

def run_calculate_AAIMON_ratios(pathdict, s, logging):
    """Runs calculate_AAIMON_ratios for each protein, using multiprocessing Pool.

    Parameters
    ----------
    pathdict : dict
        Dictionary of the key paths and files associated with that List number.
    s : dict
        Settings dictionary extracted from excel settings file.
    logging : logging.Logger
        Logger for printing to console and/or logfile.
        If multiprocessing == True, logging.info etc will only print to console.
    """
    logging.info('~~~~~~~~~~~~      starting run_calculate_AAIMON_ratios        ~~~~~~~~~~~~')
    # if multiprocessing is used, log only to the console
    p_dict_logging = logging if s["use_multiprocessing"] != True else utils.Log_Only_To_Console()
    # set current working directory as the data_dir/homol, where temp files will be saved before moving to zip
    os.chdir(os.path.join(s["data_dir"], "homol"))
    # get list of accessions that could not be downloaded, and can immediately be excluded
    not_in_homol_db = utils.get_list_not_in_homol_db(pathdict)
    # create list of protein dictionaries to process
    list_p = korbinian.utils.convert_summary_csv_to_input_list(s, pathdict, p_dict_logging, list_excluded_acc=not_in_homol_db)
    # number of processes is the number the settings, or the number of proteins, whichever is smallest
    n_processes = s["multiprocessing_cores"] if s["multiprocessing_cores"] < len(list_p) else len(list_p)

    if s["use_multiprocessing"]:
        with Pool(processes=n_processes) as pool:
            calc_AAIMON_list = pool.map(korbinian.cons_ratio.cons_ratio.calculate_AAIMON_ratios, list_p)
            # log the list of protein results (e.g. acc, "simap", True) to the actual logfile, not just the console
            logging.info("calc_AAIMON_list : {}".format(calc_AAIMON_list))
    else:
        for p in list_p:
            korbinian.cons_ratio.cons_ratio.calculate_AAIMON_ratios(p)
    logging.info("~~~~~~~~~~~~     run_calculate_AAIMON_ratios is finished      ~~~~~~~~~~~~")

def calculate_AAIMON_ratios(p):
    """Calculate the AAIMON ratios for a particular protein

    Parameters
    ----------
    p : dict
        Protein Dictionary. Contains all input settings, sequences and filepaths related to a single protein.
        Protein-specific data is extracted from one row of the the list summary, e.g. List05_summary.csv, which is read as df.
        p also contains the GENERAL korbinian settings and filepaths for that list (pathdict, s, logging)

        Components of p :
            pathdict : dict
                Dictionary of the key paths and files associated with that List number.
            s : dict
                Settings dictionary extracted from excel settings file.
            logging : logging.Logger
                Logger for printing to console and/or logfile.
                If multiprocessing == True, logging.info etc will only print to console.
            p : protein-specific dictionary components
                acc, list_of_TMDs, description, TM01_seq, etc

    Saved Files and Figures
    -----------------------
    homol_cr_ratios_zip(zipfile)
        PROTEIN_NAME_cr_ratios.zip (E.g. A6BM72_MEG11_HUMAN_cr_ratios.zip)
            A6BM72_MEG11_HUMAN_AAIMON_hist_0.png : png
                Histograms of AAIMON ratios for homologues of each TMD.
            A6BM72_MEG11_HUMAN_cr_mean.csv : csv
                Summary file for that protein. Contains conservation ratios means.
                Will be gathered for all proteins by the gather_AAIMON_ratios function.
            A6BM72_MEG11_HUMAN_nonTMD_cr_df.pickle : pickled pd.DataFrame
                Dataframe containing the percentage_identity etc for sliced nonTMD region.
            A6BM72_MEG11_HUMAN_SP01_cr_df.pickle : pickled pd.DataFrame
                Dataframe containing the percentage_identity etc for that particular TMD/region (in this case, the signal peptide).
            A6BM72_MEG11_HUMAN_TM01_cr_df.pickle
                Dataframe containing the percentage_identity etc for that particular TMD/region (in this case, TM01).

    Returns
    -------
    In all cases, a tuple (str, bool, str) is returned.

    if successful:
        return acc, True, "0"
    if not successful:
        return acc, False, "specific warning or reason why protein failed"
    """
    pathdict, s, logging = p["pathdict"], p["s"], p["logging"]
    acc = p["acc"]
    protein_name = p["protein_name"]
    if not os.path.exists(p['homol_df_orig_zip']):
        message = "{} Protein skipped. File does not exist".format(p['homol_df_orig_zip'])
        logging.info(message)
        return acc, False, message

    dfh = utils.open_df_from_pickle_zip(p['homol_df_orig_zip'],filename=os.path.basename(p['homol_df_orig_pickle']), delete_corrupt=True)
    if dfh.empty:
        message = "{} Protein skipped, file deleted as it is possibly corrupt.".format(p['homol_df_orig_zip'])
        logging.info(message)
        return acc, False, message

    # create an output series for that protein, containing mean AAIMON values, etc.
    mean_ser = pd.Series()
    mean_ser["acc"] = acc
    mean_ser["protein_name"] = p['protein_name']
    mean_ser["organism"] = p['organism']
    mean_ser["prot_descr"] = p['prot_descr']
    mean_ser['list_of_TMDs'] = p['list_of_TMDs']
    list_of_TMDs = ast.literal_eval(p['list_of_TMDs'])

    ########################################################################################
    #                                                                                      #
    #                        Calculate AAIMON conservation ratios                          #
    #                                                                                      #
    ########################################################################################

    homol_cr_ratios_zip = p['homol_cr_ratios_zip']
    mean_ser_filename = "{}_cr_mean.csv".format(protein_name)

    #assume af first that there is no previous data, and that the calculations can be re-run
    if s["overwrite_prev_calculated_AAIMON_ratios"] == False:
        if os.path.isfile(homol_cr_ratios_zip):
            with zipfile.ZipFile(homol_cr_ratios_zip, mode="r", compression=zipfile.ZIP_DEFLATED) as cr_zip:
                if mean_ser_filename in cr_zip.namelist():
                    # if the means are saved as a csv in the homol cr_ratios zipfile, skip this protein
                    message = '{} AAIMON_ratios skipped, file with mean AAIMON ratios exists (in settings, overwrite_prev_calculated_AAIMON_ratios = True)'.format(acc)
                    logging.info(message)
                    return acc, False, message

    ########################################################################################
    #                                                                                      #
    #           Filter based on homol hit properties (non-TMD-specific)                    #
    #                                                                                      #
    ########################################################################################
    cr_homol_query_str = 'FASTA_gapped_identity > {min_ident} & ' \
                        'FASTA_gapped_identity < {max_ident} & ' \
                        'hit_contains_SW_node == True & ' \
                        'disallowed_words_not_in_descr == True &' \
                        'X_in_match_seq == False'.format(min_ident=s["cr_min_identity_of_full_protein"],
                                                         max_ident=s["cr_max_identity_of_full_protein"])

    # filter based on the query string
    dfh.query(cr_homol_query_str, inplace=True)

    #Calculate fasta values, add to mean_ser output series for each protein.
    # The mean fasta identity of homologues should almost always be around the centre of the min and max in the settings.
    mean_ser['FASTA_ident_mean'] = float('%0.2f' % dfh['FASTA_identity'].mean())
    # number of identical residues in each FASTA alignment can be calculated from identity and overlap
    dfh['FASTA_num_ident_res'] = dfh['FASTA_identity'] * dfh['FASTA_overlap']
    mean_ser['FASTA_num_ident_res_mean'] = float('%0.2f' % dfh['FASTA_num_ident_res'].mean())

    if not os.path.exists(p['fa_cr_sliced_TMDs_zip']):
        message = "{} Protein skipped. File does not exist".format(p['fa_cr_sliced_TMDs_zip'])
        logging.info(message)
        return acc, False, message

    nonTMD_pickle_name = "{}_nonTMD_sliced_df.pickle".format(protein_name)
    df_nonTMD = utils.open_df_from_pickle_zip(p['fa_cr_sliced_TMDs_zip'], filename=nonTMD_pickle_name, delete_corrupt=True)
    if df_nonTMD.empty:
        message = "{} file was corrupt and was deleted".format(p['fa_cr_sliced_TMDs_zip'])
        logging.info(message)
        return acc, False, message

    # filter based on dfh above, for general homologue settings (e.g. % identity of full protein)
    try:
        df_nonTMD = df_nonTMD.loc[dfh.index, :]
    except KeyError:
        # in very rare cases, none of the dfh.index is actually found in df_nonTMD, and therefore the protein should be skipped
        # since the df_nonTMD depends on ALL TMDs being found, this occurs when none of the homologues contain all TMDs
        message = "{} none of the dfh.index is actually found in df_nonTMD".format(acc)
        logging.info(message)
        return acc, False, message

    ########################################################################################
    #                                                                                      #
    #                 Calculate the nonTMD percentage identity and gaps                    #
    #                                                                                      #
    ########################################################################################
    # extract the length of the nonTMD region in the original query here, which should be in the original list summary file, created during uniprot_parse or OMPdb_get_TM_indices_and_slice
    len_nonTMD_orig_q = p['len_nonTMD']

    df_nonTMD, mean_ser = korbinian.cons_ratio.calc.calc_nonTMD_perc_ident_and_gaps(df_nonTMD, mean_ser, len_nonTMD_orig_q)

    ########################################################################################
    #                                                                                      #
    #                Calculation of normalization factor for each homologue                #
    #                                                                                      #
    ########################################################################################
    dfh['norm_factor'] = dfh['FASTA_gapped_identity'].apply(korbinian.cons_ratio.norm.calc_AAIMON_aa_prop_norm_factor, args=(0.136, 0.055)) # random factors according to Oberai 2009

    with zipfile.ZipFile(homol_cr_ratios_zip, mode="w", compression=zipfile.ZIP_DEFLATED) as zipout:

        # save the nonTMD dataframe
        nonTMD_cr_outfile_pickle = "{}_nonTMD_cr_df.pickle".format(protein_name)
        with open(nonTMD_cr_outfile_pickle, "wb") as f:
            pickle.dump(df_nonTMD, f, protocol=pickle.HIGHEST_PROTOCOL)
        zipout.write(nonTMD_cr_outfile_pickle, arcname=nonTMD_cr_outfile_pickle)
        os.remove(nonTMD_cr_outfile_pickle)

        # filter nonTMD dataframe to only contain entries where nonTMD_perc_ident is not zero
        # filter to remove short nonTMD regions
        # note this filtering is AFTER the full dataframe has been saved to file, preventing loss of data
        nonTMD_query_str = "nonTMD_perc_ident != 0 & " \
                           "nonTMD_len >= {min_nonTMD_len}".format(min_nonTMD_len=s["cr_min_len_nonTMD"])

        df_nonTMD.query(nonTMD_query_str, inplace=True)

        linspace_binlist = np.linspace(s["1p_smallest_bin"],
                                       s["1p_largest_bin"],
                                       s["1p_number_of_bins"])
        # add 30 as the last bin, to make sure 100% of the data is added to the histogram, including major outliers
        binarray = np.append(linspace_binlist, s["1p_final_highest_bin"])

        # se default font size for text in the plot
        fontsize = 4
        # use a dictionary to organise the saving of multiple plots in multiple figures, with a certain number of plots per figure
        n_plots_per_fig = 4
        nrows_in_each_fig = 2
        ncols_in_each_fig = 2
        dict_organising_subplots = utils.create_dict_organising_subplots(
            n_plots_per_fig=n_plots_per_fig,
            n_rows=nrows_in_each_fig,
            n_cols=ncols_in_each_fig)
        # make the IDE happy
        fig, axarr = None, None

        list_of_AAIMON_all_TMD = {}

        for TMD_Nr, TMD in enumerate(list_of_TMDs):
            # find the TMD number (starting from 1)
            TMD_Nr = list_of_TMDs.index(TMD) + 1
            ########################################################################################
            #                                                                                      #
            #                    Add nonTMD info to df_cr for each TMD.                            #
            #                                                                                      #
            ########################################################################################
            # open the dataframe containing the sequences, gap counts, etc for that TMD only
            df_cr = utils.open_df_from_pickle_zip(p['fa_cr_sliced_TMDs_zip'], filename="{}_{}_sliced_df.pickle".format(protein_name, TMD), delete_corrupt=True)
            # add the nonTMD percentage identity, etc. NOTE THAT THE INDEX SHOULD STILL BE COMPATIBLE, as the hit_num!
            df_cr['nonTMD_perc_ident'] = df_nonTMD['nonTMD_perc_ident']
            df_cr['nonTMD_perc_sim_plus_ident'] = df_nonTMD['nonTMD_perc_sim_plus_ident']
            df_cr['nonTMD_SW_align_len_excl_gaps'] = df_nonTMD['nonTMD_SW_align_len_excl_gaps']
            df_cr['len_full_match_seq'] = dfh['len_full_match_seq']
            df_cr['truncation_ratio_nonTMD'] = df_nonTMD['truncation_ratio_nonTMD']
            df_cr['len_nonTMD_orig_q_minus_nonTMD_SW_align_len_excl_gaps'] = df_nonTMD['len_nonTMD_orig_q_minus_nonTMD_SW_align_len_excl_gaps']
            # filter based on dfh above, for general homologue settings (e.g. % identity of full protein), and df_nonTMD (for nonTMD_perc_ident is not zero, etc)
            df_cr = df_cr.loc[df_nonTMD.index,:]
            """FILTERING BY TMD IDENTITY HERE SHOULD NO LONGER BE NECESSARY. DIVIDE BY 0 AVOIDED WITH THE REPLACE FUNCTION LATER"""
            ## following the general filters, filter to only analyse sequences with TMD identity above cutoff, and a nonTMD_perc_ident above zero ,to avoid a divide by zero error
            #df_cr = df_cr.loc[df_cr['%s_perc_ident' % TMD] >= s['cr_min_identity_of_TMD_initial_filter']]

            ########################################################################################
            #                                                                                      #
            #                       Calculate AAIMON, AASMON for each TMD                          #
            #                                                                                      #
            ########################################################################################
            len_query_TMD = p["%s_end"%TMD] - p["%s_start"%TMD]
            df_cr = korbinian.cons_ratio.calc.calc_AAIMON(TMD, df_cr, len_query_TMD)
            df_cr['norm_factor'] = dfh['norm_factor']
            df_cr['%s_AAIMON_ratio_n'%TMD] = df_cr['%s_AAIMON_ratio'%TMD] / df_cr['norm_factor']
            df_cr['FASTA_gapped_identity'] = dfh['FASTA_gapped_identity']

            list_of_AAIMON_all_TMD['%s_AAIMON_ratio'%TMD]= df_cr['%s_AAIMON_ratio'%TMD].dropna()

            ########################################################################################
            #                                                                                      #
            #                       Calculate AAIMON_slope, AAIMON_n_slope                         #
            #                                                                                      #
            ########################################################################################
            # drop every row (hit) in df_cr that contains NaN in column TMxy_AAIMON_ratio - important for line fit that can't handle NAN
            df_cr = df_cr[np.isfinite(df_cr['%s_AAIMON_ratio'%TMD])]

            # set up variables
            FASTA_gapped_identity = df_cr['FASTA_gapped_identity']  # x-axis of plot
            AAIMON_ratios = df_cr['%s_AAIMON_ratio'%TMD]            # y-axis
            AAIMON_ratios_n = df_cr['%s_AAIMON_ratio_n'%TMD]        # y-axis

            if len(FASTA_gapped_identity) == 0 or len(AAIMON_ratios) == 0:
                # There is no gapped identity for these homologues, skip to next TMD
                continue
            # linear regression for non-norm. and norm. AAIMON with fixed 100% identity at AAIMON 1.0
            AAIMON_slope, x_data, y_data = curve_fitting_fixed_100(FASTA_gapped_identity, AAIMON_ratios)
            mean_ser['%s_AAIMON_slope' % TMD] = AAIMON_slope
            AAIMON_n_slope, x_data_n, y_data_n = curve_fitting_fixed_100(FASTA_gapped_identity, AAIMON_ratios_n)
            mean_ser['%s_AAIMON_n_slope' % TMD] = AAIMON_n_slope



            # # linear regression for non-normalised AAIMON
            # linear_regression_AAIMON = np.polyfit(FASTA_gapped_identity, AAIMON_ratios, 1)
            # fit_fn_AAIMON = np.poly1d(linear_regression_AAIMON)
            # fitted_data_AAIMON = fit_fn_AAIMON(FASTA_gapped_identity)
            # mean_ser['%s_AAIMON_slope' %TMD] = linear_regression_AAIMON[0]
            #
            # # linear regression for normalised AAIMON
            # linear_regression_AAIMON_n = np.polyfit(FASTA_gapped_identity, AAIMON_ratios_n, 1)
            # fit_fn_AAIMON_n = np.poly1d(linear_regression_AAIMON_n)
            # fitted_data_AAIMON_n = fit_fn_AAIMON_n(FASTA_gapped_identity)
            # mean_ser['%s_AAIMON_n_slope' %TMD] = linear_regression_AAIMON_n[0]

            if '{TMD}_SW_match_seq_hydro'.format(TMD=TMD) not in df_cr.columns:
                message = "{} {}_SW_match_seq_hydro not found in columns. Slice file is out of date and will be deleted.".format(acc, TMD)
                os.remove(p['fa_cr_sliced_TMDs_zip'])
                return acc, False, message

            ###########################################################################
            #                                                                         #
            #       important if you want to filter for truncated alignments          #
            #                             DO NOT DELETE!!!                            #
            #                                                                         #
            ###########################################################################
            """
            # save the unfiltered dataframe for that TMD
            TM_cr_outfile_pickle = "{}_{}_cr_df_RAW.pickle".format(protein_name, TMD)
            with open(TM_cr_outfile_pickle, "wb") as f:
                pickle.dump(df_cr, f, protocol=pickle.HIGHEST_PROTOCOL)
            zipout.write(TM_cr_outfile_pickle, arcname=TM_cr_outfile_pickle)
            os.remove(TM_cr_outfile_pickle)
            """

            TM_cr_outfile_pickle = "{}_{}_cr_df.pickle".format(protein_name, TMD)
            with open(TM_cr_outfile_pickle, "wb") as f:
                pickle.dump(df_cr, f, protocol=pickle.HIGHEST_PROTOCOL)
            zipout.write(TM_cr_outfile_pickle, arcname=TM_cr_outfile_pickle)
            os.remove(TM_cr_outfile_pickle)


            max_gaps  = s["cr_max_n_gaps_in_TMD"]
            max_hydro = s["cr_max_hydrophilicity_Hessa"]
            min_ident = s["cr_min_identity_of_TMD"]
            # filter by TMD-specific values (e.g. max_gaps_in_TMD and then calculate all the mean values for AAIMON, etc)
            # note that this is done AFTER the full df_cr is saved, so df_cr can be filtered and reduced directly without losing data
            mean_ser = korbinian.cons_ratio.calc.filt_and_save_AAIMON_mean(TMD, df_cr, mean_ser, max_gaps, max_hydro, min_ident)

            if TMD == "TM01":
                # number of homologues for TM01. since ALL TMDs have to be in each homologue before AAIMON is calculated, this number is the same for all TMDs
                mean_ser['TM01_AAIMON_n_homol'] = df_cr['TM01_AAIMON_ratio'].dropna().shape[0]

            logging.info('%s AAIMON_mean %s: %0.2f' % (acc, TMD, mean_ser['%s_AAIMON_ratio_mean' % TMD]))
            logging.info('%s AAIMON_n_mean %s: %0.2f' % (acc, TMD, mean_ser['%s_AAIMON_ratio_mean_n' % TMD]))
            logging.info('%s AAIMON_slope %s: %0.5f' % (acc, TMD, mean_ser['%s_AAIMON_slope' % TMD]))
            logging.info('%s AAIMON_n_slope %s: %0.5f' % (acc, TMD, mean_ser['%s_AAIMON_n_slope' % TMD]))
            # logging.info('%s AASMON MEAN %s: %0.2f' % (acc, TMD, mean_ser['%s_AASMON_ratio_mean'%TMD]))

            # use the dictionary to obtain the figure number, plot number in figure, plot indices, etc
            newfig, savefig, fig_nr, plot_nr_in_fig, row_nr, col_nr = dict_organising_subplots[TMD_Nr]
            # if the TMD is the last one, the figure should be saved
            if TMD_Nr == len(list_of_TMDs):
                savefig = True
            # if a new figure should be created (either because the orig is full, or the last TMD is analysed)
            if newfig:
                # create a new figure
                fig, axarr = plt.subplots(nrows=nrows_in_each_fig, ncols=ncols_in_each_fig)  # sharex=True

            #" NOT STABLE! NEED TO CHANGE save_hist_AAIMON_ratio_single_protein SO THAT IT RUNS WITHIN THE FOR LOOP ABOVE, AND TAKES A SINGLE TMD AS INPUT, RATHER THAN LIST OF TMDS" / 4
            AAIMON_hist_path_prefix = p['AAIMON_hist_path_prefix']
            ########################################################################################
            #                                                                                      #
            #       Save histograms for each TMD of that protein, with relative conservation       #
            #                                                                                      #
            ########################################################################################
            if axarr is None:
                # avoid bug where the axarr is still not created, as newfig was not True for first figure?
                continue
            korbinian.cons_ratio.histogram.save_hist_AAIMON_ratio_single_protein(fig_nr, fig, axarr, df_cr, s, TMD, binarray, zipout, row_nr, col_nr, fontsize, savefig, AAIMON_hist_path_prefix)

        ########################################################################################
        #                                                                                      #
        #               AAIMON normalization and save fig for each protein                     #
        #                                                                                      #
        ########################################################################################
        df_AAIMON_all_TMD = pd.DataFrame(list_of_AAIMON_all_TMD)
        df_AAIMON_all_TMD['AAIMON_ratio_mean_all_TMDs_1_homol'] = df_AAIMON_all_TMD.mean(axis=1)
        df_AAIMON_all_TMD['gapped_ident'] = dfh['FASTA_gapped_identity'].loc[df_AAIMON_all_TMD.index]
        df_AAIMON_all_TMD['norm_factor'] = dfh['norm_factor'].loc[df_AAIMON_all_TMD.index]
        df_AAIMON_all_TMD['AAIMON_ratio_mean_all_TMDs_1_homol_n'] = df_AAIMON_all_TMD['AAIMON_ratio_mean_all_TMDs_1_homol'] / df_AAIMON_all_TMD['norm_factor']
#        print(df_AAIMON_all_TMD)
        korbinian.cons_ratio.norm.save_graph_for_normalized_AAIMON(acc,  df_AAIMON_all_TMD['AAIMON_ratio_mean_all_TMDs_1_homol'],
                                                                     df_AAIMON_all_TMD['AAIMON_ratio_mean_all_TMDs_1_homol_n'],
                                                                     df_AAIMON_all_TMD['gapped_ident'], zipout, protein_name)
        # save the dataframe containing normalisation factor and normalised AAIMON to zipout
        df_AAIMON_all_TMD.to_csv(protein_name + '_AAIMON_all_TMD.csv')
        zipout.write(protein_name + '_AAIMON_all_TMD.csv', arcname=protein_name + '_AAIMON_all_TMD.csv')
        os.remove(protein_name + '_AAIMON_all_TMD.csv')

        value_counts_hit_contains_SW_node = dfh['hit_contains_SW_node'].value_counts()
        if True in value_counts_hit_contains_SW_node:
            mean_ser['num_hits_with_SW_align_node'] = value_counts_hit_contains_SW_node[True]
        else:
            logging.warning("{} num_hits_with_SW_align_node = 0".format(protein_name))
            mean_ser['num_hits_with_SW_align_node'] = 0
        # save the pandas series with the means to a csv in the cr_ratios zip file
        mean_ser.to_csv(mean_ser_filename)
        zipout.write(mean_ser_filename, arcname=mean_ser_filename)
        os.remove(mean_ser_filename)
        return acc, True, "0"


def throw_out_truncated_sequences(pathdict, s, logging):

    logging.info('~~~~~~~~~~~~      starting run_filter_truncated_alignments        ~~~~~~~~~~~~')
    # if multiprocessing is used, log only to the console
    p_dict_logging = logging if s["use_multiprocessing"] != True else utils.Log_Only_To_Console()
    # set current working directory as the data_dir/homol, where temp files will be saved before moving to zip
    os.chdir(os.path.join(s["data_dir"], "homol"))
    # get list of accessions that could not be downloaded, and can immediately be excluded
    not_in_homol_db = utils.get_list_not_in_homol_db(pathdict)
    # create list of protein dictionaries to process
    list_p = korbinian.utils.convert_summary_csv_to_input_list(s, pathdict, p_dict_logging, list_excluded_acc=not_in_homol_db)
    # number of processes is the number the settings, or the number of proteins, whichever is smallest
    n_processes = s["multiprocessing_cores"] if s["multiprocessing_cores"] < len(list_p) else len(list_p)

    if s["use_multiprocessing"]:
        with Pool(processes=n_processes) as pool:
            calc_AAIMON_list = pool.map(korbinian.cons_ratio.cons_ratio.truncation_filter, list_p)
            # log the list of protein results (e.g. acc, "simap", True) to the actual logfile, not just the console
            logging.info("calc_AAIMON_list : {}".format(calc_AAIMON_list))
    else:
        for p in list_p:
            korbinian.cons_ratio.cons_ratio.truncation_filter(p)
    logging.info("~~~~~~~~~~~~     run_filter_truncated_alignments is finished      ~~~~~~~~~~~~")


def truncation_filter(p):
    pathdict, s, logging = p["pathdict"], p["s"], p["logging"]
    acc = p["acc"]
    protein_name = p["protein_name"]
    uniprot_acc = p['uniprot_acc']
    homol_cr_ratios_zip = p['homol_cr_ratios_zip']
    list_of_TMDs = ast.literal_eval(p['list_of_TMDs'])
    nonTMD_truncation_cutoff = s['truncation_cutoff']

    if not os.path.exists(homol_cr_ratios_zip):
        message = "{} Protein skipped. File does not exist".format(homol_cr_ratios_zip)
        logging.info(message)
        return acc, False, message

    # read data from disk
    in_zipfile = homol_cr_ratios_zip
    # open every single original TMD dataframe in zip
    for TMD in list_of_TMDs:
        in_file = "{}_{}_cr_df_RAW.pickle".format(protein_name, TMD)
        # with zipfile.ZipFile(in_zipfile, "r", zipfile.ZIP_DEFLATED) as openzip:
        try:
            df_cr = pickle.load(zipfile.ZipFile(in_zipfile, "r", zipfile.ZIP_DEFLATED).open(in_file, "r"))
            sys.stdout.write('{}, {}: ' .format(uniprot_acc, TMD))
            # filtering step
            df_cr = utils.filter_for_truncated_sequences(nonTMD_truncation_cutoff, df_cr)
            # save filtered dataframe to pickle
            out_file = "{}_{}_cr_df.pickle".format(protein_name, TMD)
            with zipfile.ZipFile(in_zipfile, mode="a", compression=zipfile.ZIP_DEFLATED) as zipout:
                pickle.dump(df_cr, open(out_file, "wb"), protocol=pickle.HIGHEST_PROTOCOL)
                zipout.write(out_file, arcname=out_file)
            os.remove(out_file)
        except:
            logging.info('pickle {} not found in zipfile - excluded'.format(in_file))
            pass



##############################################################
#                                                            #
#       define functions for AAIMON_slope curve fitting      #
#                                                            #
##############################################################

def residuals(constants, function, x, y):
    """
    Function used to optimise the fit of the curve to the data.
    It calculates the distance between y-value from real data and y-value from the function (sigmoid/sine/etc).
    """
    return y - function(constants, x)

def lin_AAIMON_slope_eq(a, x):
    y = a * x + 1 - 100 * a
    return y

def get_line_data_to_plot(a_constant, x_low=40, x_high=100):
    y_low = lin_AAIMON_slope_eq(a_constant, x_low)
    y_high = lin_AAIMON_slope_eq(a_constant, x_high)
    x_data = np.array([x_low, x_high])
    y_data = np.array([y_low, y_high])
    return x_data, y_data

def curve_fitting_fixed_100(x_array, y_array, a_constant_guess = -0.1):
    a_constant = leastsq(residuals, a_constant_guess, args = (lin_AAIMON_slope_eq, x_array, y_array))
    x_data, y_data = get_line_data_to_plot(a_constant[0])
    return float(a_constant[0]), x_data, y_data