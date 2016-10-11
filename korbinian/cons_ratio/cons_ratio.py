import ast
import os
import pickle
import zipfile
import korbinian
import korbinian.cons_ratio.calc
import korbinian.utils as utils
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from multiprocessing import Pool

def run_calculate_AAIMON_ratios(pathdict, s, logging):
    logging.info('~~~~~~~~~~~~      starting run_calculate_AAIMON_ratios        ~~~~~~~~~~~~')
    # if multiprocessing is used, log only to the console
    p_dict_logging = logging if s["use_multiprocessing"] != True else utils.Log_Only_To_Console()
    # set current working directory as the data_dir/homol, where temp files will be saved before moving to zip
    os.chdir(os.path.join(s["data_dir"], "homol"))
    # create list of protein dictionaries to process
    list_p = korbinian.utils.convert_summary_csv_to_input_list(s, pathdict, p_dict_logging)
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

        Contains:
            pathdict : dict
                Dictionary of the key paths and files associated with that List number.
            s : dict
                Settings dictionary extracted from excel settings file.
            logging : logging.Logger
                Logger for printing to console and logfile.

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

    fa_X_filt_full_str = " and X_in_match_seq == False" if s["fa_X_allowed_in_full_seq"] == False else ""

    fa_homol_query_str = 'FASTA_gapped_identity > {min_ident} and ' \
                        'FASTA_gapped_identity < {max_ident} and ' \
                        'hit_contains_SW_node == True and ' \
                        'disallowed_words_not_in_descr == True' \
                        '{Xfull}'.format(Xfull=fa_X_filt_full_str, min_ident=s["cr_min_identity_of_full_protein"], max_ident=s["cr_max_identity_of_full_protein"])

    # filter based on the query string
    dfh.query(fa_homol_query_str, inplace=True)

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
    df_nonTMD, mean_ser = korbinian.cons_ratio.calc.calc_nonTMD_perc_ident_and_gaps(df_nonTMD, mean_ser)

    with zipfile.ZipFile(homol_cr_ratios_zip, mode="w", compression=zipfile.ZIP_DEFLATED) as zipout:

        # save the nonTMD dataframe
        nonTMD_cr_outfile_pickle = "{}_nonTMD_cr_df.pickle".format(protein_name)
        with open(nonTMD_cr_outfile_pickle, "wb") as f:
            pickle.dump(df_nonTMD, f, protocol=pickle.HIGHEST_PROTOCOL)
        zipout.write(nonTMD_cr_outfile_pickle, arcname=nonTMD_cr_outfile_pickle)
        os.remove(nonTMD_cr_outfile_pickle)

        # filter nonTMD dataframe to only contain entries where nonTMD_perc_ident is not zero
        df_nonTMD = df_nonTMD.loc[df_nonTMD['nonTMD_perc_ident'] != 0]

        linspace_binlist = np.linspace(s["1p_smallest_bin"],
                                       s["1p_largest_bin"],
                                       s["1p_number_of_bins"])
        # add 30 as the last bin, to make sure 100% of the data is added to the histogram, including major outliers
        binlist = np.append(linspace_binlist,
                            s["1p_final_highest_bin"])

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
            df_cr, mean_ser = korbinian.cons_ratio.calc.calc_AAIMON(TMD, df_cr, len_query_TMD, mean_ser)

            logging.info('%s AAIMON MEAN %s: %0.2f' % (acc, TMD, mean_ser['%s_AAIMON_ratio_mean' % TMD]))
            # logging.info('%s AASMON MEAN %s: %0.2f' % (acc, TMD, mean_ser['%s_AASMON_ratio_mean'%TMD]))

            # save the dataframe for that TMD
            TM_cr_outfile_pickle = "{}_{}_cr_df.pickle".format(protein_name, TMD)
            with open(TM_cr_outfile_pickle, "wb") as f:
                pickle.dump(df_cr, f, protocol=pickle.HIGHEST_PROTOCOL)
            zipout.write(TM_cr_outfile_pickle, arcname=TM_cr_outfile_pickle)
            os.remove(TM_cr_outfile_pickle)

            # use the dictionary to obtain the figure number, plot number in figure, plot indices, etc
            newfig, savefig, fig_nr, plot_nr_in_fig, row_nr, col_nr = dict_organising_subplots[TMD_Nr]
            # if the TMD is the last one, the figure should be saved
            if TMD_Nr == len(list_of_TMDs):
                savefig = True
            # if a new figure should be created (either because the orig is full, or the last TMD is analysed)
            if newfig:
                # create a new figure
                fig, axarr = plt.subplots(nrows=nrows_in_each_fig,
                                          ncols=ncols_in_each_fig)  # sharex=True

            #" NOT STABLE! NEED TO CHANGE save_hist_AAIMON_ratio_single_protein SO THAT IT RUNS WITHIN THE FOR LOOP ABOVE, AND TAKES A SINGLE TMD AS INPUT, RATHER THAN LIST OF TMDS" / 4
            AAIMON_hist_path_prefix = p['AAIMON_hist_path_prefix']
            ########################################################################################
            #                                                                                      #
            #       Save histograms for each TMD of that protein, with relative conservation       #
            #                                                                                      #
            ########################################################################################
            korbinian.cons_ratio.histogram.save_hist_AAIMON_ratio_single_protein(fig_nr, fig, axarr, df_cr, s, TMD, binlist, zipout, row_nr, col_nr, fontsize, savefig, AAIMON_hist_path_prefix)

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

