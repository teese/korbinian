
#Switch matplotlib library for possible commandline execution
import matplotlib

# commented out as it causes errors in dependencies.
#matplotlib.use('Agg')

from multiprocessing import Pool
from scipy.optimize import leastsq
import ast
import itertools
import korbinian
import korbinian.utils as utils
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import pickle
from scipy.stats import sem
import sys
import zipfile
import io
# import debugging tools
from korbinian.utils import pr, pc, pn, aaa

def run_calculate_AAIMONs(pathdict, s, logging):
    """Runs calculate_AAIMONs for each protein, using multiprocessing Pool.

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
    logging.info('~~~~~~~~~~~~                      starting run_calculate_AAIMONs                    ~~~~~~~~~~~~')
    # if multiprocessing is used, log only to the console
    p_dict_logging = logging if s["use_multiprocessing"] != True else utils.Log_Only_To_Console()
    # set current working directory as the data_dir/homol, where temp files will be saved before moving to zip
    os.chdir(os.path.join(s["data_dir"], "homol"))
    # get list of accessions that could not be downloaded, and can immediately be excluded
    not_in_homol_db = utils.get_acc_list_from_txt(pathdict["acc_not_in_homol_db_txt"])

    ########################################################################################
    #                                                                                      #
    #   Extract accurate randTM and rand_nonTM values from csv files, where available      #
    #                                                                                      #
    ########################################################################################
    rand_ident_TM_csv = pathdict["rand_ident_TM_csv"]
    rand_ident_nonTM_csv = pathdict["rand_ident_nonTM_csv"]
    if os.path.isfile(rand_ident_TM_csv) and os.path.isfile(rand_ident_nonTM_csv):
        # open csv as series
        TM_ser = pd.Series.from_csv(rand_ident_TM_csv, sep="\t")
        nonTM_ser = pd.Series.from_csv(rand_ident_nonTM_csv, sep="\t")
        # extract the random values
        # REPLACES THE ORIGINAL VALUES FROM THE LISTS TAB OF THE SETTINGS FILE
        s["rand_TM"] = TM_ser["random_sequence_identity_output"]
        s["rand_nonTM"] = nonTM_ser["random_sequence_identity_output"]
    else:
        raise FileNotFoundError("CSV files with calculated random identity were not found.\nRe-run prepare lists.\n{}\n{}".format(rand_ident_TM_csv, rand_ident_nonTM_csv))

    # create list of protein dictionaries to process
    list_p = korbinian.utils.convert_summary_csv_to_input_list(s, pathdict, p_dict_logging, list_excluded_acc=not_in_homol_db)
    # number of processes is the number the settings, or the number of proteins, whichever is smallest
    n_processes = s["multiprocessing_cores"] if s["multiprocessing_cores"] < len(list_p) else len(list_p)

    if s["use_multiprocessing"]:
        with Pool(processes=n_processes) as pool:
            calc_AAIMON_list = pool.map(korbinian.cons_ratio.cons_ratio.calculate_AAIMONs, list_p)
            # log the list of protein results (e.g. acc, "simap", True) to the actual logfile, not just the console
            logging.info("\ncalc_AAIMON_list : {}".format(calc_AAIMON_list))
    else:
        for p in list_p:
            korbinian.cons_ratio.cons_ratio.calculate_AAIMONs(p)
    logging.info('~~~~~~~~~~~~                      finished run_calculate_AAIMONs                    ~~~~~~~~~~~~')

def calculate_AAIMONs(p):
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
            A6BM72_AAIMON_hist_0.png : png
                Histograms of AAIMON ratios for homologues of each TMD.
            A6BM72_cr_mean.csv : csv
                Summary file for that protein. Contains conservation ratios means.
                Will be gathered for all proteins by the gather_AAIMONs function.
            A6BM72_nonTMD_cr_df.pickle : pickled pd.DataFrame
                Dataframe containing the percentage_identity etc for sliced nonTMD region.
            A6BM72_SP01_cr_df.pickle : pickled pd.DataFrame
                Dataframe containing the percentage_identity etc for that particular TMD/region (in this case, the signal peptide).
            A6BM72_TM01_cr_df.pickle
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
    sys.stdout.write("|{}|".format(acc)), sys.stdout.flush()
    protein_name = p["protein_name"]
    fraction_TM_residues = p["perc_TMD"] / 100
    max_gaps = s["maxgaps_TMD"]
    # get the lipophilicity cutoffs for the homologue TMDs
    max_lipo_homol_settings_file = s["max_lipo_homol"]
    lipo_buffer = s["lipo_buffer"]

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
    mean_ser["number_of_TMDs"] = p['number_of_TMDs']
    mean_ser["number_of_TMDs_excl_SP"] = p['number_of_TMDs_excl_SP']

    list_of_TMDs = ast.literal_eval(p['list_of_TMDs'])
    mean_ser['list_of_TMDs'] = list_of_TMDs
    list_of_TMDs_excl_SP = ast.literal_eval(p['list_of_TMDs_excl_SP'])
    mean_ser['list_of_TMDs_excl_SP'] = list_of_TMDs_excl_SP

    ########################################################################################
    #                                                                                      #
    #                        Calculate AAIMON conservation ratios                          #
    #                                                                                      #
    ########################################################################################

    homol_cr_ratios_zip = p['homol_cr_ratios_zip']
    mean_ser_filename = "{}_cr_mean.csv".format(acc)



    #assume af first that there is no previous data, and that the calculations can be re-run
    if s["overwrite_prev_calculated_AAIMON_ratios"] == False:
        if os.path.isfile(homol_cr_ratios_zip):
            with zipfile.ZipFile(homol_cr_ratios_zip, mode="r", compression=zipfile.ZIP_DEFLATED) as cr_zip:
                if mean_ser_filename in cr_zip.namelist():
                    # if the means are saved as a csv in the homol cr_ratios zipfile, skip this protein
                    message = '{} AAIMONs skipped, file with mean AAIMON ratios exists (in settings, overwrite_prev_calculated_AAIMON_ratios = True)'.format(acc)
                    logging.info(message)
                    return acc, False, message

    disallowed_strings = dfh.loc[dfh["list_disallowed_words_in_descr"] != "[]"]["list_disallowed_words_in_descr"].dropna()
    disallowed_lists = list(disallowed_strings.apply(ast.literal_eval))
    mean_ser['disallowed_KW'] = list(set(itertools.chain(*disallowed_lists)))

    ########################################################################################
    #                                                                                      #
    #           Filter based on homol hit properties (% identity of full protein, etc)     #
    #                                                                                      #
    ########################################################################################
    cr_homol_query_str ='FASTA_gapped_identity > {min_ident} & ' \
                        'FASTA_gapped_identity < {max_ident} & ' \
                        'hit_contains_SW_node == True & ' \
                        'disallowed_words_not_in_descr == True &' \
                        'X_in_match_seq == False'.format(min_ident=s["min_ident"],
                                                         max_ident=s["max_ident"])

    # filter based on the query string
    dfh.query(cr_homol_query_str, inplace=True)

    #Calculate fasta values, add to mean_ser output series for each protein.
    # The mean fasta identity of homologues should almost always be around the centre of the min and max in the settings.
    mean_ser['FASTA_ident_mean'] = float('%0.2f' % dfh['FASTA_identity'].mean())
    # add mean observed changes as a percentage value to mean_ser
    if "obs_changes" not in dfh.columns:
        os.remove(p['homol_df_orig_zip'])
        message = "{} parsed simap file is out of date, will be removed".format(acc)
        logging.info(message)
        return acc, False, message
    mean_ser['obs_changes_mean'] = float('%0.2f' % dfh['obs_changes'].mean())

    # number of identical residues in each FASTA alignment can be calculated from identity and overlap
    dfh['FASTA_num_ident_res'] = dfh['FASTA_identity'] * dfh['FASTA_overlap']
    mean_ser['FASTA_num_ident_res_mean'] = float('%0.2f' % dfh['FASTA_num_ident_res'].mean())

    if not os.path.exists(p['fa_cr_sliced_TMDs_zip']):
        message = "{} Protein skipped. File does not exist".format(p['fa_cr_sliced_TMDs_zip'])
        logging.info(message)
        return acc, False, message

    if utils.file_is_old(p['fa_cr_sliced_TMDs_zip'], s["oldest_acceptable_file_date"]):
        os.remove(p['fa_cr_sliced_TMDs_zip']),
        message = "{} skipped, file is old and has been deleted".format(acc)
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
    rand_nonTM_csv = os.path.normpath(os.path.join(s["data_dir"], "summaries/{ln:02d}/List{ln:02d}_rand/List{ln:02d}_rand_nonTM.csv".format(ln=s["list_number"])))
    aa_prop_nonTM = pd.Series.from_csv(rand_nonTM_csv, sep="\t")
    # the first line is the random identity. Extract and delete.
    rand_nonTM = aa_prop_nonTM["random_sequence_identity_output"]

    # extract the length of the nonTMD region in the original query here, which should be in the original list summary file, created during uniprot_parse or OMPdb_get_TM_indices_and_slice
    df_nonTMD, mean_ser = korbinian.cons_ratio.calc.calc_nonTMD_perc_ident_and_gaps(df_nonTMD, mean_ser, p['len_nonTMD'], rand_nonTM)

    ########################################################################################
    #                                                                                      #
    #                Calculation of normalization factor for each homologue                #
    #                                                                                      #
    ########################################################################################
    dfh['norm_factor'] = dfh['FASTA_gapped_identity'].apply(korbinian.cons_ratio.norm.calc_aa_prop_norm_factor, args=(s["rand_TM"], s["rand_nonTM"], fraction_TM_residues))
    mean_ser["norm_factor_max"] = dfh['norm_factor'].max()
    mean_ser["norm_factor_min"] = dfh['norm_factor'].min()
    mean_ser["norm_factor_mean"] = dfh['norm_factor'].mean()
    # filter so that dfh is the same size as df_nonTMD
    dfh = dfh.loc[df_nonTMD.index, :]
    with zipfile.ZipFile(homol_cr_ratios_zip, mode="w", compression=zipfile.ZIP_DEFLATED) as zipout:

        # save the nonTMD dataframe BEFORE filtering
        nonTMD_cr_outfile_pickle = "{}_nonTMD_cr_df_unfiltered.pickle".format(protein_name)
        with open(nonTMD_cr_outfile_pickle, "wb") as f:
            pickle.dump(df_nonTMD, f, protocol=pickle.HIGHEST_PROTOCOL)
        zipout.write(nonTMD_cr_outfile_pickle, arcname=nonTMD_cr_outfile_pickle)
        os.remove(nonTMD_cr_outfile_pickle)

        # filter nonTMD dataframe to only contain entries where nonTMD_perc_ident is not zero
        # filter to remove short nonTMD regions
        # note this filtering is AFTER the full dataframe has been saved to file, preventing loss of data
        nonTMD_query_str = "nonTMD_perc_ident != 0 & " \
                           "perc_nonTMD_coverage > {min_perc_nonTMD} & "\
                           "nonTMD_len >= {min_nonTMD_len}".format(min_nonTMD_len=s["cr_min_len_nonTMD"],
                                                                   min_perc_nonTMD=s["min_perc_nonTMD"])

        n_homol_before_nonTMD_query = df_nonTMD.shape[0]
        df_nonTMD.query(nonTMD_query_str, inplace=True)
        n_homol_after_nonTMD_query = df_nonTMD.shape[0]
        mean_ser["n_homol_excluded_after_nonTMD_query"] = "__ {}/{} __".format(n_homol_before_nonTMD_query - n_homol_after_nonTMD_query, n_homol_before_nonTMD_query)

        # save the nonTMD dataframe AFTER filtering
        nonTMD_cr_outfile_pickle = "{}_nonTMD_cr_df.pickle".format(protein_name)
        with open(nonTMD_cr_outfile_pickle, "wb") as f:
            pickle.dump(df_nonTMD, f, protocol=pickle.HIGHEST_PROTOCOL)
        zipout.write(nonTMD_cr_outfile_pickle, arcname=nonTMD_cr_outfile_pickle)
        os.remove(nonTMD_cr_outfile_pickle)

        # Save the percentage nonTMD for all homologues
        # note that this might not match the final homologues with available AAIMON ratios
        mean_ser['perc_nonTMD_coverage_mean'] = df_nonTMD.perc_nonTMD_coverage.mean()

        #sys.stdout.write("\nn_homol_excluded_after_nonTMD_query {}".format(mean_ser["n_homol_excluded_after_nonTMD_query"]))

        ########################################################################################
        #                                                                                      #
        #                       Prepare some variables for the histograms                      #
        #                                                                                      #
        ########################################################################################

        linspace_binlist = np.linspace(s["1p_smallest_bin"],
                                       s["1p_largest_bin"],
                                       s["1p_number_of_bins"])
        # add 30 as the last bin, to make sure 100% of the data is added to the histogram, including major outliers
        binarray = np.append(linspace_binlist, s["1p_final_highest_bin"])

        # se default font size for text in the plot
        fontsize = 6
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

        AAIMON_all_TMD_dict = {}
        SW_num_ident_res_dict = {}
        list_homol_excluded_in_TMD_filter = []

        for TMD_Nr, TMD in enumerate(list_of_TMDs):
            # open the dataframe containing the sequences, gap counts, etc for that TMD only
            df_cr = utils.open_df_from_pickle_zip(p['fa_cr_sliced_TMDs_zip'], filename="{}_{}_sliced_df.pickle".format(protein_name, TMD), delete_corrupt=True)

            ########################################################################################
            #                                                                                      #
            #           Add columns from dfh and df_nonTMD to the df_cr for that TMD               #
            #                                                                                      #
            ########################################################################################
            # Add several columns from the original dataframe with homologues, directly parsed from SIMAP XML
            # SOME OF THIS MIGHT NOT BE NECESSARY ANY MORE
            df_cr['obs_changes'] = dfh['obs_changes']
            df_cr['norm_factor'] = dfh['norm_factor']
            # df_cr['FASTA_gapped_identity'] = dfh['FASTA_gapped_identity']
            # df_cr['norm_factor'] = dfh['norm_factor']
            # df_cr['len_full_match_seq'] = dfh['len_full_match_seq']
            #
            # # add several columns from the df_nonTMD, such as the nonTMD percentage identity
            # # SOME OF THIS MIGHT NOT BE NECESSARY ANY MORE
            # df_cr['nonTMD_perc_ident'] = df_nonTMD['nonTMD_perc_ident']
            # df_cr['nonTMD_perc_sim_plus_ident'] = df_nonTMD['nonTMD_perc_sim_plus_ident']
            # df_cr['nonTMD_SW_align_len_excl_gaps'] = df_nonTMD['nonTMD_SW_align_len_excl_gaps']
            # df_cr['perc_nonTMD_coverage'] = df_nonTMD['perc_nonTMD_coverage']
            # df_cr['len_nonTMD_orig_q_minus_nonTMD_SW_align_len_excl_gaps'] = df_nonTMD['len_nonTMD_orig_q_minus_nonTMD_SW_align_len_excl_gaps']

            # FILTER BASED ON df_nonTMD, which already has excluded homologues that are unsuitable (e.g. % identity of full protein, nonTMD_perc_ident is not zero, etc)
            df_cr = df_cr.loc[df_nonTMD.index,:]

            ########################################################################################
            ########################################################################################
            #                                                                                      #
            #                       Calculate AAIMON, AASMON for each TMD                          #
            #                                                                                      #
            ########################################################################################
            ########################################################################################

            len_query_TMD = p["%s_end"%TMD] - p["%s_start"%TMD]
            #df_cr = korbinian.cons_ratio.calc.calc_AAIMON(TMD, df_cr, len_query_TMD)


            ########################################################################################
            #                                                                                      #
            #              count number of identical (and/or similar) residues in                  #
            #                  the TMD region of interest, for each homologue                      #
            #                                                                                      #
            ########################################################################################
            # count identical residues between query and match TMDs by counting the number of pipes in the markup string
            # NOTE THAT df_cr['%s_SW_markup_seq'%TMD].str.count('|') DOES NOT WORK, as "|" has a function in regex and needs to be escaped
            df_cr['%s_SW_num_ident_res' % TMD] = df_cr['%s_SW_markup_seq' % TMD].str.count('\|')
            df_cr['%s_SW_num_sim_res' % TMD] = df_cr['%s_SW_markup_seq' % TMD].str.count(':')

            ########################################################################################
            #                                                                                      #
            #     calculate the length of the alignment of the TMD region, including gaps          #
            #          (should equal the len_query_TMD + n_gaps_in_q + n_gaps_in_m)                #
            #                                                                                      #
            ########################################################################################
            df_cr['%s_SW_align_len' % TMD] = df_cr['%s_SW_match_seq' % TMD].str.len()

            ##############################################################################################################
            #                                                                                                            #
            #            calculate the length of the alignment of the TMD region, EXCLUDING GAPS                         #
            #    TM01_SW_align_len_excl_gaps = TM01_SW_align_len - TM01_SW_query_num_gaps - TM01_SW_match_num_gaps       #
            #                                                                                                            #
            ##############################################################################################################
            df_cr['%s_SW_align_len_excl_gaps' % TMD] = df_cr['%s_SW_align_len' % TMD] - df_cr['%s_SW_query_num_gaps' % TMD] - df_cr['%s_SW_match_num_gaps' % TMD]

            ########################################################################################
            #                                                                                      #
            #              calculate the percentage identity of the TMD region                     #
            #       TM01_perc_ident = TM01_SW_num_ident_res / TM01_SW_align_len_excl_gaps          #
            #                                                                                      #
            ########################################################################################
            # the percentage identity of that TMD is defined as the number of identical residues (pipes in markup)
            # divided by the length of the the aligned residues (excluding gaps)
            # note that the nonTMD percentage identity is calculated the same way
            df_cr['%s_perc_ident' % TMD] = df_cr['%s_SW_num_ident_res' % TMD] / df_cr['%s_SW_align_len_excl_gaps' % TMD]
            # calculate percentage similar residues
            df_cr['%s_perc_sim' % TMD] = df_cr['%s_SW_num_sim_res' % TMD] / df_cr['%s_SW_align_len' % TMD]
            # add together to obtain the percentage similar + identical residues
            df_cr['%s_perc_sim_plus_ident' % TMD] = df_cr['%s_perc_ident' % TMD] + df_cr['%s_perc_sim' % TMD]

            ########################################################################################
            #                                                                                      #
            #          calculate Amino Acid Identity : Membranous Over Nonmembranous               #
            #             (number of gaps)/(length of sequence excluding gaps)                     #
            #                                                                                      #
            ########################################################################################
            df_cr['%s_AAIMON' % TMD] = df_cr['%s_perc_ident' % TMD] / df_nonTMD['nonTMD_perc_ident']
            # calculate the Amino Acid Similarity : Membranous Over Nonmembranous (AASMON) (includes similarity + identity based on the matrix used in the SW alignment of SIMAP)
            df_cr['%s_AASMON' % TMD] = df_cr['%s_perc_sim_plus_ident' % TMD] / df_nonTMD['nonTMD_perc_sim_plus_ident']
            # calculate the AAIMON normalised by the random_AA_identity, to exclude identity due to lipophilicity
            df_cr['%s_AAIMON_n' % TMD] = df_cr['%s_AAIMON' % TMD] / dfh['norm_factor']

            ########################################################################################
            #                                                                                      #
            #  calculate ratio of length of TMD to length of nonTMD excl gaps & full match seq     #
            #                                                                                      #
            ########################################################################################
            df_cr['%s_ratio_len_TMD_to_len_nonTMD' % TMD] = len_query_TMD / df_nonTMD['nonTMD_SW_align_len_excl_gaps']
            df_cr['%s_ratio_len_TMD_to_len_full_match_seq' % TMD] = len_query_TMD / dfh['len_full_match_seq']

            ########################################################################################
            #                                                                                      #
            #             FILTERING TMD SEQS based on lipophilicity, etc                           #
            #                                                                                      #
            ########################################################################################

            # Set the lipophilicity cutoff as either the absolute value in the settings file,
            # or the lipophilicity of the original TM plus a buffer. Whichever is highest.
            orig_lipo = p["{}_lipo".format(TMD)]
            orig_TM_lipo_plus_buffer = orig_lipo + lipo_buffer
            # take whichever value is highest
            max_lipo_homol = np.array([max_lipo_homol_settings_file, orig_TM_lipo_plus_buffer]).max()

            """ Filtering of the homologues for this particular TMD
             - max number of gaps in query
             - max number of gaps in match
             - max polarity (lipophilicity), which mostly excludes frameshifts
             # OPTIONAL
             - percentage identity of TMD does not equal zero (gives a TM/nonTM ratio of 0)
            """
            ## count the number of gaps in the query and match sequences
            cr_TMD_query_str = '{TMD}_SW_query_num_gaps <= {max_gaps} & ' \
                               '{TMD}_SW_match_num_gaps <= {max_gaps} & ' \
                               '{TMD}_SW_match_lipo <= {max_lipo_homol}'.format(TMD=TMD, max_gaps=max_gaps, max_lipo_homol=max_lipo_homol)
            n_homol_before_TMD_filter = df_cr.shape[0]
            # filter by the above query
            df_cr.query(cr_TMD_query_str, inplace=True)

            # filter to exclude TMDs where the % identity of that TMD is zero (no identical residues)
            # this may remove some valid homologues, especially when far homologues are included in the analysis
            # Would be better to replace with number_of_ident_residues_all_TMDs / total_len_TMD_in_query, but this would disturb analyses of individual TM conservation.
            if s["exclude_homol_with_0_ident_residues_in_any_TMD"]:
                df_cr.query('{TMD}_perc_ident != 0'.format(TMD=TMD), inplace=True)

            n_homol_after_TMD_filter = df_cr.shape[0]

            # FOR SOME REASON -- 280/280 -- in the nonTMD query (e.g. no homologues with min % ident) gives -- 280/280 280/280 280/280 280/280 -- for all the TMDs.
            # If the first filter removes all homologues, then there ARE NO homologues to remove for the remaining filters
            if n_homol_after_nonTMD_query != 0:
                excluded_string = " {}/{} ".format(n_homol_before_TMD_filter - n_homol_after_TMD_filter, n_homol_before_TMD_filter)
                list_homol_excluded_in_TMD_filter.append(excluded_string)
            else:
                list_homol_excluded_in_TMD_filter = "__ 0/0 __"
                sys.stdout.write("\n{} {}:{}/{} valid homol.".format(acc, TMD, n_homol_after_TMD_filter, n_homol_before_TMD_filter))
                sys.stdout.flush()
            #sys.stdout.write("\nn_homol_excluded_after_TMD_filter", mean_ser["n_homol_excluded_after_TMD_filter"]), sys.stdout.flush()

            ########################################################################################
            #                                                                                      #
            #           Add the AAIMON values for this TMD to AAIMON_all_TMD_dict.                 #
            #  This will be used as a filter to select homologues that have good data for all TMDs #
            #   Also create an array of the number of identical residues
            #                                                                                      #
            ########################################################################################
            AAIMON_all_TMD_dict['%s_AAIMON' % TMD] = df_cr['%s_AAIMON' % TMD].dropna()
            SW_num_ident_res_dict['%s_SW_num_ident_res' % TMD] = df_cr['%s_SW_num_ident_res' % TMD].dropna()

            TM_cr_pickle = "{}_{}_cr_df.pickle".format(acc, TMD)
            with open(TM_cr_pickle, "wb") as f:
                pickle.dump(df_cr, f, protocol=pickle.HIGHEST_PROTOCOL)
            zipout.write(TM_cr_pickle, arcname=TM_cr_pickle)
            # hopefully the file is actually removed here
            #os.remove(TM_cr_pickle)

        # make a readable list of the number of homologues excluded for each TMD, based on filtering for gaps and lipophilicity, etc
        mean_ser["n_homol_excluded_after_TMD_filter"] = "__ {} __".format("".join(list_homol_excluded_in_TMD_filter))

        # save how many hits have a Smith Waterman alignment node in the XML
        # shows the status of the XML files, which can be corrupt, and also number of homologues
        value_counts_hit_contains_SW_node = dfh['hit_contains_SW_node'].value_counts()
        if True in value_counts_hit_contains_SW_node:
            mean_ser['num_hits_with_SW_align_node'] = value_counts_hit_contains_SW_node[True]
        else:
            logging.warning("{} num_hits_with_SW_align_node = 0".format(protein_name))
            mean_ser['num_hits_with_SW_align_node'] = 0

        ########################################################################################
        #                                                                                      #
        #               AAIMON normalization and save fig for each protein                     #
        #  	df_AAIMON_all_TMD:                                                                 #
        #       index = hit_num                                                                #
        #       columns : TM01_AAIMON 	TM02_AAIMON 	TM03_AAIMON                            #
        #                                                                                      #
        ########################################################################################
        df_AAIMON_all_TMD = pd.DataFrame(AAIMON_all_TMD_dict)
        if df_AAIMON_all_TMD.empty:
            # if this returns an empty dataframe, it means that no homologues with calculable AAIMONs were available
            # skip to the next protein
            zipout.close()
            os.remove(homol_cr_ratios_zip)
            message = "{} skipped, no homologues with calculable AAIMONs were available. homol_cr_ratios_zip will be deleted".format(acc)
            logging.info(message)
            return acc, False, message

        # get the original column names (M01_AAIMON, TM02_AAIMON, etc)
        AAIMON_cols = ["{}_AAIMON".format(TMD) for TMD in list_of_TMDs_excl_SP]
        # Get the number of TMDs whose AAIMON ratio was calculable
        # For truncated alignments, C and N-term TMDs may be missing
        # steps : 0) select just AAIMON data 1) convert nan to "", 2) count numbers with np.isreal,
        #         3) sum TRUE FALSE TRUE etc to get the number of TMDs with calculable AAIMON for that homologue
        df_AAIMON_all_TMD["n_TMDs_with_measurable_AAIMON"] = df_AAIMON_all_TMD.loc[:, AAIMON_cols].fillna("").applymap(np.isreal).sum(axis=1)
        # create a bool if all TMDs have calculable AAIMON ratios
        df_AAIMON_all_TMD["all_TMDs_have_AAIMON"] = df_AAIMON_all_TMD["n_TMDs_with_measurable_AAIMON"] == p["number_of_TMDs_excl_SP"]


        df_AAIMON_all_TMD['gapped_ident'] = dfh.loc[df_AAIMON_all_TMD.index, 'FASTA_gapped_identity']
        df_AAIMON_all_TMD['norm_factor'] = dfh.loc[df_AAIMON_all_TMD.index, 'norm_factor']
        df_AAIMON_all_TMD['perc_nonTMD_coverage'] = df_nonTMD.loc[df_AAIMON_all_TMD.index, 'perc_nonTMD_coverage']


        ########################################################################################
        #                                                                                      #
        #                ADD AN ARRAY OF THE NUMBER OF IDENTICAL RESIDUES                      #
        #    CALC % ident of TM region for all residues combined
        #     This is important for datasets with far homologues, where TMDs may have
        #      homologues without any identical residues (TM/nonTM = 0/x = 0)
        #                                                                                      #
        ########################################################################################
        # convert dict to dataframe. index = hit number, columns = Index(['TM01_SW_num_ident_res', 'TM02_SW_num_ident_res', etc
        df_SW_num_ident_res = pd.DataFrame(SW_num_ident_res_dict)
        # get the original column names (M01_AAIMON, TM02_AAIMON, etc)
        num_ident_res_orig_cols = ['{}_SW_num_ident_res'.format(TMD) for TMD in list_of_TMDs_excl_SP]

        # save original cols with data (other columns will be added to same dataframe)
        num_ident_res_orig_cols = df_SW_num_ident_res.columns
        # add the evolutionary distance of each homologue (% identity of full protein)
        df_SW_num_ident_res["obs_changes"] = dfh["obs_changes"]
        df_SW_num_ident_res["norm_factor"] = dfh["norm_factor"]
        df_SW_num_ident_res["nonTMD_perc_ident"] = df_nonTMD["nonTMD_perc_ident"]

        # keep only homologues that gave valid AAIMONs in the preceding dataframe
        df_SW_num_ident_res = df_SW_num_ident_res.loc[df_AAIMON_all_TMD.index, :]

        TMD_seq_joined_len = p["TMD_seq_joined_len"]

        df_SW_num_ident_res["SW_num_ident_res_all_TM"] = df_SW_num_ident_res.loc[:, num_ident_res_orig_cols].sum(axis=1)
        df_SW_num_ident_res["perc_ident_all_TM_res"] = df_SW_num_ident_res["SW_num_ident_res_all_TM"] / TMD_seq_joined_len

        # for each protein, save the mean % identity of TM and nonTM regions, to plot separately in scatter graph
        mean_ser['TMD_perc_ident_mean'] = df_SW_num_ident_res["perc_ident_all_TM_res"].mean()
        mean_ser['nonTMD_perc_ident_mean'] = df_SW_num_ident_res["nonTMD_perc_ident"].mean()

        df_SW_num_ident_res["AAIMON_all_TM_res"] = df_SW_num_ident_res["perc_ident_all_TM_res"] / df_SW_num_ident_res["nonTMD_perc_ident"]
        df_SW_num_ident_res["AAIMON_n_all_TM_res"] = df_SW_num_ident_res["AAIMON_all_TM_res"] / df_SW_num_ident_res["norm_factor"]

        # linear regression. X-axis = evolutionary distance (observed_changes), y-axis = AAIMON. Slope calculated with fixed 100% identity at AAIMON 1.0
        AAIMON_slope, x_data, y_data = fit_data_to_linear_function(df_SW_num_ident_res["obs_changes"], df_SW_num_ident_res["AAIMON_all_TM_res"])
        AAIMON_n_slope, x_data_n, y_data_n = fit_data_to_linear_function(df_SW_num_ident_res["obs_changes"], df_SW_num_ident_res["AAIMON_n_all_TM_res"])

        mean_ser['AAIMON_mean_all_TM_res'] = df_SW_num_ident_res["AAIMON_all_TM_res"].mean()
        mean_ser['AAIMON_std_all_TM_res'] = df_SW_num_ident_res["AAIMON_all_TM_res"].std()

        mean_ser['AAIMON_n_mean_all_TM_res'] = df_SW_num_ident_res["AAIMON_n_all_TM_res"].mean()
        mean_ser['AAIMON_n_std_all_TM_res'] = df_SW_num_ident_res["AAIMON_n_all_TM_res"].std()

        mean_ser['AAIMON_slope_all_TM_res'] = AAIMON_slope
        mean_ser['AAIMON_n_slope_all_TM_res'] = AAIMON_n_slope

        mean_ser['perc_ident_mean'] = 100 - df_SW_num_ident_res["obs_changes"].mean()

        ########################################################################################
        #                                                                                      #
        #          Plot the evol distance against AAIMON for all residues in protein           #
        #          (Note that this excludes signal peptides where examined!)                   #
        #                                                                                      #
        ########################################################################################
        fig, ax = plt.subplots()
        # define data to plot
        datapointsize = 0.5
        x_data_obs_changes = df_SW_num_ident_res["obs_changes"]
        scatter_data_AAIMON = df_SW_num_ident_res["AAIMON_all_TM_res"]
        scatter_data_AAIMON_n = df_SW_num_ident_res["AAIMON_n_all_TM_res"]

        # calculate angle between AAIMON_slope and AAIMON_n_slope
        angle = korbinian.cons_ratio.histogram.angle_between_slopes((x_data[0], y_data[0]), (x_data[1], y_data[1]))
        mean_ser['angle_between_slopes_mean_all_TM_res'] = angle

        xlim_min = 0
        xlim_max = 60

        ax.scatter(x_data_obs_changes, scatter_data_AAIMON, color="k", alpha=0.3, s=datapointsize)
        ax.scatter(x_data_obs_changes, scatter_data_AAIMON_n, color=(0.843, 0.098, 0.1098), marker='^', alpha=0.3, s=datapointsize)
        ax.plot(x_data, y_data, color="k", alpha=0.3)
        ax.plot(x_data, y_data_n, color=(0.843, 0.098, 0.1098), alpha=0.3)

        ax.set_ylabel('%s AAIMON' % TMD, rotation='vertical', fontsize=fontsize)
        # ax.set_xlabel('% identity')
        ax.set_ylim(0.0, 3)
        ax.annotate(s='AAIMON slope: {a:0.3f}, AAIMON_n slope: {b:0.3f}, angle {c:.2f}°'.format(a=AAIMON_slope, b=AAIMON_n_slope, c=angle),
                                       xy=(0.01, 1.01), xytext=None, xycoords='axes fraction', alpha=0.75, fontsize=fontsize)
        # ax.set_xticks(range(xlim_min,xlim_max+1,10))

        # set x-axis min
        ax.set_xlim(xlim_min, xlim_max)
        # set x-axis ticks
        # use the slide selection to select every second item in the list as an xtick(axis label)
        ax.set_xticks(range(xlim_min, xlim_max + 1, 10))
        ax.set_xlabel('% observed aa substitutions', fontsize=fontsize)
        # change axis font size
        ax.tick_params(labelsize=fontsize)
        # create legend?#http://stackoverflow.com/questions/9834452/how-do-i-make-a-single-legend-for-many-subplots-with-matplotlib
        ax.legend(['AAIMON slope', 'AAIMON_n slope', 'AAIMON', 'AAIMON_n'], loc='upper right', fontsize=fontsize)
        # add background grid
        ax.grid(True, color='0.75', alpha=0.5)
        # automatically tighten the layout of plots in the figure
        fig.tight_layout()
        figpath = p['norm_scatter_path_prefix'] + '_all_TM_res.png'
        # save files
        fig.savefig(figpath, format='png', dpi=200)
        # close figure
        plt.close('all')
        # add to zipfile
        zipout.write(figpath, arcname=os.path.basename(figpath))
        # delete temporory files
        os.remove(figpath)
        ########################################################################################

        mean_ser['%s_angle_between_slopes' % TMD] = angle

        df_AAIMON_all_TMD['AAIMON_mean_all_TMDs_1_homol'] = df_AAIMON_all_TMD.loc[:, AAIMON_cols].mean(axis=1)
        df_AAIMON_all_TMD['AAIMON_n_mean_all_TMDs_1_homol'] = df_AAIMON_all_TMD['AAIMON_mean_all_TMDs_1_homol'] / df_AAIMON_all_TMD['norm_factor']

        ########################################################################################
        #                                                                                      #
        #                Create an index of homologues with data for all TMDs                  #
        #           (Important filter to prevent truncation and ensure balance)                #
        #                                                                                      #
        ########################################################################################
        # first get a list of all the homologues that have AAIMON ratios for all TMDs
        df_AAIMON_all_TMD["AAIMON_avail_all_TMDs"] = df_AAIMON_all_TMD.n_TMDs_with_measurable_AAIMON == p["number_of_TMDs_excl_SP"]
        filt_index = df_AAIMON_all_TMD[df_AAIMON_all_TMD["AAIMON_avail_all_TMDs"] == True].index.tolist()

        # add the number of homologues with AAIMON for all TMDs
        mean_ser['AAIMON_n_homol'] = len(filt_index)

        if mean_ser['AAIMON_n_homol'] == 0:
            # save the pandas series with the means to a csv in the cr_ratios zip file
            mean_ser.to_csv(mean_ser_filename)
            zipout.write(mean_ser_filename, arcname=mean_ser_filename)
            os.remove(mean_ser_filename)
            message = "{} has no valid homologues.".format(acc)
            logging.info(message)
            return acc, False, message

        # save the dataframe containing normalisation factor and normalised AAIMON to zipout
        df_AAIMON_all_TMD.to_csv(protein_name + '_AAIMON_all_TMD.csv')
        zipout.write(protein_name + '_AAIMON_all_TMD.csv', arcname=protein_name + '_AAIMON_all_TMD.csv')
        os.remove(protein_name + '_AAIMON_all_TMD.csv')

        # count how many TMDs actually have data
        n_TMDs_with_data = 0
        # since the data is saved in separate files, iterate through them again
        for TMD in list_of_TMDs:
            # find the TMD number (starting from 1)
            TMD_Nr = list_of_TMDs.index(TMD) + 1
            # re-open dataframe from zipped pickle
            # keeping the temp pickle file from above and deleting later was discontinued, as there were troubles deleting the file for some reason

            TM_cr_pickle = "{}_{}_cr_df.pickle".format(acc, TMD)
            # df_cr = utils.open_df_from_pickle_zip(homol_cr_ratios_zip, filename=TM_cr_pickle, delete_corrupt=False)

            if utils.file_is_old(TM_cr_pickle, s["oldest_acceptable_file_date"]):
                os.remove(TM_cr_pickle),
                message = "{} skipped, file is old and has been deleted".format(acc)
                logging.info(message)
                return acc, False, message

            df_cr = pd.read_pickle(TM_cr_pickle)

            # try:
            #     logging.info(os.path.getmtime(f))
            #     #df_cr = pd.read_pickle(f)
            #     nonTMD_cr_outfile_pickle
            # except (EOFError, ModuleNotFoundError, io.UnsupportedOperation) as e:
            #     # pickle files created by an older version of python have compatibility issues
            #     message = "{} is created by an older version of pandas and has been deleted".format(TM_cr_pickle)
            #     logging.info(message)
            #     os.remove(TM_cr_pickle)
            #     return acc, False, message
            #with open(TM_cr_pickle, 'rb') as f:
               #df_cr = pickle.load(f)

            # check that it's definitely a dataframe
            assert isinstance(df_cr, (pd.DataFrame))

            # # keep only the homologues that had measurable AAIMON ratios for all proteins
            # df_cr = df_cr.loc[filt_index, :]

            #DEPRECATED. SEEMS TO BE SOLVED BY RE-SLICING (NOPE!)
            # There is a bug that caused a keyerror for an unknown protein which had only one valid homologue : KeyError: 'None of [[2044.0]] are in the [index]' (D4H7V8, J1FM30, J7QHH8?)
            # re-slicing did not fix this (M2BWP2?, S9SJA4, E6W3I6, A9L0Y0)
            try:
                df_cr = df_cr.loc[filt_index, :]
            except KeyError:
                #df_cr = pd.DataFrame(columns=df_cr.columns)
                message = "{} has {} valid homologues({}), and gives a KeyError when filtering, suggesting they are not in the current" \
                          " homologue list.".format(acc, mean_ser['AAIMON_n_homol'], filt_index)
                return acc, False, message

            ########################################################################################
            #                                                                                      #
            #                       Calculate AAIMON_slope, AAIMON_n_slope                         #
            #                                                                                      #
            ########################################################################################
            # OLD DROPNA SHOULD BE OBSOLETE, AS ONLY HOMOLOGUES WITH AAIMON FOR ALL TMDs ARE KEPT
            # drop every row (hit) in df_cr that contains NaN in column TMxy_AAIMON - important for line fit that can't handle NAN
            #df_cr = df_cr[np.isfinite(df_cr['%s_AAIMON'%TMD])]
            #df_cr.dropna(subset=['%s_AAIMON'%TMD], inplace=True)

            obs_changes = df_cr['obs_changes']
            AAIMON = df_cr['%s_AAIMON'%TMD]            # y-axis
            AAIMON_n = df_cr['%s_AAIMON_n'%TMD]        # y-axis

            if len(obs_changes) == 0 or len(AAIMON) == 0:
                # There is no gapped identity or AAIMON data for these homologues, skip to next TMD
                logging.info("{} {} No observed changes or AAIMON data for this TMD".format(acc, TMD))
                continue
            else:
                n_TMDs_with_data += 1
            # linear regression for non-norm. and norm. AAIMON with fixed 100% identity at AAIMON 1.0
            AAIMON_slope, x_data, y_data = fit_data_to_linear_function(obs_changes, AAIMON)
            mean_ser['%s_AAIMON_slope' % TMD] = AAIMON_slope
            AAIMON_n_slope, x_data_n, y_data_n = fit_data_to_linear_function(obs_changes, AAIMON_n)
            mean_ser['%s_AAIMON_n_slope' % TMD] = AAIMON_n_slope

            if '{TMD}_SW_match_lipo'.format(TMD=TMD) not in df_cr.columns:
                message = "{} {}_SW_match_lipo not found in columns. Slice file is out of date and will be deleted.".format(acc, TMD)
                os.remove(p['fa_cr_sliced_TMDs_zip'])
                return acc, False, message

            ###########################################################################
            #                                                                         #
            #            DEPRECATED: to filter for truncated alignments               #
            #                             DO NOT DELETE!!!                            #
            #                                                                         #
            ###########################################################################
            """
            # save the unfiltered dataframe for that TMD
            TM_cr_pickle = "{}_{}_cr_df_RAW.pickle".format(protein_name, TMD)
            with open(TM_cr_pickle, "wb") as f:
                pickle.dump(df_cr, f, protocol=pickle.HIGHEST_PROTOCOL)
            zipout.write(TM_cr_pickle, arcname=TM_cr_pickle)
            os.remove(TM_cr_pickle)
            """

            # DEPRECATED
            # number of homologues with valid AAIMON ratios for that TMD.
            # since ALL TMDs have to be in each homologue before AAIMON is calculated, this number is often the same for all TMDs
            # this will be different for each TMD when some are excluded due to low lipophilicity, for example
            #mean_ser['{}_AAIMON_n_homol'.format(TMD)] = df_cr.shape[0]
            #sys.stdout.write("_AAIMON_n_homol", mean_ser['{}_AAIMON_n_homol']), sys.stdout.flush()

            logging.info('%s AAIMON_slope %s: %0.5f' % (acc, TMD, mean_ser['%s_AAIMON_slope' % TMD]))
            # OTHER VALUES TO LOG, IF DESIRED
            #logging.info('%s AAIMON_mean %s: %0.2f' % (acc, TMD, mean_ser['%s_AAIMON_mean' % TMD]))
            #logging.info('%s AAIMON_n_mean %s: %0.2f' % (acc, TMD, mean_ser['%s_AAIMON_n_mean' % TMD]))
            #logging.info('%s AAIMON_n_slope %s: %0.5f' % (acc, TMD, mean_ser['%s_AAIMON_n_slope' % TMD]))
            # logging.info('%s AASMON MEAN %s: %0.2f' % (acc, TMD, mean_ser['%s_AASMON_mean'%TMD]))

            # use the dictionary to obtain the figure number, plot number in figure, plot indices, etc
            newfig, savefig, fig_nr, plot_nr_in_fig, row_nr, col_nr = dict_organising_subplots[TMD_Nr]
            # if the TMD is the last one, the figure should be saved
            if TMD_Nr == len(list_of_TMDs):
                savefig = True
            # if a new figure should be created (either because the orig is full, or the last TMD is analysed)
            if newfig:
                # create a new figure for histograms
                fig, axarr = plt.subplots(nrows=nrows_in_each_fig, ncols=ncols_in_each_fig)  # sharex=True
                # create a new figure for scatter plot
                fig2, axarr2 = plt.subplots(nrows=nrows_in_each_fig, ncols=ncols_in_each_fig)  # sharex=True

            #" NOT STABLE! NEED TO CHANGE save_hist_AAIMON_single_protein SO THAT IT RUNS WITHIN THE FOR LOOP ABOVE, AND TAKES A SINGLE TMD AS INPUT, RATHER THAN LIST OF TMDS" / 4
            AAIMON_hist_path_prefix = p['AAIMON_hist_path_prefix']
            norm_scatter_path_prefix = p['norm_scatter_path_prefix']
            ########################################################################################
            #                                                                                      #
            #       Save histograms for each TMD of that protein, with relative conservation       #
            #                                                                                      #
            ########################################################################################
            # NECESSARY???
            # if axarr is None:
            #     # avoid bug where the axarr is still not created, as newfig was not True for first figure?
            #     continue

            # create histograms for this protein
            korbinian.cons_ratio.histogram.save_hist_AAIMON_single_protein(fig_nr, fig, axarr, df_cr, s, TMD, binarray, zipout, row_nr, col_nr, fontsize, savefig, AAIMON_hist_path_prefix)
            # create scatterplots for this protein
            angle = korbinian.cons_ratio.histogram.save_scatter_AAIMON_norm_and_AAIMON_slope_single_protein(fig_nr, fig2, axarr2, df_cr, x_data, y_data, y_data_n, AAIMON_slope, AAIMON_n_slope, TMD, zipout, row_nr, col_nr, fontsize, savefig, norm_scatter_path_prefix)
            mean_ser['%s_angle_between_slopes' %TMD] = angle

            ########################################################################################
            ########################################################################################
            #                                                                                      #
            #                SAVING MEAN VALUES TO SERIES FOR LATER GATHER SCRIPT                  #
            #                                                                                      #
            ########################################################################################
            ########################################################################################
            # percentage identities
            mean_ser['%s_perc_ident_mean' % TMD] = df_cr['%s_perc_ident' % TMD].mean()
            mean_ser['%s_perc_sim_mean' % TMD] = df_cr['%s_perc_sim' % TMD].mean()
            mean_ser['%s_perc_sim_plus_ident_mean' % TMD] = df_cr['%s_perc_sim_plus_ident' % TMD].mean()
            # AAIMON ratios
            mean_ser['%s_AAIMON_mean' % TMD] = float(df_cr['%s_AAIMON' % TMD].mean())
            mean_ser['%s_AAIMON_n_mean' % TMD] = float(df_cr['%s_AAIMON_n' % TMD].mean())
            mean_ser['%s_AAIMON_std' % TMD] = df_cr['%s_AAIMON' % TMD].std()
            mean_ser['%s_AASMON_mean' % TMD] = df_cr['%s_AASMON' % TMD].mean()
            mean_ser['%s_AASMON_std' % TMD] = df_cr['%s_AASMON' % TMD].std()
            # ratios for length of TMDs
            mean_ser['%s_ratio_len_TMD_to_len_nonTMD_mean' % TMD] = float('%0.2f' % df_cr['%s_ratio_len_TMD_to_len_nonTMD' % TMD].dropna().mean())
            mean_ser['%s_ratio_len_TMD_to_len_full_match_seq_mean' % TMD] = float('%0.2f' % df_cr['%s_ratio_len_TMD_to_len_full_match_seq' % TMD].dropna().mean())
            # gaps per residue
            mean_ser['%s_SW_q_gaps_per_q_residue_mean' % TMD] = df_cr['%s_SW_q_gaps_per_q_residue' % TMD].dropna().mean()

        if n_TMDs_with_data == 0:
            # save the pandas series with the means to a csv in the cr_ratios zip file
            mean_ser.to_csv(mean_ser_filename)
            zipout.write(mean_ser_filename, arcname=mean_ser_filename)
            os.remove(mean_ser_filename)
            # skip protein, as no data is available
            return acc, False, "{} skipped. No TMDs contained homologue data.".format(acc)

        # create a dict to cycle through the functions : mean, std, SEM, etc
        function_dict = {"mean" : np.mean, "std" : np.std, "SEM": sem}
        # ignore std and SEM if there are only one TMD in the list
        fxn_list = ["mean"] if len(list_of_TMDs_excl_SP) == 1 else function_dict.keys()
        # choose the values to apply the statistical functions
        colnames = ["{}_AAIMON_mean", "{}_AAIMON_n_mean", "{}_AAIMON_slope", "{}_AAIMON_n_slope", "{}_angle_between_slopes", "{}_AASMON_mean"]
        for colname in colnames:
            sel_cols = [colname.format(TMD) for TMD in list_of_TMDs_excl_SP]
            # if all the selected cols are in the index
            if len(set(sel_cols).intersection(set(mean_ser.index))) == len(sel_cols):

                for function_name in fxn_list:
                    # get function (e.g. sem)
                    function = function_dict[function_name]
                    # make new variable name (eg. AAIMON_mean_all_TMDs_mean)
                    new_var_name = colname[3:] + "_all_TMDs_" + function_name
                    # add the mean, std or SEM to the output series
                    mean_ser[new_var_name] = function(mean_ser.loc[sel_cols])

        # add AAIMON_slope of the last TMD
        last_TMD = p['last_TMD']
        mean_ser["AAIMON_slope_last_TMD"] = mean_ser['%s_AAIMON_slope' % last_TMD]
        mean_ser["AAIMON_n_slope_last_TMD"] = mean_ser['%s_AAIMON_n_slope' %last_TMD]

        # iterate through each TMD, and calculate mean AAIMON slopes for central TMDs
        if mean_ser['number_of_TMDs'] >= 3:
            #for n, acc in enumerate(mean_ser['number_of_TMDs']):
            list_of_central_TMDs = mean_ser['list_of_TMDs'][1:-1]
            list_mean_slope_central_TMDs = []
            for TMD in list_of_central_TMDs:
                list_mean_slope_central_TMDs.append(pd.to_numeric(mean_ser['%s_AAIMON_slope' % TMD]))
            mean_ser['AAIMON_slope_central_TMDs'] = np.mean(list_mean_slope_central_TMDs)
        else:
            mean_ser['AAIMON_slope_central_TMDs'] = np.nan

        # for interest sake, check the difference between AAIMON slopes for all TM residues combined, and the TMDs separately
        mean_ser["AAIMON_slope_all_TM_res_minus_all_TMDs"] = mean_ser["AAIMON_slope_all_TM_res"] - mean_ser["AAIMON_slope_all_TMDs_mean"]

        logging.info('%s AAIMON_slope all TM residues: %0.5f' % (acc, mean_ser["AAIMON_slope_all_TM_res"]))

        # save the pandas series with the means to a csv in the cr_ratios zip file
        mean_ser.to_csv(mean_ser_filename)
        zipout.write(mean_ser_filename, arcname=mean_ser_filename)
        os.remove(mean_ser_filename)

    # the TM_cr_pickle files remained stuck in the homol directory. Hopefully this removes the offending files.
    for TMD in list_of_TMDs:
        TM_cr_pickle = "{}_{}_cr_df.pickle".format(acc, TMD)
        try:
            os.remove(TM_cr_pickle)
        except:
            sys.stdout.write("{} could not be deleted.".format(TM_cr_pickle))


    return acc, True, "0"


def throw_out_truncated_sequences(pathdict, s, logging):
    """DEPRECATED TRUNCATION FUNCTION.
    A cutoff for min % nonTMD is now in the main script.
    :param pathdict: dict
        Dictionary of the key paths and files associated with that List number.
    :param s:
        Settings dictionary extracted from excel settings file.
    :param logging: logging.Logger
        Logger for printing to console and/or logfile.
    :return:
        returns an altered dataframe that does not contain the homologues that do not match the requirements for validated homologues
        due to truncatiion of nonTMD sequence during local aligning of query and match sequence
    """
    logging.info('~~~~~~~~~~~~      starting run_filter_truncated_alignments        ~~~~~~~~~~~~')
    # if multiprocessing is used, log only to the console
    p_dict_logging = logging if s["use_multiprocessing"] != True else utils.Log_Only_To_Console()
    # set current working directory as the data_dir/homol, where temp files will be saved before moving to zip
    os.chdir(os.path.join(s["data_dir"], "homol"))
    # get list of accessions that could not be downloaded, and can immediately be excluded
    not_in_homol_db = utils.get_acc_list_from_txt(pathdict["acc_not_in_homol_db_txt"])
    # create list of protein dictionaries to process
    list_p = korbinian.utils.convert_summary_csv_to_input_list(s, pathdict, p_dict_logging, list_excluded_acc=not_in_homol_db)
    # number of processes is the number the settings, or the number of proteins, whichever is smallest
    n_processes = s["multiprocessing_cores"] if s["multiprocessing_cores"] < len(list_p) else len(list_p)

    if s["use_multiprocessing"]:
        with Pool(processes=n_processes) as pool:
            calc_AAIMON_list = pool.map(korbinian.cons_ratio.cons_ratio.truncation_filter, list_p)
            # log the list of protein results (e.g. acc, "simap", True) to the actual logfile, not just the console
            logging.info("\ncalc_AAIMON_list : {}".format(calc_AAIMON_list))
    else:
        for p in list_p:
            korbinian.cons_ratio.cons_ratio.truncation_filter(p)
    logging.info("~~~~~~~~~~~~     run_filter_truncated_alignments is finished      ~~~~~~~~~~~~")


def truncation_filter(p):
    """DEPRECATED TRUNCATION FUNCTION.
    A cutoff for min % nonTMD is now in the main script.

    Parameters
    ----------
    p

    Returns
    -------

    """
    pathdict, s, logging = p["pathdict"], p["s"], p["logging"]
    acc = p["acc"]
    protein_name = p["protein_name"]
    uniprot_acc = p['uniprot_acc']
    homol_cr_ratios_zip = p['homol_cr_ratios_zip']
    list_of_TMDs = ast.literal_eval(p['list_of_TMDs'])
    min_perc_nonTMD = s['min_perc_nonTMD']

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
            #df_cr = pickle.load(zipfile.ZipFile(in_zipfile, "r", zipfile.ZIP_DEFLATED).open(in_file, "r"))
            df_cr = pd.read_pickle(zipfile.ZipFile(in_zipfile, "r", zipfile.ZIP_DEFLATED))
            sys.stdout.write('{}, {}: ' .format(uniprot_acc, TMD))
            # filtering step
            df_cr = utils.filter_for_truncated_sequences(min_perc_nonTMD, df_cr)
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
    """Function for linear slope equation

    Parameters
    ----------
    a : float
        Value for the slope
    x : float
        Value on the x-axis

    Returns
    -------
    y : float
        Y value for linear slope, based on x, where y-axis intercept is 1.0.
    """
    y = a * x + 1
    return y

def get_line_data_to_plot(a_constant, x_low=0, x_high=60):
    """Uses fitted function to create a line to plot

    Parameters
    ----------
    a_constant : float
        Constant A, fitted to data.
    x_low : float
        Lowest value of x to plot.
    x_high : float
        Highest value of x to plot.

    Returns
    -------
    x_data : np.ndarray
        array of lowest and highest x values
    y_data : np.ndarray
        array of lowest and highest y values
    """
    y_low = lin_AAIMON_slope_eq(a_constant, x_low)
    y_high = lin_AAIMON_slope_eq(a_constant, x_high)
    x_data = np.array([x_low, x_high])
    y_data = np.array([y_low, y_high])
    #AAIMON_at_80 = lin_AAIMON_slope_eq(a_constant, 80)
    return x_data, y_data

def fit_data_to_linear_function(x_array, y_array, a_constant_guess = 0.01):
    """Fits an array of x and y values to linear function with fixed y-intercept = 1.0

    Parameters
    ----------
    x_array : np.ndarray
        array of x values
    y_array : np.ndarray
        array of y values
    a_constant_guess : float
        Initial guess of A, used for leastsq function.

    Returns
    -------
    a_constant : float
        Fitted value of a to data (slope)
    x_data : np.ndarray
        array of lowest and highest x values
    y_data : np.ndarray
        array of lowest and highest y values
    """
    a_constant = leastsq(residuals, a_constant_guess, args = (lin_AAIMON_slope_eq, x_array, y_array))[0][0]
    x_data, y_data = get_line_data_to_plot(a_constant)
    a_constant = float(a_constant)
    return a_constant, x_data, y_data
