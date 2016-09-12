import ast
import numpy as np
import csv
import os
import pickle
import korbinian
import korbinian.mtutils as utils
import pandas as pd
import zipfile

def calculate_AAIMON_ratios(pathdict, set_, logging):
    logging.info('~~~~~~~~~~~~      starting run_calculate_AAIMON_ratios        ~~~~~~~~~~~~')
    df = pd.read_csv(pathdict["list_summary_csv"], sep=",", quoting=csv.QUOTE_NONNUMERIC, index_col=0)
    # set current working directory as the homol folder, where temp files will be saved before moving to zip
    os.chdir(set_["simap_database_dir"])
    # iterate over the dataframe for proteins with an existing list_of_TMDs. acc = uniprot accession.
    for acc in df.loc[df['list_of_TMDs'].notnull()].loc[df['list_of_TMDs'] != 'nan'].index:
        logging.info(df.loc[acc, 'protein_name'])
        if not os.path.exists(df.loc[acc, 'homol_df_orig_zip']):
            logging.info("{} Protein skipped. File does not exist".format(df.loc[acc, 'homol_df_orig_zip']))
            continue
        dfh = utils.open_df_from_pickle_zip(df.loc[acc, 'homol_df_orig_zip'],filename=os.path.basename(df.loc[acc, 'homol_df_orig_pickle']), delete_corrupt=True)
        if dfh.empty:
            logging.info("{} Protein skipped, file deleted as it is possibly corrupt.".format(df.loc[acc, 'homol_df_orig_zip']))
            continue

        # create an output series containing the mean values
        mean_ser = pd.Series()
        mean_ser["acc"] = acc
        mean_ser["protein_name"] = df.loc[acc, 'protein_name']
        mean_ser["organism"] = df.loc[acc, 'organism']
        mean_ser["prot_descr"] = df.loc[acc, 'prot_descr']
        mean_ser['list_of_TMDs'] = df.loc[acc, 'list_of_TMDs']
        list_of_TMDs = df.loc[acc, 'list_of_TMDs'].strip("[']").split(", ")

        ########################################################################################
        #                                                                                      #
        #                        Calculate AAIMON conservation ratios                          #
        #                                                                                      #
        ########################################################################################

        homol_cr_ratios_zip = df.loc[acc, 'homol_cr_ratios_zip']
        mean_ser_filename = "{}_cr_mean.csv".format(acc)

        #assume af first that there is no previous data, and that the calculations can be re-run
        prev_calc_AAIMON_ratio_for_this_protein_exists = False
        if set_["overwrite_prev_calculated_AAIMON_ratios"] == False:
            if os.path.isfile(homol_cr_ratios_zip):
                with zipfile.ZipFile(homol_cr_ratios_zip, mode="r", compression=zipfile.ZIP_DEFLATED) as zip:
                    if mean_ser_filename in zip.namelist():
                        logging.info('{} AAIMON_ratios skipped, file with mean AAIMON ratios exists (in settings, overwrite_prev_calculated_AAIMON_ratios = True)'.format(df.loc[acc, 'uniprot_acc']))
                        # if the means are saved as a csv in the homol cr_ratios zipfile, skip to the next protein
                        continue

        ########################################################################################
        #                                                                                      #
        #           Filter based on homol hit properties (non-TMD-specific)                    #
        #                                                                                      #
        ########################################################################################

        fa_X_filt_full_str = " and X_in_match_seq == False" if set_["fa_X_allowed_in_full_seq"] == False else ""

        fa_homol_query_str = 'FASTA_gapped_identity > {min_ident} and ' \
                            'FASTA_gapped_identity < {max_ident} and ' \
                            'hit_contains_SW_node == True and ' \
                            'disallowed_words_not_in_descr == True' \
                            '{Xfull}'.format(Xfull=fa_X_filt_full_str, min_ident=set_["cr_min_identity_of_full_protein"], max_ident=set_["cr_max_identity_of_full_protein"])

        # filter based on the query string
        dfh.query(fa_homol_query_str, inplace=True)

        '''Calculate average values, add to original dataframe.
           1) values associated with the FASTA output of SIMAP
        '''
        # fasta identity
        mean_ser['FASTA_ident_mean'] = float('%0.2f' % dfh['FASTA_identity'].mean())
        # number of identical residues in FASTA alignment
        dfh['FASTA_num_ident_res'] = dfh['FASTA_identity'] / 100 * dfh['FASTA_overlap']
        mean_ser['FASTA_num_ident_res'] = float('%0.2f' % dfh['FASTA_identity'].mean())

        nonTMD_pickle_name = "{}_nonTMD_sliced_df.pickle".format(acc)
        df_nonTMD = utils.open_df_from_pickle_zip(df.loc[acc, 'fa_cr_sliced_TMDs_zip'], filename=nonTMD_pickle_name, delete_corrupt=True)
        if df_nonTMD.empty:
            logging.info("{}file was corrupt and was deleted".format(df.loc[acc, 'fa_cr_sliced_TMDs_zip']))
            continue

        # filter based on dfh above, for general homologue settings (e.g. % identity of full protein)
        df_nonTMD = df_nonTMD.loc[dfh.index, :]

        # calculate the nonTMD percentage identity and gaps
        mean_ser, df_nonTMD = korbinian.cons_ratio.calc_nonTMD_perc_ident_and_gaps(acc, df_nonTMD, mean_ser, logging)

        with zipfile.ZipFile(homol_cr_ratios_zip, mode="w", compression=zipfile.ZIP_DEFLATED) as zipout:

            # save the nonTMD dataframe
            nonTMD_cr_outfile_pickle = "{}_nonTMD_cr_df.pickle".format(acc)
            with open(nonTMD_cr_outfile_pickle, "wb") as p:
                pickle.dump(df_nonTMD, p, protocol=pickle.HIGHEST_PROTOCOL)
            zipout.write(nonTMD_cr_outfile_pickle, arcname=nonTMD_cr_outfile_pickle)
            os.remove(nonTMD_cr_outfile_pickle)

            # filter nonTMD dataframe to only contain entries where nonTMD_perc_ident is not zero
            df_nonTMD = df_nonTMD.loc[df_nonTMD['nonTMD_perc_ident'] != 0]

            for TMD in list_of_TMDs:
                # open the dataframe containing the sequences, gap counts, etc for that TMD only
                df_cr = utils.open_df_from_pickle_zip(df.loc[acc, 'fa_cr_sliced_TMDs_zip'], filename="{}_{}_sliced_df.pickle".format(acc, TMD), delete_corrupt=True)
                # add the nonTMD percentage identity, etc. NOTE THAT THE INDEX SHOULD STILL BE COMPATIBLE, as the hit_num!
                df_cr['nonTMD_perc_ident'] = df_nonTMD['nonTMD_perc_ident']
                df_cr['nonTMD_perc_sim_plus_ident'] = df_nonTMD['nonTMD_perc_sim_plus_ident']
                df_cr['FASTA_overlap'] = dfh['FASTA_overlap']
                df_cr['len_full_match_seq'] = dfh[ 'len_full_match_seq']
                # filter based on dfh above, for general homologue settings (e.g. % identity of full protein), and df_nonTMD (for nonTMD_perc_ident is not zero, etc)
                df_cr = df_cr.loc[df_nonTMD.index,:]
                # following the general filters, filter to only analyse sequences with TMD identity above cutoff, and a nonTMD_perc_ident above zero ,to avoid a divide by zero error
                df_cr = df_cr.loc[df_cr['%s_perc_ident' % TMD] >= set_['cr_min_identity_of_TMD_initial_filter']]
                # calculate AAISMON etc for each TMD
                mean_ser, df_cr = korbinian.cons_ratio.calc_AAIMON(TMD, df_cr, mean_ser, logging)

                AAIMON_hist_path_prefix = df.loc[acc, 'AAIMON_hist_path_prefix']
                # save histograms for each TMD of that protein, with relative conservation
                korbinian.cons_ratio.save_hist_AAIMON_ratio_single_protein(df_cr, set_, list_of_TMDs, zipout, AAIMON_hist_path_prefix)

                logging.info('AAIMON MEAN %s: %0.2f' % (TMD, mean_ser['%s_AAIMON_ratio_mean' % TMD]))
                # logging.info('AASMON MEAN %s: %0.2f' % (TMD, mean_ser['%s_AASMON_ratio_mean'%TMD]))

                # save the dataframe for that TMD
                TM_cr_outfile_pickle = "{}_{}_cr_df.pickle".format(acc, TMD)
                with open(TM_cr_outfile_pickle, "wb") as p:
                    pickle.dump(df_cr, p, protocol=pickle.HIGHEST_PROTOCOL)
                zipout.write(TM_cr_outfile_pickle, arcname=TM_cr_outfile_pickle)
                os.remove(TM_cr_outfile_pickle)

            # # remove columns to make output csv smaller
            # if set_['drop_columns_to_reduce_csv_filesize']:
            #     list_cols_to_drop = ['match_align_seq', 'query_align_seq', 'align_markup_seq',
            #                          'nonTMD_seq_query', 'nonTMD_markup']
            #     for col in list_cols_to_drop:
            #         if col in df_cr.columns:
            #             df_cr.drop(col, axis=1, inplace=True)
            #df_cr.to_csv(df.loc[acc, 'SIMAP_csv_analysed_path'], sep=",", quoting=csv.QUOTE_NONNUMERIC)
            # save df_cr with homologues for a single protein, as a single zipped csv
            #utils.save_df_to_csv_zip(df_cr, mean_ser['homol_df_orig_zip'], open_method="w")

            mean_ser['num_hits_with_SW_align_node'] = dfh['hit_contains_SW_node'].value_counts()[True]
            #logging.info('num_hits_with_SW_align_node: %s' % mean_ser['num_hits_with_SW_align_node'])

            # save to csv after each protein is analysed, incrementally adding the extra data
            #df.to_csv(pathdict["list_summary_csv"], sep=",", quoting=csv.QUOTE_NONNUMERIC)
            # save the pandas series with the means to a csv in the cr_ratios zip file
            mean_ser.to_csv(mean_ser_filename)
            zipout.write(mean_ser_filename, arcname=mean_ser_filename)
            os.remove(mean_ser_filename)

    logging.info("~~~~~~~~~~~~     run_calculate_AAIMON_ratios is finished      ~~~~~~~~~~~~")

