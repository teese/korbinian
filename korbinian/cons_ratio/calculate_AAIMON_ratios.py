import ast
import numpy as np
import csv
import os
import pickle
import korbinian
import korbinian.mtutils as utils
import zipfile
import pandas as pd

def run_calculate_AAIMON_ratios(pathdict, set_, logging):
    logging.info('~~~~~~~~~~~~starting run_calculate_AAIMON_ratios~~~~~~~~~~~~')
    df = pd.read_csv(pathdict["list_summary_csv"], sep=",", quoting=csv.QUOTE_NONNUMERIC, index_col=0)
    # iterate over the dataframe for proteins with an existing list_of_TMDs. acc = uniprot accession.
    for acc in df.loc[df['list_of_TMDs'].notnull()].loc[df['list_of_TMDs'] != 'nan'].index:
        if not os.path.exists(df.loc[acc, 'homol_df_orig_zip']):
            logging.info("{} Protein skipped. File does not exist".format(df.loc[acc, 'homol_df_orig_zip']))
            continue
        dfh = utils.open_df_from_pickle_zip(df.loc[acc, 'homol_df_orig_zip'],filename=os.path.basename(df.loc[acc, 'homol_df_orig_pickle']), delete_corrupt=True)
        if dfh.empty:
            logging.info("{} Protein skipped, file deleted as it is possibly corrupt.".format(df.loc[acc, 'homol_df_orig_zip']))
            continue

        list_of_TMDs = df.loc[acc, 'list_of_TMDs'].strip("[']").split(", ")

        ########################################################################################
        #                                                                                      #
        #                        Calculate AAIMON conservation ratios                          #
        #                                                                                      #
        ########################################################################################
        # determine if the dataframe already contains some previously analysed data
        dataframe_contains_prev_calc_AAIMON_ratios = True if 'TM01_AAIMON_ratio_mean' in df.columns else False
        logging.info('dataframe_contains_prev_calc_AAIMON_ratios = %s' % dataframe_contains_prev_calc_AAIMON_ratios)

        #assume af first that there is no previous data, and that the calculations can be re-run
        prev_calc_AAIMON_ratio_for_this_protein_exists = False
        if set_["overwrite_prev_calculated_AAIMON_ratios"] == False:
            if dataframe_contains_prev_calc_AAIMON_ratios == True:
                cell_is_empty = np.isnan(df.loc[acc, 'TM01_AAIMON_ratio_mean'])
                #if the cell is empty (np.nan), the AAIMON ratios haven't been calculated for that protein yet
                if cell_is_empty:
                    prev_calc_AAIMON_ratio_for_this_protein_exists = False
                else:
                    #if the data is there, it should be a float and not a np.nan, therefore there is no need to repeat all the calculations
                    prev_calc_AAIMON_ratio_for_this_protein_exists = True
                    logging.info('calculate_AAIMON_ratios skipped, prev_calc_AAIMON_ratio_for_this_protein_exists = True for %s' % df.loc[acc, 'uniprot_acc'])
                    # skip the rest of the analysis for this protein
                    continue
        # else:
        #     #if the settings says to overwrite the data, ignore its existence
        #     prev_calc_AAIMON_ratio_for_this_protein_exists = False
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
        df.loc[acc, 'FASTA_ident_mean'] = float('%0.2f' % dfh['FASTA_identity'].mean())
        # number of identical residues in FASTA alignment
        dfh['FASTA_num_ident_res'] = dfh['FASTA_identity'] / 100 * dfh['FASTA_overlap']
        df.loc[acc, 'FASTA_num_ident_res'] = float('%0.2f' % dfh['FASTA_identity'].mean())

        nonTMD_pickle_name = "{}_nonTMD_sliced.pickle".format(acc)
        df_nonTMD = utils.open_df_from_pickle_zip(df.loc[acc, 'fa_cr_sliced_TMDs_zip'], filename=nonTMD_pickle_name, delete_corrupt=True)
        if df_nonTMD.empty:
            logging.info("{}file was corrupt and was deleted".format(df.loc[acc, 'fa_cr_sliced_TMDs_zip']))
            continue

        # filter based on dfh above, for general homologue settings (e.g. % identity of full protein)
        df_nonTMD = df_nonTMD.loc[dfh.index, :]

        # calculate the nonTMD percentage identity and gaps
        df, df_nonTMD = korbinian.cons_ratio.calc_nonTMD_perc_ident_and_gaps(acc, df_nonTMD, df, logging)

        # save the nonTMD dataframe

        # filter nonTMD dataframe to only contain entries where nonTMD_perc_ident is not zero
        df_nonTMD = df_nonTMD.loc[df_nonTMD['nonTMD_perc_ident'] != 0]

        for TMD in list_of_TMDs:
            # open the dataframe containing the sequences, gap counts, etc for that TMD only
            df_cr = utils.open_df_from_pickle_zip(df.loc[acc, 'fa_cr_sliced_TMDs_zip'], filename="{}_{}_sliced.pickle".format(acc, TMD), delete_corrupt=True)
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
            df, df_cr = korbinian.cons_ratio.calc_AAIMON(acc, TMD, df_cr, df, logging)

            homol_cr_ratios_zip = df.loc[acc, 'homol_cr_ratios_zip']
            AAIMON_hist_path_prefix = df.loc[acc, 'AAIMON_hist_path_prefix']
            # save histograms for each TMD of that protein, with relative conservation
            korbinian.cons_ratio.save_hist_AAIMON_ratio_single_protein(df_cr, set_, list_of_TMDs, homol_cr_ratios_zip, AAIMON_hist_path_prefix)

        # remove columns to make output csv smaller
        if set_['drop_columns_to_reduce_csv_filesize']:
            list_cols_to_drop = ['match_align_seq', 'query_align_seq', 'align_markup_seq',
                                 'nonTMD_seq_query', 'nonTMD_markup']
            for col in list_cols_to_drop:
                if col in df_cr.columns:
                    df_cr.drop(col, axis=1, inplace=True)
        #df_cr.to_csv(df.loc[acc, 'SIMAP_csv_analysed_path'], sep=",", quoting=csv.QUOTE_NONNUMERIC)
        # save df_cr with homologues for a single protein, as a single zipped csv
        #utils.save_df_to_csv_zip(df_cr, df.loc[acc, 'homol_df_orig_zip'], open_method="w")
        print("important, work out a save method here" / 4)

        df.loc[acc, 'num_hits_with_SW_align_node'] = df_cr['hit_contains_SW_node'].value_counts()[True]
        logging.info('num_hits_with_SW_align_node: %s' % df.loc[acc, 'num_hits_with_SW_align_node'])

        df.loc[acc, "create_fasta_or_cons_ratio_single_protein"] = True
        # save to csv after each protein is analysed, incrementally adding the extra data
        df.to_csv(pathdict["list_summary_csv"], sep=",", quoting=csv.QUOTE_NONNUMERIC)

        logging.info("~~~~~run_calculate_AAIMON_ratios is finished ~~~~~")

