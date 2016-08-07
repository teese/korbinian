import ast
import pandas as pd
import numpy as np
import csv
import os
import korbinian
import korbinian.mtutils as utils
import zipfile

def create_fasta_or_calculate_AAIMON_ratios(pathdict, set_, logging):
    logging.info('~~~~~~~~~~~~starting create_fasta_or_calculate_AAIMON_ratios~~~~~~~~~~~~')
    overwrite_prev_calculated_AAIMON_ratios = set_["overwrite_prev_calculated_AAIMON_ratios"]
    df = pd.read_csv(pathdict["list_summary_csv"])

    #filter to remove sequences where no TMDs are found,
    df = df.loc[df['list_of_TMDs'].notnull()]
    # #filter to remove sequences where no TMDs are found (if string)
    # df.loc[df['list_of_TMDs'] != 'nan']
    # #filter to remove sequences where no TMDs are found (if string)
    # df = df.loc[df['list_of_TMDs'] != 'nan']
    #determine if the dataframe already contains some previously analysed data
    dataframe_contains_prev_calc_AAIMON_ratios = True if 'TM01_AAIMON_ratio_mean' in df.columns else False
    logging.info('dataframe_contains_prev_calc_AAIMON_ratios = %s' % dataframe_contains_prev_calc_AAIMON_ratios)
    #iterate over the dataframe for proteins with an existing list_of_TMDs. Note that acc = uniprot accession here.
    for acc in df.loc[df['list_of_TMDs'] != 'nan'].index:
        #assume af first that there is no previous data, and that the calculations can be re-run
        prev_calc_AAIMON_ratio_for_this_protein_exists = False
        if overwrite_prev_calculated_AAIMON_ratios == False:
            if dataframe_contains_prev_calc_AAIMON_ratios == True:
                cell_is_empty = np.isnan(df.loc[acc, 'TM01_AAIMON_ratio_mean'])
                #if the cell is empty (np.nan), the AAIMON ratios haven't been calculated for that protein yet
                if cell_is_empty:
                    prev_calc_AAIMON_ratio_for_this_protein_exists = False
                else:
                    #if the data is there, it should be a float and not a np.nan, therefore there is no need to repeat all the calculations
                    prev_calc_AAIMON_ratio_for_this_protein_exists = True
                    logging.info('calculate_AAIMON_ratios skipped, prev_calc_AAIMON_ratio_for_this_protein_exists = True for %s' % df.loc[acc, 'uniprot_acc'])
        else:
            #if the settings says to overwrite the data, ignore its existence
            prev_calc_AAIMON_ratio_for_this_protein_exists = False
        if overwrite_prev_calculated_AAIMON_ratios == True or prev_calc_AAIMON_ratio_for_this_protein_exists == False:
            protein_name = df.loc[acc, 'protein_name']
            logging.info('%s' % protein_name)
            # SIMAP_csv_from_XML_tarfile = df.loc[acc, 'SIMAP_csv_from_XML_tarfile']
            # if os.path.exists(SIMAP_csv_from_XML_tarfile):
            #     try:
            #         with tarfile.open(df.loc[acc, 'SIMAP_csv_from_XML_tarfile'], mode='r:gz') as tar:
            #             if df.loc[acc, 'SIMAP_csv_from_XML'] in [tarinfo.name for tarinfo in tar]:
            #                 SIMAP_csv_from_XML_exists = True
            #             else:
            #                 SIMAP_csv_from_XML_exists = False
            #     except (EOFError, tarfile.ReadError, OSError):
            #         #file may be corrupted, if script stopped unexpectedly before compression was finished
            #         logging.info('%s seems to be corrupted. File will be deleted.'
            #                      'Will need to be parsed again next time program is run' %
            #                       df.loc[acc, 'SIMAP_csv_from_XML_tarfile'])
            #         SIMAP_csv_from_XML_exists = False
            #         os.remove(df.loc[acc, 'SIMAP_csv_from_XML_tarfile'])
            # else:
            #     SIMAP_csv_from_XML_exists = False

            SIMAP_csv_from_XML_exists = False
            if os.path.exists(df.loc[acc, 'homol_csv_zip']):
                try:
                    with zipfile.ZipFile(df.loc[acc, 'homol_csv_zip'], "r", zipfile.ZIP_DEFLATED) as openzip:
                        if openzip.namelist()[0][-4:] == ".csv":
                            SIMAP_csv_from_XML_exists = True
                except:
                    #file may be corrupted, if script stopped unexpectedly before compression was finished
                    logging.info('%s seems to be corrupted. File will be deleted.'
                                 'Will need to be parsed again next time program is run' % df.loc[acc, 'homol_csv_zip'])
            if SIMAP_csv_from_XML_exists:
                # run the analysis function, which slices homologues, calculates identity, and saves output files
                analyse_homologues_single_protein(df.loc[acc, 'homol_csv_zip'], acc, protein_name, set_, df, pathdict, logging)

    logging.info('calculate_AAIMON_ratios is finished.')

def juxta_function_1(dfs, TMD):
    dfs['start_juxta_after_%s'%TMD] = np.where(utils.isNaN(dfs['TM%.2d_start_in_SW_alignment'%(int(TMD[2:])+1)])==True,np.nan,dfs['%s_end_in_SW_alignment'%TMD])
    dfs['end_juxta_before_%s'%TMD] = np.where(dfs["%s_start_in_SW_alignment"%TMD]!=0,dfs["%s_start_in_SW_alignment"%TMD],np.nan)
    dfs['end_juxta_after_%s'%TMD] = dfs["%s_end_in_SW_alignment"%TMD]+((dfs["TM%.2d_start_in_SW_alignment"%(int(TMD[2:])+1)]-dfs["%s_end_in_SW_alignment"%TMD])/2).apply(lambda x :int(x) if not np.isnan(x) else np.nan)
    dfs['start_juxta_before_%s'%TMD] = np.where(dfs["end_juxta_after_TM%.2d"%(int(TMD[2:])-1)] == dfs['end_juxta_before_%s'%TMD] ,dfs["end_juxta_after_TM%.2d"%(int(TMD[2:])-1)],dfs["end_juxta_after_TM%.2d"%(int(TMD[2:])-1)])
    return dfs

def analyse_homologues_single_protein(homol_csv_zip, acc, protein_name, set_, df, pathdict, logging):
    # with tarfile.open(SIMAP_csv_from_XML_tarfile, 'r:gz') as tar_in:
    #SIMAP_csv_from_XML_extracted = tar_in.extractfile(df.loc[acc, 'SIMAP_csv_from_XML'])
    dfs = utils.open_df_from_csv_zip(homol_csv_zip)
    # reopen the csv file containing all the homologue data for that particular protein as a pandas dataframe (labelled Data Frame SIMAP, or dfs)
    #dfs = pd.read_csv(SIMAP_csv_from_XML_extracted, sep=",", index_col=0, quoting=csv.QUOTE_NONNUMERIC)
    # create a column with the length of the Smith Waterman alignment (same length for query, markup and match)
    if True in dfs['hit_contains_SW_node']:
        at_least_one_hit_contains_SW_node = True
    else:
        at_least_one_hit_contains_SW_node = False
    if at_least_one_hit_contains_SW_node:
        try:
            dfs['len_query_alignment_sequence'] = dfs['query_alignment_sequence'].str.len()
        except KeyError:
            # dataframe does not contain query_alignment_sequence,
            # which means that the XML file is probably damaged somehow
            logging.warning('SIMAP_csv_from_XML seems to be damaged for %s' % protein_name)
            dfs['query_alignment_sequence'] = ''
            dfs['len_query_alignment_sequence'] = 0
    else:
        dfs['query_alignment_sequence'] = ''
        dfs['len_query_alignment_sequence'] = 0

    list_of_TMDs = ast.literal_eval(df.loc[acc, 'list_of_TMDs'])
    number_of_TMDs_containing_some_homologue_data = 0

    for TMD in list_of_TMDs:
        ########################################################################################
        #                                                                                      #
        #                Slice out TMD regions [fasta and AAIMON]                              #
        #                                                                                      #
        ########################################################################################
        '''slice the TMD regions out of the alignment markup and match sequences
        the method first finds the indices of the TMD region in the query, and uses these indices to slice
        the filters are not applied yet, so the sequences can be viewed in the csv file for analysis
        '''
        # create regex string for each TMD
        query_TMD_sequence = df.loc[acc, '%s_seq'%TMD]
        # create TMD regex search string (e.g. L.*L.*M.*L.*L.*M.*L.*L.*M.*L.*L.*M.*L.*L.*M.*)
        TMD_regex_ss = utils.create_regex_string(query_TMD_sequence)
        # select to only include data where the XML contained a SW node, and then apply function for a regex search
        query_seqs_ser = dfs.query('hit_contains_SW_node == True')['query_alignment_sequence']

        dfs['%s_start_end_tuple_in_SW_alignment'%TMD] = query_seqs_ser.apply(utils.get_start_and_end_of_TMD_in_query,args=(TMD_regex_ss,))
        '''the output of the regex search is a tuple with three components.
        1) a match (e.g. True),
        2) the start of the match (e.g. 124),
        3) the end of the match (e.g. 144)
        '''
        # for convenience, separate the components into separate columns
        columns_from_regex_output = ['%s_in_SW_alignment'%TMD, '%s_start_in_SW_alignment'%TMD,'%s_end_in_SW_alignment'%TMD]
        # n = the index (1,2,3) in the tuple
        # col = column item in the list (e.g. 'TM09_in_SW_alignment')
        for n, col in enumerate(columns_from_regex_output):
            # first filter to analyse only columns that contain a SW node
            df_match = dfs.query('hit_contains_SW_node == True')
            # add a new column which is named TM01_start, etc, and insert the appropriate integer (start or stop) or bool from the tuple
            dfs[col] = df_match['%s_start_end_tuple_in_SW_alignment'%TMD].apply(lambda x: x[n])

        ########################################################################################
        #                                                                                      #
        #          Apply function to slice homologues and count gaps [fasta and AAIMON]        #
        #                                                                                      #
        ########################################################################################

        # in some cases, there is no data to obtain as the hit_contains_SW_node = False for too many sequences, giving no start_in_SW_alignment
        number_of_rows_containing_data = dfs[dfs['%s_start_in_SW_alignment' % TMD].notnull()].shape[0]
        if number_of_rows_containing_data == 0:
            logging.info('%s does not have any valid homologues for %s. '
                         'Re-downloading simap homologue XML may be necessary.' % (protein_name, TMD))
        if number_of_rows_containing_data != 0:
            number_of_TMDs_containing_some_homologue_data += 1
            # apply the slicing function to the homologues
            dfs = korbinian.cons_ratio.slice_homologues_and_count_gaps(acc, TMD, df, dfs, set_)


    if set_["run_create_fasta"]:
        ########################################################################################
        #                                                                                      #
        #           Filter and create FastA files with TMDs from homologues                    #
        #                                                                                      #
        ########################################################################################
        dfs = korbinian.fasta.filter_and_save_fasta(df, dfs, acc, TMD, set_, logging)

    if set_["run_calculate_AAIMON_ratios"]:
        ########################################################################################
        #                                                                                      #
        #                        Calculate AAIMON conservation ratios                          #
        #                                                                                      #
        ########################################################################################

        """Find disallowed words (e.g. synthetic, patent) [AAIMON] """

        # add the list of words to the globals, to be accessed by utils.find_disallowed_words
        cr_words_not_allowed_in_description = ast.literal_eval(set_["cr_words_not_allowed_in_description"])
        # collect disallowed words in hit protein description (patent, synthetic, etc)
        dfs['cr_list_disallowed_words_in_descr'] = dfs['uniprot_description'].dropna().apply(utils.find_disallowed_words, args=(cr_words_not_allowed_in_description,))
        # create a boolean column to select hits that do not contain these words in the description
        dfs['cr_disallowed_words_not_in_descr'] = dfs['cr_list_disallowed_words_in_descr'] == '[]'

        # check if there are non-IUPAC amino acids in the sequence (frequently large gaps from NG sequencing data)
        dfs['X_in_match_seq'] = 'X' in dfs['match_alignment_sequence']

        for TMD in list_of_TMDs:
            ########################################################################################
            #                                                                                      #
            #        Define juxtamembrane regions associated with each TMD  [AAIMON]               #
            #                                                                                      #
            ########################################################################################
            # convert the tuple of (True, 32, 53) into separate dataframes.
            # http://stackoverflow.com/questions/29550414/how-to-split-column-of-tuples-in-pandas-dataframe

            dfs_with_match = dfs['%s_start_end_tuple_in_SW_alignment' % TMD].dropna()
            utils.aaa(dfs_with_match)
            df_match_start_end = pd.DataFrame(dfs_with_match.values.tolist(), low_memory=False)
            df_match_start_end.index = dfs_with_match.index
            dfs["%s_start_in_SW_alignment" % TMD] = df_match_start_end[1]
            dfs["%s_end_in_SW_alignment" % TMD] = df_match_start_end[2]

            last_TMD_of_acc = list_of_TMDs[-1]

            if TMD == "TM01":
                dfs['start_juxta_before_TM01'] = np.where(dfs['TM01_start_in_SW_alignment'] > 0, 0, np.nan)
                dfs['end_juxta_before_TM01'] = np.where(dfs['TM01_start_in_SW_alignment'] == 0, np.nan,
                                                        dfs['TM01_start_in_SW_alignment'])
                dfs['start_juxta_after_TM01'] = dfs['TM01_end_in_SW_alignment']
                if len(list_of_TMDs) == 1:
                    # if there is only one TMD, TM01 == last_TMD_of_acc
                    dfs['end_juxta_after_TM01'] = np.where(utils.isNaN(dfs['start_juxta_after_TM01']) == True, np.nan,
                                                           dfs['len_query_alignment_sequence'])
                elif len(list_of_TMDs) > 1:
                    #problem('dfs["TM02_start_in_SW_alignment"] cannot exist yet, because the script iterates through the TMDs one at a time')
                    dfs['end_juxta_after_TM01'] = dfs["TM01_end_in_SW_alignment"] + (
                    (dfs["TM02_start_in_SW_alignment"] - dfs["TM01_end_in_SW_alignment"]) / 2).apply(lambda x: int(x) if not np.isnan(x) else np.nan)
                    # dfs['seq_juxta_after_TM01_in_query'] = dfs[dfs['start_juxta_after_TM01'].notnull()].apply(utils.slice_juxta_after_TMD_in_query, args = (TMD,), axis=1)
                    # dfs['seq_juxta_after_TM01_in_match'] = dfs[dfs['end_juxta_after_TM01'].notnull()].apply(utils.slice_juxta_after_TMD_in_match, args = (TMD,), axis=1)

            # the analysis is slow, so don't repeat TM01 if there is only one TM helix in the protein
            if len(list_of_TMDs) > 1:
                if not TMD == "TM01" and not TMD == last_TMD_of_acc:
                    dfs = juxta_function_1(dfs, TMD)
                    # dfs['start_juxta_after_%s'%TMD] = np.where(utils.isNaN(dfs['TM%.2d_start_in_SW_alignment'%(int(TMD[2:])+1)])==True,np.nan,dfs['%s_end_in_SW_alignment'%TMD])
                    # dfs['end_juxta_before_%s'%TMD] = np.where(dfs["%s_start_in_SW_alignment"%TMD]!=0,dfs["%s_start_in_SW_alignment"%TMD],np.nan)
                    # dfs['end_juxta_after_%s'%TMD] = dfs["%s_end_in_SW_alignment"%TMD]+((dfs["TM%.2d_start_in_SW_alignment"%(int(TMD[2:])+1)]-dfs["%s_end_in_SW_alignment"%TMD])/2).apply(lambda x :int(x) if not np.isnan(x) else np.nan)
                    # dfs['start_juxta_before_%s'%TMD] = np.where(dfs["end_juxta_after_TM%.2d"%(int(TMD[2:])-1)] == dfs['end_juxta_before_%s'%TMD] ,dfs["end_juxta_after_TM%.2d"%(int(TMD[2:])-1)],dfs["end_juxta_after_TM%.2d"%(int(TMD[2:])-1)])
                    # dfs['seq_juxta_after_%s_in_query'%TMD] = dfs[dfs['start_juxta_after_%s'%TMD].notnull()].apply(utils.slice_juxta_after_TMD_in_query, args = (TMD,), axis=1)
                    # dfs['seq_juxta_after_%s_in_match'%TMD] = dfs[dfs['end_juxta_after_%s'%TMD].notnull()].apply(utils.slice_juxta_after_TMD_in_match, args = (TMD,), axis=1)

                if TMD == last_TMD_of_acc:
                    dfs['start_juxta_before_%s'%TMD] = dfs['end_juxta_after_TM%.2d' % (int(TMD[2:]) - 1)]
                    dfs['end_juxta_before_%s'%TMD] = dfs['%s_start_in_SW_alignment'%TMD]
                    dfs['start_juxta_after_%s'%TMD] = np.where(
                        dfs['%s_end_in_SW_alignment'%TMD] == dfs['len_query_alignment_sequence'], np.nan,
                        dfs['%s_end_in_SW_alignment'%TMD])
                    dfs['end_juxta_after_%s'%TMD] = np.where(utils.isNaN(dfs['start_juxta_after_%s'%TMD]) == True, np.nan,
                                                               dfs['len_query_alignment_sequence'])
                    # dfs['seq_juxta_after_%s_in_query'%TMD] = dfs[dfs['start_juxta_after_%s'%TMD].notnull()].apply(utils.slice_juxta_after_TMD_in_query, args = (TMD,), axis=1)
                    # dfs['seq_juxta_after_%s_in_query'%TMD] = dfs.query_alignment_sequence[int(dfs['start_juxta_after_TM10']):int(dfs['end_juxta_after_TM10'])]
                    # dfs['seq_juxta_after_%s_in_match'%TMD] =


            last_TMD_of_acc = list_of_TMDs[-1]
            dfs['seq_juxta_before_%s_in_query'%TMD] = dfs[dfs['start_juxta_before_%s'%TMD].notnull()].apply(
                utils.slice_juxta_before_TMD_in_query, args=(TMD,), axis=1)
            dfs['seq_juxta_before_%s_in_match'%TMD] = dfs[dfs['start_juxta_before_%s'%TMD].notnull()].apply(
                utils.slice_juxta_before_TMD_in_match, args=(TMD,), axis=1)
            if not TMD == last_TMD_of_acc:
                dfs['seq_juxta_after_%s_in_query'%TMD] = dfs[dfs['end_juxta_after_%s'%TMD].notnull()].apply(
                    utils.slice_juxta_after_TMD_in_query, args=(TMD,), axis=1)
                dfs['seq_juxta_after_%s_in_match'%TMD] = dfs[dfs['end_juxta_after_%s'%TMD].notnull()].apply(
                    utils.slice_juxta_after_TMD_in_match, args=(TMD,), axis=1)
            else:
                dfs['seq_juxta_after_%s_in_query'%TMD] = np.nan
                dfs['seq_juxta_after_%s_in_match'%TMD] = np.nan
                for hit in dfs.index:
                    if not utils.isNaN(dfs['start_juxta_after_%s'%TMD])[hit]:
                        # altered to .loc rather than ['seq_juxta_after_%s_in_match'%TMD][hit] after SettingWithCopyWarning
                        dfs.loc[hit, 'seq_juxta_after_%s_in_match'%TMD] = dfs.match_alignment_sequence[hit][
                                                                            int(dfs.loc[hit, "start_juxta_after_%s"%TMD]):int(
                                                                                dfs.loc[hit, "end_juxta_after_%s"%TMD])]
                        dfs.loc[hit, 'seq_juxta_after_%s_in_query'%TMD] = dfs.query_alignment_sequence[hit][
                                                                            int(dfs.loc[hit, "start_juxta_after_%s"%TMD]):int(
                                                                                dfs.loc[hit, "end_juxta_after_%s"%TMD])]

        if number_of_TMDs_containing_some_homologue_data == 0:
            # there was no data obtained from the csv file, which probably means the original XML file was not properly downloaded
            # there is no need to continue the script for this protein
            logging.warning('%s does not contain any valid homologues at all. '
                            'CSV and/or XML file is damaged. Re-downloading simap homologue XML may be necessary. '
                            'AAIMON ratios will not be calculated for this protein.' % protein_name)

        if number_of_TMDs_containing_some_homologue_data > 0:
            # create a boolean column that describel whether the sequence is above the minimum gapped identity
            cr_minimum_identity_of_full_protein = set_["cr_minimum_identity_of_full_protein"]
            dfs['cr_gapped_ident_above_cutoff'] = dfs['FASTA_gapped_identity'] > cr_minimum_identity_of_full_protein

            '''
            3) values associated each TMD, such as average AAIMON ratio
            '''
            # calculate the nonTMD percentage identity and gaps
            korbinian.cons_ratio.calc_nonTMD_perc_ident_and_gaps(acc, dfs, set_, df, list_of_TMDs, logging)

            # calculate AAISMON etc for each TMD
            for TMD in list_of_TMDs:
                df, dfs = korbinian.cons_ratio.calc_AAIMON(acc, TMD, dfs, set_, df, logging)

            homol_cons_ratio_zip = df.loc[acc, 'homol_cons_ratio_zip']
            AAIMON_hist_path_prefix = df.loc[acc, 'AAIMON_hist_path_prefix']
            # save histograms for each TMD of that protein, with relative conservation
            korbinian.cons_ratio.save_hist_AAIMON_ratio_single_protein(dfs, set_, list_of_TMDs, homol_cons_ratio_zip, AAIMON_hist_path_prefix)

        # remove columns to make output csv smaller
        if set_['drop_columns_to_reduce_csv_filesize']:
            list_cols_to_drop = ['match_alignment_sequence', 'query_alignment_sequence', 'alignment_markup',
                                 'nonTMD_seq_query', 'nonTMD_markup']
            for col in list_cols_to_drop:
                if col in dfs.columns:
                    dfs.drop(col, axis=1, inplace=True)
        #dfs.to_csv(df.loc[acc, 'SIMAP_csv_analysed_path'], sep=",", quoting=csv.QUOTE_NONNUMERIC)
        # save dfs with homologues for a single protein, as a single zipped csv
        utils.save_df_to_csv_zip(dfs, df.loc[acc, 'homol_csv_zip'], open_method="w")

        df.loc[acc, 'num_hits_with_SW_align_node'] = dfs['hit_contains_SW_node'].value_counts()[True]
        logging.info('num_hits_with_SW_align_node: %s' % df.loc[acc, 'num_hits_with_SW_align_node'])

    df.loc[acc, "analyse_homologues_single_protein"] = True
    # save to csv after each protein is analysed, incrementally adding the extra data
    df.to_csv(pathdict["list_summary_csv"], sep=",", quoting=csv.QUOTE_NONNUMERIC)

