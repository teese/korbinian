import ast
import numpy as np
import csv
import os
import pickle
import korbinian
import korbinian.mtutils as utils
import zipfile
import pandas as pd

def slice_TMDs_from_homologues(pathdict, set_, logging):
    logging.info('~~~~~~~~~~~~       starting slice_TMDs_from_homologues        ~~~~~~~~~~~~')
    df = pd.read_csv(pathdict["list_summary_csv"], sep=",", quoting=csv.QUOTE_NONNUMERIC, index_col=0)
    #iterate over the dataframe for proteins with an existing list_of_TMDs. Note that acc = uniprot accession here.
    for acc in df.loc[df['list_of_TMDs'].notnull()].loc[df['list_of_TMDs'] != 'nan'].index:
        logging.info(df.loc[acc, 'protein_name'])

        if not os.path.exists(df.loc[acc, 'homol_df_orig_zip']):
            logging.info("{} Protein skipped. File does not exist".format(df.loc[acc, 'homol_df_orig_zip']))
            continue

        dfs = utils.open_df_from_pickle_zip(df.loc[acc, 'homol_df_orig_zip'], delete_corrupt=True)
        if dfs.empty:
            logging.info("{} Protein skipped, file deleted as it is possibly corrupt.".format(df.loc[acc, 'homol_df_orig_zip']))
            continue

        if set_["slice_juxtamembrane_regions"]:
            ########################################################################################
            #                                                                                      #
            #        Define juxtamembrane regions associated with each TMD  [AAIMON]               #
            #                                                                                      #
            ########################################################################################
            # convert the tuple of (True, 32, 53) into separate dataframes.
            # http://stackoverflow.com/questions/29550414/how-to-split-column-of-tuples-in-pandas-dataframe

            # dfs_with_match = dfs['%s_start_end_list_in_SW_alignment' % TMD].dropna()
            # df_match_start_end = pd.DataFrame(dfs_with_match.values.tolist())
            # df_match_start_end.index = dfs_with_match.index
            # dfs["%s_start_in_SW_alignment" % TMD] = df_match_start_end[1]
            # dfs["%s_end_in_SW_alignment" % TMD] = df_match_start_end[2]

            last_TMD_of_acc = list_of_TMDs[-1]

            if set_["slice_juxtamembrane_regions"] == True:
                if TMD == "TM01":
                    dfs['start_juxta_before_TM01'] = np.where(dfs['TM01_start_in_SW_alignment'] > 0, 0, np.nan)
                    dfs['end_juxta_before_TM01'] = np.where(dfs['TM01_start_in_SW_alignment'] == 0, np.nan,
                                                            dfs['TM01_start_in_SW_alignment'])
                    dfs['start_juxta_after_TM01'] = dfs['TM01_end_in_SW_alignment']
                    if len(list_of_TMDs) == 1:
                        # if there is only one TMD, TM01 == last_TMD_of_acc
                        dfs['end_juxta_after_TM01'] = np.where(utils.isNaN(dfs['start_juxta_after_TM01']) == True, np.nan,
                                                               dfs['len_query_align_seq'])
                    elif len(list_of_TMDs) > 1:
                        # problem('dfs["TM02_start_in_SW_alignment"] cannot exist yet, because the script iterates through the TMDs one at a time')
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
                        dfs['start_juxta_before_%s' % TMD] = dfs['end_juxta_after_TM%.2d' % (int(TMD[2:]) - 1)]
                        dfs['end_juxta_before_%s' % TMD] = dfs['%s_start_in_SW_alignment' % TMD]
                        dfs['start_juxta_after_%s' % TMD] = np.where(
                            dfs['%s_end_in_SW_alignment' % TMD] == dfs['len_query_align_seq'], np.nan,
                            dfs['%s_end_in_SW_alignment' % TMD])
                        dfs['end_juxta_after_%s' % TMD] = np.where(utils.isNaN(dfs['start_juxta_after_%s' % TMD]) == True, np.nan,
                                                                   dfs['len_query_align_seq'])
                        # dfs['seq_juxta_after_%s_in_query'%TMD] = dfs[dfs['start_juxta_after_%s'%TMD].notnull()].apply(utils.slice_juxta_after_TMD_in_query, args = (TMD,), axis=1)
                        # dfs['seq_juxta_after_%s_in_query'%TMD] = dfs.query_align_seq[int(dfs['start_juxta_after_TM10']):int(dfs['end_juxta_after_TM10'])]
                        # dfs['seq_juxta_after_%s_in_match'%TMD] =

                last_TMD_of_acc = list_of_TMDs[-1]
                dfs['seq_juxta_before_%s_in_query' % TMD] = dfs[dfs['start_juxta_before_%s' % TMD].notnull()].apply(
                    utils.slice_juxta_before_TMD_in_query, args=(TMD,), axis=1)
                dfs['seq_juxta_before_%s_in_match' % TMD] = dfs[dfs['start_juxta_before_%s' % TMD].notnull()].apply(
                    utils.slice_juxta_before_TMD_in_match, args=(TMD,), axis=1)
                if not TMD == last_TMD_of_acc:
                    dfs['seq_juxta_after_%s_in_query' % TMD] = dfs[dfs['end_juxta_after_%s' % TMD].notnull()].apply(
                        utils.slice_juxta_after_TMD_in_query, args=(TMD,), axis=1)
                    dfs['seq_juxta_after_%s_in_match' % TMD] = dfs[dfs['end_juxta_after_%s' % TMD].notnull()].apply(
                        utils.slice_juxta_after_TMD_in_match, args=(TMD,), axis=1)
                else:
                    dfs['seq_juxta_after_%s_in_query' % TMD] = np.nan
                    dfs['seq_juxta_after_%s_in_match' % TMD] = np.nan
                    for hit in dfs.index:
                        if not utils.isNaN(dfs['start_juxta_after_%s' % TMD])[hit]:
                            # altered to .loc rather than ['seq_juxta_after_%s_in_match'%TMD][hit] after SettingWithCopyWarning
                            dfs.loc[hit, 'seq_juxta_after_%s_in_match' % TMD] = dfs.match_align_seq[hit][
                                                                                int(dfs.loc[hit, "start_juxta_after_%s" % TMD]):int(
                                                                                    dfs.loc[hit, "end_juxta_after_%s" % TMD])]
                            dfs.loc[hit, 'seq_juxta_after_%s_in_query' % TMD] = dfs.query_align_seq[hit][
                                                                                int(dfs.loc[hit, "start_juxta_after_%s" % TMD]):int(
                                                                                    dfs.loc[hit, "end_juxta_after_%s" % TMD])]




        # # get length of seq. Previously this was a lambda function.
        # if True in dfs['hit_contains_SW_node'].tolist():
        #     dfs['len_query_align_seq'] = dfs['query_align_seq'].str.len()
        # else:
        #     dfs['query_align_seq'] = np.nan
        #     dfs['len_query_align_seq'] = np.nan

        #     try:
        #         dfs['len_query_align_seq'] = dfs['query_align_seq'].str.len()
        #     except KeyError:
        #         pass
                # dataframe does not contain query_align_seq,
                # which means that the XML file is probably damaged somehow
        #         logging.warning('SIMAP_csv_from_XML seems to be damaged for %s' % protein_name)
        #         dfs['query_align_seq'] = ''
        #         dfs['len_query_align_seq'] = 0

        list_of_TMDs = ast.literal_eval(df.loc[acc, 'list_of_TMDs'])
        # create counter for number of TMDs with some homologue data
        n_TMDs_w_homol = 0

        fa_cr_sliced_TMDs_zip = df.loc[acc, 'fa_cr_sliced_TMDs_zip']
        # delete any existing sliced zipfile
        if os.path.isfile(fa_cr_sliced_TMDs_zip):
            os.remove(fa_cr_sliced_TMDs_zip)
        # open new zipfile (NOTE, it must be closed later!!)
        with zipfile.ZipFile(df.loc[acc, 'fa_cr_sliced_TMDs_zip'], mode="a", compression=zipfile.ZIP_DEFLATED) as homol_sliced_zip:

            # get directory for zip (and other temp files to be transferred)
            homol_dir = os.path.dirname(fa_cr_sliced_TMDs_zip)
            # create a specific dataframe to hold the nonTMD region, including indices (True, start, end) of all the TMD segments
            # add the FASTA_gapped_identity and length of the alignment sequence from dfs, to act as the "end" of all the nonTMD regions
            df_nonTMD_sliced = dfs[['len_query_align_seq']].copy()
            for TMD in list_of_TMDs:
                df_TMD = korbinian.cons_ratio.slice_TMD_homol_and_count_gaps(acc, TMD, df, dfs, set_, logging, n_TMDs_w_homol)

                # transfer the columns with indices across to the df_nonTMD_sliced
                cols = ['%s_in_SW_alignment' % TMD, '%s_start_in_SW_alignment' % TMD, '%s_end_in_SW_alignment' % TMD]
                for col in cols:
                    df_nonTMD_sliced[col] = df_TMD[col]

                TM_temp_pickle = os.path.join(homol_dir, "{}_{}_sliced_df.pickle".format(acc, TMD))
                with open(TM_temp_pickle, "wb") as p:
                    pickle.dump(df_TMD, p, protocol=pickle.HIGHEST_PROTOCOL)
                homol_sliced_zip.write(TM_temp_pickle, arcname=os.path.basename(TM_temp_pickle))
                os.remove(TM_temp_pickle)
                #korbinian.cons_ratio.slice_nonTMD_seqs(dfs, df_nonTMD_sliced, list_of_TMDs)

            df_nonTMD_sliced = korbinian.cons_ratio.slice_nonTMD_seqs(dfs, df_nonTMD_sliced, list_of_TMDs)

            df_nonTMD_temp_pickle = os.path.join(homol_dir, "{}_nonTMD_sliced_df.pickle".format(acc))
            with open(df_nonTMD_temp_pickle, "wb") as p:
                pickle.dump(df_nonTMD_sliced, p, protocol=pickle.HIGHEST_PROTOCOL)
            homol_sliced_zip.write(df_nonTMD_temp_pickle, arcname=os.path.basename(df_nonTMD_temp_pickle))
            os.remove(df_nonTMD_temp_pickle)

    logging.info("~~~~~~~~~~~~     slice_TMDs_from_homologues is finished       ~~~~~~~~~~~~")

        # ########################################################################################
        # #                                                                                      #
        # #                Slice out TMD regions [fasta and AAIMON]                              #
        # #                                                                                      #
        # ########################################################################################
        #
        # # create a new dataframe to hold the sliced homol sequences for that TMD, and number of gaps, etc
        # df_TMD = pd.DataFrame()
        #
        # '''slice the TMD regions out of the alignment markup and match sequences
        # the method first finds the indices of the TMD region in the query, and uses these indices to slice
        # the filters are not applied yet, so the sequences can be viewed in the csv file for analysis
        # '''
        # # create regex string for each TMD
        # query_TMD_sequence = df.loc[acc, '%s_seq'%TMD]
        # # create TMD regex search string (e.g. L.*L.*M.*L.*L.*M.*L.*L.*M.*L.*L.*M.*L.*L.*M.*)
        # TMD_regex_ss = utils.create_regex_string(query_TMD_sequence)
        # # select to only include data where the XML contained a SW node, and then apply function for a regex search
        # dfs_with_SW_node = dfs.query('hit_contains_SW_node == True')
        # # create series of query_alignment_seqs
        # query_seqs_ser = dfs_with_SW_node['query_align_seq']
        #
        # # obtain the bool, start, end of TMD seqs in the match sequences. Add to the new TMD-specific dataframe.
        # df_TMD['%s_start_end_list_in_SW_alignment'%TMD] = query_seqs_ser.apply(utils.get_start_and_end_of_TMD_in_query,args=(TMD_regex_ss,))
        # '''the output of the regex search is a list with three components.
        # 1) a match (e.g. True),
        # 2) the start of the match (e.g. 124),
        # 3) the end of the match (e.g. 144)
        # '''
        # # for convenience, separate the components into separate columns
        # columns_from_regex_output = ['%s_in_SW_alignment'%TMD, '%s_start_in_SW_alignment'%TMD,'%s_end_in_SW_alignment'%TMD]
        # # n = the index (1,2,3) in the tuple
        # # col = column item in the list (e.g. 'TM09_in_SW_alignment')
        # for n, col in enumerate(columns_from_regex_output):
        #     # add a new column which is named TM01_start, etc, and insert the appropriate integer (start or stop) or bool from the tuple
        #     df_TMD[col] = df_TMD['%s_start_end_list_in_SW_alignment'%TMD].dropna().apply(lambda x: x[n])
        # # drop the original listlike from the regex search
        # df_TMD.drop('%s_start_end_list_in_SW_alignment' % TMD, inplace=True)
        #
        # ########################################################################################
        # #                                                                                      #
        # #          Apply function to slice homologues and count gaps [fasta and AAIMON]        #
        # #                                                                                      #
        # ########################################################################################
        #
        # # in some cases, there is no data to obtain as the hit_contains_SW_node = False for too many sequences, giving no start_in_SW_alignment
        # number_of_rows_containing_data = df_TMD[df_TMD['%s_start_end_list_in_SW_alignment' % TMD].notnull()].shape[0]
        # if number_of_rows_containing_data == 0:
        #     logging.info('%s does not have any valid homologues for %s. '
        #                  'Re-downloading simap homologue XML may be necessary.' % (protein_name, TMD))
        # if number_of_rows_containing_data != 0:
        #     n_TMDs_w_homol += 1
        #     len_query_TMD = len(df.loc[acc, '%s_seq' % TMD])
        #     # apply the slicing function to the homologues
        #     df_TMD = korbinian.cons_ratio.slice_TMD_homol_and_count_gaps(TMD, len_query_TMD, df_TMD, set_["cr_max_n_gaps_in_query_TMD"], set_["cr_max_n_gaps_in_match_TMD"])



def juxta_function_1(dfs, TMD):
    dfs['start_juxta_after_%s'%TMD] = np.where(utils.isNaN(dfs['TM%.2d_start_in_SW_alignment'%(int(TMD[2:])+1)])==True,np.nan,dfs['%s_end_in_SW_alignment'%TMD])
    dfs['end_juxta_before_%s'%TMD] = np.where(dfs["%s_start_in_SW_alignment"%TMD]!=0,dfs["%s_start_in_SW_alignment"%TMD],np.nan)
    dfs['end_juxta_after_%s'%TMD] = dfs["%s_end_in_SW_alignment"%TMD]+((dfs["TM%.2d_start_in_SW_alignment"%(int(TMD[2:])+1)]-dfs["%s_end_in_SW_alignment"%TMD])/2).apply(lambda x :int(x) if not np.isnan(x) else np.nan)
    dfs['start_juxta_before_%s'%TMD] = np.where(dfs["end_juxta_after_TM%.2d"%(int(TMD[2:])-1)] == dfs['end_juxta_before_%s'%TMD] ,dfs["end_juxta_after_TM%.2d"%(int(TMD[2:])-1)],dfs["end_juxta_after_TM%.2d"%(int(TMD[2:])-1)])
    return dfs