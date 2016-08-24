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
    logging.info('~~~~~~~~~~~~starting slice_TMDs_from_homologues~~~~~~~~~~~~')
    df = pd.read_csv(pathdict["list_summary_csv"], sep=",", quoting=csv.QUOTE_NONNUMERIC, index_col=0)
    # if "uniprot_acc" in df.columns:
    #     df.set_index("uniprot_acc", drop=False, inplace=True)
    # else:
    #     df["uniprot_acc"] = df.index
    #iterate over the dataframe for proteins with an existing list_of_TMDs. Note that acc = uniprot accession here.
    for acc in df.loc[df['list_of_TMDs'].notnull()].loc[df['list_of_TMDs'] != 'nan'].index:
        # if set_["overwrite_prev_calculated_AAIMON_ratios"] == True or prev_calc_AAIMON_ratio_for_this_protein_exists == False:
        protein_name = df.loc[acc, 'protein_name']
        logging.info('%s' % protein_name)

        homo_orig_table_zip_exists = False
        #FIX! problem with corrupted files. Best to check and delete existing if overwrite is true.
        if os.path.exists(df.loc[acc, 'homol_df_orig_zip']):
            #try:
                # with zipfile.ZipFile(df.loc[acc, 'homol_df_orig_zip'], "r", zipfile.ZIP_DEFLATED) as openzip:
                #     if openzip.namelist()[0][-4:] == ".csv":
                #
            dfs = utils.open_df_from_pickle_zip(df.loc[acc, 'homol_df_orig_zip'], delete_corrupt=True)
            if dfs is not None:
                homo_orig_table_zip_exists = True
            # except:
            #     #file may be corrupted, if script stopped unexpectedly before compression was finished
            #     logging.info('%s seems to be corrupted. File will be deleted.'
            #                  'Will need to be parsed again next time program is run' % df.loc[acc, 'homol_df_orig_zip'])
            #     #os.remove(df.loc[acc, 'homol_df_orig_zip'])
        else:
            logging.info("file does not exist : df.loc[acc, 'homol_df_orig_zip'] : {}".format(df.loc[acc, 'homol_df_orig_zip']))

        if homo_orig_table_zip_exists:
            # run the analysis function, which slices homologues, calculates identity, and saves output files
            #create_fasta_or_cons_ratio_single_protein(df.loc[acc, 'homol_df_orig_zip'], acc, protein_name, set_, df, pathdict, logging)

            #def create_fasta_or_cons_ratio_single_protein(homol_df_orig_zip, acc, protein_name, set_, df, pathdict, logging):
            # get the dataframe pickled and zipped. Note that if there are two pickle files, it should be named!
            # the index should still be the hit_num
            #dfs = utils.open_df_from_pickle_zip(df.loc[acc, 'homol_df_orig_zip'])

            # get length of seq. Previously this was a lambda function.
            if True in dfs['hit_contains_SW_node'].tolist():
                dfs['len_query_alignment_sequence'] = dfs['query_alignment_sequence'].str.len()
            else:
                dfs['query_alignment_sequence'] = np.nan
                dfs['len_query_alignment_sequence'] = np.nan

            #     try:
            #         dfs['len_query_alignment_sequence'] = dfs['query_alignment_sequence'].str.len()
            #     except KeyError:
            #         pass
                    # dataframe does not contain query_alignment_sequence,
                    # which means that the XML file is probably damaged somehow
            #         logging.warning('SIMAP_csv_from_XML seems to be damaged for %s' % protein_name)
            #         dfs['query_alignment_sequence'] = ''
            #         dfs['len_query_alignment_sequence'] = 0

            list_of_TMDs = ast.literal_eval(df.loc[acc, 'list_of_TMDs'])
            # create counter for number of TMDs with some homologue data
            n_TMDs_w_homol = 0

            fa_cr_sliced_TMDs_zip = df.loc[acc, 'fa_cr_sliced_TMDs_zip']
            # delete any existing sliced zipfile
            if os.path.isfile(fa_cr_sliced_TMDs_zip):
                os.remove(fa_cr_sliced_TMDs_zip)
            # open new zipfile (NOTE, it must be closed later!!)
            homol_sliced_zip = zipfile.ZipFile(fa_cr_sliced_TMDs_zip, mode="a", compression=zipfile.ZIP_DEFLATED)
            if set_["run_create_fasta"]:
                if os.path.isfile(df.loc[acc, 'fa_TMDs_zip']):
                    os.remove(df.loc[acc, 'fa_TMDs_zip'])
                # open new zipfile (NOTE, it must be closed later!!)
                #fa_TMDs_zip = zipfile.ZipFile(set_["fa_TMDs_zip"], mode="a", compression=zipfile.ZIP_DEFLATED)

            # get directory for zip (and other temp files to be transferred)
            homol_dir = os.path.dirname(fa_cr_sliced_TMDs_zip)
            for TMD in list_of_TMDs:
                df_TMD = korbinian.cons_ratio.slice_homologues_and_count_gaps(acc, TMD, df, dfs, set_, logging, n_TMDs_w_homol)
                temp_pickle = os.path.join(homol_dir, "{}_{}_sliced.pickle".format(acc, TMD))
                with open(temp_pickle, "wb") as p:
                    pickle.dump(df_TMD, p)
                homol_sliced_zip.write(temp_pickle, arcname=os.path.basename(temp_pickle))
                os.remove(temp_pickle)
            homol_sliced_zip.close()

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
            # query_seqs_ser = dfs_with_SW_node['query_alignment_sequence']
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
            #     df_TMD = korbinian.cons_ratio.slice_homologues_and_count_gaps(TMD, len_query_TMD, df_TMD, set_["cr_max_n_gaps_in_query_TMD"], set_["cr_max_n_gaps_in_match_TMD"])

def juxta_function_1(dfs, TMD):
    dfs['start_juxta_after_%s'%TMD] = np.where(utils.isNaN(dfs['TM%.2d_start_in_SW_alignment'%(int(TMD[2:])+1)])==True,np.nan,dfs['%s_end_in_SW_alignment'%TMD])
    dfs['end_juxta_before_%s'%TMD] = np.where(dfs["%s_start_in_SW_alignment"%TMD]!=0,dfs["%s_start_in_SW_alignment"%TMD],np.nan)
    dfs['end_juxta_after_%s'%TMD] = dfs["%s_end_in_SW_alignment"%TMD]+((dfs["TM%.2d_start_in_SW_alignment"%(int(TMD[2:])+1)]-dfs["%s_end_in_SW_alignment"%TMD])/2).apply(lambda x :int(x) if not np.isnan(x) else np.nan)
    dfs['start_juxta_before_%s'%TMD] = np.where(dfs["end_juxta_after_TM%.2d"%(int(TMD[2:])-1)] == dfs['end_juxta_before_%s'%TMD] ,dfs["end_juxta_after_TM%.2d"%(int(TMD[2:])-1)],dfs["end_juxta_after_TM%.2d"%(int(TMD[2:])-1)])
    return dfs
