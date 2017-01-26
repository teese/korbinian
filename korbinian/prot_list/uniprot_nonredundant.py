import os
import pandas as pd
import korbinian
import numpy as np
from Bio import SeqIO

def create_nonred_uniprot_flatfile_via_uniref(s, uniprot_dir_sel, list_number, selected_uniprot_records_flatfile, logging):
    """ Creates a non-redundant UniProt flatfile from redundant redundant UniProt tab files, redundant flatfiles and UniRef cluster tab file.

    The final output is the selected list of flatfiles, in the uniprot/selected folder (E.g. List08_selected_uniprot_records_flatfile.txt)

    Parameters
    ----------
    s : dict
        Settings dictionary extracted from excel settings file.
    uniprot_dir_sel : str
        Path to uniprot/selected folder.
    list_number : int
        List number (e.g. 8), determining the input and output files.
    selected_uniprot_records_flatfile : str
        Path to output UniProt flatfile of selected records. E.g. List08_selected_uniprot_records_flatfile.txt
    logging : logging.Logger
        Logger for printing to console and logfile.
    """
    logging.info('~~~~~~~~~~~~starting create_nonred_uniprot_flatfile_via_uniref~~~~~~~~~~~~')
    # load uniref cutoff used (typically 50, for UniRef50)
    uniref_cutoff = s["uniref_cluster_cutoff"]
    # define path to csv file containing the list of redundant uniprot accessions, e.g. List08_redundant_list_uniprot_acc.tab
    redundant_uniprot_acc_tab = os.path.join(uniprot_dir_sel, "List%02d_redundant_list_uniprot_acc.tab" % list_number)
    # define path to uniprot flatfile containing the redundant protein records, e.g. List08_redundant_uniprot_flatfile.txt
    redundant_uniprot_flatfile = os.path.join(uniprot_dir_sel, "List%02d_redundant_uniprot_flatfile.txt" % list_number)
    # define path to the csv file containing the relevant uniref clusters applicable to this list of proteins, e.g. List08_UniRef50_clusters.tab
    uniref_clusters_tab = os.path.join(uniprot_dir_sel, "List%02d_UniRef%02d_clusters.tab" % (list_number, uniref_cutoff))
    # output uniprot list with redundancy determined
    nonred_uniprot_acc_csv = os.path.join(uniprot_dir_sel, "List%02d_nonred_list_uniprot_acc.csv" % list_number)

    if os.path.isfile(selected_uniprot_records_flatfile) == False or s["overwrite_selected_ff"] == True:
        korbinian.prot_list.uniprot_nonredundant.match_list_uniprot_acc_to_uniref_clusters(redundant_uniprot_acc_tab, uniref_clusters_tab, nonred_uniprot_acc_csv, uniref_cutoff, logging)
        # reopen output file
        dfu = pd.read_csv(nonred_uniprot_acc_csv, index_col=0)
        # create a list of uniprot accessions that are nonredundant
        list_nonred_acc = list(dfu.loc[dfu['nonred'] == True].index)
        # create a uniprot flatfile containing only the desired nonredundant accessions
        korbinian.prot_list.uniprot_nonredundant.retrieve_selected_uniprot_records_from_flatfile(list_nonred_acc, redundant_uniprot_flatfile, selected_uniprot_records_flatfile, logging)
    logging.info('~~~~~~~~~~~~create_nonred_uniprot_flatfile_via_uniref is finished~~~~~~~~~~~~')

def match_list_uniprot_acc_to_uniref_clusters(redundant_uniprot_acc_tab, uniref_clusters_tab, nonred_uniprot_acc_csv, uniref_cutoff, logging):
    """ Assigns UniRef clusters to each acc in a list of UniProt accessions.

    Takes an input list of uniprot accessions (downloaded in .tab format)
    For each accession, finds the respective UniRef cluster number, within a UniRef cluster list (downloaded in .tab format)
    Creates a boolean for each acc, determining whether it is "nonredundant" according to the UniRef clusters.
         - preferentially labels the UniRef cluster representative as nonredundant
         - if the acc is not the cluster representative,

     for each finds the equivalent uniref cluster or each acc,
    and labels acc as redundant or non-redundant.

    Parameters:
    -----------
    redundant_uniprot_acc_tab : str
        Path to redundant uniprot list of accessions (csv downloaded in .tab format from UniProt website)
    uniref_clusters_tab : str
        Path to list of UniRef clusters (csv downloaded in .tab format from UniProt website)
    nonred_uniprot_acc_csv : str
        Path for output file, the uniprot list of accessions with annotated UniRef cluster, and redundancy.
    uniref_cutoff : int
        UniRef % identity cutoff used for that cluster (either 50, 90 or 100).
    logging : logging.Logger
        Logger for printing to console and logfile.
    """
    # create a new dataframe with the uniref csv file, containing the accessions of the reference sequences
    dfr = pd.read_table(uniref_clusters_tab)
    # to simplify, keep only the columns with the uniprot accessions
    # for backwards compatibility, replace old header 'Cluster member(s)' with a consistent column name
    dfr.rename(columns={'Cluster member(s)': 'cluster_members', 'Cluster members': 'cluster_members', 'Cluster ID': 'cluster_ID'}, inplace=True)
    # to simplify, keep only the columns with the uniprot accessions
    dfr = dfr[['cluster_ID', 'cluster_members']]# 'Organisms'
    # remove the isoforms, which for some reason are given separate UniRef numbers (SP|P45880, SP|P45880-1, SP|P45880-2 etc)
    dfr['cluster_ID'] = dfr['cluster_ID'].apply(lambda x: x[:-2] if "-" in x[-2:] else x)
    dfr = dfr.set_index('cluster_ID', drop=False)
    # change the cluster_ID to the index
    dfr.index.names = ['cluster_ID']
    # extract the representative uniprot accession from each cluster_ID
    dfr['cluster_rep'] = dfr['cluster_ID'].str[9:]
    # convert the list of uniprot accessions in each cluster to a python list format
    #dfr['all_acc_in_cluster'] = dfr['cluster_members'].apply(lambda x : [x.strip() for x in x.split(';')])
    # delete the original list of clusters
    #dfr.drop("cluster_members", axis=1, inplace=True)
    if os.path.isfile(redundant_uniprot_acc_tab) == False:
        logging.warning('warning, file with nonredundant uniprot acc does not exist : %s' % redundant_uniprot_acc_tab)
        raise FileNotFoundError()
    # open up large csv containing the uniprot acc of all single-pass proteins in uniprot
    dfu = pd.read_table(redundant_uniprot_acc_tab)
    # set the uniprot acc as the index
    dfu = dfu.set_index('Entry', drop=False)

    ##################################################################################################################
    #                                                                                                                #
    #           Use intersection to find accessions that are cluster representatives                                 #
    #                                                                                                                #
    ##################################################################################################################

    cluster_reps_set = set(dfr['cluster_rep'])
    accs_set = set(dfu.index)
    common_accs = cluster_reps_set.intersection(accs_set)
    dfu.loc[common_accs, "nonred"] = True
    dfu.loc[common_accs, "cluster_ID"] = dfu.loc[common_accs,:].Entry.apply(lambda x : "UniRef{}_{}".format(uniref_cutoff, x))

    common_refs_set = set(dfu.loc[common_accs, "cluster_ID"])
    # create a new dataframe with the uniref csv file, containing the accessions of the reference sequences
    # find the appropriate uniref cluster containing each acc
    # acc_counter = 0
    # for acc in dfu.index:
    #     # check if the accession is already a uniref representative (saves searching time!) e.g. UniRef50_P06129
    #     uniref_equivalent = 'UniRef%i_%s' % (uniref_cutoff, acc)
    #     if uniref_equivalent in dfr.index:
    #         dfu.loc[acc, 'cluster_ID'] = 'UniRef%i_%s' % (uniref_cutoff ,acc)
    #         dfu.loc[acc, 'nonred'] = True
    #     else:
    #         # if it is not a uniref representative, go through each uniref cluster, checking to see if it is a member
    #         for acc_cluster in list(dfr['all_acc_in_cluster']):
    #             if acc in acc_cluster:
    #                 # if the acc is a member of a uniref cluster, add the cluster name to the original dataframe
    #                 dfu.loc[acc, 'cluster_ID'] = 'UniRef%i_%s' % (uniref_cutoff, acc_cluster[0])
    #                 # stop the loop, as the cluster_ID has been found
    #                 break
    #     acc_counter += 1
        # write a dot on the screen for each protein, so that it is easy to see that the script is still working
        # sys.stdout.write(". ")
        # if acc_counter % 100 == 0:
        #     logging.info('%i records checked for redundancy' % acc_counter)

    ##################################################################################################################
    #                                                                                                                #
    #       For non representatives, search each cluster separately to see if it contains the accession              #
    #                                                                                                                #
    ##################################################################################################################
    for n, acc in enumerate(dfu.loc[dfu.nonred.isnull()].index):
        cluster_found = False
        if n % 50 == 0:
            if n != 0:
                logging.info('%i records checked for redundancy' % n)
        # if it is not a uniref representative, go through each uniref cluster, checking to see if it is a member
        for cluster_ID in dfr.index:
            cluster_members = dfr.loc[cluster_ID, 'cluster_members']
            if acc in cluster_members:
                # if the acc is a member of a uniref cluster, add the cluster name to the original dataframe
                dfu.loc[acc, 'cluster_ID'] = cluster_ID
                cluster_found = True
                # stop the loop, as the cluster_ID has been found
                break
        if not cluster_found:
            # if the accession is not found in any of the cluster members, mark it as "not_found"
            dfu.loc[acc, 'cluster_ID'] = "not_found"
            dfu.loc[acc, 'nonred'] = False
    logging.info('%i records checked for redundancy' % dfu.loc[dfu.nonred.isnull()].shape[0])

    # sort the dataframe based cluster_ID
    dfu.sort_values(["cluster_ID"], inplace=True)

    ##################################################################################################################
    #                                                                                                                #
    #          Create a set of cluster_ID that are unique, and have not yet been assigned to any acc                 #
    #                                                                                                                #
    ##################################################################################################################
    # cluster IDs marked as redundant (so far)
    cluster_IDs_marked_redu_array = np.array(dfu[dfu["nonred"].isnull()]['cluster_ID'])
    # set of unique cluster IDs marked as redundant (so far)
    cluster_IDs_marked_redu_unique_array = set(cluster_IDs_marked_redu_array)
    # cluster IDs marked as redundant (so far), minus those already found due to presence of reference acc in list
    cluster_IDs_marked_redu_unique_array_excl_common_refs_set = cluster_IDs_marked_redu_unique_array - common_refs_set
    # create a list of the acc still assumed to be redundant
    dfu_red = dfu.loc[dfu.nonred.isnull()]
    # set the cluster ID as an index
    dfu_red.set_index("cluster_ID", inplace=True)
    ##################################################################################################################
    #                                                                                                                #
    #                  Get the first accession for the list of acc matching each cluster_ID                          #
    #                                                                                                                #
    ##################################################################################################################
    acc_nonred_list = []
    # selecting by the unique cluster ID will give either a single row (Series), or several rows (DataFrame)
    for uniq_CID in cluster_IDs_marked_redu_unique_array_excl_common_refs_set:
        series_or_df = dfu_red.loc[uniq_CID, :]
        if "DataFrame" in str(type(series_or_df)):
            acc = series_or_df.Entry.iloc[0]
        elif "Series" in str(type(series_or_df)):
            acc = series_or_df.Entry
        else:
            raise TypeError()
    acc_nonred_list.append(acc)
    # label all the accessions in the nonredundant list as True
    dfu.loc[acc_nonred_list, "nonred"] = True
    # convert all NaNs to "False"
    dfu["nonred"] = dfu["nonred"].fillna(False)

    ##################################################################################################################
    #                                                                                                                #
    #                                      save list and print record numbers                                        #
    #                                                                                                                #
    ##################################################################################################################
    # save list after redundancy check
    dfu.to_csv(nonred_uniprot_acc_csv)
    # log the number of redundant an nonredundant accessions
    nonred_value_counts = dfu['nonred'].value_counts()
    number_nonredundant_records = nonred_value_counts[True]
    number_total_records = dfu.shape[0]
    number_redundant_records = number_total_records - number_nonredundant_records
    if "not_found" in dfu["cluster_ID"].tolist():
        number_acc_not_found_in_UniRef = dfu["cluster_ID"].value_counts()["not_found"]
        logging.info("{} acc not found in UniRef lists".format(number_acc_not_found_in_UniRef))

    logging.info('number_total_records = {0}, number_redundant_records removed = {1}, '
                 'final number_nonredundant_records = {2}'.format(number_total_records,
                                                                  number_redundant_records,
                                                                  number_nonredundant_records))
    logging.info("~~~~~~~~~~~~match_list_uniprot_acc_to_uniref_clusters is finished~~~~~~~~~~~~")

def retrieve_selected_uniprot_records_from_flatfile(input_accession_list, large_input_uniprot_flatfile, output_flatfile, logging):
    '''Function to select records from a large uniprot flatfile, and save them as a smaller flatfile of selected records.

    Parameters:
    -----------
    input_accession_list : list
        List of UniProt accessions.
    large_input_uniprot_flatfile : str
        Path to large UniProt flatfile, from which the UniProt flatfiles will be extracted.
    output_flatfile : str
        Path to output UniProt flatfile of selected records.
    logging : logging.Logger
        Logger for printing to console and logfile.

    Notes:
    ------
    This script is quite slow for large lists, and large initial flatfiles.
    '''
    #accession_list = [line.strip() for line in open(input_accession_list, "r")]
    #open input flatfile
    uniprot_index_handle = SeqIO.index(large_input_uniprot_flatfile, "swiss")
    #create an empty list to hold all the accessions that are not in the flatfile
    list_acc_not_in_flatfile = []
    with open(output_flatfile, "wb") as output:
        for n, acc in enumerate(input_accession_list):
            try:
                #get the raw uniprot record, and write to the output file
                output.write(uniprot_index_handle.get_raw(acc))
            #if the record is not in the file
            except KeyError:
                list_acc_not_in_flatfile.append(acc)
            if n % 50 == 0:
                if n!= 0:
                    logging.info('%i records retrieved from large flatfile' % n)
    logging.info('%i-%i records retrieved from large flatfile' % (n, len(list_acc_not_in_flatfile)))
    if len(list_acc_not_in_flatfile) > 0:
        logging.info("SwissProt records not found in %s:\n%s." % (large_input_uniprot_flatfile, list_acc_not_in_flatfile))
    logging.info("~~~~~~~~~~~~retrieve_selected_uniprot_records_from_flatfile is finished~~~~~~~~~~~~")