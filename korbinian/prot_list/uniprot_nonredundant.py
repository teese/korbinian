import os
import pandas as pd
import korbinian
from Bio import SeqIO
# import debugging tools
from korbinian.utils import pr, pc, pn, aaa

def create_nonred_uniprot_flatfile_via_uniref(s, uniprot_dir, selected_uniprot_records_flatfile, logging):
    """ Creates a non-redundant UniProt flatfile from redundant redundant UniProt tab files, redundant flatfiles and UniRef cluster tab file.

    The final output is the selected list of flatfiles, in the uniprot/selected folder (E.g. List08_selected_uniprot_records_flatfile.txt)

    Parameters
    ----------
    s : dict
        Settings dictionary extracted from excel settings file.
    uniprot_dir : str
        Path to databases/uniprot folder.
    list_number : int
        List number (e.g. 8), determining the input and output files.
    selected_uniprot_records_flatfile : str
        Path to output UniProt flatfile of selected records. E.g. List08_selected_uniprot_records_flatfile.txt
    logging : logging.Logger
        Logger for printing to console and logfile.

    Saved Files and Figures
    -----------------------
    selected_uniprot_records_flatfile : flatfile
        UniProt flatfile of non-redundant sequences.

    """
    logging.info('~~~~~~~~~~~~starting create_nonred_uniprot_flatfile_via_uniref~~~~~~~~~~~~')
    # load uniref cutoff used (typically 50, for UniRef50)
    uniref_cutoff = s["uniref_cluster_cutoff"]
    # define path to csv file containing the list of redundant uniprot accessions, e.g. List08_redundant_list_uniprot_acc.tab
    redundant_uniprot_acc_tab = os.path.join(uniprot_dir, "List%02d_redundant_list_uniprot_acc.tab" % s["list_number"])
    # define path to uniprot flatfile containing the redundant protein records, e.g. List08_redundant_uniprot_flatfile.txt
    redundant_uniprot_flatfile = os.path.join(uniprot_dir, "List%02d_redundant_uniprot_flatfile.txt" % s["list_number"])
    # define path to the csv file containing the relevant uniref clusters applicable to this list of proteins, e.g. List08_UniRef50_clusters.tab
    uniref_clusters_tab = os.path.join(uniprot_dir, "List%02d_UniRef%02d_clusters.tab" % (s["list_number"], uniref_cutoff))
    # output uniprot list with redundancy determined
    nonred_uniprot_acc_csv = os.path.join(uniprot_dir, "List%02d_nonred_list_uniprot_acc.csv" % s["list_number"])

    result = korbinian.prot_list.uniprot_nonredundant.match_list_uniprot_acc_to_uniref_clusters(redundant_uniprot_acc_tab, uniref_clusters_tab, nonred_uniprot_acc_csv, uniref_cutoff, logging)
    if result is None:
        logging.warning('~~~~~~~~~~~~create_nonred_uniprot_flatfile_via_uniref NOT DONE~~~~~~~~~~~~')
        return None
    # reopen output file
    dfu = pd.read_csv(nonred_uniprot_acc_csv, index_col=0)
    # create a list of uniprot accessions that are nonredundant
    list_nonred_acc = list(dfu.loc[dfu['nonred'] == True].index)
    # create a uniprot flatfile containing only the desired nonredundant accessions
    korbinian.prot_list.uniprot_nonredundant.retrieve_selected_uniprot_records_from_flatfile(list_nonred_acc, redundant_uniprot_flatfile, selected_uniprot_records_flatfile, logging)
    logging.info('~~~~~~~~~~~~create_nonred_uniprot_flatfile_via_uniref is finished~~~~~~~~~~~~')
    return True

def match_list_uniprot_acc_to_uniref_clusters(redundant_uniprot_acc_tab, uniref_clusters_tab, nonred_uniprot_acc_csv, uniref_cutoff, logging):
    """ Assigns UniRef clusters to each acc in a list of UniProt accessions.

    Takes an input list of redundand UniProt accessions (downloaded in .tab format)
    For each accession, finds the respective UniRef cluster number, within a UniRef cluster list (downloaded in .tab format)
        [note that to speed things up, we've iteratively removed clusters as they are identified, resulting in
        redundant proteins whose UniRef representative is not found]

    Creates a boolean for each acc, determining whether it is "nonredundant" (nonred) according to the UniRef clusters.
         - preferentially labels the UniRef cluster representative as nonredundant
         - if the acc is not the cluster representative, preferentially takes the "reviewed" SwissProt recods

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

    Saved Files and Figures
    -----------------------
    nonred_uniprot_acc_csv : csv
        csv with the identified cluster_ID for each protein, and whether or not it is classified as redundant
        e.g. List01_nonred_list_uniprot_acc.csv
        index = UniProt accession
        columns = Protein names	Gene names	Organism	Length	nonred	cluster_ID
    """
    if not os.path.isfile(uniref_clusters_tab):
        logging.warning("REDUNDANCY REDUCTION NOT POSSIBLE, no uniref clusters found ({})".format(uniref_clusters_tab))
        return None

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
    # open up large csv containing the redundant list of uniprot acc
    dfu = pd.read_table(redundant_uniprot_acc_tab)
    # set the uniprot acc as the index
    dfu = dfu.set_index('Entry', drop=False)
    n_initial_records = dfu.shape[0]

    ##################################################################################################################
    #                                                                                                                #
    #           Use intersection to find accessions that are cluster representatives                                 #
    #                                                                                                                #
    ##################################################################################################################
    # create a set of cluster representatives (e.g. {'UniRef50_Q7Z7M0', 'UniRef50_Q96I36', 'UniRef50_P20333',....}
    cluster_reps_set = set(dfr['cluster_rep'])
    # create a set of all the proteins in orig redundant list
    red_protein_acc_set = set(dfu.index)
    # use intersection to find the acc that are cluster representatives
    cluster_reps_acc_set = cluster_reps_set.intersection(red_protein_acc_set)
    # annotate the cluster representatives as non-redundant, and add the cluster ID
    dfu.loc[cluster_reps_acc_set, "nonred"] = True
    dfu.loc[cluster_reps_acc_set, "cluster_ID"] = dfu.loc[cluster_reps_acc_set,:].Entry.apply(lambda x : "UniRef{}_{}".format(uniref_cutoff, x))
    # collect a set of the cluster IDs corresponding to the representatives in the list
    reps_cluster_ID_set = set(dfu.loc[cluster_reps_acc_set, "cluster_ID"])
    n_prot_that_are_ref_seqs = len(reps_cluster_ID_set)

    ##################################################################################################################
    #                                                                                                                #
    #       For non representatives, search each cluster separately to see if it contains the accession              #
    #                                                                                                                #
    ##################################################################################################################

    # OPTIONAL: to make the check faster, drop any of the UniRef clusters whose reference protein was in the non-redundant list
    # NOTE: this means that there will be less clusters to search, BUT, you won't find the cluster for proteins who belong
    # to a cluster, where the representative protein was in the original list. These will be labelled "not found".
    unassigned = set(dfr.index) - reps_cluster_ID_set
    dfr_unassigned = dfr.loc[unassigned, :].copy()
    n_clusters = dfr.shape[0]
    logging.info("Number of initial redundant proteins : {}\n"
                 "Number of uniref clusters : {}\n"
                 "Number of proteins that were ref seqs for a cluster : {}".format(n_initial_records, n_clusters, n_prot_that_are_ref_seqs))

    ##################################################################################################################
    #                                                                                                                #
    #            FOR EACH CLUSTER REPRESENTATIVE, GRAB THE INTERSECTION OF CLUSTER ACC and LIST ACC                  #
    #                                                                                                                #
    ##################################################################################################################

    # split into accessions and convert to a set
    dfr_unassigned["cluster_members"] = dfr_unassigned.cluster_members.str.split("; ")
    dfr_unassigned["cluster_members"] = dfr_unassigned["cluster_members"].apply(lambda x : set(x))

    # create a set of the accessions that still need to be assigned
    set_red_acc = set(dfu.loc[dfu.nonred.isnull()].index)

    # make a new column, showing the intersection between the list of acc in that cluster, and the list of acc of unassigned proteins
    dfr_unassigned["acc_intersection_cluster_and_unassigned"] = dfr_unassigned["cluster_members"].apply(lambda x: x.intersection(set_red_acc))
    # drop any UniRef clusters that don't contain a reference to the remaining unassigned accessions
    dfr_unassigned = dfr_unassigned.loc[dfr_unassigned["acc_intersection_cluster_and_unassigned"] != set()]


    ##################################################################################################################
    #                                                                                                                #
    #       Reverse the dataset so the ACC is the index, and the ClusterIDs are the values in a new dataframe        #
    #                                                                                                                #
    ##################################################################################################################

    # create a new dataframe to
    df_acc_to_cluster = pd.DataFrame()
    # iterate through the unassigned cluster IDs
    for cluster_ID in dfr_unassigned.index:
        # grab the set of acc for that cluster, which are also in the protein list
        acc_intersection_cluster_and_unassigned = dfr_unassigned.loc[cluster_ID, "acc_intersection_cluster_and_unassigned"]
        # for some reason, there are sometimes redundant UniRef clusters!!
        if isinstance(acc_intersection_cluster_and_unassigned, pd.Series):
            # take the first list as the correct one
            acc_intersection_cluster_and_unassigned = acc_intersection_cluster_and_unassigned.iloc[0]
        # skip the empty sets
        if acc_intersection_cluster_and_unassigned != set():
            counter = 0
            # the multiple acc all have the same cluster ID, which can be assigned in the new dataframe
            for acc in acc_intersection_cluster_and_unassigned:
                df_acc_to_cluster.loc[acc, "cluster_ID"] = cluster_ID
            # IF THERE ARE MORE THAN 2, TRY TO TAKE THE "REVIEWED" PROTEIN FROM SWISSPROT
            if len(acc_intersection_cluster_and_unassigned) >= 2:
                if "Status" in dfu.columns:
                    # first make a series whether they are "reviewed" or "unreviewed"
                    reviewed_ser = dfu.loc[acc_intersection_cluster_and_unassigned, "Status"]
                    # sort the series so that the "reviewed" proteins are at the top
                    reviewed_ser.sort_values(inplace=True)
                    for acc in reviewed_ser.index:
                        if counter == 0:
                            # take the first protein, which should be "reviewed", if available
                            df_acc_to_cluster.loc[acc, "nonred"] = True
                        else:
                            # label all following proteins as redundant, to be excluded from the final list
                            df_acc_to_cluster.loc[acc, "nonred"] = False
                        counter += 1
                else:
                    # there is no reviewed or unreviewed status in the cluster data
                    # simply take the first acc as the non-redundant, and mark the others as redundant
                    for acc in acc_intersection_cluster_and_unassigned:
                        if counter == 0:
                            # take the first protein, which should be "reviewed", if available
                            df_acc_to_cluster.loc[acc, "nonred"] = True
                        else:
                            # label all following proteins as redundant, to be excluded from the final list
                            df_acc_to_cluster.loc[acc, "nonred"] = False
                        counter += 1

    ##################################################################################################################
    #                                                                                                                #
    #    Join the "nonred" and "cluster_ID" series from the two methods:                                             #
    #        1) cluster representatives                                                                              #
    #        2) taking the first (or reviewed) from the non-cluster representatives                                  #
    #                                                                                                                #
    ##################################################################################################################
    # series 1
    nonred_from_reps = dfu["nonred"].dropna()
    # series 2
    if 'nonred' in df_acc_to_cluster.columns:
        nonred_from_cluster_search = df_acc_to_cluster["nonred"].dropna()
    else:
        nonred_from_cluster_search = df_acc_to_cluster
    # join series
    complete_series_nonred_annotation = nonred_from_reps.append(nonred_from_cluster_search)
    # add to output file
    dfu["nonred"] = complete_series_nonred_annotation

    # series 1
    cluster_IDs_from_reps_ser = dfu["cluster_ID"].dropna()
    # series 2, join series
    complete_series_cluster_IDs = cluster_IDs_from_reps_ser.append(df_acc_to_cluster["cluster_ID"].dropna())
    # add to output file
    dfu["cluster_ID"] = complete_series_cluster_IDs

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
    #number_total_records = dfu.shape[0]
    number_redundant_records = n_initial_records - number_nonredundant_records
    if "not_found" in dfu["cluster_ID"].tolist():
        number_acc_not_found_in_UniRef = dfu["cluster_ID"].value_counts()["not_found"]
        acc_not_found_list = list(dfu["cluster_ID"].isnull().index)
        logging.info("{} acc not found in UniRef clusters\n{}".format(number_acc_not_found_in_UniRef, acc_not_found_list))

    logging.info('final number of non-redundant records : {}\n({} redundant removed)'.format(number_nonredundant_records, number_redundant_records))

    logging.info("~~~~~~~~~~~~match_list_uniprot_acc_to_uniref_clusters is finished~~~~~~~~~~~~")
    return True

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
    n = 0
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
    logging.info('{} records retrieved from large flatfile ({} not found)'.format(n + 1, len(list_acc_not_in_flatfile)))
    if len(list_acc_not_in_flatfile) > 0:
        logging.info("SwissProt records not found in {}:\n{}.".format(large_input_uniprot_flatfile, list_acc_not_in_flatfile))
    logging.info("~~~~~~~~~~~~retrieve_selected_uniprot_records_from_flatfile is finished~~~~~~~~~~~~")
