import os
import pandas as pd
import korbinian.mtutils as utils

def convert_uniprot_list_to_nonred_ff_via_uniref(set_, list_number, uniprot_dir, logging, uniprot_flatfile_of_selected_records):
    """ Convert uniprot accession list to a non-redundant uniprot flatfile using UniRef clusters
    """
    # load uniref cutoff used (typically 50, for UniRef50)
    uniref_cutoff = set_["uniref_cluster_cutoff"]
    # define path to csv file containing the list of redundant uniprot accessions
    redundant_list_uniprot_acc_csv = os.path.join(uniprot_dir ,"List%02d_redundant_list_uniprot_acc.tab" % list_number)
    # define path to uniprot flatfile containing the redundant protein records
    redundant_uniprot_flatfile = os.path.join(uniprot_dir ,"List%02d_redundant_uniprot_flatfile.txt" % list_number)
    # define path to the csv file containing the relevant uniref clusters applicable to this list of proteins
    uniref_clusters_csv = os.path.join(uniprot_dir ,"List%02d_UniRef%02d_clusters.tab" % (list_number ,uniref_cutoff))
    # create a new dataframe with the uniref csv file, containing the accessions of the reference sequences
    df_uniref = pd.read_table(uniref_clusters_csv)
    # to simplify, keep only the columns with the uniprot accessions
    # for backwards compatibility, check if the old header is used
    if 'Cluster member(s)' in df_uniref.columns:
        df_uniref.rename(columns={'Cluster member(s)': 'Cluster members'})
    # to simplify, keep only the columns with the uniprot accessions
    df_uniref = df_uniref[['Cluster ID', 'Cluster members']  ]# 'Organisms'
    df_uniref = df_uniref.set_index('Cluster ID')
    # change the cluster ID to the index
    df_uniref.index.names = ['cluster_ID']
    # convert the list of uniprot accessions in each cluster to a python list format
    df_uniref['all_acc_in_cluster'] = df_uniref['Cluster members'].apply(lambda x : [x.strip() for x in x.split(';')])
    # delete the original list of clusters
    df_uniref.drop("Cluster members", axis=1, inplace=True)
    if os.path.isfile(redundant_list_uniprot_acc_csv) == False:
        logging.warning('warning, file with nonredundant uniprot acc does not exist : %s' % redundant_list_uniprot_acc_csv)
    # open up large csv containing the uniprot acc of all single-pass proteins in uniprot
    df_uniprot_acc = pd.read_table(redundant_list_uniprot_acc_csv)
    # set the uniprot acc as the index
    df_uniprot_acc = df_uniprot_acc.set_index('Entry')
    # create an empty column to contain the cluster ID
    df_uniprot_acc['Cluster ID'] = ''
    # create a new dataframe with the uniref csv file, containing the accessions of the reference sequences
    # find the appropriate uniref cluster containing each acc
    acc_counter = 0
    for acc in df_uniprot_acc.index:
        # check if the accession is already a uniref representative (saves searching time!)
        if 'UniRef%i_%s' % (uniref_cutoff, acc) in df_uniref.index:
            df_uniprot_acc.loc[acc, 'Cluster ID'] = 'UniRef%i_%s' % (uniref_cutoff ,acc)
        else:
            # if it is not a uniref representative, go through each uniref cluster, checking to see if it is a member
            for acc_cluster in list(df_uniref['all_acc_in_cluster']):
                if acc in acc_cluster:
                    # if the acc is a member of a uniref cluster, add the cluster name to the original dataframe
                    df_uniprot_acc.loc[acc, 'Cluster ID'] = 'UniRef%i_%s' % (uniref_cutoff, acc_cluster[0])
        acc_counter += 1
        # write a dot on the screen for each protein, so that it is easy to see that the script is still working
        # sys.stdout.write(". ")
        if acc_counter % 100 == 0:
            logging.info('%i records checked for redundancy' % acc_counter)
    # determine which cluster IDs are in the database only once, as these uniprot entries are already nonredundant
    series_unique_bool = df_uniprot_acc['Cluster ID'].value_counts() == 1
    list_uniref_clusters_matching_only_one_acc = series_unique_bool.loc[series_unique_bool == True].index
    # now use this list to label the original sequences that are nonredundant
    df_uniprot_acc['nonred'] = df_uniprot_acc['Cluster ID'].apply(lambda x : x in list_uniref_clusters_matching_only_one_acc)

    # create a list of uniprot accessions that are nonredundant
    df_uniprot_acc_nonred = df_uniprot_acc.loc[df_uniprot_acc['nonred'] == True]
    list_nonred_acc = list(df_uniprot_acc_nonred.index)

    # create a uniprot flatfile containing only the desired nonredundant accessions
    utils.retrieve_selected_uniprot_records_from_flatfile(list_nonred_acc,
                                                          redundant_uniprot_flatfile,
                                                          uniprot_flatfile_of_selected_records, logging)
    number_nonredundant_records = len(list_nonred_acc)
    number_total_records = df_uniprot_acc.shape[0]
    number_redundant_records = number_total_records - number_nonredundant_records

    logging.info('number_total_records = {0}, number_redundant_records removed = {1}, '
                 'final number_nonredundant_records = {2}'.format(number_total_records,
                                                                  number_redundant_records,
                                                                  number_nonredundant_records))
    logging.info("A00_convert_uniprot_list_to_nonred_ff_via_uniref is finished")