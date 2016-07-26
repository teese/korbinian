import csv
import pandas as pd
import korbinian.mtutils as utils
import platform
import os
import tarfile

def download_homologues_from_simap(pathdict, set_, logging):
    # logging.info('~~~~~~~~~~~~  starting A06_retrieve_simap_feature_table_and_homologues_from_list_in_csv   ~~~~~~~~~~~~')
    # #test if the dataframe has already been created, otherwise reopen from csv file
    # if 'df' in globals():
    #     if isinstance(df, pd.DataFrame):
    #         logging.info('first protein acc = %s, df already exists' % df.iloc[0][0])
    # else:
    #     logging.info('df loaded from %s' % pathdict["dfout04_uniprotcsv_incl_paths"])
    df = pd.read_csv(pathdict["dfout04_uniprotcsv_incl_paths"], sep=",", quoting=csv.QUOTE_NONNUMERIC,index_col=0)

    #def retrieve_simap_feature_table_and_homologues_from_list_in_csv(input_file, list_of_keys, settings):
    '''
    First prepare the csv file from the uniprot record.
    Run this to save the files based on their domain.
    '''
    #global list_of_files_with_feature_tables, list_of_files_with_homologues
    #The SIMAP download settings can be altered as desired, using the json settings file
    max_hits = set_["max_hits"]
    java_exec_str = set_["old_run_stat_analysis_sim_ratios_in_dfout05java_exec_str"]
    max_memory_allocation = set_["java_max_RAM_memory_allocated_to_simap_download"]
    taxid = set_["taxid"]  # eg.'7227' for Drosophila melanogaster

    enough_hard_drive_space = True
    try:
        byteformat = "GB"
        data_harddrive = set_["old_run_stat_analysis_sim_ratios_in_dfout05data_harddrive"]
        size = utils.get_free_space(data_harddrive, byteformat)
        logging.info('Hard disk remaining space =')
        logging.info(size)
        if size[0] < 5:
            raise utils.HardDriveSpaceException("Hard drive space limit reached, there is only %s %s space left." % (size[0], size[1]))
    except utils.HardDriveSpaceException as e:
        logging.warning(e)
    if enough_hard_drive_space:
        #iterate over each uniprot record contained in the dataframe. note that acc = uniprot accession number
        number_of_files_not_found = 0
        for acc in df.index:
            protein_name = df.loc[acc, 'A2_protein_name']
            query_sequence_length = df.loc[acc, 'uniprot_seqlen']
            input_sequence = df.loc[acc, 'uniprot_seq']
            ''' windows has a character limit in the command prompt in theory of 8191 characters,
            but the command line java command seems to cause errors with sequences above 3000 amino acids.
            Assume that this problem only applies to Windows,
            therefore in Windows systems limit the java string to proteins less than 3000 amino acids.
            The character limit can be adjusted in the settings file
            '''
            if 'Windows' in str(platform.system()):
                if query_sequence_length < set_["max_query_sequence_length"]:
                    download_homologues = True
                else:
                    download_homologues = False
                    logging.warning('%s cannot be processed into a java command in windows OS,'
                                    'as the sequence is longer than %i characters (%i). Moving to next sequence' % (
                                    protein_name, set_["max_query_sequence_length"],
                                    query_sequence_length))
            else:
                download_homologues = True
            if download_homologues == True:
                simap_data_folder = os.path.join(set_['data_folder'], 'simap')
                subfolder = os.path.join(simap_data_folder, df.loc[acc, 'first_two_letters_of_uniprot_acc'])
                if os.path.isdir(subfolder) == False:
                    os.mkdir(subfolder)
                #check which files exist. This is useful, because it is not possible to open the tarfile as 'a:gz',
                #therefore you cannot add files to an existing tarfile)
                feature_table_XML_exists, homologues_XML_exists, SIMAP_tarfile_exists = utils.check_tarfile(df, acc)
                if not SIMAP_tarfile_exists:
                    if not feature_table_XML_exists:
                        #download feature table from SIMAP
                        utils.retrieve_simap_feature_table(input_sequence,
                                                           java_exec_str=java_exec_str,
                                                           max_memory_allocation=max_memory_allocation,
                                                           output_file=df.loc[acc, 'SIMAP_feature_table_XML_file_path'],
                                                           eaSimap_path=set_["old_run_stat_analysis_sim_ratios_in_dfout05eaSimap_path"])
                    if not homologues_XML_exists:
                        #download homologue file from SIMAP
                        utils.retrieve_simap_homologues(input_sequence,
                                                        output_file=df.loc[acc, 'SIMAP_homologues_XML_file_path'],
                                                        max_hits=max_hits, java_exec_str=java_exec_str,
                                                        max_memory_allocation=max_memory_allocation, taxid=taxid,
                                                        eaSimap_path=set_["old_run_stat_analysis_sim_ratios_in_dfout05eaSimap_path"])
                        #now check again if the files exist
                    feature_table_XML_exists, homologues_XML_exists, SIMAP_tarfile_exists = utils.check_tarfile(df, acc)
                    if not homologues_XML_exists:
                        #add one to the list of consecutive failed downloads.
                        number_of_files_not_found += 1
                        #if a large number of downloads failed, then the SIMAP server is probably not working.
                        #Wait some time and try again later.
                        if number_of_files_not_found > 30:
                            utils.sleep_24_hours()
                        if number_of_files_not_found == 20:
                            utils.sleep_24_hours()
                        if number_of_files_not_found == 15:
                            utils.sleep_6_hours()
                        if number_of_files_not_found == 10:
                            utils.sleep_6_hours()
                        #utils.sleep_120_seconds()
                    else:
                        #if download is successful or file exists, the SIMAP server must be working,
                    #therefore reset the number_of_files_not_found
                        number_of_files_not_found = 0

                    #since we can't add files to the compressed tarfile, only when both the feature table
                    #and xml file are downloaded should we pack and compress them
                    if feature_table_XML_exists and homologues_XML_exists:
                        with tarfile.open(df.loc[acc, 'SIMAP_tarfile'], mode='w:gz') as tar:
                            #add the files to the compressed tarfile
                            logging.info('%s XML files will be moved into the tarball, original XML files deleted' % protein_name)
                            tar.add(df.loc[acc, 'SIMAP_feature_table_XML_file_path'],
                                    arcname=df.loc[acc, 'SIMAP_feature_table_XML_file'])
                            tar.add(df.loc[acc, 'SIMAP_homologues_XML_file_path'],
                                    arcname=df.loc[acc, 'SIMAP_homologues_XML_file'])
                        #delete the original files
                        try:
                            os.remove(df.loc[acc, 'SIMAP_feature_table_XML_file_path'])
                            os.remove(df.loc[acc, 'SIMAP_homologues_XML_file_path'])
                        except FileNotFoundError:
                            pass

                            #now add the downloaded files to the tarball, and delete the original XML files

                        #            directory = r'E:/Stephis/Projects/Programming/Python/files/learning/tarfile'
                        #            newtar = directory + r'/newtar.tar.gz'
                        #            SIMAP_feature_table_XML_file_basename = 'P29274_AA2AR_HUMAN_feature_table.xml'
                        #            SIMAP_feature_table_XML_file = '%s/%s' % (directory, SIMAP_feature_table_XML_file_basename)
                        #            SIMAP_homologues_XML_file_basename = 'P29274_AA2AR_HUMAN_homologues.xml'
                        #            SIMAP_homologues_XML_file = '%s/%s' % (directory, SIMAP_homologues_XML_file_basename)
                        #
                        #            with tarfile.open(df.loc[acc, 'SIMAP_tarfile'], 'w:gz') as tar:
                        #                tar.add(SIMAP_feature_table_XML_file, arcname = SIMAP_feature_table_XML_file_basename)
                        #                tar.add(SIMAP_homologues_XML_file, arcname = SIMAP_homologues_XML_file_basename)
                        #            with tarfile.open(newtar, 'r:gz') as tar:
                        #                for tarinfo in tar:
                        #                    print(tarinfo.isreg())
                        #                    print(tarinfo.name)
                        #                    print(tarinfo.size)
                        #
                        #            df['SIMAP_tarfile'] = df.simap_filename_base + '_SIMAP.tar.gz'
                        #            df['SIMAP_feature_table_XML_file'] = df.A2_protein_name + '_feature_table.xml'
                        #            df['SIMAP_feature_table_XML_file_path'] = df.simap_filename_base + '_feature_table.xml'
                        #            df['SIMAP_homologues_XML_file'] = df.A2_protein_name + '_homologues.xml'
                        #            df['SIMAP_homologues_XML_file_path'] = df.A2_protein_name + '_homologues.xml'
                        #            df['SIMAP_csv_from_XML'] = df.simap_filename_base + '_homologues.csv'
                        #            df['SIMAP_csv_from_XML_path'] = df.simap_filename_base + '_homologues.csv'
                        #            df['csv_SIMAP_homologues_kept_for_statistical_analysis'] = df.simap_filename_base + '_homologues_for_stat.csv'
            #print("download homologues = %s" % download_homologues)
    logging.info('retrieve_simap_feature_table_and_homologues_from_list_in_csv is finished')