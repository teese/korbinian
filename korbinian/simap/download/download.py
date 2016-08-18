from time import strftime
import csv
import pandas as pd
import korbinian.mtutils as utils
import platform
import os
import tarfile

def download_homologues_from_simap(pathdict, set_, logging):

    df = pd.read_csv(pathdict["list_summary_csv"], sep = ",", quoting = csv.QUOTE_NONNUMERIC, index_col = 0)
    #def retrieve_simap_feature_table_and_homologues_from_list_in_csv(input_file, list_of_keys, settings):
    '''
    First prepare the csv file from the uniprot record.
    Run this to save the files based on their domain.
    '''
    #global list_of_files_with_feature_tables, list_of_files_with_homologues
    #The SIMAP download settings can be altered as desired, using the json settings file
    max_hits = set_["max_hits"]
    java_exec_str = set_["java_exec_str"]
    max_memory_allocation = set_["java_max_RAM_memory_allocated_to_simap_download"]
    taxid = set_["taxid"]  # eg.'7227' for Drosophila melanogaster
    enough_hard_drive_space = True
    try:
        byteformat = "GB"
        data_harddrive = set_["data_harddrive"]
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
            protein_name = df.loc[acc, 'protein_name']
            query_sequence_length = df.loc[acc, 'uniprot_seqlen']
            input_sequence = df.loc[acc, 'full_seq']
            SIMAP_tar = df.loc[acc, 'SIMAP_tar']
            ft_xml_path = df.loc[acc, 'SIMAP_feature_table_XML_path']
            homol_xml_path = df.loc[acc, 'SIMAP_homol_XML_path']
            date_file_path = df.loc[acc, 'SIMAP_download_date_file_path']
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
                # create directories to hold file, if necessary
                utils.make_sure_path_exists(homol_xml_path, isfile=True)

                #check which files exist. This is useful, because it is not possible to open the tarfile as 'a:gz',
                #therefore you cannot add files to an existing tarfile)
                feature_table_XML_exists, homologues_XML_exists, SIMAP_tarfile_exists = utils.check_tarfile(SIMAP_tar, ft_xml_path, homol_xml_path)
                if not SIMAP_tarfile_exists or set_["overwrite_homologue_files"] == True:
                    if not feature_table_XML_exists:
                        #download feature table from SIMAP
                        utils.retrieve_simap_feature_table(input_sequence,
                                                           java_exec_str=java_exec_str,
                                                           max_memory_allocation=max_memory_allocation,
                                                           output_file=ft_xml_path,
                                                           eaSimap_path=set_["eaSimap_path"])
                    if not homologues_XML_exists:
                        #download homologue file from SIMAP
                        utils.retrieve_simap_homologues(input_sequence,
                                                        output_file=homol_xml_path,
                                                        max_hits=max_hits, java_exec_str=java_exec_str,
                                                        max_memory_allocation=max_memory_allocation, taxid=taxid,
                                                        eaSimap_path=set_["eaSimap_path"])
                    #now check again if the files exist
                    feature_table_XML_exists, homologues_XML_exists, SIMAP_tarfile_exists = utils.check_tarfile(SIMAP_tar, ft_xml_path, homol_xml_path)
                    if not homologues_XML_exists:
                        #add one to the list of consecutive failed downloads.
                        number_of_files_not_found += 1
                        #if a large number of downloads failed, then the SIMAP server is probably not working.
                        #Wait some time and try again later.
                        if number_of_files_not_found > 30:
                            utils.sleep_x_hours(24)
                        if number_of_files_not_found == 20:
                            utils.sleep_x_hours(24)
                        if number_of_files_not_found == 15:
                            utils.sleep_x_hours(6)
                        if number_of_files_not_found == 10:
                            utils.sleep_x_hours(6)
                    else:
                        #if download is successful or file exists, the SIMAP server must be working,
                        #therefore reset the number_of_files_not_found
                        number_of_files_not_found = 0

                    # create an empty text file with the download date
                    date = strftime("%Y%m%d")
                    with open(date_file_path, "w") as f:
                        f.write("{}\nEmpty file with the date.\nHave a nice day!".format(date))
                    # since we can't add files to the compressed tarfile, only when both the feature table
                    #and xml file are downloaded should we pack and compress them
                    if feature_table_XML_exists and homologues_XML_exists:
                        with tarfile.open(SIMAP_tar, mode='w:gz') as tar:
                            #add the files to the compressed tarfile
                            logging.info('%s XML files will be moved into the tarball, original XML files deleted' % protein_name)
                            tar.add(ft_xml_path,arcname=os.path.basename(ft_xml_path))
                            tar.add(homol_xml_path,arcname=os.path.basename(homol_xml_path))
                            tar.add(date_file_path, arcname=os.path.basename(date_file_path))
                        #delete the original files
                        try:
                            os.remove(ft_xml_path)
                            os.remove(homol_xml_path)
                            os.remove(date_file_path)
                        except FileNotFoundError:
                            pass

    logging.info('retrieve_simap_feature_table_and_homologues_from_list_in_csv is finished')