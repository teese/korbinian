import csv
import logging
import os
import platform
import tarfile
from time import strftime
import korbinian
import korbinian.utils as utils
import pandas as pd
import sys
# import debugging tools
from korbinian.utils import pr, pc, pn, aaa

def download_homologues_from_simap(pathdict, s, logging):
    """From the list of proteins in csv format, begins downloading homologues from the SIMAP database.

     - opens the csv file containing the list of proteins
     - opens or creates a text file with the list of failed downloads
     - checks if there is enough hard-drive space
     - checks what files currently exist (feature table, homologue, zip)
     - tries to download feature table (if download_feature_tables listed as TRUE in settings)
     - tries to download homologues
     - if both feature table and homologues exist, compresses both into a tarball and deletes original files
     - counts the number of failed downloads. Assumes most failed downloads are due to server errors on the SIMAP side.
     With more and more failed downloads, sleeps for longer and longer.

    Parameters
    ----------
    pathdict : dict
        Dictionary of the key paths and files associated with that List number.
    s : dict
        Settings dictionary extracted from excel settings file.
    logging : logging.Logger
        Logger for printing to console and logfile.

    Saved Files and Figures
    -----------------------
    PROTEIN_NAME_SIMAP.tar.gz : gzip file
        (e.g. A1A5B4_ANO9_HUMAN_SIMAP.tar.gz)

        Contains
        --------
        PROTEIN_NAME_feature_table.xml (e.g. A1A5B4_ANO9_HUMAN_feature_table.xml)
            XML feature table from SIMAP, with information regarding each protein.
        PROTEIN_NAME_homologues.xml (e.g. A1A5B4_ANO9_HUMAN_homologues.xml)
            homologues from SIMAP in SIMAP-XML (rather than BLAST-XML) format
        PROTEIN_NAME--DATE--RESEARCHERNAME.txt (e.g. A1DT13_A1DT13_HUMAN--20160811--Mark Teese.txt)
            [only in later versions] Text file showing the download date and researcher name.

    pathdict["failed_downloads_txt"] : txt
        File containing a list of accessions that could not be downloaded. At each run, the program checks
        if this file exists. If it doesn't exist, it will be created. If it exists, the settings file
        determines whether the previously failed downloads will be re-attempted.
    """
    logging.info("~~~~~~~~~~~~                 starting download_homologues_from_simap                ~~~~~~~~~~~~")
    df = pd.read_csv(pathdict["list_csv"], sep = ",", quoting = csv.QUOTE_NONNUMERIC, index_col = 0)

    if s["attempt_prev_failed_downloads"] == False:
        # get list of accessions that could not be downloaded, and can immediately be excluded
        acc_list_failed_downloads = utils.get_acc_list_from_txt(pathdict["failed_downloads_txt"])
        not_in_homol_db = utils.get_acc_list_from_txt(pathdict["acc_not_in_homol_db_txt"])
        acc_excluded_list = acc_list_failed_downloads + not_in_homol_db
        # the list of desired proteins = total_list - excluded
        acc_not_excluded = list(set(df.index) - set(acc_excluded_list))
        # filter dataframe to only contain the desired proteins, which have not been excluded
        df = df.loc[acc_not_excluded, :]

    max_hits = s["max_hits"]
    java_exec_str = s["java_exec_str"]
    max_memory_allocation = s["java_max_RAM_memory_allocated_to_simap_download"]
    taxid = s["taxid"]  # eg.'7227' for Drosophila melanogaster

    # if "Linux" in platform.system() or "Windows" in platform.system():
    #     # if Linux or Windows
    #     byteformat = "GB"
    #     data_harddrive = os.path.splitdrive(s["data_dir"])[0]
    #     # print initial hard-drive space
    #     size = utils.get_free_space(data_harddrive, byteformat)
    #     logging.info('Hard disk remaining space = {}'.format(size))
    # else:
    #     # assume the system is a mac
    #     # code works only on mac!!! reverted to previous version
    #     statvfs = os.statvfs(s["simap_dir"])
    #     available_space = statvfs.f_frsize * statvfs.f_bavail
    #     size = available_space / 1073741824
    #     # print initial hard-drive space
    #     logging.info('Hard disk remaining space = {:.2f} GB'.format(size))

    #iterate over each uniprot record contained in the dataframe. note that acc = uniprot accession number
    number_of_files_not_found = 0
    for acc in df.index:
        # check hand-drive space before each download
        # try:
        #     if "Linux" in platform.system() or "Windows" in platform.system():
        #         size = utils.get_free_space(data_harddrive, byteformat)
        #         if size[0] < 5:
        #             raise utils.HardDriveSpaceException("Hard drive space limit reached, there is only %s %s space left." % (size[0], size[1]))
        #     else:
        #         # MAC only stuff...
        #         statvfs = os.statvfs(s["simap_dir"])
        #         available_space = statvfs.f_frsize * statvfs.f_bavail
        #         size = available_space / 1073741824
        #         if size < 5:
        #             raise utils.HardDriveSpaceException("Hard drive space limit reached, there is only %s %s space left." % (size[0], size[1]))
        # except utils.HardDriveSpaceException as e:
        #     logging.warning(e)
        protein_name = df.loc[acc, 'protein_name']
        seqlen = df.loc[acc, 'seqlen']
        input_sequence = df.loc[acc, 'full_seq']
        SIMAP_tar = df.loc[acc, 'SIMAP_tar']
        ft_xml_path = df.loc[acc, 'SIMAP_feature_table_XML_path']
        homol_xml_path = df.loc[acc, 'SIMAP_homol_XML_path']
        date_file_path = df.loc[acc, 'SIMAP_download_date_file_path']

        # create directories to hold file, if necessary
        utils.make_sure_path_exists(homol_xml_path, isfile=True)

        #check which files exist and delete corrupt tarballs
        ft_XML_exists, homol_XML_exists, SIMAP_tar_exists, ff, hh = utils.check_SIMAP_tarfile(SIMAP_tar, ft_xml_path, homol_xml_path,
                                                                                                               acc, logging, delete_corrupt=True)
        ''' Windows command prompt accepts only 8191 characters.
            Limit protein length according to settings (typically max length = 3000)
        '''
        if 'Windows' in str(platform.system()):
            too_large_to_download_list = utils.get_acc_list_from_txt(pathdict["too_large_to_download_txt"])
            if seqlen > s["max_query_sequence_length"]:
                logging.warning('%s homologue download will be skipped. It cannot be processed into a java command in windows OS,'
                                'as the sequence is longer than %i characters (%i). Moving to next sequence' % (protein_name, s["max_query_sequence_length"],seqlen))
                # if the accession is not already in the text file, add it
                if acc not in too_large_to_download_list:
                    # add accession number to the list of failed downloads
                    with open(pathdict["too_large_to_download_txt"], "a") as source:
                        source.write("\n{}".format(acc))
                # skip this protein
                continue

        if SIMAP_tar_exists and s["overwrite_homologue_files"] == False:
            # skip this protein
            logging.info("{} SIMAP_tar_exists, download skipped.".format(acc))
            continue
        eaSimap_path = os.path.join(s["data_dir"], "programs", "eaSimap.jar")
        # NOTE: DOWNLOADING FEATURE TABLES IS NO LONGER CONSIDERED NECESSARY.
        if not ft_XML_exists and s["download_feature_tables"] == True:
            #download feature table from SIMAP
            korbinian.simap_download.retrieve_simap_feature_table(input_sequence,
                                                                  java_exec_str=java_exec_str,
                                                                  max_memory_allocation=500,
                                                                  output_file=ft_xml_path,
                                                                  eaSimap_path=eaSimap_path)
            utils.sleep_x_seconds(60)
        if not homol_XML_exists:
            #download homologue file from SIMAP
            korbinian.simap_download.retrieve_simap_homologues(input_sequence,
                                                               output_file=homol_xml_path,
                                                               max_hits=max_hits, java_exec_str=java_exec_str,
                                                               max_memory_allocation=max_memory_allocation, taxid=taxid,
                                                               eaSimap_path=eaSimap_path)
            # sometimes the SIMAP server seems to like a little rest in between downloads?
            utils.sleep_x_seconds(30)
        #now check again if the files exist
        ft_XML_exists, homol_XML_exists, SIMAP_tar_exists, ff, hh = utils.check_SIMAP_tarfile(SIMAP_tar, ft_xml_path, homol_xml_path,
                                                                                                               acc, logging)
        if not homol_XML_exists:
            # add accession number to the list of failed downloads
            with open(pathdict["failed_downloads_txt"], "a") as source:
                source.write("\n{}".format(acc))
            #add one to the list of consecutive failed downloads.
            number_of_files_not_found += 1
            if s["sleep_if_downloads_unsuccessful"]:
                # if a large number of downloads failed, then the SIMAP server is probably not working.
                # Wait some time and try again later.
                if number_of_files_not_found > 30:
                    sys.stdout.write("\nnumber_of_files_not_found = {}, sleeping for 24 h".format(number_of_files_not_found))
                    utils.sleep_x_hours(24)
                if number_of_files_not_found == 20:
                    sys.stdout.write("\nnumber_of_files_not_found = {}, sleeping for 6 h".format(number_of_files_not_found))
                    utils.sleep_x_hours(6)
                if number_of_files_not_found == 15:
                    sys.stdout.write("\nnumber_of_files_not_found = {}, sleeping for 1 h".format(number_of_files_not_found))
                    utils.sleep_x_hours(1)
        else:
            # if download is successful or file exists, the SIMAP server must be working,
            # therefore reset the number_of_files_not_found
            number_of_files_not_found = 0
            # create an empty text file with the download date
            date = strftime("%Y%m%d")
            with open(date_file_path, "w") as f:
                f.write("{}\nEmpty file with the date.\nHave a nice day!".format(date))
            with tarfile.open(SIMAP_tar, mode='w:gz') as tar:
                #add the files to the compressed tarfile
                logging.info('%s XML files will be moved into the tarball, original XML files deleted' % protein_name)
                tar.add(homol_xml_path, arcname=os.path.basename(homol_xml_path))
                tar.add(date_file_path, arcname=os.path.basename(date_file_path))
                if ft_XML_exists:
                    tar.add(ft_xml_path, arcname=os.path.basename(ft_xml_path))
            #delete the original files
            try:
                os.remove(homol_xml_path)
                os.remove(date_file_path)
                if ft_XML_exists:
                    os.remove(ft_xml_path)
            except FileNotFoundError:
                pass
    logging.info("~~~~~~~~~~~~                 finished download_homologues_from_simap                ~~~~~~~~~~~~")

def retrieve_simap_feature_table(input_sequence, java_exec_str, max_memory_allocation, output_file, eaSimap_path):
    """ Runs eaSimap.jar from the command line, to download the feature table XML from SIMAP.

    Downloads the feature table, which contains various information on that protein (TM regions, names, etc).

    Parameters
    ----------
    input_sequence : str
        Protein sequence.
    java_exec_str : str
        String used to run Java Runtime Environment. E.g. "java" or "%JAVA_JRE%"
    max_memory_allocation : int
        Maximum memory allocated to running the jar file. E.g. 3000, for 3GB ram.
    output_file : str
        Path for output XML file.
    eaSimap_path : str
        Path to the eaSimap.jar file in this operating system, to be run with java runtime environment.

    Saved Files and Figures
    -----------------------
    output_file : XML
        Feature table XML file.

    """
    #prepare input sequence and settings as a "command_str", and run command
    command_str = '{jes} -Xmx{mma:d}m -jar {esp} -s {s} -o {o} -f'.format(jes=java_exec_str,
                                      mma=max_memory_allocation, esp=eaSimap_path, s=input_sequence,o=output_file)
    logging.info(command_str)
    command = utils.Command(command_str)
    command.run(timeout=1500)
    logging.info("Output file:     %s\n" % output_file)
    if not os.path.exists(output_file):
        logging.info('********************SIMAP download failed for : %s***************' % output_file)


def retrieve_simap_homologues(input_sequence, output_file, max_hits, java_exec_str, max_memory_allocation, taxid, eaSimap_path):
    """Runs eaSimap.jar from the command line, to download the homologues XML from SIMAP.

    Downloads the homologues, which contains BLAST-like output in SIMAP format.

    Parameters
    ----------
    input_sequence : str
        Protein sequence.
    output_file : str
        Path for output XML file.
    max_hits : int
        Maximum number of hits from the BLAST-like output (usually 5000, the maximum).
    java_exec_str : str
        String used to run Java Runtime Environment. E.g. "java" or "%JAVA_JRE%"
    max_memory_allocation : int
        Maximum memory allocated to running the jar file. E.g. 3000, for 3GB ram.
    taxid : int
        Taxonomy ID (usually blank)
    eaSimap_path : str
        Path to the eaSimap.jar file in this operating system, to be run with java runtime environment.

    Saved Files and Figures
    -----------------------
    output_file : XML
        Homologue XML file.

    """
    '''
    Uses the java program to access the simap database and download the large file containing all homologues of that protein.
    '''
    # database selection is currently not working for download. Can be filtered later from all results.
    #database_dictionary = {313: 'uniprot_swissprot', 314: 'uniprot_trembl', 595: 'refseq', 721: 'Escherichia coli', 1296: 'Homo sapiens', 4250: 'Hot springs metagenome'}
    taxid_search_string = '' if taxid == '""' else '-i %s' % taxid
    command_str = '{jes} -Xmx{mma:d}m -jar {esp} -s {s} -m {m} -o {o} -x{tss}'.format(jes=java_exec_str, mma=max_memory_allocation, esp=eaSimap_path, s=input_sequence,
                                                                    m=max_hits, o=output_file, tss=taxid_search_string)
    logging.info(command_str)
    command = utils.Command(command_str)
    #timeout = max_hits/5 if max_hits > 500 else 100
    # set the timeout for 1000 seconds. This has never been reached, usually a Java error is returned much earlier.
    timeout = 1000
    command.run(timeout=timeout)
    if not os.path.exists(output_file):
        logging.info('********************SIMAP download failed for : %s***************' % output_file)