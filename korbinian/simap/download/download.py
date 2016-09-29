import logging
from time import strftime
import csv
import pandas as pd
import korbinian.mtutils as utils
import platform
import os
import tarfile
import korbinian

def download_homologues_from_simap(pathdict, set_, logging):
    """From the list of proteins in csv format, begins downloading homologues from the SIMAP database.

     - opens the csv file containing the list of proteins
     - opens or creates a text file with the list of failed downloads
     - checks if there is enough hard-drive space
     - checks what files currently exist (feature table, homologue, zip)
     - tries to download feature table
     - tries to download homologues
     - if both feature table and homologues exist, compresses both into a tarball and deletes original files
     - counts the number of failed downloads. Assumes most failed downloads are due to server errors on the SIMAP side.
     With more and more failed downloads, sleeps for longer and longer.


    Parameters
    ----------
    pathdict : dict
        Dictionary of the key paths and files associated with that List number.
    set_ : dict
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
    df = pd.read_csv(pathdict["list_summary_csv"], sep = ",", quoting = csv.QUOTE_NONNUMERIC, index_col = 0)
    # if "uniprot_acc" in df.columns:
    #     df.set_index("uniprot_acc", drop=False, inplace=True)
    # else:
    #     df["uniprot_acc"] = df.index
    acc_list_failed_downloads = []
    if os.path.isfile(pathdict["failed_downloads_txt"]):
        # Extracts accession numbers out of file
        with open(pathdict["failed_downloads_txt"], "r") as source:
            for line in source:
                line = line.strip()
                acc_list_failed_downloads.append(line)

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
    byteformat = "GB"
    data_harddrive = os.path.splitdrive(set_["data_dir"])[0]
    # print initial hard-drive space
    size = utils.get_free_space(data_harddrive, byteformat)
    logging.info('Hard disk remaining space = {}'.format(size))

    #iterate over each uniprot record contained in the dataframe. note that acc = uniprot accession number
    number_of_files_not_found = 0
    for acc in df.index:
        # check hand-drive space before each download
        try:
            size = utils.get_free_space(data_harddrive, byteformat)
            if size[0] < 5:
                raise utils.HardDriveSpaceException("Hard drive space limit reached, there is only %s %s space left." % (size[0], size[1]))
        except utils.HardDriveSpaceException as e:
            logging.warning(e)
        protein_name = df.loc[acc, 'protein_name']
        seqlen = df.loc[acc, 'seqlen']
        input_sequence = df.loc[acc, 'full_seq']
        SIMAP_tar = df.loc[acc, 'SIMAP_tar']
        ft_xml_path = df.loc[acc, 'SIMAP_feature_table_XML_path']
        homol_xml_path = df.loc[acc, 'SIMAP_homol_XML_path']
        date_file_path = df.loc[acc, 'SIMAP_download_date_file_path']
        if set_["attempt_prev_failed_downloads"] == False:
            if acc in acc_list_failed_downloads:
                logging.info("{} is in list of previously failed downloads. Will be skipped.".format(protein_name))
                continue

        # create directories to hold file, if necessary
        utils.make_sure_path_exists(homol_xml_path, isfile=True)

        #check which files exist. This is useful, because it is not possible to open the tarfile as 'a:gz',
        #therefore you cannot add files to an existing tarfile)
        feature_table_XML_exists, homologues_XML_exists, SIMAP_tarfile_exists = utils.check_tarfile(SIMAP_tar, ft_xml_path, homol_xml_path)

        ''' windows has a character limit in the command prompt in theory of 8191 characters,
        but the command line java command seems to cause errors with sequences above 3000 amino acids.
        Assume that this problem only applies to Windows,
        therefore in Windows systems limit the java string to proteins less than 3000 amino acids.
        The character limit can be adjusted in the settings file
        '''
        if 'Windows' in str(platform.system()):
            if seqlen > set_["max_query_sequence_length"]:
                logging.warning('%s homologue download will be skipped. It cannot be processed into a java command in windows OS,'
                                'as the sequence is longer than %i characters (%i). Moving to next sequence' % (protein_name, set_["max_query_sequence_length"],seqlen))
                # skip this protein
                continue

        if SIMAP_tarfile_exists and set_["overwrite_homologue_files"] == False:
            # skip this protein
            logging.info("{} SIMAP_tarfile_exists, download skipped.".format(acc))
            continue
        eaSimap_path = os.path.join(set_["data_dir"], "programs", "eaSimap.jar")
        if not feature_table_XML_exists:
            #download feature table from SIMAP
            korbinian.simap.download.download.retrieve_simap_feature_table(input_sequence,
                                                                           java_exec_str=java_exec_str,
                                                                           max_memory_allocation=max_memory_allocation,
                                                                           output_file=ft_xml_path,
                                                                           eaSimap_path=eaSimap_path)
        if not homologues_XML_exists:
            #download homologue file from SIMAP
            korbinian.simap.download.download.retrieve_simap_homologues(input_sequence,
                                                                        output_file=homol_xml_path,
                                                                        max_hits=max_hits, java_exec_str=java_exec_str,
                                                                        max_memory_allocation=max_memory_allocation, taxid=taxid,
                                                                        eaSimap_path=eaSimap_path)
        #now check again if the files exist
        feature_table_XML_exists, homologues_XML_exists, SIMAP_tarfile_exists = utils.check_tarfile(SIMAP_tar, ft_xml_path, homol_xml_path)
        if not homologues_XML_exists or not feature_table_XML_exists:
            # add accession number to the list of failed downloads
            with open(pathdict["failed_downloads_txt"], "a") as source:
                source.write("\n{}".format(acc))
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
        # sometimes the SIMAP server seems to like a little rest in between downloads?
        utils.sleep_x_seconds(30)

    logging.info('retrieve_simap_feature_table_and_homologues_from_list_in_csv is finished')


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
    # command_str = '%s -Xmx%im -jar %s -s %s -o %s -f' % (java_exec_str, max_memory_allocation, eaSimap_path, input_sequence, output_file)
    command_str = '{jes} -Xmx{mma:d}m -jar {esp} -s {s} -o {o} -f'.format(jes=java_exec_str,
                                      mma=max_memory_allocation, esp=eaSimap_path, s=input_sequence,o=output_file)
    logging.info(command_str)
    command = utils.Command(command_str)
    command.run(timeout=500)
    logging.info("Output file:     %s\n" % output_file),
    utils.sleep_x_seconds(5)
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
    #note that windows has a character limit in the command prompt in theory of 8191 characters, but the command line java command seems to cause errors with sequences above 3000 amino acids.
    #the 3000 character limit is currently applied in the main_simap script, rather than here
    #run command
    # command_str = '%s -Xmx%im -jar %s -s %s -m %s -o %s -x %s%s' % (java_exec_str, max_memory_allocation, eaSimap_path, input_sequence,
    #                                                                 max_hits, output_file, taxid_search_string)
    command_str = '{jes} -Xmx{mma:d}m -jar {esp} -s {s} -m {m} -o {o} -x{tss}'.format(jes=java_exec_str, mma=max_memory_allocation, esp=eaSimap_path, s=input_sequence,
                                                                    m=max_hits, o=output_file, tss=taxid_search_string)
    logging.info(command_str)
    command = utils.Command(command_str)
    #timeout = max_hits/5 if max_hits > 500 else 100
    timeout = 3000
    command.run(timeout=timeout) #give 1000 for 5000 hits to download?
    logging.info("Output file:     %s\n'file saved'" % output_file)
    #sleep_x_seconds(30)
    if not os.path.exists(output_file):
        logging.info('********************SIMAP download failed for : %s***************' % output_file)

#def retrieve_simap_from_multiple_fasta(input_file):
#    records = SeqIO.parse_uniprot(input_file, "fasta")
#    global list_of_files_with_feature_tables, list_of_files_with_homologues
#    list_of_files_with_feature_tables = []
#    list_of_files_with_homologues = []
#    recordcounter = 0
#    for record in records:
#        name = record.name.replace('|', '_')[:30]
#        accession = 'Acc'
#        label = '%s_%s' % (accession, name)
#        SIMAP_feature_table_XML_file = r"E:\\Stephis\\Projects\\Programming\\Python\\files\\learning\\simap\\%s_simap_feature_table.xml" % label
#        list_of_files_with_feature_tables.append(SIMAP_feature_table_XML_file)
#        SIMAP_homologues_XML_file = r"E:\\Stephis\\Projects\\Programming\\Python\\files\\learning\\simap\\%s_simap_homologues.xml" % label
#        list_of_files_with_homologues.append(SIMAP_homologues_XML_file)
#        retrieve_simap_feature_table(input_sequence=record.seq, output_file=SIMAP_feature_table_XML_file)
#        retrieve_simap_homologues(input_sequence=record.seq, output_file=SIMAP_homologues_XML_file, database='', max_hits='10', taxid='7227', extra_search_string='')
#        recordcounter += 1
#    logging.info('Download complete, %s SIMAP records saved.' % recordcounter)
#input_seqs_mult_fasta = r'E:\Stephis\Projects\Programming\Python\files\learning\simap\multiple_protein_seqs_in_fasta_format.txt'
#retrieve_simap_from_multiple_fasta(input_seqs_mult_fasta)
#throwaway functions, currently kept in main
#slice_TMD_seq = lambda x: x['full_seq'][int(x['%s_start'%TMD_name]-1):int(x['%s_end'%TMD_name])]
#slice_TMD_plus_surrounding_seq = lambda x: x['full_seq'][int(x['%s_start_plus_surr'%TMD_name]-1):int(x['%s_end_plus_surr'%TMD_name])]