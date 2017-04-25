from time import strftime
import glob
import korbinian
import korbinian.utils as utils
import os
import tarfile
import warnings
warnings.filterwarnings('ignore')

def run_psiblast_on_fasta_queries_in_folder(query_dir, databases_dir, psiblast_exec_str, db, retry_failed=False, retry_successful=False):
    """Runs standalone PSIBLAST on every query fasta file in a folder.

    What you need:
         - query folder with I3L0P3.fasta, P42532.fasta etc
         - standalone BLAST must be working
         - databases folder with databases\\uniref\\vertebra90\\vertebra90.fasta, or databases\\uniref\\uniref90\\uniref90.fasta etc
    Where are the output files saved?
        - databases\\BLAST\\PSI\\vertebra90\\P4\\P42532_PSIBLAST.tar.gz

    Parameters
    ----------
    query_dir : str
        Folder containing protein sequences in fasta format.
    databases_dir : str
        Databases directory, e.g. "D:\Databases"
    psiblast_exec_str : str
        Path to psiblast executable
        if you are using linux or your Windows environmental variables are working, can simply be "psiblast"
    db : str
        Database for PSI-BLAST
        e.g. "metazoa90"
        Determines the filepath for the .fasta containing the search database.
    retry_failed : bool
        If True, proteins in the list of failed acc will be re-attempted
    retry_successful : bool
        If True, previous files will be overwritten

    Usage
    -----
    %capture
    from korbinian.psiblast import run_psiblast_on_fasta_queries_in_folder

    # SET YOUR GENERAL DATABASES DIRECTORY
    databases_dir = r"D:\\Databases"
    # uniref files must be in databases\\uniref\\database

    # SET THE DIRECTORY CONTAINING THE FASTA SEQUENCES YOU WANT TO USE AS A QUERY
    query_dir = r"D:\\Databases\\xtestproteins\\BLAST_small"
    # SET YOUR BLAST EXECUTABLE (windows)
    psiblast_exec_str = r"C:\\Program Files\\NCBI\\blast-2.6.0+\\bin\\psiblast.exe"
    run_psiblast_on_fasta_queries_in_folder(query_dir, databases_dir, psiblast_exec_str)

    """
    # set location of logfile
    logfile = os.path.join(query_dir,"PSI-BLAST_logfile.txt")
    logging = korbinian.common.setup_error_logging(logfile)
    # set location of txt file containing the failed sequences
    failed_psiblast_list_txt = os.path.join(query_dir,"failed_PSIBLAST_list.txt")
    ########################################################################################
    #                                                                                      #
    #       Create a list of all FASTA files in a particular folder for analysis           #
    #                                                                                      #
    ########################################################################################
    query_fasta_list = glob.glob(query_dir + "\*.fasta")
    logging.info("query_fasta_list[0:5] : {}".format(query_fasta_list[0:5]))

    ########################################################################################
    #                                                                                      #
    #                        Get list of previously failed sequences                       #
    #                                                                                      #
    ########################################################################################
    if os.path.isfile(failed_psiblast_list_txt):
        failed_psiblast_list = utils.get_acc_list_from_txt(failed_psiblast_list_txt)
    else:
        failed_psiblast_list = []
    logging.info("failed_psiblast_list[0:5] : {}".format(failed_psiblast_list[0:5]))

    ########################################################################################
    #                                                                                      #
    #                create a dictionary, s, with various parameters                       #
    #               can be converted to a korbinian dictionary later                       #
    #                                                                                      #
    ########################################################################################
    s = {}
    s["psiblast_exec_str"] = psiblast_exec_str
    s["evalue"] = "1e-3"
    s["inclusion_ethresh"] = "1e-3"
    s["num_threads"] = 10
    # s["db"] = "metazoa90"
    s["num_descriptions"] = 3000
    s["num_alignments"] = 3000
    command_str = '"{psiblast_exec_str}" -query {query} -db {db} -out_pssm {out_pssm} -out_ascii_pssm {out_ascii_pssm} '\
    '-out {out_BLAST_xml} -evalue {evalue} -inclusion_ethresh {inclusion_ethresh} -num_iterations 3 '\
    '-use_sw_tback -seg no -num_threads {num_threads} -num_descriptions {num_descriptions} -num_alignments {num_alignments}'
    logging.info("Example of command str, before inserting variables".format(command_str))

    ########################################################################################
    #                                                                                      #
    #         Run PSI-BLAST for each query sequence:                                       #
    #             input: query.fas, database.fasta (after makeblastdb)                     #
    #             output: query.pssm, query_ascii.pssm, query_BLAST.xml, query_date.txt    #
    #                     (compressed into a tarball, query_PSIBLAST.tar.gz)               #
    #                                                                                      #
    ########################################################################################
    # define the BLAST database. Note that you should specify the fasta file.
    db_path = os.path.join(databases_dir, "uniref\{db}\{db}.fasta".format(db=db))

    for query in query_fasta_list:
        acc = os.path.basename(query).split(".")[0]
        logging.info("\n{}".format(acc))
        # get first two letters of acc, used as a subfolder
        first2 = acc[:2]
        # create a basename, e.g. "D:\Databases\BLAST\PSI\vertebra90\P4\P42532" from which files are created
        basename = r"D:\Databases\BLAST\PSI\{db}\{first2}\{acc}".format(db=db, first2=first2, acc=acc)
        # create path for output files
        out_pssm = basename + ".pssm"
        out_ascii_pssm = basename + "_ascii.pssm"
        print("out_ascii_pssm", out_ascii_pssm)
        out_BLAST_xml = basename + "_BLAST.xml"
        date_file_path = basename + "_BLAST_date.txt"
        PSIBLAST_tar = basename + "_PSIBLAST.tar.gz"

        # if the tar exists or the accession has previously failed, skip to the next protein
        if os.path.exists(PSIBLAST_tar) and retry_successful == False:
            message = "{} PSIBLAST_tar exists, file skipped".format(acc)
            logging.info(message)
            continue
        if acc in failed_psiblast_list and retry_failed == False:
            message = "{} acc is in failed_psiblast_list, file skipped".format(acc)
            logging.info(message)
            continue
        # create folders if necessary
        utils.make_sure_path_exists(out_ascii_pssm, isfile=True)

        ########################################################################################
        #                                                                                      #
        #                        run the PSI-BLAST command-line argument                       #
        #                                                                                      #
        ########################################################################################

        # create full command string to be run, as if in the console
        c = command_str.format(psiblast_exec_str=s["psiblast_exec_str"], query=query, db=db_path, out_pssm=out_pssm,
                               out_ascii_pssm=out_ascii_pssm, out_BLAST_xml=out_BLAST_xml, evalue=s["evalue"],
                               inclusion_ethresh=s["inclusion_ethresh"], num_threads=s["num_threads"],
                               num_descriptions=s["num_descriptions"], num_alignments=s["num_alignments"])
        logging.info("\n{}\n".format(c))

        command = utils.Command(c)
        timeout_minutes = 120
        command.run(timeout=timeout_minutes * 60)
        # wait 1 second. In some cases, the files are not immediately recognised as existing?
        utils.sleep_x_seconds(1, print_stuff=False)

        ########################################################################################
        #                                                                                      #
        #               if successful, move output files into a tarball                        #
        #                                                                                      #
        ########################################################################################
        # check which output files exist (e.g. [True, True, False])
        output_file_exists_list = [os.path.exists(out_pssm), os.path.exists(out_ascii_pssm), os.path.exists(out_BLAST_xml)]
        logging.info("pssm, ascii_pssm, xml exists : {}".format(output_file_exists_list))
        # if all output files exist, create a date file and move all to a tarball
        if False not in output_file_exists_list:
            # create date file
            date = strftime("%Y%m%d")
            with open(date_file_path, "w") as f:
                f.write("Acc\t{}\nDate\t{}\nDatabase\t{}\nGreeting\tHave a nice day!".format(acc, date, db))
            # move all files into the tarball
            file_list = [out_pssm, out_ascii_pssm, out_BLAST_xml, date_file_path]
            with tarfile.open(PSIBLAST_tar, mode='w:gz') as tar:
                # add the files to the compressed tarfile
                logging.info('{} files will be moved into the tarball, original files deleted'.format(acc))
                for file in file_list:
                    tar.add(file, arcname=os.path.basename(file))
            # delete the original files
            for file in file_list:
                try:
                    os.remove(file)
                except FileNotFoundError:
                    pass
        else:
            if acc not in failed_psiblast_list:
                # add accession number to the list of failed blast sequences
                with open(failed_psiblast_list_txt, "a") as source:
                    source.write("\n{}".format(acc))