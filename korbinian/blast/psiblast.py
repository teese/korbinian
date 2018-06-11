"""
Author:         Dominik Müller
Created:        May 7 00:06 2018
Dependencies:   Python >=3.3
                pandas
                Bio
                BLAST+
Purpose:        Protein Data Science
                Analysis of evolutionary sequences of transmembrane proteins
Credits:        Mark Teese
                Martin Ortner
                Shenger Wang
                Rimma Jenske
                Dominik Müller
License         Released under the permissive MIT license.
"""
import csv
import logging
import os
import gzip
import korbinian
import pandas as pd
import sys
import gzip
import shutil
import subprocess
from multiprocessing import Pool
from Bio.Blast.Applications import NcbipsiblastCommandline
from korbinian.utils import downloaderInterface
from tqdm import tqdm

def create_BLAST_database(pathdict, s, logging):
    """Automatically downloads UniRef50 fasta from the UniProt FTP server and create a BLAST database out of it.

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
    Uniref50.fasta
    BLAST database of the UniRef50.fasta
    """
    logging.info("~~~~~~~~~~~~                 starting running PSI-BLAST database creation                 ~~~~~~~~~~~~")

    #IF psiblast directory doesn't exist -> create psiblast directory in the data(base) directory
    blast_dir = os.path.join(s["data_dir"], "psiblast")
    if not os.path.exists(blast_dir):
        os.makedirs(blast_dir)
    #IF psiblast/database directory doesn't exist -> create psiblast/database directory in the data(base) directory
    db_dir = os.path.join(s["data_dir"], "psiblast", "database")
    if not os.path.exists(db_dir):
        os.makedirs(db_dir)

    #Download UniRef50 from the UniProt FTP server
    uniref50_FTP = "ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/uniref50.fasta.gz"
    uniref50_path_compressed = os.path.join(db_dir, "uniref50.fasta.gz")
    logging.info("PSI-BLAST database creation: Start downloading UniRef50 fasta")
    downloaderInterface(uniref50_FTP, uniref50_path_compressed, True)
    logging.info("PSI-BLAST database creation: Finished downloading UniRef50 fasta")

    #Decompress the FASTA file
    uniref50_path = os.path.join(db_dir, "uniref50.fasta")
    logging.info("PSI-BLAST database creation: Start decompressing the fasta file")
    with gzip.open(uniref50_path_compressed, 'rb') as file_compressed:
        with open(uniref50_path, 'wb') as file_decompressed:
            shutil.copyfileobj(file_compressed, file_decompressed)
    logging.info("PSI-BLAST database creation: Finished decompressing the fasta file")

    #Create the BLAST database
    logging.info("PSI-BLAST database creation: Start creating the BLAST database")
    uniref50_DB = os.path.join(db_dir, "uniref50")
    subprocess.call(["makeblastdb", "-in", uniref50_path, "-out", uniref50_DB, "-dbtype", "prot"])
    logging.info("PSI-BLAST database creation: Finished creating the BLAST database")

    #Delete the decompressed fasta file again
    os.remove(uniref50_path)

    logging.info("~~~~~~~~~~~~                 finished PSI-BLAST database creation                 ~~~~~~~~~~~~")

def run_BLAST(pathdict, s, logging):
    """From the list of proteins in csv format, begins PSI-BLAST searches for homologous proteins

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
    PROTEIN_NAME.blast_result.xml : xml file containing PSI-BLAST hits
        (e.g. A2A2V5.blast_result.xml)
    """
    logging.info("~~~~~~~~~~~~                 starting running PSI-BLAST                 ~~~~~~~~~~~~")

    #IF blast directory doesn't exist -> create blast directory in the data(base) directory
    blast_dir = os.path.join(s["data_dir"], "psiblast")
    if not os.path.exists(blast_dir):
        os.makedirs(blast_dir)

    #Obtain protein data frame
    df = pd.read_csv(pathdict["list_csv"], sep = ",", quoting = csv.QUOTE_NONNUMERIC, index_col = 0)

    #Obtain blast settings from the settings file
    evalue = s["blast_Evalue"];
    hitsize = s["blast_max_hits"]
    overwrite_results = s["blast_overwrite_existing_results"]
    compress_results = s["blast_compress_results"]
    #The remote option for the PSI-BLAST is deactivated due to the newest BLAST+ version have a different
    #pssm format which TMSEQ can't handle/read
    #remote_tag = s["TMSEQ_PSI-BLAST_remote"]
    remote_tag = False

    #Identify BLAST database
    database = None
    if s["TMSEQ_PSI-BLAST_database"] == "DEFAULT":
        database = os.path.join(s["data_dir"], "psiblast", "database", "uniref50")
    else:
        database = s["TMSEQ_PSI-BLAST_database"]

    #preprocessing PSI-BLAST execution by iterate over each protein
    task_list = []
    for acc in df.index:
        #Variable initializations for current protein
        protein_name = df.loc[acc, 'protein_name']
        input_sequence = df.loc[acc, 'full_seq']
        query = ">" + protein_name + "\n" + input_sequence + "\n"

        #Initialize file system in database/blast directory
        blast_proteinID_dir = os.path.join(s["data_dir"], "psiblast", protein_name[:2])
        if not os.path.exists(blast_proteinID_dir):
            os.makedirs(blast_proteinID_dir)
        output_hit_file = os.path.join(s["data_dir"], "psiblast", protein_name[:2], protein_name + ".blast_result.xml")
        output_pssm_file = os.path.join(s["data_dir"], "psiblast", protein_name[:2], protein_name + ".blast_result.pssm")

        #Save query data into task_list
        task_list.append([protein_name, query, output_hit_file, output_pssm_file,
            evalue, hitsize, overwrite_results, database, remote_tag, compress_results])

    #Execute BLASTp online searches
    if s["use_multiprocessing"]:
        # number of processes is the number the settings, or the number of proteins, whichever is smallest
        n_processes = s["multiprocessing_cores"] if s["multiprocessing_cores"] < len(task_list) else len(task_list)
        #TODO adviseable to just use 8 threads or something if the calculation will be done online and client is just waiting
        #TODO check NCBI BLAST API service how many tasks should be (are allowed to be) submitted in a second
        #TODO side note: NCBI entrez has 5 per second without an account and 10 with an account for example

        #Multithreaded execution
        with Pool(processes=n_processes) as thread_controller:
            list(tqdm(thread_controller.imap(BLAST_submission, task_list), total=len(task_list)))
    else:
        #Sequential execution
        for task in task_list:
            BLAST_submission(task)

    logging.info("~~~~~~~~~~~~                 finished PSI-BLAST search                 ~~~~~~~~~~~~")

def BLAST_submission(task):
    """Run a single PSI-BLAST search

    Parameters
    ----------
    task : array
        List/array which contains the query data and the parameters for the BLAST search

    Saved Files and Figures
    -----------------------
    PROTEIN_NAME.blast_result.xml : xml file containing PSI-BLAST hits
        (e.g. A2A2V5.blast_result.xml)
    PROTIEN_NAME.blast_result.pssm : pssm file containing PSI-BLAST similiarity scoring matrix
        (e.g. A2A2V5.blast_result.pssm)
    """
    #Obtain query and parameters
    protein_name, query, output_hit_file, output_pssm_file, evalue, hitsize, overwrite_results, database,\
        remote_tag, compress_results = task
    #IF overwrite existing results is FALSE
    if not overwrite_results:
        #check if compressed file exist and is not empty: IF yes return and don't execute PSI-BLAST
        if  os.path.exists(output_hit_file + ".gz") and os.path.getsize(output_hit_file + ".gz") > 0:
            #IF compress results is TURE -> return and don't execute PSI-BLAST
            if compress_results:
                #logging.warning("Skipped PSI-BLAST search for protein:" + "\t" + protein_name)
                return
            #ELSE decompress file, return and don't execute PSI-BLAST
            else:
                with gzip.open(output_hit_file + ".gz", 'rb') as blast_result_in, open(output_hit_file, 'wb') as blast_result_out:
                    blast_result_out.write(blast_result_in.read())
                os.remove(output_hit_file + ".gz")
                with gzip.open(output_pssm_file + ".gz", 'rb') as blast_result_in, open(output_pssm_file, 'wb') as blast_result_out:
                    blast_result_out.write(blast_result_in.read())
                os.remove(output_pssm_file + ".gz")
                #logging.warning("Skipped PSI-BLAST search for protein:" + "\t" + protein_name)
                return
        #check if non compressed file exist and is not empty:
        if  os.path.exists(output_hit_file) and os.path.getsize(output_hit_file) > 0:
            #IF compress results is TURE -> compress the existing result file, return and don't execute PSI-BLAST
            if compress_results:
                with open(output_hit_file, "rb") as blast_result_in, gzip.open(output_hit_file + ".gz", 'wb') as blast_result_out:
                    blast_result_out.writelines(blast_result_in)
                os.remove(output_hit_file)
                with open(output_pssm_file, "rb") as blast_result_in, gzip.open(output_pssm_file + ".gz", 'wb') as blast_result_out:
                    blast_result_out.writelines(blast_result_in)
                os.remove(output_pssm_file)
                #logging.warning("Skipped PSI-BLAST search for protein:" + "\t" + protein_name)
                return
            #ELSE: return and don't execute PSI-BLAST
            else:
                #logging.warning("Skipped PSI-BLAST search for protein:" + "\t" + protein_name)
                return
    #ELSE remove existing results
    else:
        #remove raw results
        if os.path.exists(output_hit_file):
            os.remove(output_file)
        if os.path.exists(output_pssm_file):
            os.remove(output_pssm_file)
        #remove compressed results
        if os.path.exists(output_hit_file + ".gz"):
            os.remove(output_hit_file + ".gz")
        if os.path.exists(output_pssm_file + ".gz"):
            os.remove(output_pssm_file + ".gz")


    #Run PSI-BLAST search
    #logging.info("Run PSI-BLAST search for protein:" + "\t" + protein_name)

    #TODO: Add -out parameter for the hits and out_ascii_pssm for the pssm file
    #TODO: Both compressed if required
    psiblast_cline = NcbipsiblastCommandline('psiblast', db=database, evalue=evalue, max_target_seqs=hitsize, outfmt=5,
                                            out=output_hit_file, out_ascii_pssm=output_pssm_file, remote=remote_tag,
                                            inclusion_ethresh=evalue, num_iterations=3, use_sw_tback=True, seg="no")
    out, err = psiblast_cline(stdin=query)

    #Exception handling of BLAST execution
    if out or err:
        logging.warning(out)
        logging.warning(err)
        return

    #IF compress results is TURE -> compress the blast xml result file
    if compress_results:
        if os.path.exists(output_hit_file):
            with open(output_hit_file, "rb") as blast_result_in, gzip.open(output_hit_file + ".gz", 'wb') as blast_result_out:
                blast_result_out.writelines(blast_result_in)
            os.remove(output_hit_file)
        if os.path.exists(output_pssm_file):
            with open(output_pssm_file, "rb") as blast_result_in, gzip.open(output_pssm_file + ".gz", 'wb') as blast_result_out:
                blast_result_out.writelines(blast_result_in)
            os.remove(output_pssm_file)
