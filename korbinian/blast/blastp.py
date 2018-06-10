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
from multiprocessing import Pool
from Bio.Blast import NCBIWWW
from Bio.Blast.Applications import NcbiblastpCommandline
# import debugging tools
from korbinian.utils import pr, pc, pn, aaa

def run_BLAST_online(pathdict, s, logging):
    """From the list of proteins in csv format, begins online BLASTp searches for homologous proteins

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
    PROTEIN_NAME.blast_result.xml : xml file containing BLASTp hits
        (e.g. A2A2V5.blast_result.xml)
    """
    logging.info("~~~~~~~~~~~~                 starting running BLASTp via online API                 ~~~~~~~~~~~~")

    #IF blast directory doesn't exist -> create blast directory in the data(base) directory
    blast_dir = os.path.join(s["data_dir"], "blast")
    if not os.path.exists(blast_dir):
        os.makedirs(blast_dir)

    #Obtain protein data frame
    df = pd.read_csv(pathdict["list_csv"], sep = ",", quoting = csv.QUOTE_NONNUMERIC, index_col = 0)

    #Obtain blast settings from the settings file
    evalue = s["blast_Evalue"];
    hitsize = s["blast_max_hits"]
    overwrite_results = s["blast_overwrite_existing_results"]
    compress_results = s["blast_compress_results"]

    #preprocessing BLASTp execution by iterate over each protein
    task_list = []
    for acc in df.index:
        #Variable initializations for current protein
        protein_name = df.loc[acc, 'protein_name']
        input_sequence = df.loc[acc, 'full_seq']
        query = ">" + protein_name + "\n" + input_sequence + "\n"

        #Initialize file system in database/blast directory
        blast_proteinID_dir = os.path.join(s["data_dir"], "blast", protein_name[:2])
        if not os.path.exists(blast_proteinID_dir):
            os.makedirs(blast_proteinID_dir)
        output_file = os.path.join(s["data_dir"], "blast", protein_name[:2], protein_name + ".blast_result.xml")

        #Save query data into task_list
        task_list.append([protein_name, query, output_file, evalue, hitsize, overwrite_results, compress_results])

    #Execute BLASTp online searches
    if s["use_multiprocessing"]:
        # number of processes is the number the settings, or the number of proteins, whichever is smallest
        n_processes = s["multiprocessing_cores"] if s["multiprocessing_cores"] < len(task_list) else len(task_list)
        #TODO adviseable to just use 8 threads or something due to calculation will be done online and client is just waiting
        #TODO check NCBI BLAST API service how many tasks should be (are allowed to be) submitted in a second
        #TODO side note: NCBI entrez has 5 per second without an account and 10 with an account for example

        #Multithreaded execution
        with Pool(processes=n_processes) as thread_controller:
            thread_controller.map(BLAST_online_submission, task_list)
    else:
        #Sequential execution
        for task in task_list:
            BLAST_online_submission(task)

    logging.info("~~~~~~~~~~~~                 finished BLASTp search                 ~~~~~~~~~~~~")

def BLAST_online_submission(task):
    """Run a single BLASTp online search

    Parameters
    ----------
    task : array
        List/array which contains the query data and the parameters for the BLAST search

    Saved Files and Figures
    -----------------------
    PROTEIN_NAME.blast_result.xml : xml file containing BLASTp hits
        (e.g. A2A2V5.blast_result.xml)
    """
    #Obtain query and parameters
    protein_name, query, output_file, evalue, hitsize, overwrite_results, compress_results = task

    #IF overwrite existing results is FALSE
    if not overwrite_results:
        #check if compressed file exist and is not empty: IF yes return and don't execute BLASTp
        if  os.path.exists(output_file + ".gz") and os.path.getsize(output_file + ".gz") > 0:
            #IF compress results is TURE -> return and don't execute BLASTp
            if compress_results:
                logging.info("Skipped BLASTp online search for protein:" + "\t" + protein_name)
                return
            #ELSE decompress file, return and don't execute BLASTp
            else:
                with gzip.open(output_file + ".gz", 'rb') as blast_result_in, open(output_file, 'wb') as blast_result_out:
                    blast_result_out.write(blast_result_in.read())
                os.remove(output_file + ".gz")
                logging.info("Skipped BLASTp online search for protein:" + "\t" + protein_name)
                return
        #check if non compressed file exist and is not empty:
        if  os.path.exists(output_file) and os.path.getsize(output_file) > 0:
            #IF compress results is TURE -> compress the existing result file, return and don't execute BLASTp
            if compress_results:
                with open(output_file, "rb") as blast_result_in, gzip.open(output_file + ".gz", 'wb') as blast_result_out:
                    blast_result_out.writelines(blast_result_in)
                os.remove(output_file)
                logging.info("Skipped BLASTp online search for protein:" + "\t" + protein_name)
                return
            #ELSE: return and don't execute BLASTp
            else:
                logging.info("Skipped BLASTp online search for protein:" + "\t" + protein_name)
                return
    #ELSE remove existing results
    else:
        #remove raw results
        if os.path.exists(output_file):
            os.remove(output_file)
        #remove compressed results
        if os.path.exists(output_file + ".gz"):
            os.remove(output_file + ".gz")


    #Run BLASTp search
    logging.info("Run BLASTp online search for protein:" + "\t" + protein_name)
    blast_result = NCBIWWW.qblast("blastp", "nr", query, expect=evalue,  hitlist_size=hitsize)

    #Output results
    #Write BLASTp result with a gzip compression into file
    if compress_results:
        with gzip.open(output_file + ".gz", 'wb') as blast_result_writer:
            blast_result_writer.write(blast_result.read().encode('utf-8'))
    #Write BLASTp result directly into file
    else:
        with open(output_file, "w") as blast_result_writer:
            blast_result_writer.write(blast_result.read())








def run_BLAST_local(pathdict, s, logging):
    """From the list of proteins in csv format, begins local BLASTp searches for homologous proteins

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
    PROTEIN_NAME.blast_result.xml : xml file containing BLASTp hits
        (e.g. A2A2V5.blast_result.xml)
    """
    logging.info("~~~~~~~~~~~~                 starting running BLASTp locally                 ~~~~~~~~~~~~")

    #IF blast directory doesn't exist -> create blast directory in the data(base) directory
    blast_dir = os.path.join(s["data_dir"], "blast")
    if not os.path.exists(blast_dir):
        os.makedirs(blast_dir)

    #Obtain protein data frame
    df = pd.read_csv(pathdict["list_csv"], sep = ",", quoting = csv.QUOTE_NONNUMERIC, index_col = 0)

    #Obtain blast settings from the settings file
    evalue = s["blast_Evalue"];
    hitsize = s["blast_max_hits"]
    database = s["BLAST_local_DB"]
    overwrite_results = s["blast_overwrite_existing_results"]
    compress_results = s["blast_compress_results"]

    #preprocessing BLASTp execution by iterate over each protein
    task_list = []
    for acc in df.index:
        #Variable initializations for current protein
        protein_name = df.loc[acc, 'protein_name']
        input_sequence = df.loc[acc, 'full_seq']
        query = ">" + protein_name + "\n" + input_sequence + "\n"

        #Initialize file system in database/blast directory
        blast_proteinID_dir = os.path.join(s["data_dir"], "blast", protein_name[:2])
        if not os.path.exists(blast_proteinID_dir):
            os.makedirs(blast_proteinID_dir)
        output_file = os.path.join(s["data_dir"], "blast", protein_name[:2], protein_name + ".blast_result.xml")

        #Save query data into task_list
        task_list.append([protein_name, query, output_file, evalue, hitsize, database, overwrite_results, compress_results])

    #Execute BLASTp local searches
    if s["use_multiprocessing"]:
        # number of processes is the number the settings, or the number of proteins, whichever is smallest
        n_processes = s["multiprocessing_cores"] if s["multiprocessing_cores"] < len(task_list) else len(task_list)

        #Multithreaded execution
        with Pool(processes=n_processes) as thread_controller:
            thread_controller.map(BLAST_local_submission, task_list)
    else:
        #Sequential execution
        for task in task_list:
            BLAST_local_submission(task)

    logging.info("~~~~~~~~~~~~                 finished BLASTp search                 ~~~~~~~~~~~~")

def BLAST_local_submission(task):
    """Run a single BLASTp local search

    Parameters
    ----------
    task : array
        List/array which contains the query data and the parameters for the BLAST search

    Saved Files and Figures
    -----------------------
    PROTEIN_NAME.blast_result.xml : xml file containing BLASTp hits
        (e.g. A2A2V5.blast_result.xml)
    """
    #Obtain query and parameters
    protein_name, query, output_file, evalue, hitsize, database, overwrite_results, compress_results = task

    #IF overwrite existing results is FALSE
    if not overwrite_results:
        #check if compressed file exist and is not empty: IF yes return and don't execute BLASTp
        if  os.path.exists(output_file + ".gz") and os.path.getsize(output_file + ".gz") > 0:
            #IF compress results is TURE -> return and don't execute BLASTp
            if compress_results:
                logging.info("Skipped BLASTp local search for protein:" + "\t" + protein_name)
                return
            #ELSE decompress file, return and don't execute BLASTp
            else:
                with gzip.open(output_file + ".gz", 'rb') as blast_result_in, open(output_file, 'wb') as blast_result_out:
                    blast_result_out.write(blast_result_in.read())
                os.remove(output_file + ".gz")
                logging.info("Skipped BLASTp local search for protein:" + "\t" + protein_name)
                return
        #check if non compressed file exist and is not empty:
        if  os.path.exists(output_file) and os.path.getsize(output_file) > 0:
            #IF compress results is TURE -> compress the existing result file, return and don't execute BLASTp
            if compress_results:
                with open(output_file, "rb") as blast_result_in, gzip.open(output_file + ".gz", 'wb') as blast_result_out:
                    blast_result_out.writelines(blast_result_in)
                os.remove(output_file)
                logging.info("Skipped BLASTp local search for protein:" + "\t" + protein_name)
                return
            #ELSE: return and don't execute BLASTp
            else:
                logging.info("Skipped BLASTp local search for protein:" + "\t" + protein_name)
                return
    #ELSE remove existing results
    else:
        #remove raw results
        if os.path.exists(output_file):
            os.remove(output_file)
        #remove compressed results
        if os.path.exists(output_file + ".gz"):
            os.remove(output_file + ".gz")

    #Run BLASTp search
    logging.info("Run BLASTp local search for protein:" + "\t" + protein_name)
    blastp_cline = NcbiblastpCommandline(db=database, evalue=evalue, max_target_seqs=hitsize, outfmt=5, out=output_file)
    out, err = blastp_cline(stdin=query)

    #Exception handling of BLAST execution
    if out or err:
        logging.warning(out)
        logging.warning(err)

    #IF compress results is TURE -> compress the blast xml result file
    if compress_results:
        with open(output_file, "rb") as blast_result_in, gzip.open(output_file + ".gz", 'wb') as blast_result_out:
            blast_result_out.writelines(blast_result_in)
        os.remove(output_file)
