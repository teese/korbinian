"""
Author:         Dominik Müller
Created:        May 7 00:06 2018
Dependencies:   Python >=3.3
                pandas
                tmseg
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
import korbinian
import sys
import subprocess
import shutil
import gzip
import tarfile
import pandas as pd
# import debugging tools
from korbinian.utils import pr, pc, pn, aaa

def run_filtering(pathdict, s, logging):
    """From the list of proteins in csv format, execeute tmseg transmembrane protein prediction
    and exclude all proteins which aren't predicted as a membrane protein.

    Parameters
    ----------
    pathdict : dict
        Dictionary of the key paths and files associated with that List number.
    s : dict
        Settings dictionary extracted from excel settings file.
    logging : logging.Logger
        Logger for printing to console and logfile.
    """
    logging.info("~~~~~~~~~~~~                 starting running filtering: tmseg                 ~~~~~~~~~~~~")

    #create tmseg output directory and create temporary FASTA & PSSM directories
    logging.info("TMSEG filtering: Preprocessing")
    tmseg_out_dir = os.path.join(pathdict['TMSEG_dir'], "output")
    if not os.path.exists(tmseg_out_dir):
        os.makedirs(tmseg_out_dir)
    ##create/link the pssm files for all TM proteins
    tmseg_pssm_dir = create_PSSM(pathdict, s["data_dir"])
    ##create fasta files for all TM protein sequences
    tmseg_fasta_dir = create_FASTA(pathdict)

    #execute tmseg
    logging.info("TMSEG filtering: Run TMSEG")
    out = subprocess.check_output(["java", "-jar", s["TMSEG_local"], "-i", tmseg_fasta_dir, "-p", tmseg_pssm_dir,
                                    "-o", tmseg_out_dir, "-m"])

    #gather predictions out of the tmseg results directory
    logging.info("TMSEG filtering: Gather results")
    collect_predictions(pathdict, tmseg_out_dir)

    logging.info("TMSEG filtering: Clean up temporary files and compress TMSEG result files")
    #delete temporary FASTA & PSSM directories
    shutil.rmtree(tmseg_pssm_dir)
    shutil.rmtree(tmseg_fasta_dir)
    #compress tmseg results directory
    with tarfile.open(tmseg_out_dir + ".tar.gz", "w:gz") as tar:
        tar.add(tmseg_out_dir, arcname=os.path.basename(tmseg_out_dir))
    shutil.rmtree(tmseg_out_dir)

    logging.info("~~~~~~~~~~~~                 finished filtering: tmseg                 ~~~~~~~~~~~~")


#create the fasta directory for TMSEG (for each protein a fasta file)
def create_FASTA(pathdict):
    # load list parsed from uniprot
    df = pd.read_csv(pathdict["list_parsed_csv"], sep=",", quoting=csv.QUOTE_NONNUMERIC, index_col=0, low_memory=False)
    #make folder for tmseg output
    tmseg_dir = pathdict['TMSEG_dir']
    if not os.path.exists(tmseg_dir):
        os.makedirs(tmseg_dir)
    #make folder for fasta output
    tmseg_fasta_dir = os.path.join(tmseg_dir, "fasta")
    if not os.path.exists(tmseg_fasta_dir):
        os.makedirs(tmseg_fasta_dir)

    #iterate over each protein
    for acc in df.index:
        #check if pssm file exist -> else skip
        if not os.path.exists(os.path.join(tmseg_dir, "pssm", acc + ".pssm")):
            continue
        #specify outfile path
        outfile = os.path.join(tmseg_fasta_dir, acc + ".fasta")
        with open(outfile, 'w') as fasta_writer:
            fasta_writer.write('>{}\n{}\n'.format(acc, df.loc[acc, 'full_seq']))
    return tmseg_fasta_dir

#create the pssm directory for TMSEG (for each protein a pssm file)
def create_PSSM(pathdict, data_dir):
    # load list parsed from uniprot
    df = pd.read_csv(pathdict["list_parsed_csv"], sep=",", quoting=csv.QUOTE_NONNUMERIC, index_col=0, low_memory=False)
    #make folder for tmseg output
    tmseg_dir = pathdict['TMSEG_dir']
    if not os.path.exists(tmseg_dir):
        os.makedirs(tmseg_dir)
    #make folder for pssm output
    tmseg_pssm_dir = os.path.join(tmseg_dir, "pssm")
    if not os.path.exists(tmseg_pssm_dir):
        os.makedirs(tmseg_pssm_dir)

    #iterate over each protein
    for acc in df.index:
        compressed = None
        #define path to the pssm file obtained from PSI-BLAST
        pssmFile_original = os.path.join(data_dir, "psiblast", acc[:2], acc + ".blast_result.pssm")
        #check if non compressed file exist and is not empty:
        if os.path.exists(pssmFile_original) and os.path.getsize(pssmFile_original) > 0:
            compressed = False
        #check if compressed file exist and is not empty:
        elif os.path.exists(pssmFile_original + ".gz") and os.path.getsize(pssmFile_original + ".gz") > 0:
            compressed = True
        #if none pssm file exist for this protein -> skip
        else:
            continue

        #define output path for pssm
        pssmFile_out = os.path.join(tmseg_pssm_dir, acc + ".pssm")
        #IF pssm is compressed -> decompress it at correct location
        if compressed:
            with gzip.open(pssmFile_original + ".gz", 'rb') as file_compressed:
                with open(pssmFile_out, 'wb') as file_decompressed:
                    shutil.copyfileobj(file_compressed, file_decompressed)
        #IF pssm is not compressed -> create a system link at correct location
        else:
            if not os.path.exists(pssmFile_out):
                os.symlink(pssmFile_original, pssmFile_out)
    return tmseg_pssm_dir


#collect all predictions from tmseg result direcotry
def collect_predictions(pathdict, tmseg_out_dir):
    #Iterate over each TMSEG result file
    result_collection = {}
    result_file_list = os.listdir(tmseg_out_dir)
    for result_file in result_file_list:
        #obtain protein accession ID
        protein_id = result_file.split(".")[0]
        #define path to result file
        rf = os.path.join(tmseg_out_dir, result_file)
        #Read TMSEG result file
        with open(rf, 'r') as rf_reader:
            isTM = False
            #Iterate over each line in the file
            for line in rf_reader:
                #Remove new line character
                if "\n" in line:
                    line = line.rstrip("\n")
                #IF file starts with '#' (header of result file) and contains 'TRANSMEM' prediction -> TM protein
                if line.startswith("#") and "TRANSMEM" in line:
                    isTM = True
                    break
        #save collected information in result_collection hash
        result_collection[protein_id] = isTM

    #load protein list from parsed uniprot csv
    df = pd.read_csv(pathdict["list_parsed_csv"], sep=",", quoting=csv.QUOTE_NONNUMERIC, index_col=0, low_memory=False)
    #iterate over each protein from the parsed csv list
    output_list = []
    pred_counter = [0,0,0]
    for acc in df.index:
        prediction = None
        #check if a TMSEG prediction for the protein exist -> get prediction
        if acc in result_collection:
            if result_collection[acc]:
                prediction = "TM"
                pred_counter[0] += 1
            else:
                prediction = "noTM"
                pred_counter[1] += 1
        #ELSE assign X as prediction (probably no PSI-BLAST hits for the protein to obtain a PSSM)
        else:
            prediction = "-"
            pred_counter[2] += 1
        #save output of collected results
        output_list.append(acc + "\t" + prediction + "\n")
    #initialize a writer to output the collected results
    out_collected_results = os.path.join(pathdict["TMSEG_dir"], "tmseg.results.tsv")
    with open(out_collected_results, 'w') as ocr_writer:
        #output header
        ocr_writer.write("#Proteins as TM predicted (TM):" + "\t\t\t" + str(pred_counter[0]) + "\n")
        ocr_writer.write("#Proteins as non TM predicted (noTM):" + "\t\t" + str(pred_counter[1]) + "\n")
        ocr_writer.write("#Proteins with no prediction at all (-):" + "\t" + str(pred_counter[2]) + "\n")
        ocr_writer.write("#id" + "\t" + "prediction" + "\n")
        #output all collected results
        for result in output_list:
            ocr_writer.write(result)
