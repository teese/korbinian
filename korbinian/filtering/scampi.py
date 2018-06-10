"""
Author:         Dominik Müller
Created:        May 7 00:06 2018
Dependencies:   Python >=3.3
                pandas
                SCAMPI2
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
import pandas as pd
# import debugging tools
from korbinian.utils import pr, pc, pn, aaa

def run_filtering(pathdict, s, logging):
    """From the list of proteins in csv format, execeute SCAMPI2 transmembrane protein prediction
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
    logging.info("~~~~~~~~~~~~                 starting running filtering: SCAMPI2                 ~~~~~~~~~~~~")

    #create fasta file for all TM protein sequences
    path_fasta = create_FASTA(s["list_number"], pathdict)

    #execute SCAMPI2
    out = subprocess.call([s["SCAMPI_local"], path_fasta, pathdict["SCAMPI_top"]])

    #read & parse SCAMPI results and create query.nonTM_list.txt and query.TM_list.txt summary result files
    with open(pathdict["SCAMPI_top"], 'r') as scampi_result:
        #initialize lists
        noTM_list = []
        TM_list = []
        protein_id = None
        prediction = None
        #Iterate over each line in the SCAMPI results file
        for line in scampi_result:
            #Remove new line character
            if "\n" in line:
                line = line.rstrip("\n")
            #IF line is a entry header -> save protein ID
            if line.startswith(">"):
                protein_id = line[1:]
            #ELSE line is a prediction
            else:
                prediction = line
                #IF the prediction result contains a "M" for membrane prediction
                #-> add it to the TM_list
                if 'M' in prediction:
                    TM_list.append(protein_id)
                #ELSE -> add it to the nonTM_list
                else:
                    noTM_list.append(protein_id)

    #output noTM_list
    noTM_list_path = pathdict["SCAMPI_nonTM"]
    with open(noTM_list_path, 'w') as SCAMPI_noTM_list:
        for item in noTM_list:
            SCAMPI_noTM_list.write(item + "\n")
    #output TM_list
    TM_list_path = os.path.join(pathdict["SCAMPI_dir"], "query.TM_list.txt")
    with open(TM_list_path, 'w') as SCAMPI_TM_list:
        for item in TM_list:
            SCAMPI_TM_list.write(item + "\n")

    logging.info("~~~~~~~~~~~~                 finished filtering: SCAMPI2                 ~~~~~~~~~~~~")


#copied from korbinian.cons_ratio.SCAMPI.generate_scampi_input_files
def create_FASTA(list_number, pathdict):
    # load list parsed from uniprot
    df = pd.read_csv(pathdict["list_parsed_csv"], sep=",", quoting=csv.QUOTE_NONNUMERIC, index_col=0, low_memory=False)
    # specify outpath
    outpath = pathdict['SCAMPI_dir']
    # make folder for output
    if not os.path.exists(outpath):
        os.makedirs(outpath)
    # specify outfile path
    outfile = os.path.join(outpath, 'List{:02d}_fasta_for_scampi.txt'.format(list_number))
    # open new .txt file and write accession and full sequence from df into file
    file = open(outfile, 'w')
    for acc in df.index:
        file.write('>{}\n{}\n'.format(acc, df.loc[acc, 'full_seq']))
    file.close()
    #return resulting fasta file path
    return outfile
