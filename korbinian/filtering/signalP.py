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
import re
import korbinian
import sys
import subprocess
import pandas as pd
# import debugging tools
from korbinian.utils import pr, pc, pn, aaa

def run_filtering(pathdict, s, logging):
    """From the list of proteins in csv format, execeute SignalP transmembrane protein prediction
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
    logging.info("~~~~~~~~~~~~                 starting running filtering: SignalP                 ~~~~~~~~~~~~")

    #create fasta file for all TM protein sequences
    path_fasta = create_FASTA(s["list_number"], pathdict)

    #parse provided SignalP parameters
    organism = s["SignalP_organism"]
    cutoff_noTM = s["SignalP_cutoff_noTM_networks"]
    cutoff_TM = s["SignalP_cutoff_TM_networks"]

    #execute SignalP
    out = subprocess.check_output([s["SignalP_local"], "-t", organism, "-u", cutoff_noTM, "-U", cutoff_TM, "-n",
                                    pathdict["SignalP_SiPe_acc"], path_fasta])

    #save results into file
    with open(os.path.join(pathdict["SignalP_dir"], "output.txt"), "wb") as signalP_out:
        signalP_out.write(out)

    #save summary of results into a summary file
    out_string = out.decode("utf-8")
    with open(os.path.join(pathdict["SignalP_dir"], "SignalP.results.tsv"), "w") as signalP_summary:
        #output header to summary file
        signalP_summary.write("#id" + "\t" + "signal" + "\t" + "used_network" + "\n")
        #iterate over the SignalP output
        for line in out_string.split("\n"):
            #IF line starts with '#' -> header -> skip line
            if line.startswith("#"):
                continue
            #parse format
            #INFO: I know this is an ugly regex :D I have to look someday if you can capture repeating groups in python somehow
            parser = re.search("^([A-Za-z0-9\_]+)\s+([0-9\.]+)\s+([0-9\.]+)\s+([0-9\.]+)\s+([0-9\.]+)\s+([0-9\.]+)\s+"\
                                + "([0-9\.]+)\s+([0-9\.]+)\s+([0-9\.]+)\s+([YN])\s+([0-9\.]+)\s+([A-Za-z\-]+)", line)
            if parser:
                protein_id =  parser.group(1)
                signal = parser.group(10)
                network = parser.group(12).split("-")[1]
                signal_boolean = None
                if signal == "Y":
                    signal_boolean = True
                else:
                    signal_boolean = False
                signalP_summary.write(protein_id + "\t" + str(signal_boolean) + "\t" + network + "\n")
            elif line != "":
                logging.warning("SignalP filtering: Parsing error of the output file! Check the output.txt file.")

    logging.info("~~~~~~~~~~~~                 finished filtering: SignalP                 ~~~~~~~~~~~~")


#copied from korbinian.cons_ratio.SCAMPI.generate_scampi_input_files
def create_FASTA(list_number, pathdict):
    # load list parsed from uniprot
    df = pd.read_csv(pathdict["list_parsed_csv"], sep=",", quoting=csv.QUOTE_NONNUMERIC, index_col=0, low_memory=False)
    # specify outpath
    outpath = pathdict['SignalP_dir']
    # make folder for output
    if not os.path.exists(outpath):
        os.makedirs(outpath)
    # specify outfile path
    outfile = os.path.join(outpath, 'List{:02d}_fasta_for_SignalP.txt'.format(list_number))
    # open new .txt file and write accession and full sequence from df into file
    file = open(outfile, 'w')
    for acc in df.index:
        file.write('>{}\n{}\n'.format(acc, df.loc[acc, 'full_seq']))
    file.close()
    #return resulting fasta file path
    return outfile
