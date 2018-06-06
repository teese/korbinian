"""
Author:         Dominik Müller
Created:        May 7 14:28 2018
Dependencies:   Python >=3.3
                pandas
                Bio
Purpose:        Protein Data Science
                Analysis of evolutionary sequences of transmembrane proteins
Credits:        Mark Teese
                Martin Ortner
                Shenger Wang
                Rimma Jenske
                Dominik Müller
License         Released under the permissive MIT license.
"""
import pandas as pd
from Bio.Blast import NCBIXML
import ast
import pickle
import zipfile
import csv
import gzip
import logging
import re
import os
import sys
from multiprocessing import Pool
import korbinian
import korbinian.utils as utils
# import debugging tools
from korbinian.utils import pr, pc, pn, aaa

def run(pathdict, s, logging):
    """ Runner of the BLAST xml results parser.

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
    homol_df_orig_zip : zipfile
        Zipfile containing the following:
            homol_df_orig_pickle : pickled pd.DataFrame
                Dataframe containing all sequence extracted from the XML file.
                This can be large, as it contains the full query, markup and match sequences
    """
    #Initialize protein list
    acc_not_in_homol_db = []
    p_dict_logging = logging if s["use_multiprocessing"] != True else utils.Log_Only_To_Console()
    list_p = korbinian.utils.convert_summary_csv_to_input_list(s, pathdict, p_dict_logging, list_excluded_acc=acc_not_in_homol_db)

    logging.info("~~~~~~~~~~~~                 starting parsing BLAST results                 ~~~~~~~~~~~~")

    #Multithreading
    if s["use_multiprocessing"]:
        # number of processes is the number the settings, or the number of proteins, whichever is smallest
        n_processes = s["multiprocessing_cores"] if s["multiprocessing_cores"] < len(list_p) else len(list_p)

        #Scatter list_p to threads and start parallel execution
        with Pool(processes=n_processes) as thread_controller:
            thread_controller.map(parse_blast_result, list_p)
    #Sequential execution
    else:
        for p in list_p:
            parse_blast_result(p)

    logging.info("\n" + "~~~~~~~~~~~~                 finished parsing BLAST results                 ~~~~~~~~~~~~")

def parse_blast_result(p):
    """ Parses the BLSAT XML file to csv for a single protein.

    Parameters
    ----------
    p : dict
        Protein Dictionary. Contains all input settings, sequences and filepaths related to a single protein.
        Protein-specific data is extracted from one row of the the list summary, e.g. List05_summary.csv, which is read as df.
        p also contains the GENERAL korbinian settings and filepaths for that list (pathdict, s, logging)

        Components
        ----------
        pathdict : dict
            Dictionary of the key paths and files associated with that List number.
        s : dict
            Settings dictionary extracted from excel settings file.
        logging : logging.Logger
            Logger for printing to console and/or logfile.
        p : protein-specific dictionary components
            acc, list_of_TMDs, description, TM01_seq, etc

    Saved Files and Figures
    -----------------------
    homol_df_orig_zip : zipfile
        Zipfile containing the following:
            homol_df_orig_pickle : pickled pd.DataFrame
                Dataframe containing all sequence extracted from the XML file.
                This can be large, as it contains the full query, markup and match sequences

    Returns
    -------
    In all cases, a tuple (str, bool, str) is returned.
    if sucsessful:
        return acc, True, "0"
    if not successful:
        return acc, False, "specific warning or reason why protein failed"
    """
    #Variable intializations for the current protein
    s, logging = p["s"], p["logging"]
    acc = p["acc"]
    protein_name = p['protein_name']
    blast_xml_path = os.path.join(s["data_dir"], "blast", protein_name[:2], protein_name + ".blast_result.xml")

    #if BLAST_xml file is empty or don't exist -> give warning and return
    if  not (os.path.exists(blast_xml_path) and os.path.getsize(blast_xml_path) > 0) \
        and not (os.path.exists(blast_xml_path + ".gz") and os.path.getsize(blast_xml_path + ".gz") > 0):
        logging.info("BLAST_parser:" + "\t" + protein_name + " BLAST results file don't exist or is empty")
        return

    #Identify compression state of blast xml results file
    blast_xml_compressed = False
    if os.path.exists(blast_xml_path + ".gz"):
        blast_xml_compressed = True

    #Initialize BLAST xml result reader
    #IF compressed -> open gzip reader
    if blast_xml_compressed:
        openBLAST = gzip.open
        blast_xml_path = blast_xml_path + ".gz"
    #ELSE -> open normal reader
    else:
        openBLAST = open

    #Read BLAST result xml file and parse required information into the match_details_dict Dictionary
    #BLAST parser structure was copied from the THOIPApy software
    #Github source: https://github.com/bojigu/thoipapy/blob/master/thoipapy/homologues/NCBI_parser.py
    BLAST_xml_file = blast_xml_path
    BLAST_csv_file = os.path.join(s["data_dir"], "blast", protein_name[:2], protein_name + ".blast_result.TMP.csv")
    match_details_dict = {}
    with open(BLAST_csv_file, 'w') as homo_out_csv_file_handle:
        with openBLAST(BLAST_xml_file) as xml_result_handle:
            xml_record = NCBIXML.read(xml_result_handle)
            hit_num = 0
            #Iterate over each alignment in a blast result file
            for alignment in xml_record.alignments:
                #Iterate over each high-scoring segment pair of an alignment
                for hsp in alignment.hsps:
                    match_details_dict['hit_num'] = hit_num
                    #Description of protein
                    if (hit_num) == 0:
                        description = "%s_NCBI_query_sequence" % acc
                    else:
                        description = alignment.title
                    match_details_dict["description"] = description
                    #Identify homolog database ID
                    datebaseID_search = re.search('\|ref\|(.*?)\|', alignment.title)
                    if datebaseID_search:
                        match_details_dict["databaseId"] = datebaseID_search.group(1)
                    else:
                        match_details_dict["databaseId"] = "-"
                    #Identify organism of homolog
                    taxonomy = re.search('\[(.*?)\]', alignment.title)
                    if taxonomy:
                        taxonomyNode = taxonomy.group(1)
                        match_details_dict["organism"] = taxonomyNode
                    else:
                        match_details_dict["organism"] = "no_organism"

                    #Subject/Match sequence length
                    match_details_dict['len_full_match_seq'] = alignment.length

                    #Alignment quality scores
                    match_details_dict["blast_e_value"] = hsp.expect
                    match_details_dict["blast_identity"] = hsp.identities / hsp.align_length
                    match_details_dict["blast_positives"] = hsp.positives / hsp.align_length

                    #Aligned sequences
                    match_details_dict['query_align_seq'] = hsp.query
                    match_details_dict['match_align_seq'] = hsp.sbjct

                    #Alignment markup (modificated to match up with the korbinian pipeline)
                    markup = re.sub("\+", ':', hsp.match)
                    markup = re.sub("[A-Za-z]", '|', markup)
                    match_details_dict['align_markup_seq'] = markup

                    #Calculate coverage for query and subject
                    match_details_dict['blast_query_coverage'] = len(re.sub('-', '', hsp.query)) / p['seqlen']
                    match_details_dict['blast_match_coverage'] = len(re.sub('-', '', hsp.sbjct)) / alignment.length

                    #Necessary korbinian pipeline variables
                    match_details_dict['hit_contains_SW_node'] = True
                    match_details_dict['SW_query_coverage'] = match_details_dict['blast_query_coverage']
                    match_details_dict['SW_match_coverage'] = match_details_dict['blast_match_coverage']

                    #Calculate gapped identity to get the number of  "observed_changes"
                    gapped_identity = hsp.identities / (hsp.align_length - hsp.gaps)
                    match_details_dict['obs_changes'] = 100 - (gapped_identity*100)

                    #Assign necessary FASTA alignment variables for the korbinian pipeline with BLAST alignment results
                    match_details_dict['FASTA_identity'] = match_details_dict["blast_identity"]
                    match_details_dict['FASTA_gapped_identity'] = gapped_identity
                    match_details_dict['FASTA_expectation'] = match_details_dict["blast_e_value"]
                    match_details_dict['FASTA_overlap'] = alignment.length
                    match_details_dict['FASTA_query_coverage'] = match_details_dict['blast_query_coverage']
                    match_details_dict['FASTA_match_coverage'] = match_details_dict['blast_match_coverage']
                    match_details_dict['FASTA_query_start'] = hsp.query_start
                    match_details_dict['FASTA_query_end'] = hsp.query_end
                    match_details_dict['FASTA_match_start'] = hsp.sbjct_start
                    match_details_dict['FASTA_match_end'] = hsp.sbjct_end

                    #write the header to the header of the temporary csv file which is required for later processing
                    if hit_num == 0:
                        csv_header_for_ncbi_homologues_file = sorted(list(match_details_dict.keys()))
                        writer = csv.writer(homo_out_csv_file_handle, delimiter=',', quotechar='"', lineterminator='\n', quoting=csv.QUOTE_NONNUMERIC, doublequote=True)
                        writer.writerow(csv_header_for_ncbi_homologues_file)
                    #save the math_details_dict into the csv file
                    writer = csv.DictWriter(homo_out_csv_file_handle, fieldnames=csv_header_for_ncbi_homologues_file, extrasaction='ignore', delimiter=',', quotechar='"', lineterminator='\n',
                                            quoting=csv.QUOTE_MINIMAL, doublequote=True)
                    writer.writerow(match_details_dict)
                    hit_num += 1

    #if csv file is empty (no BLAST hits) -> delete temporary files, give warning and return
    if os.stat(BLAST_csv_file).st_size == 0:
        os.remove(BLAST_csv_file)
        logging.info("BLAST_parser:" + "\t" + protein_name + " hadn't any BLAST hits")
        return

    #save file into pickle format for the later korbinian pipeline usage
    df_homol = pd.read_csv(BLAST_csv_file, sep=",", quoting=csv.QUOTE_NONNUMERIC, index_col="hit_num")
    if "query_align_seq" not in df_homol.columns:
        # this is a serious error in the XML file. None of the hits had a protein node. The file should probably be downloaded.
        warning = 'The homologue XML file likely has a serious error, "query_align_seq" is not in dataframe. ' \
                  'XML should probably be re-downloaded.\n' \
                  'df_homol["hit_contains_SW_node"].value_counts()\n{}'.format(df_homol["hit_contains_SW_node"].value_counts())
        logging.warning(warning)
        # skip this protein
        return acc, False, warning
    # get length of seq. Previously this was a lambda function that needed more filtering
    df_homol['len_query_align_seq'] = df_homol['query_align_seq'].str.len()

    # conduct the text searching for disallowed words
    words_not_allowed_in_description = ast.literal_eval(s["words_not_allowed_in_description"])
    # collect disallowed words in hit protein description (patent, synthetic, etc)
    df_homol['list_disallowed_words_in_descr'] = df_homol['description'].dropna().apply(utils.find_disallowed_words, args=(words_not_allowed_in_description,))
    # create a boolean column to select hits that do not contain these words in the description
    df_homol['disallowed_words_not_in_descr'] = df_homol['list_disallowed_words_in_descr'] == '[]'
    # check if there are non-IUPAC amino acids in the sequence (frequently large gaps from NG sequencing data)
    df_homol['X_in_match_seq'] = df_homol['match_align_seq'].str.contains("X")

    #Initialize homol directory file system
    homol_proteinID_dir = os.path.join(s["data_dir"], "homol", "parsed", protein_name[:2])
    if not os.path.exists(homol_proteinID_dir):
        os.makedirs(homol_proteinID_dir)

    #outputs data.table df_homol (homolog alignments) into pickle in directory homol/parsed/
    with open(p['homol_df_orig_pickle'], "wb") as pick:
        pickle.dump(df_homol, pick, protocol=pickle.HIGHEST_PROTOCOL)

    #Create new zip and add the just created df_homol data.table pickle file
    with zipfile.ZipFile(p['homol_df_orig_zip'], mode="w", compression=zipfile.ZIP_DEFLATED) as zipout:
        zipout.write(p['homol_df_orig_pickle'], arcname=os.path.basename(p['homol_df_orig_pickle']))

    #Delete temporary uncompressed files
    os.remove(p['homol_df_orig_pickle'])
    os.remove(BLAST_csv_file)

    #output current id
    logging.info("BLAST_parser:" + "\t" + protein_name + " successfully parsed")
