import ast
import csv
import korbinian.utils as utils
import pandas as pd
import numpy as np
import re
import sys
# import debugging tools
from korbinian.utils import pr, pc, pn, aaa

def extract_omp_IDs_from_nr_fasta(ListXX_OMPdb_nr_fasta, ListXX_OMPdb_nr_acc, logging):
    """Takes the OMP non-redundant list of fasta sequences, and extracts the protein IDs (fasta names).

    Parameters
    ----------
    ListXX_OMPdb_nr_fasta : str
        Path to OMPdb non-redundant list of fasta sequences
    ListXX_OMPdb_nr_acc : str
        Path to output file, with list of OMPdb accessions.
    logging : logging.Logger
        Logger for printing to console and logfile.

    Saved Files
    -----------
    ListXX_OMPdb_nr_acc : txt
        Contains list of non-redundant omp IDs, each on a new line.
    """
    with open(ListXX_OMPdb_nr_acc, "w") as f_out:
        with open(ListXX_OMPdb_nr_fasta) as f:
            # extract the fasta header (ID)
            ID_list = [lines.strip()[1:] for lines in f.readlines() if ">" in lines]
            # write to text file
            for ID in ID_list:
                f_out.write("%s\n" % ID)
    logging.info("extract_omp_IDs_from_nr_fasta is completed")

def parse_OMPdb_all_selected_to_csv(ListXX_OMPdb_nr_acc, ListXX_OMPdb_redundant_flatfile, OMPdb_list_csv, logging, s):
    """ Extracts ID, seq and topology data from the full OMPdb flatfile, saves to csv.

    Note that instead of parsing line-by-line and saving to a csv, this method store every single value into a huge dictionary, which for 3
    proteins looks like this:

    BB_SiPe                                        [True, True, False]
    Coverage(%)                                  [99.82, 96.61, 92.19]
    Description      [Pilin outer membrane usher protein SafC, Oute...
    NCBI_TaxID                                     [59201, 470, 59201]
    Organism         [Salmonella enterica I, Acinetobacter baumanni...
    SP01_start                                          [1, 1, np.nan]
    Sequence         [MKFKQPALLLFIAGVVHCANAHTYTFDASMLGDAAKGVDMSLFNQ...
    Topology         [IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII...
    Topology_Reli                                [84.13, 94.00, 93.07]
    Uniprot                       [A0A0F7J6A4, A0A090B0L0, A0A0F7J6D5]
    seqlen                                             [836, 356, 721]

    The size of this dictionary with >7000 entries may cause memory problems on a regular PC.
    Currently this method is functional, though, and is not a priority to be fixed.

    Parameters
    ----------
    ListXX_OMPdb_nr_acc : str
        Path to OMPdb list of non-redundant IDs, textfile.
    ListXX_OMPdb_redundant_flatfile : str
        Path to full OMPdb flatfile, all proteins, unzipped.
    OMPdb_list_csv : str
        Path to output csv file.

    Saved Files
    -----------
    OMPdb_list_csv : csv
        csv file derived from dfKW
        contains a row for each protein
        contains indices for TM regions

    Notes
    -----
    This script parses the original text file, rather than XML. Version on server was not working due to " ' error.
    The script works, but goes through an entire flatfile, so is unbelievably slow. Use at your own risk.
    """

    # Creating dictionary keywords
    keywords = {"Uniprot": [], "Family": [], "Gene_Name": [], "Organism": [], "NCBI_TaxID": [], "Coverage(%)": [],
                "Sequence": [], "Topology_Reli": [], "Topology": [], "Description": [], "Pfam_ID": [], "BB_SiPe": [], "seqlen": [], "fragment": []}

    # check if signal peptides should be added to the list_of_TMDs and analysed
    # signal peptides will still be detected, via "True" in BB_SiPe. This is useful for excluding potential TM01 mis-labelled as SP.
    analyse_SiPe = False
    if 'SiPe' in s['regions']:
        analyse_SiPe = True
        keywords.update({"SP01_start": [], "SP01_end": [], "SP01_seq": [], "SiPe_source": []})

    logging.info('analyse_SiPe: {}'.format(analyse_SiPe))

    # Start settings which are changed during the for loop
    take_next_seq = False
    take_next_topo = False
    take_ID = False

    # Empty lists which are filled during the for loop
    Raw_Sequences = []
    Raw_Topos = []
    ID_list = []

    # Extracts IDs out of file
    with open(ListXX_OMPdb_nr_acc) as source:
        for line in source:
            line = line.strip()
            ID_list.append(line)

    # Checking ListXX_OMPdb_redundant_flatfile(complete OMPdb in very unfriendly formatting)for IDs(which are stored in list of Potential IDs) and extracting information
    with open(ListXX_OMPdb_redundant_flatfile) as data_file:
        counter = 0
        db_cross_ref = {}
        save_db_cross_ref = False
        for line in data_file:
            line_list = line.strip().split(" ")
            # Further settings which are changed every loop
            sequence_header = False
            topo_header = False
            # If-conditions make sure, that the ID is in the list of Potential IDs and check for keywords in each line
            if "UNIPROT" in line_list and line_list[-1] in ID_list:
                keywords["Uniprot"].append(line_list[-1])
                take_ID = True
                counter += 1
                if counter % 100 == 0:
                    sys.stdout.write(". ")
                    sys.stdout.flush()
            if "FAMILY" in line_list and take_ID == True:
                keywords["Family"].append(" ".join(line_list[9:]))
            if "DESCRIPTION" in line_list and take_ID == True:
                keywords["Description"].append(" ".join(line_list[4:]))
            if "GENE_NAME" in line_list and take_ID == True:
                keywords["Gene_Name"].append(" ".join(line_list[6:]))
            if "ORGANISM" in line_list and take_ID == True:
                keywords["Organism"].append(" ".join(line_list[7:]))
            if "NCBI_TAXID" in line_list and take_ID == True:
                keywords["NCBI_TaxID"].append(line_list[-1])
            if "DB_REF" in line_list and take_ID == True:
                # add database cross references to special dict
                db_cross_ref.update({line_list[9][:-1]: line_list[10].split('|')})
            if "SIGNAL_PEPTIDE" in line_list and take_ID == True:
                if ' '.join(line_list[1:]) == 'No information available':
                    keywords["BB_SiPe"].append(False)
                    if analyse_SiPe == True:
                        keywords["SP01_start"].append(np.nan)
                        keywords["SP01_end"].append(np.nan)
                        keywords["SP01_seq"].append(np.nan)
                        keywords["SiPe_source"].append('No information available')
                else:
                    # assume there is a signal peptide that starts at 1 (not very optimum code!!!)
                    keywords["BB_SiPe"].append(True)
                    if analyse_SiPe == True:
                        keywords["SP01_start"].append(line_list[1][0])
                        keywords["SP01_end"].append(line_list[1][2:-1])
                        keywords["SP01_seq"].append(line_list[2][:-1])
                        keywords["SiPe_source"].append(' '.join(line_list[-2:]))

            if "COVERAGE(%)" in line_list and take_ID == True:
                keywords["Coverage(%)"].append(line_list[-1])
            if "SEQUENCE" in line_list and take_ID == True:
                # after the "SEQUENCE" statement in a line_list, all db cross references are collected and can be saved
                save_db_cross_ref = True
                keywords["seqlen"].append(line_list[7])
                take_next_seq = True
                sequence_header = True
                # some of the OMPdb entries are labeled as fragments. These should be removed.
                if "Fragment:" in line:
                    searchstring = ".*Fragment\:([NC])"
                    match = re.match(searchstring, line)
                    if match:
                        N_or_C = match.group(1)
                    else:
                        N_or_C = "undefined"
                    keywords["fragment"].append("{}-term".format(N_or_C))
                else:
                    keywords["fragment"].append("undefined")
            # add db cross references from previous protein to keywords dict
            if save_db_cross_ref == True:
                if "Pfam" in db_cross_ref.keys():
                    keywords["Pfam_ID"].append(db_cross_ref["Pfam"])
                else:
                    keywords["Pfam_ID"].append(np.nan)
                # reset db_cross_ref for next cycle
                save_db_cross_ref = False
                db_cross_ref = {}
            if "TOPOLOGY" in line_list and take_ID == True:
                Raw_Sequences.extend(";")
                keywords["Topology_Reli"].append(line_list[-1].strip('"').strip("%"))
                take_next_seq = False
                take_next_topo = True
                topo_header = True
            if take_next_seq == True and sequence_header != True and take_ID == True:
                Raw_Sequences.extend(line_list)
            if "//" in line_list and take_ID == True:
                Raw_Topos.extend(";")
                topo_header = False
                take_next_topo = False
                take_ID = False
            if take_next_topo == True and topo_header != True:
                Raw_Topos.extend(line_list)
        Sequences = "".join(Raw_Sequences).split(";")
        Sequences.remove("")
        keywords["Sequence"] = Sequences
        Topos = "".join(Raw_Topos).split(";")
        Topos.remove("")
        keywords["Topology"] = Topos

    # Creating Dataframe and saving it as csv
    dfKW = pd.DataFrame(keywords)
    # set the uniprot_acc as the index
    dfKW.set_index("Uniprot", inplace=True, drop=False)
    dfKW.index.name = "acc"

    # DEPRECATED. OMPdb seems to label everything as a fragment?
    # n_prot_before_dropping_fragments = dfKW.shape[0]
    # dfKW = dfKW.loc[dfKW.fragment == "no fragment annotation"]
    # n_prot_after_dropping_fragments = dfKW.shape[0]
    # n_prot_fragments_dropped = n_prot_before_dropping_fragments - n_prot_after_dropping_fragments

    n_fragments = dfKW.loc[dfKW.fragment != "N"].shape[0]
    logging.info("{}/{} proteins labeled as 'Fragment:N' in flatfile.".format(n_fragments, dfKW.shape[0]))

    utils.make_sure_path_exists(OMPdb_list_csv, isfile=True)
    dfKW.to_csv(OMPdb_list_csv)
    logging.info("parse_OMPdb_all_selected_to_csv is completed.\n"
                 "Final number of proteins = {}".format(dfKW.shape[0]))


def get_omp_TM_indices_and_slice_from_summary_table(OMPdb_list_csv, list_parsed_csv, OMPdb_topology_reliability_cutoff, logging, s):
    """ Take a csv parsed from OMPdb, get the TM indices and slice the TMDs for each protein

    Parameters:
    -----------
    OMPdb_list_csv : str
        Path to input csv with OMP sequences and membrane annotation
    list_summary_csv : str
        Path to output csv with the sliced TM sequences
    logging : logging.Logger
        Logger for printing to console and logfile.
    """
    logging.info('~~~~~~~~~starting get_omp_TM_indices_and_slice_from_summary_table~~~~~~~~~')
    df_KW = pd.read_csv(OMPdb_list_csv, sep=",", quoting=csv.QUOTE_NONNUMERIC, index_col=0)

    # check if signal peptides should be extracted, modify keywords dict
    analyse_SiPe = False
    if 'SiPe' in s['regions']:
        analyse_SiPe = True

    # get sequence length
    df_KW["seqlen"] = df_KW["Sequence"].str.len()

    # Creating new column M_indices, which contains the indices of Ms
    df_KW["M_indices"] = df_KW.Topology.apply(getting_membrane_indices)

    # Converting empty entries to NaN
    df_KW["M_indices"] = df_KW.M_indices.apply(lambda x: np.nan if x == [] else x)

    num_proteins_BEFORE_dropping_those_without_mem_indices = df_KW.shape[0]

    # Extracting entries to a new Dataframe
    df_KW = df_KW[df_KW.M_indices.notnull()]

    num_proteins_AFTER_dropping_those_without_mem_indices = df_KW.shape[0]

    # Filter, cutting of Coverages under 85% & Creating new Index
    df_KW = df_KW.loc[df_KW["Coverage(%)"] >= 85]

    # df_KW.index = range(1,len(df_KW["Uniprot"])+1)
    num_proteins_AFTER_dropping_those_with_coverage_below_85 = df_KW.shape[0]

    # Creating new list (nested list)
    nested_list_of_membrane_borders = []

    # Filling nest with lists of start and end-points
    for n in df_KW.M_indices:
        m_borders = []
        m_borders.append(n[0])
        m_borders = check_for_border(n, m_borders)
        m_borders.append(n[-1])
        nested_list_of_membrane_borders.append(m_borders)

    array_membrane_borders = np.array(nested_list_of_membrane_borders)
    array_membrane_borders_corrected = []
    for subarray in array_membrane_borders:
        # logging.info(subarray[::2] = subarray[::2]*10)
        subarray = np.array(subarray)
        subarray[1::2] = subarray[1::2] + 1
        array_membrane_borders_corrected.append(list(subarray))

    nested_list_of_membrane_borders_python_indexstyle = array_membrane_borders_corrected

    # Creating new column, which contains start and end-points
    df_KW["Membrane_Borders"] = nested_list_of_membrane_borders_python_indexstyle

    # Creating new column, which contains the Amoung of TMDS
    df_KW["number_of_TMDs"] = df_KW.Membrane_Borders.apply(lambda x: len(x) / 2)

    # Filter, filters out, if less than 8 or more than 24 TMDs
    # REMOVED. FILTERING BY NUMBER OF TMDS IS NOW DONE LATER, in PROT_LIST
    #df_KW["number_of_TMDs"] = df_KW["number_of_TMDs"].apply(lambda x: int(x) if 5 <= x <= 36 else np.nan)

    # Creating new dataframe without nan
    df_KW = df_KW[df_KW["number_of_TMDs"].notnull()]

    num_proteins_AFTER_dropping_those_without_TMDs = df_KW.shape[0]

    df_KW = df_KW[df_KW["Topology_Reli"] > OMPdb_topology_reliability_cutoff]

    num_proteins_AFTER_dropping_those_with_topology_reliability_below_cutoff = df_KW.shape[0]

    df_KW["TM_indices"] = df_KW["Membrane_Borders"].apply(lambda x: tuple(zip(x[::2], x[1::2])))

    # create a list of [TM01, TM02, TM03, etc.
    long_list_of_TMDs = []
    for i in range(1, 50):
        long_list_of_TMDs.append("TM{:02d}".format(i))

    # for the .set_value function, set dtype as object
    df_KW["list_of_TMDs"] = ""
    df_KW["list_of_TMDs"].astype(object)

    sys.stdout.write('slicing TMD and nonTMD sequences:\n')

    for row_nr, row in enumerate(df_KW.index):
        # get nested tuple of TMDs
        nested_tup_TMs = df_KW.loc[row, "TM_indices"]
        # slice long list of TMD names to get an appropriate list for that protein [TM01, TM02, TM03, etc.
        len_nested_tup_TMs = len(nested_tup_TMs)
        list_of_TMDs = long_list_of_TMDs[:len_nested_tup_TMs]
        # add that list to the dataframe (could also be added as a stringlist, but that's irritating somehow)
        df_KW.set_value(row, "list_of_TMDs", list_of_TMDs)
        # set seq for slicing
        full_seq = df_KW.loc[row, "Sequence"]
        # topology = df_KW.loc[row, "Topology"]
        # iterate through all the TMDs of that protein, slicing out the sequences
        for i in range(len(list_of_TMDs)):
            TMD = list_of_TMDs[i]
            tup = nested_tup_TMs[i]
            df_KW.loc[row, TMD + "_start"] = tup[0]
            df_KW.loc[row, TMD + "_end"] = tup[1]
            df_KW.loc[row, TMD + "_seq"] = utils.slice_with_listlike(full_seq, tup)
            # df_KW.loc[row, TMD + "_top"] = utils.slice_with_listlike(topology, tup)
        if row_nr % 50 == 0:
            sys.stdout.write(". ")
            sys.stdout.flush()
            if row_nr % 500 == 0:
                sys.stdout.write("\n")
                sys.stdout.flush()


        ''' ~~   SLICE nonTMD sequence  ~~ '''
        #list_of_TMDs = df_KW.loc[row, 'list_of_TMDs'].copy()
        if 'SP01' in list_of_TMDs:
            list_of_TMDs.remove('SP01')
        # sequence from N-term. to first TMD
        nonTMD_first = df_KW.loc[row, 'Sequence'][0: (df_KW.loc[row, 'TM01_start'] - 1).astype('int64')]
        sequence = nonTMD_first
        # only for multipass proteins, generate sequences between TMDs
        if len(list_of_TMDs) == 0:
            # no TMDs are annotated, skip to next protein
            continue
        elif len(list_of_TMDs) > 1:
            for TM_Nr in range(len(list_of_TMDs) - 1):
                # the TMD is the equivalent item in the list
                TMD = list_of_TMDs[TM_Nr]
                # the next TMD, which contains the end index, is the next item in the list
                next_TMD = list_of_TMDs[TM_Nr + 1]
                between_TM_and_TMplus1 = df_KW.loc[row, 'Sequence'][df_KW.loc[row, '%s_end' % TMD].astype('int64'): df_KW.loc[row, '%s_start' % next_TMD].astype('int64') - 1]
                sequence += between_TM_and_TMplus1
        last_TMD = list_of_TMDs[-1]
        # sequence from last TMD to C-term.
        nonTMD_last = df_KW.loc[row, 'Sequence'][df_KW.loc[row, '%s_end' % last_TMD].astype('int64'):df_KW.loc[row, 'seqlen']]
        sequence += nonTMD_last
        df_KW.loc[row, 'nonTMD_seq'] = sequence
        df_KW.loc[row, 'len_nonTMD'] = len(sequence)

        if analyse_SiPe == True:
            if pd.notnull(df_KW.loc[row, 'SP01_start']):
                list_of_TMDs.append('SP01')
                df_KW.set_value(row, "list_of_TMDs", list_of_TMDs)

    ########################################################################################
    #                                                                                      #
    #                slicing out TMD_seq_plus_surr shifted to prot_list.py                 #
    #                                                                                      #
    ########################################################################################
    # max_num_TMDs = df_KW["number_of_TMDs"].max()
    #
    # # n_aa_before_tmd = s["n_aa_before_tmd"]
    # # n_aa_after_tmd = s["n_aa_after_tmd"]
    # n_aa_before_tmd = 10
    # n_aa_after_tmd = 10
    #
    # # currently the loop is run for each TMD, based on the sequence with the most TMDs
    # for i in range(1, int(max_num_TMDs) + 1):
    #     TMD = 'TM%02d' % i
    #     # get the indices for TMD plus surrounding sequence
    #     df_KW = korbinian.prot_list.prot_list.get_indices_TMD_plus_surr_for_summary_file(df_KW, TMD, n_aa_before_tmd, n_aa_after_tmd)
    #     # slice out the TMD_seq_plus_surr for each TMD
    #     df_KW['%s_seq_plus_surr' % TMD] = df_KW[df_KW['%s_start' % TMD].notnull()].apply(utils.slice_uniprot_TMD_plus_surr_seq, args=(TMD,), axis=1)

    # rename columns to match protein lists from uniprot (Note that Family is currently translated as prot_descr)
    dict_ = {"Sequence": "full_seq", "Organism": "organism", "Uniprot": "uniprot_acc", "Gene_Name": "gene_name",
             "Topology_Reli": "topology_reliability", "Family" : "prot_descr"}
    df_KW["betabarrel"] = True
    df_KW["multipass"] = True
    df_KW["singlepass"] = False
    # since all beta-barrel proteins have the N-terminus in the periplasm, "N-term is Extracellular" is False
    # you could make 100% sure of this by checking that the first letter of "Topology" is "I", but it is not really necessary
    df_KW["n_term_ec"] = False
    df_KW.rename(columns=dict_, inplace=True)
    df_KW["acc"] = df_KW["uniprot_acc"]
    df_KW["protein_name"] = df_KW["uniprot_acc"]
    num_proteins_AFTER_get_omp_TM_indices_and_slice_from_summary_table = df_KW.shape[0]

    # save to csv (presumably in summaries folder as a list number, so it is accessible by the rest of the scripts)
    utils.make_sure_path_exists(list_parsed_csv, isfile=True)
    df_KW.to_csv(list_parsed_csv, sep=",", quoting=csv.QUOTE_NONNUMERIC)

    logging.info("\nnum_proteins_BEFORE_dropping_those_without_mem_indices : {}".format(num_proteins_BEFORE_dropping_those_without_mem_indices))
    logging.info("num_proteins_AFTER_dropping_those_without_mem_indices : {}".format(num_proteins_AFTER_dropping_those_without_mem_indices))
    logging.info("num_proteins_AFTER_dropping_those_with_coverage_below_85 : {}".format(num_proteins_AFTER_dropping_those_with_coverage_below_85))
    logging.info("num_proteins_AFTER_dropping_those_without_TMDs : {}".format(num_proteins_AFTER_dropping_those_without_TMDs))
    logging.info("num_proteins_AFTER_dropping_those_with_topology_reliability_below_cutoff : {}".format(num_proteins_AFTER_dropping_those_with_topology_reliability_below_cutoff))
    logging.info("num_proteins_AFTER_get_omp_TM_indices_and_slice_from_summary_table : {}".format(num_proteins_AFTER_get_omp_TM_indices_and_slice_from_summary_table))
    logging.info('~~~~~~~~~finished get_omp_TM_indices_and_slice_from_summary_table~~~~~~~~~')

# Function which returns list of all M-indices
def getting_membrane_indices(Topo_data):
    m_list = [i for i, topology in enumerate(Topo_data) if topology == "M"]  # find(Topo_data)
    return m_list

# Function which filters out start and end-points
def check_for_border(m_index_list, m_borders):
    """ Checks for the borders of membrane regions, from M indices.

    Parameters
    ----------
    m_indices : list
        List of membrane region indices, extracted from the flatfile, (e.g. IIIMMMMMMMMMOOOO for Inside, Membrane and Outside),
        giving the list of M indices (e.g. [28, 29, 30, 31, 32, 33, 34, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 72, 73, 74, 75, 76]
    m_borders : list
        List of borders, where the membrane regions start and stop.

    Returns
    -------
    m_borders : list
        Updated list of borders
    """
    # iterate through the m_index_list (e.g. [13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33] for a singlepass protein)
    for n in range(0, len(m_index_list) - 1):
        # if the two indices next to each other are not consecutive numbers
        if m_index_list[n] + 1 != m_index_list[n + 1]:
            # add the first one to the list
            m_borders.append(m_index_list[n])
            # add the second one to the list
            m_borders.append(m_index_list[n + 1])
    return m_borders