import csv
import korbinian.utils as utils
import korbinian.prot_list.prot_list
import pandas as pd
import numpy as np

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
    with open(ListXX_OMPdb_nr_acc, "w") as s:
        with open(ListXX_OMPdb_nr_fasta) as f:
            # extract the fasta header (ID)
            ID_list = [lines.strip()[1:] for lines in f.readlines() if ">" in lines]
            # write to text file
            for ID in ID_list:
                s.write("%s\n" % ID)
    logging.info("extract_omp_IDs_from_nr_fasta is completed")

def parse_OMPdb_all_selected_to_csv(ListXX_OMPdb_nr_acc, ListXX_OMPdb_redundant_flatfile, OMPdb_list_summary_csv, logging):
    """ Extracts ID, seq and topology data from the full OMPdb flatfile, saves to csv.

    Parameters
    ----------
    ListXX_OMPdb_nr_acc : str
        Path to OMPdb list of non-redundant IDs, textfile.
    ListXX_OMPdb_redundant_flatfile : str
        Path to full OMPdb flatfile, all proteins, unzipped.
    OMPdb_list_summary_csv : str
        Path to output csv file.

    Saved Files
    -----------
    OMPdb_list_summary_csv : csv
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
                "Sequence": [], "len_Sequence": [], "Topology_Reli": [], "Topology": []}

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
        for line in data_file:
            line = line.strip().split(" ")
            # Further settings which are changed every loop
            sequence_header = False
            topo_header = False
            # If-conditions make sure, that the ID is in the list of Potential IDs and check for keywords in each line
            if "UNIPROT" in line and line[-1] in ID_list:
                keywords["Uniprot"].append(line[-1])
                take_ID = True
                counter += 1
                if counter % 100 == 0:
                    print(". ", end="", flush=True)
            if "FAMILY" in line and take_ID == True:
                keywords["Family"].append(" ".join(line[9:]))
            if "GENE_NAME" in line and take_ID == True:
                keywords["Gene_Name"].append(" ".join(line[6:]))
            if "ORGANISM" in line and take_ID == True:
                keywords["Organism"].append(" ".join(line[7:]))
            if "NCBI_TAXID" in line and take_ID == True:
                keywords["NCBI_TaxID"].append(line[-1])
            if "COVERAGE(%)" in line and take_ID == True:
                keywords["Coverage(%)"].append(line[-1])
            if "SEQUENCE" in line and take_ID == True:
                keywords["len_Sequence"].append(line[7])
                take_next_seq = True
                sequence_header = True
            if "TOPOLOGY" in line and take_ID == True:
                Raw_Sequences.extend(";")
                keywords["Topology_Reli"].append(line[-1].strip('"').strip("%"))
                take_next_seq = False
                take_next_topo = True
                topo_header = True
            if take_next_seq == True and sequence_header != True and take_ID == True:
                Raw_Sequences.extend(line)
            if "//" in line and take_ID == True:
                Raw_Topos.extend(";")
                topo_header = False
                take_next_topo = False
                take_ID = False
            if take_next_topo == True and topo_header != True:
                Raw_Topos.extend(line)
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
    utils.make_sure_path_exists(OMPdb_list_summary_csv, isfile=True)
    dfKW.to_csv(OMPdb_list_summary_csv)
    logging.info("parse_OMPdb_all_selected_to_csv is completed. Dataframe shape = {}".format(dfKW.shape))


def get_omp_TM_indices_and_slice_from_summary_table(OMPdb_list_summary_csv, list_summary_csv, logging):
    """ Take a csv parsed from OMPdb, get the TM indices and slice the TMDs for each protein

    Parameters:
    -----------
    OMPdb_list_summary_csv : str
        Path to input csv with OMP sequences and membrane annotation
    list_summary_csv : str
        Path to output csv with the sliced TM sequences
    logging : logging.Logger
        Logger for printing to console and logfile.
    """
    df_KW = pd.read_csv(OMPdb_list_summary_csv, sep=",", quoting=csv.QUOTE_NONNUMERIC, index_col=0)

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
    df_KW["number_of_TMDs"] = df_KW["number_of_TMDs"].apply(lambda x: int(x) if 8 <= x <= 24 else np.nan)

    num_proteins_AFTER_dropping_those_without_TMs_between_8_and_24 = df_KW.shape[0]

    df_KW = df_KW[df_KW["Topology_Reli"] > 90]

    num_proteins_AFTER_dropping_those_with_topology_reliability_below_90 = df_KW.shape[0]

    # Creating new dataframe without nan
    df_KW = df_KW[df_KW["number_of_TMDs"].notnull()]

    df_KW["TM_indices"] = df_KW["Membrane_Borders"].apply(lambda x: tuple(zip(x[::2], x[1::2])))

    # create a list of [TM01, TM02, TM03, etc.
    long_list_of_TMDs = []
    for i in range(1, 50):
        long_list_of_TMDs.append("TM{:02d}".format(i))

    # for the .set_value function, set dtype as object
    df_KW["list_of_TMDs"] = ""
    df_KW["list_of_TMDs"].astype(object)

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
            print(".", end="", flush=True)
            if row_nr % 1000 == 0:
                print("", flush=True)

    max_num_TMDs = df_KW["number_of_TMDs"].max()

    # n_aa_before_tmd = set_["n_aa_before_tmd"]
    # n_aa_after_tmd = set_["n_aa_after_tmd"]
    n_aa_before_tmd = 10
    n_aa_after_tmd = 10

    # currently the loop is run for each TMD, based on the sequence with the most TMDs
    for i in range(1, int(max_num_TMDs) + 1):
        TMD = 'TM%02d' % i
        df_KW = korbinian.prot_list.prot_list.get_indices_TMD_plus_surr_for_summary_file(df_KW, TMD, n_aa_before_tmd, n_aa_after_tmd)

    # rename columns to match protein lists from uniprot (Note that Family is currently translated as prot_descr)
    dict_ = {"Sequence": "full_seq", "Organism": "organism", "Uniprot": "uniprot_acc", "Gene_Name": "gene_name",
             "Topology_Reli": "topology_reliability", "Family" : "prot_descr"}
    df_KW["betabarrel"] = True
    df_KW["multipass"] = True
    df_KW["singlepass"] = False
    df_KW.rename(columns=dict_, inplace=True)
    df_KW["acc"] = df_KW["uniprot_acc"]
    df_KW["protein_name"] = df_KW["uniprot_acc"]

    # save to csv (presumably in summaries folder as a list number, so it is accessible by the rest of the scripts)
    utils.make_sure_path_exists(list_summary_csv, isfile=True)
    df_KW.to_csv(list_summary_csv, sep=",", quoting=csv.QUOTE_NONNUMERIC)

    logging.info("\nnum_proteins_BEFORE_dropping_those_without_mem_indices : {}".format(num_proteins_BEFORE_dropping_those_without_mem_indices))
    logging.info("num_proteins_AFTER_dropping_those_without_mem_indices : {}".format(num_proteins_AFTER_dropping_those_without_mem_indices))
    logging.info("num_proteins_AFTER_dropping_those_with_coverage_below_85 : {}".format(num_proteins_AFTER_dropping_those_with_coverage_below_85))
    logging.info("num_proteins_AFTER_dropping_those_without_TMs_between_8_and_24 : {}".format(num_proteins_AFTER_dropping_those_without_TMs_between_8_and_24))
    logging.info("num_proteins_AFTER_dropping_those_with_topology_reliability_below_90 : {}".format(num_proteins_AFTER_dropping_those_with_topology_reliability_below_90))


# Function which returns list of all M-indices
def getting_membrane_indices(Topo_data):
    m_list = [i for i, topology in enumerate(Topo_data) if topology == "M"]  # find(Topo_data)
    return m_list

# Function which filters out start and end-points
def check_for_border(m_indices, m_borders):
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
    for n in range(0, len(m_indices) - 1):
        if m_indices[n] + 1 != m_indices[n + 1]:
            m_borders.append(m_indices[n])
            m_borders.append(m_indices[n + 1])
    return m_borders