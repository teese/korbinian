import pandas as pd

# location of unzipped fasta seqs from OMPdb, nonredundant(e.g. to 30% aa identity)
omp_nr_fasta = r"D:\Databases\OMPdb\OMPdb.30"
omp_ID_nr_txt = omp_nr_fasta + "_IDS.txt"
#omp_ID_nr_txt = r"D:\Databases\OMPdb\test_list_IDs.txt"
OMPdb_all_flatfile = r"D:\Databases\OMPdb\OMPdb.flat"
#OMPdb_all_flatfile = r"D:\Databases\OMPdb\OMPdbsmaller.txt"
OMPdb_summary_nr_csv = omp_nr_fasta + "_flatfiles_nr.csv"



























def extract_omp_IDs_from_nr_fasta(omp_nr_fasta):
    with open(omp_ID_nr_txt, "w") as s:
        with open(omp_nr_fasta) as f:
            # extract the fasta header (ID)
            ID_list = [lines.strip()[1:] for lines in f.readlines() if ">" in lines]
            # write to text file
            for ID in ID_list:
                s.write("%s\n" % ID)

def parse_OMPdb_all_selected_to_csv(omp_ID_nr_txt, OMPdb_all_flatfile, OMPdb_summary_nr_csv):
    """ Extracts ID, seq and topology data from the full OMPdb flatfile, saves to csv.

    Parameters
    ----------
    omp_ID_nr_txt : str
    Path to OMPdb list of non-redundant IDs, textfile.
    OMPdb_all_flatfile : str
    Path to full OMPdb flatfile, all proteins, unzipped.
    OMPdb_summary_nr_csv : str
    Path to output csv file.

    Notes:
    ------
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
    with open(omp_ID_nr_txt) as source:
        for line in source:
            line = line.strip()
            ID_list.append(line)

    # Checking OMPdb_all_flatfile(complete OMPdb in very unfriendly formatting)for IDs(which are stored in list of Potential IDs) and extracting information
    with open(OMPdb_all_flatfile) as data_file:
        for line in data_file:
            line = line.strip().split(" ")
            # Further settings which are changed every loop
            sequence_header = False
            topo_header = False
            # If-conditions make sure, that the ID is in the list of Potential IDs and check for keywords in each line
            if "UNIPROT" in line and line[-1] in ID_list:
                keywords["Uniprot"].append(line[-1])
                take_ID = True
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

    # Checking if keywords-lists are equally long
    print(len(keywords["Uniprot"])
          , len(keywords["Family"])
          , len(keywords["Gene_Name"])
          , len(keywords["Organism"])
          , len(keywords["NCBI_TaxID"])
          , len(keywords["Coverage(%)"])
          , len(keywords["Sequence"])
          , len(keywords["len_Sequence"])
          , len(keywords["Topology_Reli"])
          , len(keywords["Topology"]))

    # Creating Dataframe and saving it as csv
    dfKW = pd.DataFrame(keywords)
    dfKW.to_csv(OMPdb_summary_nr_csv,index=False)
