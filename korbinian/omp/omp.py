import korbinian.mtutils as utils

############################################################
#
#        from Ipython notebook, IDs_of_OMPdb30
#
############################################################


save_file =open("Potential_IDs.txt", "w+")

with open (r"C:\Users\Rima\Desktop\Bachelorthesis\Databases\Download\OMPdb.30\OMPdb.30") as text_file:

    Potential_IDs = [lines.strip()[1:] for lines in text_file.readlines() if ">" in lines]
    for ID in Potential_IDs:
        save_file.write("%s\n" % ID)


######################################################################################
#
#        from Ipython notebook, OMPdb entries of Potential IDs
#
######################################################################################


import pandas as pd

text_file = r"C:\Users\Rima\Desktop\Bachelorthesis\Databases\Download\download_OMPdb.txt"
id_file = r"C:\Users\Rima\Desktop\Bachelorthesis\Lists\Potential_IDs_of_OMPdb30.txt"

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
IDs = []

# Extracts IDs out of file
with open(id_file) as source:
    for line in source:
        line = line.strip()
        IDs.append(line)
# Checking text_file(complete OMPdb in very unfriendly formatting)for IDs(which are stored in list of Potential IDs) and extracting information
with open(text_file) as data_file:
    for line in data_file:
        line = line.strip().split(" ")
        # Further settings which are changed every loop
        sequence_header = False
        topo_header = False
        # If-conditions make sure, that the ID is in the list of Potential IDs and check for keywords in each line
        if "UNIPROT" in line and line[-1] in IDs:
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
        if '"SEQUENCE' in line and take_ID == True:
            keywords["len_Sequence"].append(line[7])
            take_next_seq = True
            sequence_header = True
        if '"TOPOLOGY' in line and take_ID == True:
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


# Creating Dataframe and saving it as csv and Excelfile

dfKW = pd.DataFrame(keywords)
dfKW.to_csv("OMPdb_Selected_by_IDs.csv",index=False)

writer = pd.ExcelWriter("OMPdb_Selected_by_IDs.xlsx")
dfKW.to_excel(writer)

