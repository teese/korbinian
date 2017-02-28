import pandas as pd
import numpy as np
import csv
import korbinian
import sys
import korbinian.utils as utils

def read_scampi_data(pathdict, s, logging, df):
    #pathdict['SCAMPI'] = '/Volumes/Musik/Databases/summaries/01/List01_SCAMPI/query.top.txt'
    logging.info('~~~~~~~~~~~~                    reading scampi data                    ~~~~~~~~~~~~')

    # read text file from disk, Check SCAMPI output for accession numbers in acc list and read the toplology
    topo_list = []

    with open(pathdict['SCAMPI']) as data_file:
        for line in data_file:
            line = line.strip()
            if line[0] == '>':
                topo_list.append(line[1:])
            else:
                topo_list.append(line)

    # initialise pandas dataframe with accession as index - here every 2nd item in list is used
    dft = pd.DataFrame(index=topo_list[0::2])
    # add topology from topo_list to dataframe
    dft['SCAMPI_topo'] = topo_list[1::2]

    # get indices for TMD residues
    dft['SCAMPI_M_indices'] = dft.SCAMPI_topo.apply(korbinian.prot_list.parse_OMPdb.getting_membrane_indices)
    # remove proteins without transmembrane regions
    dft = dft[dft.astype(str)['SCAMPI_M_indices'] != '[]']

    # Creating new list (nested list)
    nested_list_of_membrane_borders = []

    # Filling nest with lists of start and end-points
    for n in dft.SCAMPI_M_indices:
        m_borders = []
        m_borders.append(n[0])
        m_borders = korbinian.prot_list.parse_OMPdb.check_for_border(n, m_borders)
        m_borders.append(n[-1])
        nested_list_of_membrane_borders.append(m_borders)

    array_membrane_borders = np.array(nested_list_of_membrane_borders)
    array_membrane_borders_corrected = []
    for subarray in array_membrane_borders:
        subarray = np.array(subarray)
        subarray[1::2] = subarray[1::2] + 1
        array_membrane_borders_corrected.append(list(subarray))

    nested_list_of_membrane_borders_python_indexstyle = array_membrane_borders_corrected

    # Creating new column, which contains start and end-points
    dft["Membrane_Borders"] = nested_list_of_membrane_borders_python_indexstyle

    # Creating new column, which contains the number of TMDS
    dft["number_of_TMDs"] = dft.Membrane_Borders.apply(lambda x: len(x) / 2)
    dft["TM_indices"] = dft["Membrane_Borders"].apply(lambda x: tuple(zip(x[::2], x[1::2])))
    # create a list of [TM01, TM02, TM03, etc.
    long_list_of_TMDs = []
    for i in range(1, 50):
        long_list_of_TMDs.append("TM{:02d}".format(i))
    ## for the .set_value function, set dtype as object
    dft["list_of_TMDs"] = ""
    #dft["list_of_TMDs"].astype(str)

    sys.stdout.write('slicing TMD and nonTMD sequences:\n')

    for row_nr, row in enumerate(dft.index):
        # get nested tuple of TMDs
        nested_tup_TMs = dft.loc[row, "TM_indices"]
        # slice long list of TMD names to get an appropriate list for that protein [TM01, TM02, TM03, etc.
        len_nested_tup_TMs = len(nested_tup_TMs)
        list_of_TMDs = long_list_of_TMDs[:len_nested_tup_TMs]
        # add that list to the dataframe (could also be added as a stringlist, but that's irritating somehow)
        dft.set_value(row, "list_of_TMDs", list_of_TMDs)
        # set seq for slicing
        full_seq = df.loc[row, "full_seq"]
        # topology = dft.loc[row, "Topology"]
        # iterate through all the TMDs of that protein, slicing out the sequences
        for i in range(len(list_of_TMDs)):
            TMD = list_of_TMDs[i]
            tup = nested_tup_TMs[i]
            dft.loc[row, "%s_start" % TMD] = tup[0]
            dft.loc[row, "%s_end" % TMD] = tup[1]
            dft.loc[row, "%s_seq" % TMD] = utils.slice_with_listlike(full_seq, tup)
            dft.loc[row, "%s_seqlen" % TMD] = len(dft.loc[row, "%s_seq" % TMD])
            # dft.loc[row, TMD + "_top"] = utils.slice_with_listlike(topology, tup)
        if row_nr % 50 == 0:
            sys.stdout.write(". ")
            sys.stdout.flush()
            if row_nr % 500 == 0:
                sys.stdout.write("\n")
                sys.stdout.flush()

    # # define columns to replace
    # columns_to_replace = ['nonTMD_seq', 'number_of_TMDs', 'list_of_TMDs']
    # for n in range(1, df.number_of_TMDs.max().astype('int') + 1):
    #     TMD = 'TM{:02d}'.format(n)
    #     drop = '{a}_description, {a}_end, {a}_seq_variants, {a}_start, {a}_end_plus_surr, {a}_start_plus_surr'.format(a=TMD)
    #     for element in drop.split(', '):
    #         columns_to_replace.append(element)
    # # define columns to keep in df
    # columns_to_keep = [x for x in list(df.columns) if x not in columns_to_replace]
    # print(columns_to_keep)
    # # select df for columns to keep
    # df = df[columns_to_keep]
    # merge SCAMPI output into the uniprot-cleaned dataframe
    df = pd.merge(df, dft, left_index=True, right_index=True, suffixes=('_uniprot', ''))

    return df

