import pandas as pd
import numpy as np
import csv
import korbinian
import sys
import korbinian.utils as utils
import os
# import debugging tools
from korbinian.utils import pr, pc, pn, aaa

def read_scampi_data(pathdict, s, logging, df):
    """ Reads the output of the topology predictor SCAMPI (http://scampi.bioinfo.se).

    Parameters
    ----------
    s : dict
        Dictionary of settings derived from settings excel file.
        Columns "Parameter" and "Value" are converted to key/value in the dictionary, respectively.
    pathdict : dict
        Dictionary of the key paths and files associated with that List number.

    Returned Files and Figures
    -----------------------
    Input pandas dataframe is overwritten, TMD regions from Uniprot are removed and replaced with SCAMPI output.

    Required Files
    -----------------------
    query.top.txt : holds all topology information of SCAMPI output.

    """
    #pathdict['SCAMPI_top'] = '/Volumes/Musik/Databases/summaries/01/List01_SCAMPI/query.top.txt'
    logging.info('\n~~~~~~~~~~~~                           using scampi data                            ~~~~~~~~~~~~')
    # define columns to replace
    columns_to_replace = ['nonTMD_seq', 'number_of_TMDs', 'list_of_TMDs']
    for n in range(1, df.number_of_TMDs.max().astype('int') + 1):
        TMD = 'TM{:02d}'.format(n)
        drop = '{a}_seq_plus_surr, {a}_seq, {a}_description, {a}_end, {a}_seq_variants, {a}_start, {a}_end_plus_surr, {a}_start_plus_surr'.format(a=TMD)
        for element in drop.split(', '):
            columns_to_replace.append(element)
    # define columns to keep in df
    columns_to_keep = [x for x in list(df.columns) if x not in columns_to_replace]
    # select df for columns to keep
    df = df[columns_to_keep]

    # read text file from disk, Check SCAMPI output for accession numbers in acc list and read the toplology
    topo_list = []

    with open(pathdict['SCAMPI_top']) as data_file:
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

    # merge dataframe containing SCAMPI data with previous, uniprot-cleaned dataframe
    df = pd.merge(df, dft, left_index=True, right_index=True, suffixes=('_uniprot', ''))

    sys.stdout.write('slicing TMD and nonTMD sequences:\n')

    for row_nr, row in enumerate(dft.index):
        # get nested tuple of TMDs
        nested_tup_TMs = dft.loc[row, "TM_indices"]
        # slice long list of TMD names to get an appropriate list for that protein [TM01, TM02, TM03, etc.
        len_nested_tup_TMs = len(nested_tup_TMs)
        list_of_TMDs = long_list_of_TMDs[:len_nested_tup_TMs]
        # add that list to the dataframe (could also be added as a stringlist, but that's irritating somehow)
        df.set_value(row, "list_of_TMDs", list_of_TMDs)
        # set seq for slicing
        full_seq = df.loc[row, "full_seq"]
        # topology = dft.loc[row, "Topology"]
        # iterate through all the TMDs of that protein, slicing out the sequences
        for i in range(len(list_of_TMDs)):
            TMD = list_of_TMDs[i]
            tup = nested_tup_TMs[i]
            df.loc[row, "%s_start" % TMD] = tup[0]
            df.loc[row, "%s_end" % TMD] = tup[1]
            df.loc[row, "%s_seq" % TMD] = utils.slice_with_listlike(full_seq, tup)
            #df.loc[row, "%s_seqlen" % TMD] = len(df.loc[row, "%s_seq" % TMD])
            # dft.loc[row, TMD + "_top"] = utils.slice_with_listlike(topology, tup)
        if row_nr % 50 == 0 and n != 0:
            sys.stdout.write(". ")
            sys.stdout.flush()
            if row_nr % 500 == 0:
                sys.stdout.write("\n")
                sys.stdout.flush()

        ''' ~~   SLICE nonTMD sequence  ~~ '''
        list_of_TMDs = df.loc[row, 'list_of_TMDs']
        if 'SP01' in list_of_TMDs:
            list_of_TMDs.remove('SP01')
        # sequence from N-term. to first TMD
        nonTMD_first = df.loc[row, 'full_seq'][0: (df.loc[row, 'TM01_start'] - 1).astype('int64')]
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
                between_TM_and_TMplus1 = df.loc[row, 'full_seq'][df.loc[row, '%s_end' % TMD].astype('int64'): df.loc[row, '%s_start' % next_TMD].astype('int64') - 1]
                sequence += between_TM_and_TMplus1
        last_TMD = list_of_TMDs[-1]
        # sequence from last TMD to C-term.
        nonTMD_last = df.loc[row, 'full_seq'][df.loc[row, '%s_end' % last_TMD].astype('int64'):df.loc[row, 'seqlen'].astype('int64')]
        sequence += nonTMD_last
        df.loc[row, 'nonTMD_seq'] = sequence
        df.loc[row, 'len_nonTMD'] = len(sequence)

    logging.info('\n~~~~~~~~~~~~                     scampi data replaced uniprot                       ~~~~~~~~~~~~\n')
    return df



def generate_scampi_input_files(pathdict, s, logging):

    logging.info('~~~~~~~~~~~~                 starting generate_scampi_input_files                   ~~~~~~~~~~~~')

    list_number = s["list_number"]
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

    logging.info('~~~~~~~~~~~~                generate_scampi_input_files is finished                 ~~~~~~~~~~~~')


def generate_SignalP_input_files(pathdict, s, logging):

    logging.info('~~~~~~~~~~~~                starting generate_SignalP_input_files                   ~~~~~~~~~~~~')

    list_number = s["list_number"]
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
        full_seq = df.loc[acc, 'full_seq']
        if len(full_seq) < 6000:
            file.write('>{}\n{}\n'.format(acc, full_seq))
        elif len(full_seq) > 6000:
            file.write('>{}\n{}\n'.format(acc, full_seq[:5999]))
        else:
            return 'MADNESS happened, please panic!'
    file.close()

    logging.info('~~~~~~~~~~~~               generate_SignalP_input_files is finished                 ~~~~~~~~~~~~')

def get_SignalP_SiPe_acc (SignalP_SiPe_path):
    acc_list = []
    with open(SignalP_SiPe_path) as source:
        for line in source:
            if line[0] != '#':
                line = line.strip()
                line = line.split('\t')
                acc_list.append(line[0])
    return acc_list

def get_PrediSi_SiPe_acc(PrediSi_SiPe_path):
    PrediSi_SiPe_list = []
    with open(PrediSi_SiPe_path) as source:
        for line in source:
            line = line.strip()
            PrediSi_SiPe_list.append(line)
    return PrediSi_SiPe_list