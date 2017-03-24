import pandas as pd
import numpy as np
import korbinian
import sys
import korbinian.utils as utils
import os
import csv
import ast

def parse_TMSEG_results(analyse_sp, pathdict, s, logging):
    logging.info("~~~~~~~~~~~~                        starting parse_TMSEG_results                    ~~~~~~~~~~~~")

    list_number = s['list_number']

    # define the uniprot directory with selected records
    uniprot_dir_sel = os.path.join(s["data_dir"], 'uniprot', 'selected')
    selected_uniprot_records_flatfile = os.path.join(uniprot_dir_sel, 'List%02d_selected_uniprot_records_flatfile.txt' % list_number)
    n_aa_before_tmd = s["n_aa_before_tmd"]
    n_aa_after_tmd = s["n_aa_after_tmd"]
    list_parsed_csv = pathdict["list_parsed_csv"]
    analyse_signal_peptides = s['SiPe']
    output = korbinian.prot_list.uniprot_parse.create_csv_from_uniprot_flatfile(selected_uniprot_records_flatfile, n_aa_before_tmd, n_aa_after_tmd, analyse_signal_peptides, logging, list_parsed_csv, slice=False)
    logging.info(output)

    TMSEG_results_filepath = pathdict['TMSEG_top']
    TMSEG_nonTM_outpath = pathdict['TMSEG_nonTM']

    list_parsed = pd.read_csv(pathdict["list_parsed_csv"], sep=",", quoting=csv.QUOTE_NONNUMERIC, index_col=0, low_memory=False)

    columns_to_keep = ['organism_domain', 'create_csv_from_uniprot_flatfile', 'uniprot_acc', 'uniprot_all_accessions', 'uniprot_entry_name', 'uniprot_features',
                       'uniprot_orgclass', 'uniprot_SiPe', 'singlepass', 'typeI', 'typeII', 'uniprot_KW', 'organism', 'prot_descr', 'membrane',
                       'multipass', 'gene_name', 'comments_subcellular_location_uniprot']
    list_indices = list(list_parsed.index)
    list_parsed = list_parsed[columns_to_keep]

    # read data from file
    input_data = []
    with open(TMSEG_results_filepath) as data_file:
        for line in data_file:
            line = line.strip()
            if line[0] == '>':
                line = line[1:]
                line = line.split(' ')
                comment = ' '.join(line[1:])
                line = line[0].split('|')
                uniprot_acc = line[0]
                uniprot_entry_name = line[1]
                input_data.append(uniprot_acc)
                input_data.append(uniprot_entry_name)
                input_data.append(comment)
            else:
                input_data.append(line)


    # initialise pandas dataframe with uniprot accession as index
    df = pd.DataFrame(index=input_data[0::5])

    # add selected columns from input_data list
    df['uniprot_entry_name'] = input_data[1::5]
    df['prot_descr'] = input_data[2::5]
    df['full_seq'] = input_data[3::5]
    df['topology'] = input_data[4::5]

    keep = []
    for acc in df.index:
        if acc in list_indices:
            keep.append(acc)
    df = df.loc[keep,:]
    if df.shape[0] == 0:
        return sys.stdout.write('no remaining proteins in list!')

    # get list of uniprot accessions of proteins where no transmembrane region was predicted
    list_nonTMD = []
    for acc in df.index:
        if 'N' in df.loc[acc, 'topology']:
            list_nonTMD.append(acc)

    # write list of nonTM proteins to file
    # outpath = '/Volumes/Musik/Databases/TMSEG/humanU90_nonTMD.txt'
    file = open(TMSEG_nonTM_outpath, 'w')
    for line in list_nonTMD:
        file.write('{}\n'.format(line))
    file.close()

    # drop proteins that do not contain TM regions
    df = df.drop(list_nonTMD)

    # add seqlen and indices for all TMD and SiPe regions
    df["seqlen"] = df.full_seq.apply(lambda x: len(x))
    df['M_indices'] = df.topology.apply(getting_membrane_indices_from_helix_symbol)
    df['SiPe_indices'] = df.topology.apply(getting_SiPe_indices_from_symbol)



    # Creating new list (nested list)
    nested_list_of_membrane_borders = []

    # Filling nest with lists of start and end-points
    for n in df.M_indices:
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
    df["Membrane_Borders"] = nested_list_of_membrane_borders_python_indexstyle

    # Creating new column, which contains the number of TMDS
    #df["number_of_TMDs"] = df.Membrane_Borders.apply(lambda x: len(x) / 2)
    df["TM_indices"] = df["Membrane_Borders"].apply(lambda x: tuple(zip(x[::2], x[1::2])))
    # create a list of [TM01, TM02, TM03, etc.
    long_list_of_TMDs = []
    for i in range(1, 50):
        long_list_of_TMDs.append("TM{:02d}".format(i))

    ## for the .set_value function, set dtype as object
    df["list_of_TMDs"] = ""
    # dft["list_of_TMDs"].astype(object)


    sys.stdout.write('slicing TMD and nonTMD sequences:\n')

    for n, acc in enumerate(df.index):
        # get nested tuple of TMDs
        nested_tup_TMs = df.loc[acc, "TM_indices"]
        # slice long list of TMD names to get an appropriate list for that protein [TM01, TM02, TM03, etc.
        len_nested_tup_TMs = len(nested_tup_TMs)
        list_of_TMDs = long_list_of_TMDs[:len_nested_tup_TMs]
        # add that list to the dataframe (could also be added as a stringlist, but that's irritating somehow)
        df.set_value(acc, "list_of_TMDs", list_of_TMDs)
        # set seq for slicing
        full_seq = df.loc[acc, "full_seq"]
        # topology = dft.loc[acc, "Topology"]
        # iterate through all the TMDs of that protein, slicing out the sequences
        for i in range(len(list_of_TMDs)):
            TMD = list_of_TMDs[i]
            tup = nested_tup_TMs[i]
            df.loc[acc, "%s_start" % TMD] = tup[0]
            df.loc[acc, "%s_end" % TMD] = tup[1]
            df.loc[acc, "%s_seq" % TMD] = utils.slice_with_listlike(full_seq, tup)
            df.loc[acc, "%s_seqlen" % TMD] = len(df.loc[acc, "%s_seq" % TMD])
            # dft.loc[acc, TMD + "_top"] = utils.slice_with_listlike(topology, tup)
        # add signal peptides and their corresponding values to list_of_TMDs
        if analyse_sp == True:
            SiPe_indices = df.loc[acc, 'SiPe_indices']
            if SiPe_indices != []:
                df.loc[acc, 'SP01_start'] = SiPe_indices[0]
                df.loc[acc, 'SP01_end'] = SiPe_indices[-1]
                df.loc[acc, 'SP01_seq'] = full_seq[SiPe_indices[0]:SiPe_indices[-1]+1]
                list_of_TMDs.append('SP01')
                df.set_value(acc, "list_of_TMDs", list_of_TMDs)

        if n % 50 == 0 and n != 0:
            sys.stdout.write(". ")
            sys.stdout.flush()
            if n % 500 == 0:
                sys.stdout.write("\n")
                sys.stdout.flush()

        df.loc[acc, "number_of_TMDs"] = len(list_of_TMDs)

        ''' ~~   SLICE nonTMD sequence  ~~ '''
        list_of_TMDs = df.loc[acc, 'list_of_TMDs'].copy()
        if 'SP01' in list_of_TMDs:
            list_of_TMDs.remove('SP01')
        # sequence from N-term. to first TMD
        nonTMD_first = df.loc[acc, 'full_seq'][0: (df.loc[acc, 'TM01_start'] - 1).astype('int64')]
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
                between_TM_and_TMplus1 = df.loc[acc, 'full_seq'][df.loc[acc, '%s_end' % TMD].astype('int64'): df.loc[acc, '%s_start' % next_TMD].astype('int64') - 1]
                sequence += between_TM_and_TMplus1
        last_TMD = list_of_TMDs[-1]
        # sequence from last TMD to C-term.
        nonTMD_last = df.loc[acc, 'full_seq'][df.loc[acc, '%s_end' % last_TMD].astype('int64'):df.loc[acc, 'seqlen'].astype('int64')]
        sequence += nonTMD_last
        df.loc[acc, 'nonTMD_seq'] = sequence
        df.loc[acc, 'len_nonTMD'] = len(sequence)




    cols_to_drop = ['topology', 'M_indices', 'SiPe_indices', 'Membrane_Borders', 'TM_indices']
    df = pd.merge(df, list_parsed, left_index=True, right_index=True, suffixes=('', '_list_parsed'))
    df.drop(cols_to_drop, axis=1, inplace=True)
    df['parse_TMSEG'] = True
    df.to_csv(pathdict["list_parsed_csv"], sep=",", quoting=csv.QUOTE_NONNUMERIC)

    logging.info("\n~~~~~~~~~~~~                       parse_TMSEG_results is finished                  ~~~~~~~~~~~~")





def getting_membrane_indices_from_helix_symbol(Topo_data):
    m_list = [i for i, topology in enumerate(Topo_data) if topology == "H"]  # find(Topo_data)
    return m_list

def getting_SiPe_indices_from_symbol(Topo_data):
    m_list = [i for i, topology in enumerate(Topo_data) if topology == "S"]  # find(Topo_data)
    return m_list