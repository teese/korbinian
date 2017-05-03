import ast
import csv
import korbinian
import korbinian.utils as utils
import numpy as np
import os
import pandas as pd
import sys
import time

def parse_TMSEG_results(analyse_sp, pathdict, s, logging):
    logging.info("~~~~~~~~~~~~                        starting parse_TMSEG_results                    ~~~~~~~~~~~~")

    list_number = s['list_number']

    # define the uniprot directory with selected records
    uniprot_dir_sel = os.path.join(s["data_dir"], 'uniprot', 'selected')
    selected_uniprot_records_flatfile = os.path.join(uniprot_dir_sel, 'List%02d_selected_uniprot_records_flatfile.txt' % list_number)
    n_aa_before_tmd = s["n_aa_before_tmd"]
    n_aa_after_tmd = s["n_aa_after_tmd"]
    list_parsed_csv = pathdict["list_parsed_csv"]
    analyse_signal_peptides = True if "SiPe" in s["regions"] else False
    output = korbinian.prot_list.uniprot_parse.create_protein_list(selected_uniprot_records_flatfile, n_aa_before_tmd, n_aa_after_tmd, analyse_signal_peptides, logging, list_parsed_csv, slice=False)
    logging.info(output)

    TMSEG_fastalike_path = pathdict['TMSEG_fastalike']
    TMSEG_top_txtoutput_path = pathdict['TMSEG_top_txtoutput']
    TMSEG_nonTM_outpath = pathdict['TMSEG_nonTM']

    df_parsed = pd.read_csv(pathdict["list_parsed_csv"], sep=",", quoting=csv.QUOTE_NONNUMERIC, index_col=0, low_memory=False)

    columns_to_keep = ['organism_domain', 'create_protein_list', 'uniprot_acc', 'uniprot_all_accessions', 'uniprot_entry_name', 'uniprot_features',
                       'uniprot_orgclass', 'uniprot_SiPe', 'singlepass', 'typeI', 'typeII', 'uniprot_KW', 'organism', 'prot_descr', 'membrane',
                       'multipass', 'gene_name', 'comments_subcellular_location_uniprot', 'uniprot_SiPe', 'full_seq']

    # for datasets without SP found, turn off analyse_sp
    if analyse_sp == True and 'SP01_start' in df_parsed.columns:
        columns_to_keep = columns_to_keep + ['SP01_start', 'SP01_end', 'SP01_seq']
    else:
        analyse_sp == False

    list_indices = list(df_parsed.index)

    df_parsed = df_parsed[columns_to_keep]

    if os.path.isfile(TMSEG_fastalike_path):
        # drop the full sequence, and get from TMSEG
        df_parsed.drop('full_seq', axis=1, inplace=True)
        # read data from file
        input_data = []
        with open(TMSEG_fastalike_path) as data_file:
            for line in data_file:
                line = line.strip()
                if line[0] == '>':
                    line = line[1:]
                    line = line.split(' ')
                    line = line[0].split('|')
                    uniprot_acc = line[0]
                    input_data.append(uniprot_acc)
                else:
                    input_data.append(line)


        # initialise pandas dataframe with uniprot accession as index
        df = pd.DataFrame(index=input_data[0::3])

        # add selected columns from input_data list
        #df['uniprot_entry_name'] = input_data[1::5]
        #df['prot_descr'] = input_data[2::5]
        df['full_seq'] = input_data[1::3]
        df['topo'] = input_data[2::3]

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
            if 'N' in df.loc[acc, 'topo']:
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
        df['M_indices'] = df.topo.apply(getting_membrane_indices_from_helix_symbol)
        df['SiPe_indices'] = df.topo.apply(getting_SiPe_indices_from_symbol)



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
        df["list_of_TMDs"].astype(object)

        sys.stdout.write('slicing TMD and nonTMD sequences:\n')

        for n, acc in enumerate(df.index):
            # get nested tuple of TMDs
            nested_tup_TMs = df.loc[acc, "TM_indices"]
            # slice long list of TMD names to get an appropriate list for that protein [TM01, TM02, TM03, etc.
            len_nested_tup_TMs = len(nested_tup_TMs)
            list_of_TMDs = long_list_of_TMDs[:len_nested_tup_TMs]
            # add that list to the dataframe (could also be added as a stringlist, but that's irritating somehow)
            #df.loc[acc, 'list_of_TMDs'] = list_of_TMDs
            df.set_value(acc, "list_of_TMDs", list_of_TMDs)
            # set seq for slicing
            full_seq = df.loc[acc, "full_seq"]
            # topo = dft.loc[acc, "Topology"]
            # iterate through all the TMDs of that protein, slicing out the sequences
            for i in range(len(list_of_TMDs)):
                TMD = list_of_TMDs[i]
                tup = nested_tup_TMs[i]
                df.loc[acc, "%s_start" % TMD] = tup[0]
                df.loc[acc, "%s_end" % TMD] = tup[1]
                df.loc[acc, "%s_seq" % TMD] = utils.slice_with_listlike(full_seq, tup)
                df.loc[acc, "%s_seqlen" % TMD] = len(df.loc[acc, "%s_seq" % TMD])
                # dft.loc[acc, TMD + "_top"] = utils.slice_with_listlike(topo, tup)
            # add signal peptides and their corresponding values to list_of_TMDs
            if analyse_sp == True:
                if type(df_parsed.loc[acc, 'SP01_seq']) == str:
                    list_of_TMDs.append('SP01')
                    df.set_value(acc, "list_of_TMDs", list_of_TMDs)

                # # code necessary for TMSEG signal peptides - depreciated by MO 20.04.2017
                # SiPe_indices = df.loc[acc, 'SiPe_indices']
                # if SiPe_indices != []:
                #     df.loc[acc, 'SP01_start'] = SiPe_indices[0]
                #     df.loc[acc, 'SP01_end'] = SiPe_indices[-1]
                #     df.loc[acc, 'SP01_seq'] = full_seq[SiPe_indices[0]:SiPe_indices[-1]+1]
                #     list_of_TMDs.append('SP01')
                #     df.set_value(acc, "list_of_TMDs", list_of_TMDs)

            if n % 50 == 0 and n != 0:
                sys.stdout.write(". ")
                sys.stdout.flush()
                if n % 500 == 0:
                    sys.stdout.write("\n")
                    sys.stdout.flush()

        # slice out the nonTM segments with a function
        # note that for some reason, this is very slow after merging the dataframes
        df = slice_nonTMD_in_prot_list(df)


        cols_to_drop = ['M_indices', 'SiPe_indices', 'Membrane_Borders', 'TM_indices']
        df = pd.merge(df, df_parsed, left_index=True, right_index=True, suffixes=('', '_list_parsed'))
        df.drop(cols_to_drop, axis=1, inplace=True)

    elif os.path.isfile(TMSEG_top_txtoutput_path):
        """ PARSE DATA WITH THE FOLLOWING FORMAT, proteins listed one after each other

        IMPORTANT : this format is sub-optimal, because the sequences come from uniprot, and the predictions from TMPRED

        Can only be trusted when they are from the same date: best to use TMPRED output which also contains the orig sequence.

        ---
        ID: A4ZUB1
        # TRANSMEM	6	18	4
        # TRANSMEM	50	67	7
        SIG: SIGNAL 1 22 {ECO:0000255}.
        TMH: TRANSMEM 53 69 Helical. {ECO:0000255}.
        ---
        """
        # if the regions column in the lists tab is "TM01" instead of the usual "TM", take only the first TM
        take_only_the_first_TM = s["regions"] == "TM01"

        # create dataframe for text topology (dftt)
        dftt = pd.DataFrame()
        with open(TMSEG_top_txtoutput_path, "r") as f:
            acc_counter = 0
            for line in f:
                if line[0:4] == "ID: ":
                    acc = line.split(" ")[1].strip("\n")
                    dftt.loc[acc_counter, "acc"] = acc
                    acc_counter += 1
                    # reset the TM_counter
                    TM_counter = 1
                if line[0:10] == "# TRANSMEM":
                    if TM_counter > 1:
                        if take_only_the_first_TM:
                            # skip to next line, as the first TM is already taken
                            continue

                    # split by tab
                    split = line.split("\t")
                    # the start is split[1] (end is not really necessary here)
                    start = split[1]
                    # note that acc_counter += 1 is already + 1 for the next protein,
                    # therefore the dftt.loc is acc_counter-1
                    dftt.loc[acc_counter - 1, "TM{:02d}_start".format(TM_counter)] = start
                    end = split[2]
                    # note that acc_counter += 1 is already + 1 for the next protein,
                    # therefore the dftt.loc is acc_counter-1
                    dftt.loc[acc_counter - 1, "TM{:02d}_end".format(TM_counter)] = end
                    TM_counter += 1
        # add an extra number_of_TMDs column, so they can be counted consistently
        dftt["number_of_TMDs"] = 0
        for row in dftt.index:
            # drop TM02_start etc if they don't contain data
            subset = dftt.loc[row, :].dropna()
            # count columns
            n_cols = subset.shape[0]
            # calculate number of columns (TM01_start, TM01_end) /2, which is the number of TMDs
            number_of_TMDs = int((n_cols - 2) / 2)
            dftt.loc[row, "number_of_TMDs"] = number_of_TMDs
            dftt.loc[row, "list_of_TMDs"] = str(["TM{:02d}".format(n) for n in range(1, number_of_TMDs + 1)])
        # set the acc as the index, so it can be merged with df_parsed
        dftt.set_index("acc", drop=False, inplace=True)
        # save temp csv with TMSEG output
        TMSEG_txtoutput_parsed_csv = TMSEG_top_txtoutput_path[:-4] + "TMSEG_txtoutput_parsed.csv"
        dftt.to_csv(TMSEG_txtoutput_parsed_csv)

        df = pd.merge(dftt, df_parsed, left_index=True, right_index=True, suffixes=('', '_list_parsed'))

        # convert from string to python list
        if isinstance(df['list_of_TMDs'][0], str):
            df['list_of_TMDs'] = df['list_of_TMDs'].dropna().apply(lambda x: ast.literal_eval(x))

        # (re)define sequence length
        df["seqlen"] = df["full_seq"].str.len()

        # slice out all the TMD sequences
        for n, acc in enumerate(df.index):
            list_of_TMDs = df.loc[acc, "list_of_TMDs"]
            # add that list to the dataframe (could also be added as a stringlist, but that's irritating somehow)
            # set seq for slicing
            full_seq = df.loc[acc, "full_seq"]
            # iterate through all the TMDs of that protein, slicing out the sequences
            for i in range(len(list_of_TMDs)):
                TMD = list_of_TMDs[i]
                tuple_slice_indices = (df.loc[acc, "%s_start" % TMD], df.loc[acc, "%s_end" % TMD])
                df.loc[acc, "%s_seq" % TMD] = utils.slice_with_listlike(full_seq, tuple_slice_indices)
                df.loc[acc, "%s_seqlen" % TMD] = len(df.loc[acc, "%s_seq" % TMD])
            # add signal peptides and their corresponding values to list_of_TMDs
            if analyse_sp == True:
                if type(df_parsed.loc[acc, 'SP01_seq']) == str:
                    list_of_TMDs.append('SP01')
                    df.set_value(acc, "list_of_TMDs", list_of_TMDs)

        start = time.clock()
        # slice out the nonTM segments with a function
        # note that for some reason, this is very slow after merging the dataframes
        df = slice_nonTMD_in_prot_list(df)
        sys.stdout.write("\ntime taken : {:0.03f} s".format(time.clock() - start))


    else:
        raise FileNotFoundError("None of the TMSEG combined output files were found.")

    # define number of TMDs (includes Signal peptides!)
    df["number_of_TMDs"] = df["list_of_TMDs"].dropna().apply(lambda x : len(x))
    df['parse_TMSEG'] = True
    df.to_csv(pathdict["list_parsed_csv"], sep=",", quoting=csv.QUOTE_NONNUMERIC)
    logging.info("\n~~~~~~~~~~~~                       parse_TMSEG_results is finished                  ~~~~~~~~~~~~")

def slice_nonTMD_in_prot_list(df):
    """Using existing indices and sequence, slices out all the TMD sequences.

    Originally from TMSEG fasta parse code.

    Parameters
    ----------
    df : pd.DataFrame

    Returns
    -------
    df : pd.Dataframe
        returns the same dataframe, with added sliced sequences
    """
    # glance at the watch
    start = time.clock()

    for n, acc in enumerate(df.index):
        ''' ~~   SLICE nonTMD sequence  ~~ '''
        list_of_TMDs = df.loc[acc, 'list_of_TMDs']#.copy()
        # if 'SP01' in list_of_TMDs:
        #     list_of_TMDs.remove('SP01')
        # sequence from N-term. to first TMD
        TM01_start = int(df.loc[acc, 'TM01_start'])
        nonTMD_first = df.loc[acc, 'full_seq'][0: TM01_start - 1]
        # start the sequence with the first segment
        sequence_list = [nonTMD_first]
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
                # define start of next TMD
                start_next = int(df.loc[acc, '%s_start' % next_TMD])
                # end of current TMD
                end = int(df.loc[acc, '%s_end' % TMD])
                # middle sequence between TMDs
                between_TM_and_TMplus1 = df.loc[acc, 'full_seq'][end: start_next - 1]
                sequence_list.append(between_TM_and_TMplus1)
        last_TMD = list_of_TMDs[-1]
        # sequence from last TMD to C-term.
        lastTM_end = int(df.loc[acc, '%s_end' % last_TMD])
        seqlen = int(df.loc[acc, 'seqlen'])
        nonTMD_last = df.loc[acc, 'full_seq'][lastTM_end:seqlen]
        sequence_list.append(nonTMD_last)
        # join all the sequences together
        sequence = "".join(sequence_list)

        df.loc[acc, 'nonTMD_seq'] = sequence
        df.loc[acc, 'len_nonTMD'] = len(sequence)

        if n % 50 == 0 and n != 0:
            sys.stdout.write(".")
            sys.stdout.flush()
            if n % 500 == 0:
                sys.stdout.write("\n")
                sys.stdout.flush()
    # glance at the watch again. Ruminate on time passed
    sys.stdout.write("\ntime taken to slice nonTMD sequences : {:0.03f} s".format(time.clock() - start))

    return df

def getting_membrane_indices_from_helix_symbol(Topo_data):
    m_list = [i for i, topo in enumerate(Topo_data) if topo == "H"]  # find(Topo_data)
    return m_list

def getting_SiPe_indices_from_symbol(Topo_data):
    m_list = [i for i, topo in enumerate(Topo_data) if topo == "S"]  # find(Topo_data)
    return m_list