import ast
import csv
import korbinian
import korbinian.utils as utils
import numpy as np
import os
import pandas as pd
import sys
import time
# import debugging tools
from korbinian.utils import pr, pc, pn, aaa

def get_TM_indices_from_TMSEG_topo_str(topo_str, TM_symbol="H"):
    """Get TM indices from TMSEG topology string.

    Code is not very elegant in comparison to a regex approach, but works fine.

    Parameters
    ----------
    topo_str : str
        Topology string output from TMSEG.
        H = TM helix
        Note that TM orientation (N-cyto or N-out) is currently not extracted.
        E.g. "11111111111111HHHHHHHHHHHHHHHHHHH222222222222222222222222222222222222222222222222222222222222222222222222HHHHHHHHHHHHHHHHHHHHHHHHH"
        "111111111111111111111111HHHHHHHHHHHHHHHHHHHHH222222222222222HHHHHHHHHHHHHHHHHHHH111111111111111111111111111111111111HHHHHHHHHHHHHHHHHHHHHHH"
        "22222222222222222222222222222222HHHHHHHHHHHHHHHHHHHHHH1111111111111111111111111HHHHHHHHHHHHHHHHHHHHHHH22222222222222222222222222222222222222"
        "2222HHHHHHHHHHHHHHHHHHHHH11111111111111111111111111111111111111"

    Returns
    -------
    TM_indices : tuple
        Nested tuple with start and end of all TM helices in topology string.
        UniProt indexing is used ("HHH111" is (1:3), not (0:3))
        E.g.
        ((15, 33),1
         (106, 130),
         (155, 175),
         (191, 210),
         (247, 269),
         (302, 323),
         (349, 371),
         (414, 434))
    """
    if TM_symbol in topo_str:
        # get indices (eg. [28, 29, 30, 31, 32, 33, 34, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 72, 73, 74, 75, 76])
        M_indices = get_list_TM_residues_from_topo_string(topo_str, TM_symbol)
        #SiPe_indices = get_signal_peptide_indices_from_TMSEG_topo(topo_str)
        # get borders to TM regions(eg. [28, 34, 58, 68, 72, 76])
        m_borders = []
        m_borders.append(M_indices[0])
        m_borders = korbinian.prot_list.parse_OMPdb.check_for_border(M_indices, m_borders)
        # add the last membrane index (e.g. 33)
        m_borders.append(M_indices[-1])
        # convert to nested tuples
        TM_indices = convert_alternating_list_to_nested_tuples(m_borders)
        return TM_indices
    else:
        return ()

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

    # set to be an empty string, which avoids the error related to inserting a python list into a cell
    # "ValueError: Must have equal len keys and value when setting with an iterable"
    df['list_of_TMDs_excl_SP'] = ""

    for n, acc in enumerate(df.index):
        ''' ~~   SLICE nonTMD sequence  ~~ '''
        # list of TMDs excluding signal peptides
        list_of_TMDs_excl_SP = df.loc[acc, 'list_of_TMDs']
        # set value to avoid errors adding a list to a cell
        df.set_value(acc, 'list_of_TMDs_excl_SP', list_of_TMDs_excl_SP)

        seqstart = 0
        # if any protein in list conatains a SP
        if 'SP01_end' in df.columns:
            # if THIS PARTICULAR PROTEIN contains a signal peptide sequence
            if isinstance(df.loc[acc, 'SP01_seq'], str):
                # change sequence start for nonTM to the end of the signal peptide
                seqstart = int(df.loc[acc, 'SP01_end'])
                # add the SP01 to the list of TMDs
                df.set_value(acc, 'list_of_TMDs', ["SP01"] + list_of_TMDs_excl_SP)

        # sequence from N-term. to first TMD
        TM01_start = int(df.loc[acc, 'TM01_start'])

        # NOTE THAT THIS USED TO BE nonTMD_first = df.loc[acc, 'full_seq'][0: TM01_start -1], but indexing missed the last nonTM residue.
        nonTMD_first = df.loc[acc, 'full_seq'][seqstart: TM01_start - 1]
        # start the sequence with the first segment
        sequence_list = [nonTMD_first]
        # only for multipass proteins, generate sequences between TMDs
        if len(list_of_TMDs_excl_SP) == 0:
            # no TMDs are annotated, skip to next protein
            continue
        # for multipass proteins
        elif len(list_of_TMDs_excl_SP) > 1:
            for TM_Nr in range(len(list_of_TMDs_excl_SP) - 1):
                # the TMD is the equivalent item in the list
                TMD = list_of_TMDs_excl_SP[TM_Nr]
                # the next TMD, which contains the end index, is the next item in the list
                next_TMD = list_of_TMDs_excl_SP[TM_Nr + 1]
                # define start of next TMD
                start_next = int(df.loc[acc, '%s_start' % next_TMD])
                # end of current TMD
                end_current = int(df.loc[acc, '%s_end' % TMD])
                # middle sequence between TMDs
                # note the "start_next - 1", used to convert uniprot indices to python indices
                between_TM_and_TMplus1 = df.loc[acc, 'full_seq'][end_current: start_next - 1]
                sequence_list.append(between_TM_and_TMplus1)
        last_TMD = list_of_TMDs_excl_SP[-1]
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

def get_list_TM_residues_from_topo_string(Topo_data, TM_symbol):
    # get list of membrane indices
    # note that this is UNIPROT indexing, not python indexing
    m_list = [i+1 for i, topo in enumerate(Topo_data) if topo == TM_symbol]  # find(Topo_data)
    return m_list

def convert_alternating_list_to_nested_tuples(x):
    return tuple(zip(x[::2], x[1::2]))

def parse_TMSEG_results_DEPRECATED(pathdict, s, logging):
    """DEPRECATED METHOD BASED ON LARGE FILE OF ALL TMSEG RESULTS

    USE METHODS BASED ON INDIVIDUAL TMSEG DATAFILES INSTEAD.

    """
    logging.info("~~~~~~~~~~~~                        starting parse_TMSEG_results_DEPRECATED                    ~~~~~~~~~~~~")
    # create or open dataframe for protein list summary
    if os.path.isfile(pathdict["prot_list_summary_csv"]):
        df_PLS = pd.read_csv(pathdict["prot_list_summary_csv"], index_col=0)
    else:
        df_PLS = pd.DataFrame(columns=["v", "date"])
    # get the timestamp for current time
    t = time.ctime(time.time())

    list_number = s['list_number']

    # define the uniprot directory with selected records
    uniprot_dir = os.path.join(s["data_dir"], 'uniprot')
    selected_uniprot_records_flatfile = os.path.join(uniprot_dir, 'selected', 'List%02d_selected_uniprot_records_flatfile.txt' % list_number)
    n_aa_before_tmd = s["n_aa_before_tmd"]
    n_aa_after_tmd = s["n_aa_after_tmd"]
    list_parsed_csv = pathdict["list_parsed_csv"]
    # check if the lists tab says to analyse the signal peptides
    analyse_sp = True if "SiPe" in s["regions"] else False
    output = korbinian.prot_list.uniprot_parse.parse_flatfile_to_csv(selected_uniprot_records_flatfile, n_aa_before_tmd, n_aa_after_tmd, analyse_sp, logging, list_parsed_csv, slice=False)
    logging.info(output)

    TMSEG_fastalike_path = pathdict['TMSEG_fastalike']
    TMSEG_top_txtoutput_path = pathdict['TMSEG_top_txtoutput']
    TMSEG_nonTM_outpath = pathdict['TMSEG_nonTM']

    df_parsed = pd.read_csv(pathdict["list_parsed_csv"], sep=",", quoting=csv.QUOTE_NONNUMERIC, index_col=0, low_memory=False)

    columns_to_keep = ['organism_domain', 'uniprot_acc', 'uniprot_all_accessions', 'uniprot_entry_name', 'uniprot_features',
                       'uniprot_orgclass', 'uniprot_SiPe', 'singlepass', 'typeI', 'typeII', 'uniprot_KW', 'organism', 'prot_descr', 'membrane',
                       'multipass', 'gene_name', 'comments_subcellular_location_uniprot', 'uniprot_SiPe', 'full_seq']

    # # for datasets without SP found, turn off analyse_sp
    # if analyse_sp == True and 'SP01_start' in df_parsed.columns:
    #     columns_to_keep = ['SP01_start', 'SP01_end', 'SP01_seq']
    # else:
    #     analyse_sp == False

    acc_list_orig = list(df_parsed.index)

    if os.path.isfile(TMSEG_fastalike_path):
        df_PLS.loc["TMSEG_fastalike_path", :] = ("exists", t)
        sys.stdout.write("Extracting topology from TMSEG_fastalike file.")
        # DEPRECATED drop the full sequence, and get from TMSEG
        #df_parsed.drop('full_seq', axis=1, inplace=True)

        # read data from file
        # list will have acc, seq, topology, acc, seq, topology etc
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
        df_TMSEG = pd.DataFrame(index=input_data[0::3])

        # add the signal peptide definitions from UniProt, to be used for slicing the nonTMD etc later
        if analyse_sp:
            for col in ['SP01_start', 'SP01_end', 'SP01_seq']:
                df_TMSEG[col] = df_parsed[col]

        # drop unnecessary columns from df_parsed, to be merged later
        df_parsed = df_parsed[columns_to_keep]

        # add selected columns from input_data list
        #df_TMSEG['uniprot_entry_name'] = input_data[1::5]
        #df_TMSEG['prot_descr'] = input_data[2::5]
        df_TMSEG['full_seq'] = input_data[1::3]
        df_TMSEG['topo'] = input_data[2::3]

        acc_list_TMSEG = df_TMSEG.index.tolist()

        TMSEG_avail_list = set(acc_list_TMSEG).intersection(set(acc_list_orig))
        TMSEG_unavail_list = list(set(acc_list_orig) - set(acc_list_TMSEG))

        df_PLS.loc["n_prot_TMSEG_file"] = (len(acc_list_TMSEG), t)

        # create a boolean whether the TMSEG topology is available
        df_parsed.loc[TMSEG_avail_list,"TMSEG_avail"] = True
        df_parsed.loc[TMSEG_unavail_list, "TMSEG_avail"] = False

        # drop proteins from df_TMSEG that are not in the listxx_parsed.csv
        df_TMSEG = df_TMSEG.loc[TMSEG_avail_list, :]

        fa_dir = pathdict['TMSEG_unavail_fa_dir']
        utils.make_sure_path_exists(fa_dir)
        for acc in TMSEG_unavail_list:
            out_fasta = os.path.join(fa_dir, "{}.fasta".format(acc))
            seq = df_parsed.loc[acc, "full_seq"]
            with open(out_fasta, "w") as f:
                f.write(">{}\n{}".format(acc, seq))

        n_prot_TMSEG_file_not_in_list = len(set(acc_list_TMSEG) - set(acc_list_orig))
        logging.info("n_prot_TMSEG_file_not_in_list as not in listxx_parsed.csv = {} ({} remaining)".format(n_prot_TMSEG_file_not_in_list, len(TMSEG_avail_list)))
        df_PLS.loc["n_prot_TMSEG_file_not_in_list"] = (n_prot_TMSEG_file_not_in_list, t)

        if df_TMSEG.shape[0] == 0:
            return sys.stdout.write('no remaining proteins in list!')

        # get list of uniprot accessions of proteins where no transmembrane region was predicted
        list_nonTMD = []
        for acc in df_TMSEG.index:
            if 'N' in df_TMSEG.loc[acc, 'topo']:
                list_nonTMD.append(acc)

        # write list of nonTM proteins to file
        # outpath = '/Volumes/Musik/Databases/TMSEG/humanU90_nonTMD.txt'
        file = open(TMSEG_nonTM_outpath, 'w')
        for line in list_nonTMD:
            file.write('{}\n'.format(line))
        file.close()

        # drop proteins that do not contain TM regions
        df_TMSEG = df_TMSEG.drop(list_nonTMD)

        # create a boolean whether the TMSEG topology is available
        TMSEG_avail_and_TM = set(TMSEG_avail_list) - set(list_nonTMD)
        TMSEG_avail_but_SOL = set(acc_list_orig).intersection(set(list_nonTMD))
        df_parsed["membrane"] = np.nan
        df_parsed.loc[TMSEG_avail_and_TM, "membrane"] = True
        df_parsed.loc[TMSEG_avail_but_SOL, "membrane"] = False

        # add seqlen and indices for all TMD and SiPe regions
        df_TMSEG["seqlen"] = df_TMSEG.full_seq.apply(lambda x: len(x))
        #df_TMSEG['M_indices'] = df_TMSEG.topo.apply(get_list_TM_residues_from_topo_string)
        #df_TMSEG['SiPe_indices'] = df_TMSEG.topo.apply(get_list_TM_residues_from_topo_string, args=("S"))

        df_TMSEG['TM_indices'] = df_TMSEG.topo.apply(get_TM_indices_from_TMSEG_topo_str)
        df_TMSEG['SiPe_indices'] = df_TMSEG.topo.apply(get_TM_indices_from_TMSEG_topo_str, args=("S"))

        # # Creating new list (nested list)
        # nested_list_of_membrane_borders = []
        #
        # ########################################################################################
        # #                                                                                      #
        # #              Extract the membrane indices in UniProt Indexing style                  #
        # #                                                                                      #
        # ########################################################################################
        # # Filling nest with lists of start and end-points
        # for m_index_list in df_TMSEG.M_indices:
        #     m_borders = []
        #     # add the first membrane index (e.g. 13)
        #     m_borders.append(m_index_list[0])
        #     m_borders = korbinian.prot_list.parse_OMPdb.check_for_border(m_index_list, m_borders)
        #     # add the last membrane index (e.g. 33)
        #     m_borders.append(m_index_list[-1])
        #     nested_list_of_membrane_borders.append(m_borders)
        #
        # # DEPRECATED
        # #FOR CONSISTENCY, LEAVE INDEXING STYLE AS UNIPROT
        # # ########################################################################################
        # # #                                                                                      #
        # # #            Convert to python indexing style (NECESSARY?? NOT COMPAT WITH UNIPROT!)   #
        # # #                                                                                      #
        # # ########################################################################################
        # # array_membrane_borders = np.array(nested_list_of_membrane_borders)
        # # nested_list_m_borders_python_indexstyle = []
        # # for subarray in array_membrane_borders:
        # #     # convert to array
        # #     subarray = np.array(subarray)
        # #     # add 1 to the second index number, to allow slicing
        # #     subarray[1::2] = subarray[1::2] + 1
        # #     # add to list with corrected values, python index style
        # #     nested_list_m_borders_python_indexstyle.append(list(subarray))
        #
        # # Creating new column, which contains start and end-points
        # #df_TMSEG["Membrane_Borders"] = nested_list_m_borders_python_indexstyle
        #
        # df_TMSEG["Membrane_Borders"] = nested_list_of_membrane_borders
        #
        # # Creating new column, which contains the number of TMDS
        # #df_TMSEG["number_of_TMDs"] = df_TMSEG.Membrane_Borders.apply(lambda x: len(x) / 2)
        #
        # df_TMSEG["TM_indices"] = df_TMSEG["Membrane_Borders"].apply(lambda x: tuple(zip(x[::2], x[1::2])))

        # create a list of [TM01, TM02, TM03, etc.
        long_list_of_TMDs = []
        for i in range(1, 50):
            long_list_of_TMDs.append("TM{:02d}".format(i))

        ## for the .set_value function, set dtype as object
        df_TMSEG["list_of_TMDs"] = ""
        df_TMSEG["list_of_TMDs"].astype(object)

        sys.stdout.write('slicing TMD and nonTMD sequences:\n')

        for n, acc in enumerate(df_TMSEG.index):
            # get nested tuple of TMDs
            nested_tup_TMs = df_TMSEG.loc[acc, "TM_indices"]
            # slice long list of TMD names to get an appropriate list for that protein [TM01, TM02, TM03, etc.
            len_nested_tup_TMs = len(nested_tup_TMs)
            list_of_TMDs = long_list_of_TMDs[:len_nested_tup_TMs]
            # add that list to the dataframe (could also be added as a stringlist, but that's irritating somehow)
            #df_TMSEG.loc[acc, 'list_of_TMDs'] = list_of_TMDs
            df_TMSEG.set_value(acc, "list_of_TMDs", list_of_TMDs)
            # set seq for slicing
            full_seq = df_TMSEG.loc[acc, "full_seq"]
            # topo = dft.loc[acc, "Topology"]
            # iterate through all the TMDs of that protein, slicing out the sequences
            for i, TMD in enumerate(list_of_TMDs):
                TMD = list_of_TMDs[i]
                start, end = nested_tup_TMs[i]
                # with UniProt indexing, need to slice with -1, not like python index style
                df_TMSEG.loc[acc, "%s_start" % TMD] = start
                df_TMSEG.loc[acc, "%s_end" % TMD] = end
                # for python indexing of the TMD rather than uniprot, the start should be minus 1
                python_indexing_tuple = (start - 1, end)
                df_TMSEG.loc[acc, "%s_seq" % TMD] = utils.slice_with_listlike(full_seq, python_indexing_tuple)
                df_TMSEG.loc[acc, "%s_seqlen" % TMD] = len(df_TMSEG.loc[acc, "%s_seq" % TMD])
                # dft.loc[acc, TMD + "_top"] = utils.slice_with_listlike(topo, tup)

            #DEPRECATED, ONLY REINSTATE IF YOU REALLY WANT TMSEG SP DEFINITIONS TO STAY
            # # add signal peptides and their corresponding values to list_of_TMDs
            # if analyse_sp == True:
            #     if type(df_parsed.loc[acc, 'SP01_seq']) == str:
            #         list_of_TMDs.append('SP01')
            #         df_TMSEG.set_value(acc, "list_of_TMDs", list_of_TMDs)

                # # code necessary for TMSEG signal peptides - depreciated by MO 20.04.2017
                # SiPe_indices = df_TMSEG.loc[acc, 'SiPe_indices']
                # if SiPe_indices != []:
                #     df_TMSEG.loc[acc, 'SP01_start'] = SiPe_indices[0]
                #     df_TMSEG.loc[acc, 'SP01_end'] = SiPe_indices[-1]
                #     df_TMSEG.loc[acc, 'SP01_seq'] = full_seq[SiPe_indices[0]:SiPe_indices[-1]+1]
                #     list_of_TMDs.append('SP01')
                #     df_TMSEG.set_value(acc, "list_of_TMDs", list_of_TMDs)

            if n % 50 == 0 and n != 0:
                sys.stdout.write(". ")
                sys.stdout.flush()
                if n % 500 == 0:
                    sys.stdout.write("\n")
                    sys.stdout.flush()

        # slice out the nonTM segments with a function
        # note that for some reason, this is very slow after merging the dataframes
        df_TMSEG = slice_nonTMD_in_prot_list(df_TMSEG)

        #df_TOP = pd.merge(df_parsed, df_TMSEG, how="left", left_on=True, suffixes=('_list_parsed', ""))# left_index=True, right_index=False,
        df_TOP = df_parsed.merge(df_TMSEG, how="left", suffixes=('_list_parsed', ""))  # left_index=True, right_index=False,

        # actually, I'd prefer to keep these for troubleshooting purposes
        # cols_to_drop = ['M_indices', 'SiPe_indices', 'Membrane_Borders', 'TM_indices']
        # df_TMSEG.drop(cols_to_drop, axis=1, inplace=True)

    elif os.path.isfile(TMSEG_top_txtoutput_path):
        df_PLS.loc["TMSEG_top_txtoutput_path", :] = ("exists", t)
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
        df_TOP = slice_nonTMD_in_prot_list(df)
        sys.stdout.write("\ntime taken : {:0.03f} s".format(time.clock() - start))

    else:
        raise FileNotFoundError("None of the TMSEG combined output files were found.")

    # define number of TMDs (includes Signal peptides!)
    df_TOP["number_of_TMDs"] = df_TOP["list_of_TMDs"].dropna().apply(lambda x : len(x))
    df_TOP['parse_TMSEG'] = True
    df_TOP.to_csv(pathdict["list_parsed_csv"], sep=",", quoting=csv.QUOTE_NONNUMERIC)
    logging.info("\n~~~~~~~~~~~~                       parse_TMSEG_results_DEPRECATED is finished                  ~~~~~~~~~~~~")

# def get_signal_peptide_indices_from_TMSEG_topo(Topo_data):
#     # as above for membrane regions
#     sp_list = [i for i, topo in enumerate(Topo_data) if topo == "S"]  # find(Topo_data)
#     return sp_list
