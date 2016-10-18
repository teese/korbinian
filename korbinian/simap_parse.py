import ast
import csv
import os
import pandas as pd
import pickle
import tarfile
import xml.etree.ElementTree as ET
import korbinian
import korbinian.utils as utils
import zipfile
from multiprocessing import Pool

def run_parse_simap_to_csv(pathdict, s, logging):
    """For a dataframe containing a list of proteins, for each protein parses the SIMAP XML file to a csv file.

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
    acc_not_in_homol_db_txt : txt
        List of uniprot accessions that ar not in the homologue database (e.g. SIMAP)
        Any XML files with "Query failed: could not find the query sequence (check your query parameters)"
        will be added to this list

    For each protein (see parse_SIMAP_to_csv() below):
        homol_df_orig_zip : zipfile
            Zipfile containing the following:
                SIMAP_align_pretty_csv : csv
                    CSV file containing the hit_number protein description and the pretty alignment for each homologue
                homol_df_orig_pickle : pickled pd.DataFrame
                    Dataframe containing all sequence extracted from the XML file.
                    This can be large, as it contains the full query, markup and match sequences

    """
    logging.info('~~~~~~~~~~~~            starting parse_SIMAP_to_csv           ~~~~~~~~~~~~')
    acc_not_in_homol_db = []
    if os.path.isfile(pathdict["acc_not_in_homol_db_txt"]):
        # Extracts accession numbers out of file
        with open(pathdict["acc_not_in_homol_db_txt"], "r") as source:
            for line in source:
                line = line.strip()
                acc_not_in_homol_db.append(line)

    # if multiprocessing is used, log only to the console
    p_dict_logging = logging if s["use_multiprocessing"] != True else utils.Log_Only_To_Console()
    # create list of protein dictionaries to process
    list_p = korbinian.utils.convert_summary_csv_to_input_list(s, pathdict, p_dict_logging, list_excluded_acc=acc_not_in_homol_db)

    # number of processes is the number the settings, or the number of proteins, whichever is smallest
    n_processes = s["multiprocessing_cores"] if s["multiprocessing_cores"] < len(list_p) else len(list_p)

    if s["use_multiprocessing"]:
        with Pool(processes=n_processes) as pool:
            parse_simap_list = pool.map(parse_SIMAP_to_csv, list_p)
        # log the list of protein results to the actual logfile, not just the console
        logging.info(parse_simap_list)
        try:
            df_parsed = pd.DataFrame(parse_simap_list)
            df_parsed.set_index(0, inplace=True)
            df_parsed.index.name = "acc"
            df_parsed.columns = ["finished", "error"]
            not_finished_df = df_parsed.loc[df_parsed.finished == False]
            finished_df = df_parsed.loc[df_parsed.finished == True]
            if not not_finished_df.empty:
                logging.info("\nparse_SIMAP_to_csv proteins not finished :\n\n{}\n".format(df_parsed.loc[df_parsed.finished == False]))
            if not finished_df.empty:
                logging.info("\nparse_SIMAP_to_csv proteins finished correctly:\n\n{}\n".format(df_parsed.loc[df_parsed.finished == True]))
            df_parsed["not_in_database"] = df_parsed.error.str.contains("not in simap database")
            new_acc_not_in_db_list = list(df_parsed.loc[df_parsed["not_in_database"]].index)
            new_acc_not_in_db_nr_set = set(new_acc_not_in_db_list) - set(acc_not_in_homol_db)
            # add accession number to the list of failed downloads
            with open(pathdict["acc_not_in_homol_db_txt"], "a") as source:
                for acc in new_acc_not_in_db_nr_set:
                    source.write("\n{}".format(acc))
        except (TypeError, IndexError):
            logging.info(parse_simap_list)
            print("TypeError, IndexError, parse_simap_list is not a list of 3-item tuples for some reason.")
    else:
        for p in list_p:
            parse_SIMAP_to_csv(p)
        logging.info('~~~~~~~~~~~~          parse_SIMAP_to_csv is finished          ~~~~~~~~~~~~')

def parse_SIMAP_to_csv(p):
    """ Parses the SIMAP XML file to csv for a single protein.

    Designed for use in multiprocessing, where logging.info will only print to the console, and the logfile will
    contain the messages in the return statements, telling if that protein was successful.

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
            If multiprocessing == True, logging.info etc will only print to console.
        p : protein-specific dictionary components
            acc, list_of_TMDs, description, TM01_seq, etc

    Saved Files and Figures
    -----------------------
    homol_df_orig_zip : zipfile
        Zipfile containing the following:
            SIMAP_align_pretty_csv : csv
                CSV file containing the hit_number protein description and the pretty alignment for each homologue
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
    pathdict, s, logging = p["pathdict"], p["s"], p["logging"]
    acc = p["acc"]
    print(acc, end=", ", flush=True)
    # if overwrite_simap_parsed_to_csv is False, skip proteins where the homol_df_orig_zip file exists
    if s["overwrite_simap_parsed_to_csv"] == False:
        if os.path.isfile(p['homol_df_orig_zip']):
            warning = "{} skipped, homologues already parsed to csv".format(p['protein_name'])
            logging.info(warning)
            return acc, False, warning
    #set up counters
    number_of_hits_missing_protein_node = 0
    num_hits_with_SW_align_node = 0
    number_of_hits_missing_smithWatermanAlignment_node = 0
    #number_of_hits_kept_for_statistical_analysis = 0  # number_of_hits
    protein_name = p['protein_name']
    ft_xml_path = p['SIMAP_feature_table_XML_path']
    homol_xml_path = p['SIMAP_homol_XML_path']
    SIMAP_tar = p['SIMAP_tar']
    homol_xml_filename = os.path.basename(homol_xml_path)
    logging.info('%s' % protein_name)
    #check which files exist
    homol_in_tar = utils.check_SIMAP_tarfile(SIMAP_tar, ft_xml_path, homol_xml_path, acc, logging, delete_corrupt=True)[-1]

    # NEW: XML is parsed if only the homol_in_tar (feature tables are not necessary)
    if homol_in_tar:
        # create subfolders, if they don't exist
        subfolder = os.path.dirname(p['homol_df_orig_zip'])
        utils.make_sure_path_exists(subfolder)

        #if the setting is "False", and you don't want to overwrite the files, skip this section
        if s["overwrite_homologue_csv_files"]:
            create_homol_csv = True
        else:
            #check if output file already exists
            if os.path.isfile(p['homol_df_orig_zip']):
                try:
                    dfs_test = utils.open_df_from_csv_zip(p['homol_df_orig_zip'])
                    description_of_first_hit = dfs_test.loc[1, 'description']
                    logging.info('Protein %s: homologues already converted to csv. (%s)' % (p["acc"], description_of_first_hit))
                    create_homol_csv = False
                except (EOFError, KeyError):
                    #file may be corrupted, if script stopped unexpectedly before compression was finished
                    logging.info('%s seems to be corrupted. File will be deleted.' % p['homol_df_orig_zip'])
                    os.remove(p['homol_df_orig_zip'])
                    create_homol_csv = True
            else:
                logging.info('%s not found, create_homol_csv = True' % p['homol_df_orig_zip'])
                create_homol_csv = True
        #if the files don't exist, or you want to overwrite them
        if create_homol_csv:
            #extract the tarfile so that it can be read as xml
            with tarfile.open(p['SIMAP_tar'], 'r:gz') as tar:
                SIMAP_homologues_XML_file_extracted = tar.extractfile(homol_xml_filename)

                #parse_uniprot the XML file with elementtree, define the 'root' of the XML file
                simap_homologue_tree = ET.parse(SIMAP_homologues_XML_file_extracted)
                simap_homologue_root = simap_homologue_tree.getroot()

                try:
                    error = simap_homologue_root[0][0][1][0].text
                    if "could not find the query sequence" in error:
                        message = "{} not in simap database".format(acc)
                        logging.info(message)
                        return acc, False, message
                except:
                    message = "{} XML could not be opened".format(acc)
                    logging.info(message)
                    return acc, False, message

                try:
                    p['SIMAP_created'] = simap_homologue_root[0][0][0][0][2][1][0].attrib["created"]

                    for parameters in simap_homologue_root[0][0][0][0].iter('parameters'):
                        p['SIMAP_input_seq_details_dict'] = str(parameters[0][0].attrib)
                        for SIMAP_filter in parameters.iter('filter'):
                            SIMAP_filter_string = SIMAP_filter.text
                        p['SIMAP_filter_string'] = str(SIMAP_filter_string)
                        for resultSpecification in parameters.iter('resultSpecification'):
                            SIMAP_resultSpecification_dict = resultSpecification.attrib
                        p['SIMAP_resultSpecification_dict'] = '"%s"' % SIMAP_resultSpecification_dict
                        for databases in parameters.iter('databases'):
                            database_details_dict = databases[0].attrib
                        p['database_details_dict'] = '"%s"' % database_details_dict
                        p['simap_version'] = simap_homologue_root[0][0][0][0][0].attrib['version']
                        p['SIMAP_total_hits'] = int(simap_homologue_root[0][0][0][1][0].attrib['total'])

                    if p['simap_version'] != '4.0':
                        logging.warning('WARNING! Your XML file is simap version %s,'
                                        'however this SIMAP parser was developed for SIMAP version 4.0.' %
                                         p['simap_version'])

                    query_sequence_node = simap_homologue_root[0][0][0][0][2][0][0]
                    ''' xxxx CURRENTLY THE df is filled with nan values,
                        but that doesn't make sense as the script seems to work
                    '''
                    p['query_md5'] = query_sequence_node.attrib['md5']
                    p['seqlen'] = int(query_sequence_node.attrib['length'])
                    p['query_selfscore'] = query_sequence_node.attrib['selfscore']
                    p['query_sequenceid'] = query_sequence_node.attrib['sequenceid']
                    p['total_number_of_simap_hits'] = query_sequence_node[0].attrib['number_hits']
                    p['query_sequence_from_homologue_XML_file'] = query_sequence_node[0][0].text
                    p['number_of_hits_in_homologue_XML_file'] = int(
                        simap_homologue_root[0][0][0][1][0].attrib['total'])
                except (IndexError, KeyError):
                    p['homol_XML_damaged'] = True
                    warning = "{} skipped, homologue XML seems to be damaged. Error in reading general query details.".format(protein_name)
                    logging.warning("{} skipped, homologue XML seems to be damaged. Error in reading general query details.".format(protein_name))
                    # skip to the next protein
                    return acc, False, warning

                #for each hit, save all the relevant data in the form of a dictionary,
                # so it can be added to a csv file or used in other calculations
                simap_homologue_hits = simap_homologue_root[0][0][0][1][0]

                #see if there are any hits at all
                try:
                    test2 = simap_homologue_root[0][0][0][1][0][0]
                    xml_contains_simap_homologue_hits = True
                except IndexError:
                    xml_contains_simap_homologue_hits = False

                if xml_contains_simap_homologue_hits:
                    """OLD AMINO ACID SUBSTITUTION CODE. THIS IS SLOW, AND GIVES NO SIGNIFICANT DIFFERENCE TO
                    AAIMON OR AASMON WITH THE SIMAP SMITH-WATERMAN MATRIX"""
                    #load the amino acid substitution matrices from the settings file
                    #list_of_aa_sub_matrices = s['aa_sub_matrices']
                    #import the amino acid substitution matrices
                    #utils.import_amino_acid_substitution_matrices()
                    #add the similarity ratios to the csv_header_for_SIMAP_homologue_file.
                    # These will depend on the individual settings
                    #                    if s['["mp_calculate_TMD_conservation_with_aa_matrices']:
                    #                        for j in range(s["gap_open_penalty_min"],
                    #                                       s["gap_open_penalty_max"],
                    #                                       s["gap_open_penalty_increment"]):
                    #                            gap_open_penalty = j
                    #                            gap_extension_penalty = j
                    #                            for matrix_name in list_of_aa_sub_matrices:
                    #                                column_name = 'sim_ratio_%s_gapo%i' % (matrix_name.replace("'", "")[0:-7], j)
                    #                                csv_header_for_SIMAP_homologue_file.append(column_name)
                    #import the necessary matrices
                    #for matrix_name in list_of_aa_sub_matrices:
                    #matrix = matrix_name[0:-7]
                    #from Bio.SubsMat.MatrixInfo import matrix as matrix_name

                    SIMAP_orig_csv = p['homol_df_orig_zip'][:-4] + ".csv"
                    #fasta_file_path = p['fasta_file_path']

                    #create an empty file
                    open(SIMAP_orig_csv, 'w').close()

                    #reopen to add match details iteratively from dictionary
                    with open(SIMAP_orig_csv, 'a') as csvfile:

                        #set up a bool to catch those files where not a single hit actually gives data
                        at_least_one_hit_contains_SW_node = False

                        for hit in simap_homologue_hits:
                            match_details_dict = {}

                            #add desired hit information to the dictionary for transfer to csv
                            hit_num = int(hit.attrib['number'])
                            match_details_dict['hit_num'] = hit_num
                            match_details_dict['md5'] = hit[1].attrib['md5']

                            #define the major nodes in the XML-file
                            try:
                                protein_node = hit[1][1]
                                hit_contains_protein_node = True
                            except IndexError:
                                hit_contains_protein_node = False
                                number_of_hits_missing_protein_node += 1
                                logging.warning('%s hit %s contains no protein node' % (protein_name, match_details_dict['md5']))
                            if hit_contains_protein_node:
                                try:
                                    smithWatermanAlignment_node = hit[0][0][14]
                                    hit_contains_SW_node = True
                                    num_hits_with_SW_align_node += 1
                                except IndexError:
                                    hit_contains_SW_node = False
                                match_details_dict['hit_contains_SW_node'] = hit_contains_SW_node
                                #add the description. Add a custom name if it is the first (query) hit
                                if hit_num == 1:
                                    description = '%s_SIMAP_query_sequence' % protein_name
                                else:
                                    description = protein_node.attrib['description']
                                match_details_dict['description'] = description
                                try:
                                    databaseId = int(protein_node[1].attrib['databaseId'])
                                    match_details_dict['databaseId'] = int(protein_node[1].attrib['databaseId'])
                                except KeyError:
                                    databaseId = 0
                                    #match_details_dict['databaseId'] = int(0)
                                #databaseId = int(protein_node[1].attrib['databaseId'])
                                databasenode = protein_node[1]
                                match_details_dict['database'] = databasenode.attrib['name']
                                try:
                                    taxonomyNode = protein_node[2]
                                    match_details_dict['organism'] = taxonomyNode.attrib['name']
                                    match_details_dict['taxonomy_node_id'] = taxonomyNode.attrib['node_id']
                                    match_details_dict['taxonomy_rank'] = taxonomyNode.attrib['rank']
                                except IndexError:
                                    #sequence is from an unknown organism, as it has no database node
                                    match_details_dict['description'] += ', no_database_node'
                                    match_details_dict['organism'] = 'no_database_node'
                                    match_details_dict['taxonomy_node_id'] = 'no_database_node'
                                    match_details_dict['taxonomy_rank'] = 'no_database_node'
                                match_details_dict['len_full_match_seq'] = len(hit[1][0][0].text)
                                #len_full_match_seq = len(full_match_seq)
                                alignment_node = hit[0][0]
                                #E-value for hit
                                match_details_dict['FASTA_expectation'] = float(alignment_node[1].text)
                                #convert identity from e.g. 80 (80%) to 0.8
                                match_details_dict['FASTA_identity'] = float(alignment_node[3].text) / 100
                                #strangely, I think gappedIdentity is the identity EXCLUDING gaps, which is a better value to base judgements on. convert identity from e.g. 80 (80%) to 0.8
                                match_details_dict['FASTA_gapped_identity'] = float(alignment_node[4].text) / 100
                                '''xxx notes on the gapped identity
                                N.B The FASTA_gapped_identity data here is from the FASTA algorithm, that precedes the SW algorithm.
                                Occasionally they don’t match!!!
                                I calculate the TMD identity manually from the SW alignment, BUT
                                currently for the calculation of membranous/nonmembranous I use the gappedIdentity from the FASTA output
                                (the SW output inly has identity including gaps)
                                -    if I counted the gaps from the SW alignment, I COULD recalculate the gappedIdentity for the SW alignment
                                -    OR: I could simply remove the data where the FASTA and SW don’t match.
                                '''
                                #FASTA overlap should be the length of the aligned region after running the FASTA algorithm (alignment is not shown by SIMAP)
                                match_details_dict['FASTA_overlap'] = int(alignment_node[5].text)
                                match_details_dict['FASTA_query_coverage'] = float(alignment_node[11].text)
                                match_details_dict['FASTA_match_coverage'] = float(alignment_node[12].text)
                                #find the start and the stop of the hsp
                                querySeq = alignment_node[6]
                                match_details_dict['FASTA_query_start'] = int(querySeq.attrib['start'])
                                match_details_dict['FASTA_query_end'] = int(querySeq.attrib['end'])
                                matchSeq = alignment_node[7]
                                match_details_dict['FASTA_match_start'] = int(matchSeq.attrib['start'])
                                match_details_dict['FASTA_match_end'] = int(matchSeq.attrib['end'])
                                """OLD CALCULATIONS THAT ARE NOW CONVERTED TO PANDAS ARRAY-WISE FUNCTIONS"""
                                #some parameters that are needed for identity calculations later
                                #FASTA_num_ident_res = FASTA_identity / 100.0 * FASTA_overlap
                                #is_start_of_TMD_in_FASTA = True if FASTA_query_start <= TMDstart else False
                                #is_end_of_TMD_in_FASTA = True if TMDend <= FASTA_query_end else False
                                #is_TMD_in_FASTA_alignment = True if all([is_start_of_TMD_in_FASTA, is_end_of_TMD_in_FASTA]) else False
                                '''***********************if the TMD region is actually covered by the hsp, then conduct some further analyses of the match TMD region*************************'''
                                if hit_contains_SW_node:
                                    query_align_seq = ''
                                    '''For the moment, there is no need to put the whole match hsp sequence into the csv file'''
                                    #for smithWatermanAlignment in alignment_node.iter('smithWatermanAlignment'):
                                    match_details_dict['SW_query_score_ratio'] = smithWatermanAlignment_node[0].text
                                    match_details_dict['SW_match_score_ratio'] = smithWatermanAlignment_node[1].text
                                    match_details_dict['SW_query_coverage'] = smithWatermanAlignment_node[2].text
                                    match_details_dict['SW_match_coverage'] = smithWatermanAlignment_node[3].text
                                    match_details_dict['SW_coverage_ratio'] = smithWatermanAlignment_node[4].text
                                    match_details_dict['align_pretty'] = smithWatermanAlignment_node[8].text
                                    match_details_dict['SW_alignment_seq1offset'] = int(smithWatermanAlignment_node.attrib['alignment-seq1offset'])
                                    match_details_dict['SW_alignment_seq2offset'] = int(smithWatermanAlignment_node.attrib['alignment-seq2offset'])
                                    match_details_dict['SW_identity'] = float(smithWatermanAlignment_node.attrib['identity'])
                                    match_details_dict['SW_similarity'] = float(smithWatermanAlignment_node.attrib['similarity'])
                                    #Get the full sequences. Note that they greatly increase the size of the csv file.
                                    match_details_dict['query_align_seq'] = smithWatermanAlignment_node[5].text
                                    match_details_dict['align_markup_seq'] = smithWatermanAlignment_node[6].text
                                    match_details_dict['match_align_seq'] = smithWatermanAlignment_node[7].text
                                else:
                                    number_of_hits_missing_smithWatermanAlignment_node += 1
                                if hit_num == 1:
                                    #sort
                                    csv_header_for_SIMAP_homologue_file = sorted(list(match_details_dict.keys()))
                                    #save the csv header to the csv file
                                    writer = csv.writer(csvfile, delimiter=',', quotechar='"', lineterminator='\n',quoting=csv.QUOTE_NONNUMERIC, doublequote=True)
                                    writer.writerow(csv_header_for_SIMAP_homologue_file)
                                #save the match_details_dict as a line in the csv file
                                writer = csv.DictWriter(csvfile, fieldnames=csv_header_for_SIMAP_homologue_file,
                                                        extrasaction='ignore', delimiter=',', quotechar='"',
                                                        lineterminator='\n', quoting=csv.QUOTE_NONNUMERIC,
                                                        doublequote=True)
                                writer.writerow(match_details_dict)

                    # open csv as a dataframe,
                    df_homol = pd.read_csv(SIMAP_orig_csv, sep=",", quoting=csv.QUOTE_NONNUMERIC, index_col="hit_num")
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

                    # restrict to just a few columns including the align_pretty that might be useful to check manually
                    df_pretty = df_homol[["FASTA_gapped_identity", "organism", "description", "align_pretty"]]
                    # save the align_pretty to csv
                    df_pretty.to_csv(p['SIMAP_align_pretty_csv'], sep=',', quoting=csv.QUOTE_NONNUMERIC)
                    # drop the align_pretty column from the orig dataframe
                    df_homol.drop('align_pretty', axis=1, inplace=True)
                    # save the whole dataframe as a pickle for faster opening later
                    with open(p['homol_df_orig_pickle'], "wb") as pick:
                        pickle.dump(df_homol, pick, protocol=pickle.HIGHEST_PROTOCOL)
                    # either create new zip and add ("w"), or open existing zip and add "a"
                    with zipfile.ZipFile(p['homol_df_orig_zip'], mode="w", compression=zipfile.ZIP_DEFLATED) as zipout:
                        #zipout.write(SIMAP_orig_csv, arcname=os.path.basename(SIMAP_orig_csv))
                        zipout.write(p['SIMAP_align_pretty_csv'], arcname=os.path.basename(p['SIMAP_align_pretty_csv']))
                        zipout.write(p['homol_df_orig_pickle'], arcname=os.path.basename(p['homol_df_orig_pickle']))
                    # delete temporary uncompressed files
                    os.remove(SIMAP_orig_csv)
                    os.remove(p['SIMAP_align_pretty_csv'])
                    os.remove(p['homol_df_orig_pickle'])
                    return acc, True, "0"

def get_phobius_TMD_region(feature_table_root):
    """Old function, no longer in use."""
    for feature in feature_table_root[0]:
        if 'PHOBIUS' and 'TMHelix' in feature.attrib.values():
            for begin in feature.iter('begin'):
                TMD_start = int(begin.attrib['position']) #same as feature[0][0].attrib['position'], but more resistant to parser breaking
            for end in feature.iter('end'):
                TMD_end = int(end.attrib['position'])
            TMD_length = TMD_end - TMD_start + 1
            #logging.info('phobius prediction: TMD start = %s, TMD end = %s' % (TMD_start, TMD_end))
                #logging.info(begin.attrib['position']) #same as feature[0][0].attrib['position'], but more resistant to parser breaking
        else:
            TMD_start, TMD_end, TMD_length = 0, 0, 0
    return TMD_start, TMD_end, TMD_length

def get_TMHMM_TMD_region(root):
    """Old function, no longer in use."""
    for feature in root[0]:
        for subfeature in feature:
            if 'TMHMM' and 'TMHelix' in subfeature.attrib.values():
                for begin in subfeature.iter('begin'):
                    TMD_start = begin.attrib['position'] #same as feature[0][0].attrib['position'], but more resistant to parser breaking
                for end in subfeature.iter('end'):
                    TMD_end = end.attrib['position']
                #logging.info('TMHMM prediction: TMD start = %s, TMD end = %s' % (TMD_start, TMD_end))
                #logging.info(begin.attrib['position']) #same as feature[0][0].attrib['position'], but more resistant to parser breaking
            else:
                TMD_start, TMD_end = 0, 0
    return TMD_start, TMD_end