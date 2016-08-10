import ast
import csv
import os
import pandas as pd
import tarfile
import xml.etree.ElementTree as ET
import korbinian.mtutils as utils
import zipfile

def parse_SIMAP_to_csv(pathdict, set_, logging):
    counter_XML_to_CSV = 0
    logging.info('~~~~~~~~~~~~  starting parse_SIMAP_to_csv   ~~~~~~~~~~~~')
    df = pd.read_csv(pathdict["list_summary_csv"])
    # #filter to remove sequences where no TMDs are found
    # df = df.loc[df['list_of_TMDs'].notnull()]
    # #filter to remove sequences where no TMDs are found (if string)
    # df = df.loc[df['list_of_TMDs'] != 'nan']
    #iterate over the dataframe, excluding any proteins that do not have a list of TMDs
    for acc in df.loc[df['list_of_TMDs'].notnull()].loc[df['list_of_TMDs'] != 'nan'].index:
        #set up counters
        number_of_hits_missing_protein_node = 0
        num_hits_with_SW_align_node = 0
        number_of_hits_missing_smithWatermanAlignment_node = 0
        #number_of_hits_kept_for_statistical_analysis = 0  # number_of_hits
        protein_name = df.loc[acc, 'protein_name']
        ft_xml_path = df.loc[acc, 'SIMAP_feature_table_XML_file_path']
        homol_xml_path = df.loc[acc, 'SIMAP_homologues_XML_file_path']
        SIMAP_tar = df.loc[acc, 'SIMAP_tar']
        ft_xml_filename = os.path.basename(ft_xml_path)
        homol_xml_filename = os.path.basename(homol_xml_path)
        logging.info('%s' % protein_name)
        #check which files exist
        feature_table_XML_exists, homologues_XML_exists, SIMAP_tarfile_exists = utils.check_tarfile(SIMAP_tar, ft_xml_path, homol_xml_path)
        #check if the feature table and homologue files actually exist
        if os.path.isfile(df.loc[acc, 'SIMAP_tar']):
            SIMAP_tarfile_exists = True
        else:
            SIMAP_tarfile_exists = False
            #at the moment, we'll only create the tarfile if both the feature table and
        #homologue XML files downloaded successfully, but this might change depending on preference
        if SIMAP_tarfile_exists:
            try:
                with tarfile.open(SIMAP_tar, mode='r:gz') as tar:
                    if ft_xml_filename in [tarinfo.name for tarinfo in tar]:
                        feature_table_in_tarfile = True
                    else:
                       feature_table_in_tarfile = False
                    if homol_xml_filename in [tarinfo.name for tarinfo in tar]:
                        homologues_XML_in_tarfile = True
                    else:
                       homologues_XML_in_tarfile = False
            except (EOFError, tarfile.ReadError):
                #file may be corrupted, if script stopped unexpectedly before compression was finished
                logging.info('%s seems to be corrupted.'
                             'File will be deleted, and will need to be re-downloaded next time program is run' %
                             df.loc[acc, 'SIMAP_tar'])
                SIMAP_tarfile_exists = False
                feature_table_in_tarfile = False
                homologues_XML_in_tarfile = False
                os.remove(df.loc[acc, 'SIMAP_tar'])
        else:
            feature_table_in_tarfile = False
            homologues_XML_in_tarfile = False

        # if only one of the files is in the tarfile, it is corrupt and should be deleted:
        if True in [feature_table_in_tarfile, homologues_XML_in_tarfile] and False in [feature_table_in_tarfile, homologues_XML_in_tarfile]:
            os.remove(df.loc[acc, 'SIMAP_tar'])
            logging.warning("Tarfile containing homologues may be corrupt, contains only one of two XML files. "
                            "File will be deleted:\n{}".format(df.loc[acc, 'SIMAP_tar']))

        if all([feature_table_in_tarfile, homologues_XML_in_tarfile]):
            '''get the Phobius and TMHMM predictions from the feature table of the query sequence
            NOT USED, PHOBIUS PRED OFTEN MISSING, in the future the TMD region taken from uniprot record
            '''
            # create subfolders, if they don't exist
            subfolder = os.path.dirname(df.loc[acc, 'homol_csv_zip'])
            utils.make_sure_path_exists(subfolder)

            #if the setting is "False", and you don't want to overwrite the files, skip this section
            if set_["overwrite_homologue_csv_files"]:
                create_homol_csv = True
            else:
                #check if output file already exists
                if os.path.isfile(df.loc[acc, 'homol_csv_zip']):
                    try:
                        dfs_test = utils.open_df_from_csv_zip(df.loc[acc, 'homol_csv_zip'])
                        description_of_first_hit = dfs_test.loc[1, 'uniprot_description']
                        logging.info('%s homologues already converted to csv. (%s)' % (acc, description_of_first_hit))
                        create_homol_csv = False
                        # with tarfile.open(df.loc[acc, 'SIMAP_csv_from_XML_tarfile'], 'r:gz') as tar:
                        #     #create a list of files
                        #     files_in_output_tarball = [t.name for t in tar]
                        #     #check that the analysed files are actually there
                        #     if df.loc[acc, 'SIMAP_csv_from_XML'] in files_in_output_tarball:
                        #         #read output csv in tarfile
                        #         dfs = pd.read_csv(tar.extractfile(df.loc[acc, 'SIMAP_csv_from_XML']), index_col = 0)
                        #         description_of_first_hit = dfs.loc[1, 'uniprot_description']
                        #         logging.info('%s homologues already converted to csv. (%s)' % (acc, description_of_first_hit))
                        #         #filter to include only desired hits
                        #         '''OLD STUFF, from when XML to CSV was not saved separately
                        #         dfs_filt = dfs.query(
                        #             'gapped_ident_above_cutoff == True and hit_contains_SW_node == True and cr_disallowed_words_not_in_descr == True')
                        #         #avoid a divide by zero error in the unlikely case that there are no_identical_residues_in_alignment
                        #         dfs_filt_AAIMON = dfs_filt.loc[dfs_filt['nonTMD_perc_ident'] != 0]
                        #         list_of_TMDs = ast.literal_eval(df.loc[acc, 'list_of_TMDs'])
                        #         for TMD in list_of_TMDs:
                        #             #following the general filters, filter to only analyse sequences with TMD identity above cutoff,
                        #             #and a nonTMD_perc_ident above zero ,to avoid a divide by zero error
                        #             dfs_filt_AAIMON = dfs_filt_AAIMON.loc[
                        #                 dfs['%s_perc_ident'%TMD] >= set_['min_identity_of_TMD_initial_filter']]
                        #             #add to original dataframe with the list of sequences
                        #             df.loc[acc, '%s_AAIMON_ratio_mean'%TMD] = dfs_filt_AAIMON['%s_AAIMON_ratio'%TMD].mean()
                        #             df.loc[acc, '%s_AAIMON_ratio_std'%TMD] = dfs_filt_AAIMON['%s_AAIMON_ratio'%TMD].std()
                        #             logging.info('AIMON MEAN %s: %0.2f' % (TMD, df.loc[acc, '%s_AAIMON_ratio_mean'%TMD]))
                        #         '''
                        #     else:
                        #         logging.info('%s not in output file tarball, tarball will be deleted' % df.loc[acc, 'SIMAP_csv_from_XML'])
                        #         os.remove(df.loc[acc, 'SIMAP_csv_from_XML_tarfile'])
                        #         create_homol_csv = True
                        # logging.info('%s already converted to csv, moving to next sequence' %
                        #             df.loc[acc, 'SIMAP_csv_from_XML'])
                        # create_homol_csv = False
                    except (EOFError, KeyError):
                        #file may be corrupted, if script stopped unexpectedly before compression was finished
                        logging.info('%s seems to be corrupted. File will be deleted.' % df.loc[acc, 'homol_csv_zip'])
                        os.remove(df.loc[acc, 'homol_csv_zip'])
                        create_homol_csv = True
                else:
                    logging.info('%s not found, create_homol_csv = True' % df.loc[acc, 'homol_csv_zip'])
                    create_homol_csv = True
            #if the files don't exist, or you want to overwrite them
            if create_homol_csv:
                #extract the tarfile so that it can be read as xml
                with tarfile.open(df.loc[acc, 'SIMAP_tar'], 'r:gz') as tar:
                    SIMAP_homologues_XML_file_extracted = tar.extractfile(homol_xml_filename)

                    #parse the XML file with elementtree, define the 'root' of the XML file
                    simap_homologue_tree = ET.parse(SIMAP_homologues_XML_file_extracted)
                    simap_homologue_root = simap_homologue_tree.getroot()

                    #print the simap search details (database, e-value cutoff etc)
                    #dict_with_query_data = print_query_details_from_homologue_XML(simap_homologue_root, dict_with_query_data)
                    for parameters in simap_homologue_root[0][0][0][0].iter('parameters'):
                        df.loc[acc, 'SIMAP_input_seq_details_dict'] = str(parameters[0][0].attrib)
                        for SIMAP_filter in parameters.iter('filter'):
                            SIMAP_filter_string = SIMAP_filter.text
                        df.loc[acc, 'SIMAP_filter_string'] = str(SIMAP_filter_string)
                        for resultSpecification in parameters.iter('resultSpecification'):
                            SIMAP_resultSpecification_dict = resultSpecification.attrib
                        df.loc[acc, 'SIMAP_resultSpecification_dict'] = '"%s"' % SIMAP_resultSpecification_dict
                        for databases in parameters.iter('databases'):
                            database_details_dict = databases[0].attrib
                        df.loc[acc, 'database_details_dict'] = '"%s"' % database_details_dict
                        df.loc[acc, 'simap_version'] = simap_homologue_root[0][0][0][0][0].attrib['version']
                        df.loc[acc, 'SIMAP_total_hits'] = int(simap_homologue_root[0][0][0][1][0].attrib['total'])
                    if df.loc[acc, 'simap_version'] != '4.0':
                        logging.warning('WARNING! Your XML file is simap version %s,'
                                        'however this SIMAP parser was developed for SIMAP version 4.0.' %
                                         df.loc[acc, 'simap_version'])
                    counter_XML_to_CSV += 1

                    query_sequence_node = simap_homologue_root[0][0][0][0][2][0][0]
                    ''' xxxx CURRENTLY THE df is filled with nan values,
                        but that doesn't make sense as the script seems to work
                    '''
                    df.loc[acc, 'query_md5'] = query_sequence_node.attrib['md5']
                    df.loc[acc, 'query_length'] = int(query_sequence_node.attrib['length'])
                    df.loc[acc, 'query_selfscore'] = query_sequence_node.attrib['selfscore']
                    df.loc[acc, 'query_sequenceid'] = query_sequence_node.attrib['sequenceid']
                    df.loc[acc, 'total_number_of_simap_hits'] = query_sequence_node[0].attrib['number_hits']
                    df.loc[acc, 'query_sequence_from_homologue_XML_file'] = query_sequence_node[0][0].text
                    df.loc[acc, 'number_of_hits_in_homologue_XML_file'] = int(
                        simap_homologue_root[0][0][0][1][0].attrib['total'])
                    '''
                    Create an updated csv_file_with_uniprot_data to include the data from SIMAP regarding the query
                    '''
                    #save current df to csv
                    with open(pathdict["list_summary_csv"], 'w') as f:
                        df.to_csv(f, sep=",", quoting=csv.QUOTE_NONNUMERIC)

                    # #create new files to store the fasta sequences, and save the query sequence (the row here is the uniprot number in th df index)
                    # for row in df.index:
                    #     #for a protein with TMDs, the list of TMD names should be ['TM01','TM02']
                    #     list_of_TMDs = ast.literal_eval(df.loc[row, 'list_of_TMDs'])

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
                        #load the amino acid substitution matrices from the settings file
                        list_of_aa_sub_matrices = set_['aa_sub_matrices']

                        #import the amino acid substitution matrices
                        utils.import_amino_acid_substitution_matrices()

                        #add the similarity ratios to the csv_header_for_SIMAP_homologue_file.
                        # These will depend on the individual settings
                        #                    if set_['["mp_calculate_TMD_conservation_with_aa_matrices']:
                        #                        for j in range(set_["gap_open_penalty_min"],
                        #                                       set_["gap_open_penalty_max"],
                        #                                       set_["gap_open_penalty_increment"]):
                        #                            gap_open_penalty = j
                        #                            gap_extension_penalty = j
                        #                            for matrix_name in list_of_aa_sub_matrices:
                        #                                column_name = 'sim_ratio_%s_gapo%i' % (matrix_name.replace("'", "")[0:-7], j)
                        #                                csv_header_for_SIMAP_homologue_file.append(column_name)

                        #import the necessary matrices
                        #for matrix_name in list_of_aa_sub_matrices:
                        #matrix = matrix_name[0:-7]
                        #print (matrix)
                        #from Bio.SubsMat.MatrixInfo import matrix as matrix_name

                        SIMAP_temp_csv_from_XML_path = df.loc[acc,'homol_csv_zip'][:-4]
                        #fasta_file_path = df.loc[acc, 'fasta_file_path']

                        #create an empty file
                        open(SIMAP_temp_csv_from_XML_path, 'w').close()

                        #reopen to add match details iteratively from dictionary
                        with open(SIMAP_temp_csv_from_XML_path, 'a') as csvfile:

                            #set up a bool to catch those files where not a single hit actually gives data
                            at_least_one_hit_contains_SW_node = False

                            for hit in simap_homologue_hits:
                                match_details_dict = {}

                                #add desired hit information to the dictionary for transfer to csv
                                A1_hit_number = int(hit.attrib['number'])
                                match_details_dict['A1_hit_number'] = A1_hit_number
                                match_details_dict['A3_md5'] = hit[1].attrib['md5']

                                #define the major nodes in the XML-file
                                try:
                                    protein_node = hit[1][1]
                                    hit_contains_protein_node = True
                                except IndexError:
                                    hit_contains_protein_node = False
                                    number_of_hits_missing_protein_node += 1
                                    logging.warning('%s hit %s contains no protein node' % (protein_name,
                                                                                            match_details_dict['A3_md5'])
                                                                                            )
                                if hit_contains_protein_node:
                                    try:
                                        smithWatermanAlignment_node = hit[0][0][14]
                                        hit_contains_SW_node = True
                                        num_hits_with_SW_align_node += 1
                                    except IndexError:
                                        hit_contains_SW_node = False
                                    match_details_dict['hit_contains_SW_node'] = hit_contains_SW_node
                                    #add the description. Make an empty string if it is the first (query) hit, preventing the presence of np.nan in the later dataframe
                                    if A1_hit_number == 1:
                                        uniprot_description = '%s_SIMAP_query_sequence' % protein_name
                                    else:
                                        uniprot_description = protein_node.attrib['description']
                                    match_details_dict['uniprot_description'] = uniprot_description
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
                                        #sequence is probably synthetic, as it has no database node
                                        match_details_dict['uniprot_description'] += ', synthetic'
                                        match_details_dict['organism'] = 'synthetic'
                                        match_details_dict['taxonomy_node_id'] = 'synthetic'
                                        match_details_dict['taxonomy_rank'] = 'synthetic'
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
                                    #some parameters that are needed for identity calculations later
                                    #FASTA_num_ident_res = FASTA_identity / 100.0 * FASTA_overlap
                                    #check if the TMD is in the smith waterman alignment. Note that start and stop in the alignment node is based on FASTA alignment,
                                    #which is not shown. Occasionally, this will not match the SW alignment.
                                    #xxx it might be better to insert a function that determines of the TMD is after the start of the SW alignment
                                    #is_start_of_TMD_in_FASTA = True if FASTA_query_start <= TMDstart else False
                                    #is_end_of_TMD_in_FASTA = True if TMDend <= FASTA_query_end else False
                                    #is_TMD_in_FASTA_alignment = True if all(
                                    #    [is_start_of_TMD_in_FASTA, is_end_of_TMD_in_FASTA]) else False
                                    #mmmmm TEMP
                                    #logging.info('gapped_identity_too_low: %s' % gapped_identity_too_low)
                                    '''***********************if the TMD region is actually covered by the hsp, then conduct some further analyses of the match TMD region*************************'''
                                    if hit_contains_SW_node:
                                        #check that at least one hit gives data
                                        at_least_one_hit_contains_SW_node = True
                                        query_alignment_sequence = ''
                                        '''For the moment, there is no need to put the whole match hsp sequence into the csv file'''
                                        #for smithWatermanAlignment in alignment_node.iter('smithWatermanAlignment'):
                                        match_details_dict['SW_query_score_ratio'] = smithWatermanAlignment_node[0].text
                                        match_details_dict['SW_match_score_ratio'] = smithWatermanAlignment_node[1].text
                                        match_details_dict['SW_query_coverage'] = smithWatermanAlignment_node[2].text
                                        match_details_dict['SW_match_coverage'] = smithWatermanAlignment_node[3].text
                                        match_details_dict['SW_coverage_ratio'] = smithWatermanAlignment_node[4].text
                                        match_details_dict['A5_alignment_pretty'] = smithWatermanAlignment_node[8].text
                                        match_details_dict['SW_alignment_seq1offset'] = int(smithWatermanAlignment_node.attrib['alignment-seq1offset'])
                                        match_details_dict['SW_alignment_seq2offset'] = int(smithWatermanAlignment_node.attrib['alignment-seq2offset'])
                                        match_details_dict['SW_identity'] = float(smithWatermanAlignment_node.attrib['identity'])
                                        match_details_dict['SW_similarity'] = float(smithWatermanAlignment_node.attrib['similarity'])
                                        #Get the full sequences. Note that they greatly increase the size of the csv file.
                                        match_details_dict['query_alignment_sequence'] = smithWatermanAlignment_node[5].text
                                        match_details_dict['alignment_markup'] = smithWatermanAlignment_node[6].text
                                        match_details_dict['match_alignment_sequence'] = smithWatermanAlignment_node[7].text
                                        #create a list of TMD names to be used in the loops below (TM01, TM02 etc)
                                        # list_of_TMDs = ast.literal_eval(df.loc[acc, 'list_of_TMDs'])
                                        # #run the search using the regular expression that will find the TMD even if it contains gaps
                                        # for TMD in list_of_TMDs:
                                        #     #if is_TMD_in_FASTA_alignment:
                                        #     query_TMD_sequence = df.loc[acc, '%s_seq'%TMD]
                                    else:
                                        number_of_hits_missing_smithWatermanAlignment_node += 1
                                    if A1_hit_number == 1:
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
                        # either create new zip and add ("w"), or open existing zip and add "a"
                        with zipfile.ZipFile(df.loc[acc,'homol_csv_zip'], mode="w", compression=zipfile.ZIP_DEFLATED) as zipout:
                            zipout.write(SIMAP_temp_csv_from_XML_path, arcname=os.path.basename(df.loc[acc,'homol_csv_zip'])[:-4])
                        # delete temporary csv file
                        os.remove(df.loc[acc,'homol_csv_zip'][:-4])
                        # with tarfile.open(df.loc[acc, 'SIMAP_csv_from_XML_tarfile'], 'w:gz') as tar_SIMAP_out:
                        #     tar_SIMAP_out.add(SIMAP_temp_csv_from_XML_path, arcname=df.loc[acc, 'SIMAP_csv_from_XML'])
                        # os.remove(SIMAP_temp_csv_from_XML_path)
                        logging.info('%s homologous sequences parsed from SIMAP XML to csv' % int(df.loc[acc, 'SIMAP_total_hits']))
    logging.info('number_of_hits_missing_smithWatermanAlignment_node: %i' % number_of_hits_missing_smithWatermanAlignment_node)
    logging.info('number_of_hits_missing_protein_node: %i' % number_of_hits_missing_protein_node)
    logging.info('****parse_SIMAP_to_csv finished!!****\n%g files parsed from SIMAP XML to csv' % counter_XML_to_CSV)