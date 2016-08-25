import csv
import os
import pandas as pd
import pickle
import tarfile
import xml.etree.ElementTree as ET
import zipfile

def parse_SIMAP_to_csv(pathdict, set_, logging):
    # @BO: you don't need to look at this stuff. skip to the next # @BO: comment :)
    counter_XML_to_CSV = 0
    logging.info('~~~~~~~~~~~~  starting parse_SIMAP_to_csv   ~~~~~~~~~~~~')
    # open dataframe with list of proteins
    df = pd.read_csv(pathdict["list_summary_csv"], sep = ",", quoting = csv.QUOTE_NONNUMERIC, index_col = 0)
    #iterate over the dataframe, excluding any proteins that do not have a list of TMDs
    for acc in df.loc[df['list_of_TMDs'].notnull()].loc[df['list_of_TMDs'] != 'nan'].index:
        protein_name = df.loc[acc, 'protein_name']
        homol_xml_path = df.loc[acc, 'SIMAP_homol_XML_path']
        homol_xml_filename = os.path.basename(homol_xml_path)
        logging.info('%s' % protein_name)
        #check if the feature table and homologue files actually exist
        if not os.path.isfile(df.loc[acc, 'SIMAP_tar']):
            # skip to next protein
            continue

        #extract the tarfile so that it can be read as xml
        with tarfile.open(df.loc[acc, 'SIMAP_tar'], 'r:gz') as tar:
            SIMAP_homologues_XML_file_extracted = tar.extractfile(homol_xml_filename)
            #parse_uniprot the XML file with elementtree, define the 'root' of the XML file
            simap_homologue_tree = ET.parse(SIMAP_homologues_XML_file_extracted)
            simap_homologue_root = simap_homologue_tree.getroot()

            SIMAP_orig_csv = df.loc[acc,'homol_df_orig_zip'][:-4]
            #create an empty file
            open(SIMAP_orig_csv, 'w').close()

            # @BO: add stuff one line at a time to CSV
            #reopen to add match details iteratively from dictionary
            with open(SIMAP_orig_csv, 'a') as csvfile:

                simap_homologue_hits = simap_homologue_root[0][0][0][1][0]

                for hit in simap_homologue_hits:
                    # @BO: necessary
                    hit_num = int(hit.attrib['number'])
                    # @BO: dict used to write to csv
                    match_details_dict = {}

                    #add desired hit information to the dictionary for transfer to csv
                    match_details_dict['hit_num'] = hit_num

                    protein_node = hit[1][1]

                    #add the description. Add a custom name if it is the first (query) hit
                    if hit_num == 1:
                        # @BO: OPTIONAL: Uniprot name used for first hit
                        description = '%s_SIMAP_query_sequence' % protein_name
                    else:
                        description = protein_node.attrib['description']
                    # @BO: necessary (protein description, or long protein name)
                    match_details_dict['description'] = description

                    try:
                        taxonomyNode = protein_node[2]
                        # @BO: necessary
                        match_details_dict['organism'] = taxonomyNode.attrib['name']
                    except IndexError:
                        #sequence has no database node
                        match_details_dict['organism'] = 'no_database_node'
                    #len_full_match_seq = len(full_match_seq)
                    alignment_node = hit[0][0]
                    # @BO: necessary (I don't use e-values, but journal reviewers might ask for it)
                    #E-value for hit
                    match_details_dict['FASTA_expectation'] = float(alignment_node[1].text)
                    # @BO: necessary.
                    #convert identity from e.g. 80 (80%) to 0.8
                    match_details_dict['FASTA_identity'] = float(alignment_node[3].text) / 100
                    # @BO: % identity excluding gaps is really useful, but I don't know if it is available in BLAST data
                    #strangely, I think gappedIdentity is the identity EXCLUDING gaps, which is a better value to base judgements on. convert identity from e.g. 80 (80%) to 0.8
                    match_details_dict['FASTA_gapped_identity'] = float(alignment_node[4].text) / 100

                    smithWatermanAlignment_node = hit[0][0][14]
                    # @BO: optional
                    match_details_dict['align_pretty'] = smithWatermanAlignment_node[8].text

                    #Get the full sequences. Note that they greatly increase the size of the csv file.
                    # @BO: necessary
                    match_details_dict['query_align_seq'] = smithWatermanAlignment_node[5].text
                    match_details_dict['align_markup_seq'] = smithWatermanAlignment_node[6].text
                    match_details_dict['match_align_seq'] = smithWatermanAlignment_node[7].text

                    # @BO: This is how I write the header to the CSV
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
                    # @BO: This is how I write the line in the csv containing the extracted data for that hit
                    writer.writerow(match_details_dict)

            # @BO: Here I re-open the csv, split it into two files (pretty.csv, and homol_df_orig.pickle), move to zip, and delete the non-zipped originals

            # open csv as a dataframe,
            df_homol = pd.read_csv(SIMAP_orig_csv, sep=",", quoting=csv.QUOTE_NONNUMERIC, index_col="hit_num")
            # restrict to just a few columns including the align_pretty that might be useful to check manually
            df_pretty = df_homol[["FASTA_gapped_identity", "organism", "description", "align_pretty"]]
            # save the align_pretty to csv
            df_pretty.to_csv(df.loc[acc,'SIMAP_align_pretty_csv'], sep=',', quoting=csv.QUOTE_NONNUMERIC)
            # drop the align_pretty column from the orig dataframe
            df_homol.drop('align_pretty', axis=1, inplace=True)
            # save the whole dataframe as a pickle for faster opening later
            with open(df.loc[acc,'homol_df_orig_pickle'], "wb") as p:
                pickle.dump(df_homol, p)
            # either create new zip and add ("w"), or open existing zip and add "a"
            with zipfile.ZipFile(df.loc[acc,'homol_df_orig_zip'], mode="w", compression=zipfile.ZIP_DEFLATED) as zipout:
                #zipout.write(SIMAP_orig_csv, arcname=os.path.basename(SIMAP_orig_csv))
                zipout.write(df.loc[acc,'SIMAP_align_pretty_csv'], arcname=os.path.basename(df.loc[acc,'SIMAP_align_pretty_csv']))
                zipout.write(df.loc[acc,'homol_df_orig_pickle'], arcname=os.path.basename(df.loc[acc,'homol_df_orig_pickle']))
            # delete temporary uncompressed files
            os.remove(SIMAP_orig_csv)
            os.remove(df.loc[acc,'SIMAP_align_pretty_csv'])
            os.remove(df.loc[acc,'homol_df_orig_pickle'])

    logging.info('****parse_SIMAP_to_csv finished!!****\n%g files parsed from SIMAP XML to csv' % counter_XML_to_CSV)

