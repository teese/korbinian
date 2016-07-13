import os
from Bio import SeqIO
import pandas as pd

def parse_large_flatfile_with_list_uniprot_accessions(uniprot_folder, list_number, logging, uniprot_flatfile_of_selected_records):
    """ Retrieves UniProt flatfiles from a large flatfile (e.g. All UniProt), based on a list of accession numbers

    """
    input_accession_list = "NEED TO ADD VARIABLE SOMEWHERE"
    logging.info('~~~~~~~~~~~~  starting A01_parse_large_flatfile_with_list_uniprot_accessions   ~~~~~~~~~~~~')
    # parse_large_flatfile_with_list_uniprot_accessions(list_of_uniprot_accessions, uniprot_flatfile_all_single_pass, uniprot_flatfile_of_selected_records)
    # def parse_large_flatfile_with_list_uniprot_accessions(input_accession_list, input_uniprot_flatfile, output_uniprot_flatfile):
    # define path to large uniprot flatfile containing the protein records to be extracted
    input_uniprot_flatfile = os.path.join(uniprot_folder, "List%02d_large_uniprot_flatfile.txt" % list_number)
    output_uniprot_flatfile = uniprot_flatfile_of_selected_records
    # from Bio import SeqIO
    # create a list of all the uniprot accessions of the proteins to be selected from the larger uniprot file
    accession_list = [line.strip() for line in open(input_accession_list, "r")]
    uniprot_index_handle = SeqIO.index(input_uniprot_flatfile, "swiss")
    with open(output_uniprot_flatfile, "wb") as output:
        for acc in accession_list:
            try:
                # add the selected records to the file, but adds a new line after each line! Doesn't affect later conversion to SeqRecord object
                output.write(uniprot_index_handle.get_raw(acc))
            except KeyError:
                logging.info("No SwissProt record found in %s for %s." % (input_uniprot_flatfile, acc))

def retrieve_uniprot_data_for_acc_list_in_xlsx_file(excelfile_with_uniprot_accessions, logging, uniprot_flatfile_of_selected_records):
    uniprot_flatfile_all_human_membrane_compressed = "NEED TO ADD VARIABLE SOMEWHERE"
    logging.info('~~~~~~~~~~~~  starting A02_retrieve_uniprot_data_for_acc_list_in_xlsx_file   ~~~~~~~~~~~~')
    # take list of acc, search in default uniprot flatfile. If missing, download from uniprot server.
    input_uniprot_flatfile = uniprot_flatfile_all_human_membrane_compressed
    df_uniprot_accessions = pd.read_excel(excelfile_with_uniprot_accessions, sheetname='uniprot_numbers')
    # remove proteins that are marked as 'not included in analysis'
    df_uniprot_accessions = df_uniprot_accessions[df_uniprot_accessions['include_in_analysis'] == True]
    # accession_list = [line.strip() for line in open(input_accession_list, "r")]
    uniprot_index_handle = SeqIO.index(input_uniprot_flatfile, "swiss")
    with open(uniprot_flatfile_of_selected_records, "wb") as output:
        for uniprot_accession in df_uniprot_accessions['A1_uniprot_accession']:
            try:
                # this adds the selected records to the file, but adds a new line after each line!
                # Doesn't affect conversion to SeqRecord object)
                assert isinstance(uniprot_index_handle, object)
                output.write(uniprot_index_handle.get_raw(uniprot_accession))
            except KeyError:
                print("No SwissProt record found in %s for %s." % (input_uniprot_flatfile, uniprot_accession))