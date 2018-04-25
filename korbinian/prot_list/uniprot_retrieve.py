import os
from Bio import SeqIO
import pandas as pd
import sys
# import debugging tools
from korbinian.utils import pr, pc, pn, aaa

def parse_large_flatfile_with_list_uniprot_accessions(s, input_accession_list, uniprot_dir, logging, selected_uniprot_records_flatfile):
    """Retrieves UniProt flatfiles from a large flatfile (e.g. All UniProt), based on a list of accession numbers.

    Parameters
    ----------
    input_accession_list : list
        List of accessions.
    uniprot_dir : str
        Path to UniProt folder, where selected_uniprot_records_flatfile will be saved
    list_number : int
        List number used to determine the output filename.
    logging : logging.Logger
        Logger for printing to console and logfile.
    selected_uniprot_records_flatfile : str
        Path to UniProt flatfile containing selected records for analysis. In this case, the output file.

    """
    """ Retrieves UniProt flatfiles from a large flatfile (e.g. All UniProt), based on a list of accession numbers

    """
    logging.info('~~~~~~~~~~~~  starting A01_parse_large_flatfile_with_list_uniprot_accessions   ~~~~~~~~~~~~')
    # parse_large_flatfile_with_list_uniprot_accessions(list_of_uniprot_accessions, uniprot_flatfile_all_single_pass, selected_uniprot_records_flatfile)
    # def parse_large_flatfile_with_list_uniprot_accessions(input_accession_list, input_uniprot_flatfile, output_uniprot_flatfile):
    # define path to large uniprot flatfile containing the protein records to be extracted
    input_uniprot_flatfile = os.path.join(uniprot_dir, "List%02d_large_uniprot_flatfile.txt" % s["list_number"])
    output_uniprot_flatfile = selected_uniprot_records_flatfile
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

def retrieve_uniprot_data_for_acc_list_in_xlsx_file(excelfile_with_uniprot_accessions, input_uniprot_flatfile, selected_uniprot_records_flatfile, logging):
    """ From a list of uniprot accessions in excel, select out desired records from a large UniProt flatfile.

    Parameters
    ----------
    excelfile_with_uniprot_accessions : str
        Path to excel input file.
    logging : logging.Logger
        Logger for printing to console and logfile.
    selected_uniprot_records_flatfile : str
        Path to output UniProt flatfile containing selected records for analysis.

    """
    logging.info('~~~~~~~~~~~~  starting retrieve_uniprot_data_for_acc_list_in_xlsx_file   ~~~~~~~~~~~~')
    # take list of acc, search in default uniprot flatfile. If missing, download from uniprot server.
    df_uniprot_accessions = pd.read_excel(excelfile_with_uniprot_accessions, sheetname='uniprot_numbers')
    # remove proteins that are marked as 'not included in analysis'
    df_uniprot_accessions = df_uniprot_accessions[df_uniprot_accessions['include_in_analysis'] == True]
    # accession_list = [line.strip() for line in open(input_accession_list, "r")]
    uniprot_index_handle = SeqIO.index(input_uniprot_flatfile, "swiss")
    with open(selected_uniprot_records_flatfile, "wb") as output:
        for uniprot_accession in df_uniprot_accessions['uniprot_acc']:
            try:
                # this adds the selected records to the file, but adds a new line after each line!
                # Doesn't affect conversion to SeqRecord object)
                assert isinstance(uniprot_index_handle, object)
                output.write(uniprot_index_handle.get_raw(uniprot_accession))
            except KeyError:
                sys.stdout.write("No SwissProt record found in %s for %s." % (input_uniprot_flatfile, uniprot_accession))