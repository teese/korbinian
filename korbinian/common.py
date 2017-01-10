from time import strftime
import json
import korbinian
import korbinian.utils as utils
import logging
import logging.config # not a duplicate, avoids logging error.
import os
import pandas as pd
import platform
import psutil
import signal
import sys

def create_settingsdict(excel_file_with_settings):
    """Creates a python dictionary from an excel file with various user settings.

    Parameters
    ----------
    excel_file_with_settings : str
        Path to excel file containing the settings for the korbinian module.

    Returns
    -------
    s : dict
        Dictionary derived from the settings file, with the "parameter" column as keys,
        and the "value" column as values. Notes are excluded.
    """
    sheetnames = ["run_settings", "file_locations", "variables", "figs_settings"]
    s = {}
    for sheetname in sheetnames:
        # open excel file as pandas dataframe
        dfset = pd.read_excel(excel_file_with_settings, sheetname=sheetname)
        if "email" in dfset.columns:
            # exclude row with notes, set parameter as index
            dfset_email = dfset[["parameter", "email"]].dropna().copy()
            dfset_email["parameter"] = dfset_email["parameter"] + "_email"
            dfset_email.set_index("parameter", inplace=True)
            # convert true-like strings to True, and false-like strings to False
            dfset_email.email = dfset_email.email.apply(utils.convert_truelike_to_bool, convert_nontrue=False)
            dfset_email.email = dfset_email.email.apply(utils.convert_falselike_to_bool)
            # convert to dictionary
            sheet_as_dict_email = dfset_email.to_dict()["email"]
            # join dictionaries together
            s.update(sheet_as_dict_email)

        # exclude row with notes, set parameter as index
        dfset = dfset[["parameter", "value"]].dropna()
        dfset.set_index("parameter", inplace=True)
        # convert true-like strings to True, and false-like strings to False
        dfset.value = dfset.value.apply(utils.convert_truelike_to_bool, convert_nontrue=False)
        dfset.value = dfset.value.apply(utils.convert_falselike_to_bool)
        # convert to dictionary
        sheet_as_dict = dfset.to_dict()["value"]
        # join dictionaries together
        s.update(sheet_as_dict)

    # normalise the paths for the data directory
    s["data_dir"] = os.path.normpath(s["data_dir"])

    return s

def setup_keyboard_interrupt_and_error_logging(s, list_number):
    """Sets up keyboard interrupt (in unavailable) and error logging.

    Parameters
    ----------
    s : dict
        Dictionary derived from the settings file, with the "parameter" column as keys,
        and the "value" column as values. Notes are excluded.
    list_number : int
        Number describing the list of proteins.
        Used to defines the input and output files, E.g. 8 for List08_selected_uniprot_records_flatfile.txt
    Returns
    -------
    logging : logging.Logger
        Logging object, for printing to console and logfile
    """

    def ctrlc(sig, frame):
        raise KeyboardInterrupt("CTRL-C!")
    signal.signal(signal.SIGINT, ctrlc)
    '''+++++++++++++++LOGGING++++++++++++++++++'''
    date_string = strftime("%Y%m%d_%H_%M_%S")

    # designate the output logfile, within a folder in tha data_dir called "logfiles"
    logfile = os.path.join(s["data_dir"],"logfiles",'List%02d_%s_logfile.log' % (list_number, date_string))

    level_console = s["logging_level_console"]
    level_logfile = s["logging_level_logfile"]

    logging = korbinian.common.setup_error_logging(logfile, level_console, level_logfile)
    return logging


def setup_error_logging(logfile, level_console="DEBUG", level_logfile="DEBUG"):
    """ Sets up error logging, and logs a number of system settings.

    Parameters
    ----------
    logfile : str
        Path to output logfile. If size exceeds limit set below in JSON settings, path.1, path.2 etc will be created.
    level_console : str
        Logging level for printing to console. DEBUG, WARNING or CRITICAL
    level_logfile : str
        Logging level for printing to logfile. DEBUG, WARNING or CRITICAL

    Returns
    -------
    logging : logging.Logger
        Logging object, for printing to console and logfile
    """
    # load the log settings in json format
    logsettings = json.dumps({
        "handlers": {
            "console": {
                "formatter": "brief",
                "class": "logging.StreamHandler",
                "stream": "ext://sys.stdout",
                "level": "DEBUG"
            },
            "file": {
                "maxBytes": 10000000,
                "formatter": "precise",
                "backupCount": 3,
                "class": "logging.handlers.RotatingFileHandler",
                "level": "DEBUG",
                "filename": "logfile.txt"
            }
        },
        "version": 1,
        "root": {
            "handlers": [
                "console",
                "file"
            ],
            "propagate": "no",
            "level": "DEBUG"
        },
        "formatters": {
            "simple": {
                "format": "format=%(asctime)s - %(name)s - %(levelname)s - %(message)s"
            },
            "precise": {
                "format": "%(asctime)s %(name)-15s %(levelname)-8s %(message)s"
            },
            "brief": {
                "format": "%(message)s"
            }
        }
    }, skipkeys=True, sort_keys=True, indent=4, separators=(',', ': '))

    config=json.loads(logsettings)
    # add user parameters to the logging settings (logfile, and logging levels)
    config['handlers']['file']['filename'] = logfile
    config['handlers']['console']['level'] = level_console
    config['handlers']['file']['level'] = level_logfile

    # create folder if necessary
    utils.make_sure_path_exists(logfile, isfile=True)
    #create a blank logging file
    with open(logfile, 'w') as f:
        pass

    #clear any previous logging handlers that might have been previously run in the console
    logging.getLogger('').handlers = []
    #load the logging settings from the modified json string
    logging.config.dictConfig(config)
    # collect a number of system settings that could be useful for troubleshooting
    system_settings_dict = {}
    system_settings_dict["system description"] = platform.uname()
    system_settings_dict["system"] = platform.system()
    system_settings_dict["architecture"] = platform.architecture()
    system_settings_dict["network_name"] = platform.node()
    system_settings_dict["release"] = platform.release()
    system_settings_dict["version"] = platform.version()
    system_settings_dict["machine"] = platform.machine()
    system_settings_dict["processor"] = platform.processor()
    system_settings_dict["python_version"] = platform.python_version()
    system_settings_dict["python_build"] = platform.python_build()
    system_settings_dict["python_compiler"] = platform.python_compiler()
    system_settings_dict["argv"] = sys.argv
    system_settings_dict["dirname(argv[0])"] = os.path.abspath(os.path.expanduser(os.path.dirname(sys.argv[0])))
    system_settings_dict["pwd"] = os.path.abspath(os.path.expanduser(os.path.curdir))
    system_settings_dict["total_ram"] = "{:0.2f} GB".format(psutil.virtual_memory()[0] / 1000000000)
    system_settings_dict["available_ram"] = "{:0.2f} GB ({}% used)".format(psutil.virtual_memory()[1] / 1000000000, psutil.virtual_memory()[2])
    # log the system settings
    logging.warning("system description : {}".format(system_settings_dict))
    #test error message reporting
    #logging.warning('LOGGING TEST:')
    #try:
    #    open('/path/to/does/not/exist', 'rb')
    #except (SystemExit, KeyboardInterrupt):
    #    raise
    #except Exception:
    #    logging.error('Failed to open file', exc_info=True)
    logging.warning('logging setup is successful (logging levels: console={}, logfile={}). \n'.format(level_console, level_logfile))
    return logging


def create_pathdict(base_filename_summaries):
    pathdict = {}
    pathdict["base_filename_summaries"] = base_filename_summaries
    # currently the protein list summary (each row is a protein, from uniprot etc) is a csv file
    pathdict["list_summary_csv"] = '%s_summary.csv' % base_filename_summaries
    # define path for csv with conservation ratios(mean AAIMON) (output from the gather_AAIMON_ratios function)
    pathdict["list_cr_summary_csv"] = '%s_cr_summary.csv' % base_filename_summaries
    # define path for csv with gap densities (from gather_gap_densities)
    pathdict["list_gap_summary_csv"] = '%s_gap_summary.csv' % base_filename_summaries
    # define path for pickle file with gap raw data
    pathdict["gap_data_pickle"] = '%s_gap_data.pickle' % base_filename_summaries
    # for SIMAP or BLAST downloads, gives a list of accessions that failed
    pathdict["failed_downloads_txt"] = '%s_failed_downloads.txt' % base_filename_summaries
    # for SIMAP or BLAST, a list of accessions not found in the homologue database
    pathdict["acc_not_in_homol_db_txt"] = '%s_acc_not_in_homol_db.txt' % base_filename_summaries
    # add the base path for the sub-sequences (SE01, SE02, etc) added by the user
    pathdict["list_user_subseqs_csv"] = '%s_user_subseqs.csv' % base_filename_summaries
    pathdict["list_user_subseqs_xlsx"] = '%s_user_subseqs.xlsx' % base_filename_summaries
    # create output png from gap analysis
    pathdict["gap_density_fig_path"] = '%s_create_graph_of_gap_density.png' % base_filename_summaries
    pathdict["gap_density_testfig_png"] = '%s_gap_test_out.png' % base_filename_summaries
    pathdict["gap_fastagap_all_pos_pickle"] = '%s_gap_all_pos.pickle' % base_filename_summaries
    pathdict["single_list_fig_path"] = '%s_figs' %base_filename_summaries
    pathdict["save_df_characterising_each_homol_TMD"] = '%s_characterising_each_homol_TMD.zip' % base_filename_summaries
    pathdict["list_keywords_csv"] = '%s_keywords.csv' % base_filename_summaries


    """OUTPUT FILES FOR OLD SCRIPTS"""
    pathdict["dfout08_simap_AAIMON"] = '%s_simap_AAIMON.csv' % base_filename_summaries
    pathdict["dfout09_simap_AAIMON_02"] = '%s_simap_AAIMON_02.csv' % base_filename_summaries
    pathdict["dfout10_uniprot_gaps"] = '%s_gap_densities.csv' % base_filename_summaries
    return pathdict
