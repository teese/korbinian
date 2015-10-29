
import pandas as pd
import numpy as np
import re as re

# Function to store Dataframes in an Excelfile; Converting lists etc. into strings
def df_to_excel (dataframe,path):
    temp_df = pd.DataFrame(dataframe, dtype = str)
    temp_df.to_excel(path)
  
# Creates subfolders, based on first to letters of files; Distributes files to these subfolders
def distribute_files_to_subfolders (path):
    list_of_filenames = os.listdir(path)
    for n in list(set(n[0:2] for n in list_of_filenames)):
        os.makedirs(path+"/%s" %n)    
    for n in list_of_filenames:
        directory = path+"/"+n
        new_directory = path+"/"+n[0:2]+"/"+n 
        os.replace(directory,new_directory)
   
# Moves files from a subfolder to the root folder; 
# Length of subfoldername = amount of letters; necessary for path-identification         
def move_files_from_subfolder_to_folder(path_of_subfolder,length_of_subfoldername):
    for n in os.listdir(path_of_subfolder):
        directory = path_of_subfolder+"/"+n
        new_directory = str(path_of_subfolder)[:-length_of_subfoldername]+"/"+n
        os.replace(directory,new_directory)

# this functions works exclusivly with dataframes; query, start, stop, and new_name refer to columns
# and as well, it does not work yet, still working on it
def slicing (columname_of_sequence, start, stop, columnname_for_spliced_sequence):
    for n in df["%s"%columname_of_sequence]:
        df["%s"%columnname_for_spliced_sequence] = n[start,stop]
        
        
def getting_list_of_indices_of_M_in_a_string(string):
    ind_list = [i for i, element in enumerate(string) if element == "M"]  #find(Topo_data)
    return ind_list

def getting_list_of_gapindices_of_string(string):
    gap_list = [i for i, element in enumerate(string) if element == "-"]  #find(Topo_data)
    return gap_list

def setup_error_logging(logfile):
    import os
    import logging
    import logging.config
    import sys
    import json
    import platform
    
    #you can either adjust the log settings from an external file, or paste them in teh script 
    #external_log_settings_file = r'E:\Stephis\Projects\Programming\Python\files\learning\json\logging_settings.json'
    
    #load the log settings in json format, so it is easy to modify
    logsettings = json.dumps({
        "handlers": {
            "console": {
                "formatter": "brief", 
                "class": "logging.StreamHandler", 
                "stream": "ext://sys.stdout", 
                "level": "DEBUG"
            }, 
            "file": {
                "maxBytes": 5000000, 
                "formatter": "precise", 
                "backupCount": 3, 
                "class": "logging.handlers.RotatingFileHandler", 
                "filename": "logfile.txt"
            }
        }, 
        "loggers": {
            "simpleExample": {
                "handlers": [
                    "console", 
                    "file"
                ], 
                "propagate": "no", 
                "level": "DEBUG"
            }
        }, 
        "version": 1, 
        "root": {
            "handlers": [
                "console", 
                "file"
            ], 
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
                "format": "%(levelname)-8s: %(name)-15s: %(message)s"
            }
        }
    }, skipkeys=True, sort_keys=True, indent=4, separators=(',', ': '))
    
    #modify the location of the log file to suit your needs
    config=json.loads(logsettings)
    config['handlers']['file']['filename'] = logfile
    #logging.info(logsettings)
    
    #save the logging settings in an external file (if desired)
    #with open(settings_file_output, 'w') as f:
    #    f.write(json.dumps(config, f, indent=4, sort_keys=True))
    
    #create a blank logging file (if desired)
    f = open(logfile, 'w')
    f.close()
    
    #clear any previous logging handlers that might have been previously run in the console
    logging.getLogger('').handlers = [] 
    #load the logging settings from the modified json string
    logging.config.dictConfig(config)
    
    #write system settings to logfile
    logging.warning('LOGGING LEVEL')
    logging.critical('Example of critical-level error. Current logging settings are level %s. At level DEBUG this logfile should also contain examples of WARNING and INFO level reports.' % config['handlers']['console']['level'])
    logging.warning('Example of warning-level error')
    logging.info('Example of info-level error\n')
    logging.info('SYSTEM INFORMATION')
    logging.info('system description: %s' % str(platform.uname()))
    logging.info('system       : %s', str(platform.system()))
    logging.info('architecture : %s', str(platform.architecture()))
    logging.info('network_name : %s', str(platform.node()))
    logging.info('release      : %s', str(platform.release()))
    logging.info('version      : %s', str(platform.version()))
    logging.info('machine      : %s', str(platform.machine()))
    logging.info('processor    : %s', str(platform.processor()))
    logging.info('python_version: %s', str(platform.python_version()))
    logging.info('python_build : %s', str(platform.python_build()))
    logging.info('python_compiler: %s\n', str(platform.python_compiler()))
    logging.info('PATH INFORMATION')
    logging.info("argv: %r"%(sys.argv,))
    logging.info("dirname(argv[0]): %s" % os.path.abspath(os.path.expanduser(os.path.dirname(sys.argv[0]))))
    logging.info("pwd: %s\n" % os.path.abspath(os.path.expanduser(os.path.curdir)))
    
    #save the logging settings in the logfile    
    logging.info('LOGGING SETTINGS FOR THIS RUN IN JSON FORMAT')    
    logging.info("%s\n" % config)
    
    #test error message reporting
    logging.warning('LOGGING TEST:')
    try:
        open('/path/to/does/not/exist', 'rb')
    except (SystemExit, KeyboardInterrupt):
        raise
    except Exception:
        logging.error('Failed to open file', exc_info=True)
    logging.warning('LOGGING SETUP IS SUCCESSFUL\n\nLOG INFORMATION STARTS HERE:\n')

def getting_list_of_gapindices_of_string(string):
    gap_list = [i for i, element in enumerate(string) if element == "-"]  #find(Topo_data)
    return gap_list
	
	
def create_border_list (index_list):
    m_borders = []
    m_borders.append(index_list[0])
    for n in range (0,len(index_list)-1):
        if index_list[n]+1 != index_list[n+1]:
            m_borders.append(index_list[n]+1)
            m_borders.append(index_list[n+1])
    m_borders.append(index_list[-1]+1)
    return m_borders
	
def create_list_of_TMDs (amount_of_TMDs):
    list_of_TMDs = []
    for n in range (1,int(amount_of_TMDs)+1):
        list_of_TMDs.append("TM%.2d"%n)
    return list_of_TMDs
    
def isEven(number):
    return number % 2 == 0
    
def isOdd(number):
    return number %2 != 0 
    
def isNaN(num):
    return num != num
    
def sum_gaps(df_column):
    sum = 0
    for n in df_column: 
        if not isNaN(n):
            sum = sum + n        
    return sum
    
def frequency_of_tmd(int_of_tmd, column_containing_tmd_amount):
    frequency = 0
    for n in column_containing_tmd_amount:
        if int_of_tmd <= n:
            frequency = frequency+1
    return frequency
    
def create_regex_string_for_juxta(inputseq):
    ''' adds '-*' between each aa or nt/aa in a DNA or protein sequence, so that a particular
    aligned sequence can be identified via a regex search, even if it contains gaps
    inputseq : 'LQQLWNA'
    output   : 'L-*Q-*Q-*L-*W-*N-*A'
    '''
    search_string = ''
    for letter in inputseq:
        letter_with_underscore = letter + '-*'
        search_string += letter_with_underscore
    return "-*" + search_string
    


#müsste passen
def get_end_juxta_before_TMD (x,input_TMD):
    TM_int = int(input_TMD[2:])
    if input_TMD == "TM01":
        x['end_juxta_before_%s_in_query'%input_TMD] = np.where(x['%s_start_in_SW_alignment'%input_TMD]==0,0,x['%s_start_in_SW_alignment'%input_TMD]-1)
    else:
        x["end_juxta_before_%s_in_query"%input_TMD]=x["%s_start_in_SW_alignment"%input_TMD]-1

    
def get_end_juxta_after_TMD (x,input_TMD):
    TM_int = int(input_TMD[2:])
    last_TMD = list_of_tmds[-1]
    if input_TMD == last_TMD:
         x["end_juxta_after_%s"%input_TMD] = np.where((x["%s_end_in_SW_alignment"]+30)<x["len_query_aligment_sequence"],x["%s_end_in_SW_alignment"]+30,x["len_query_aligment_sequence"])
    else:    
        x["end_juxta_after_%s"%input_TMD]= dfs["%s_end_in_SW_alignment"%input_TMD]+((dfs["TM%.2d_start_in_SW_alignment"%(TM_int+1)]-dfs["%s_end_in_SW_alignment"%input_TMD])/2).apply(lambda x :int(x) if not np.isnan(x) else np.nan)
           
 
 
def get_start_and_end_of_TMD_in_query(x, TMD_for_regular_expression_search):
    '''
    define function to obtain regex output (start, stop, etc) as a tuple    
    the functions were originally in MAIN, as I am not sure how to apply **args and **kwargs in functions that apply to pandas
    note that TMD_for_regular_expression_search is therefore a global variable
    '''
    #print('x  ', x, 'TMD_for_regular_expression_search     ', TMD_for_regular_expression_search)
    m = re.search(TMD_for_regular_expression_search, x)
    if m:
        #if the tmd is in the query, return True, start, stop
        return (bool(m), m.start(), m.end())
    else:
        #if the tmd is not in the query, return False, NaN, NaN
        return (bool(m), np.nan, np.nan)    
    
def slice_juxta_before_TMD_in_query(x, TMD):
    return x['query_alignment_sequence'][int(x['start_juxta_before_%s' % TMD]):int(x['end_juxta_before_%s' % TMD])]

def slice_juxta_after_TMD_in_query(x,TMD):
    return x['query_alignment_sequence'][int(x['start_juxta_after_%s' % TMD]):int(x['end_juxta_after_%s' % TMD])]

def slice_juxta_before_TMD_in_match(x,TMD):
    return x['match_alignment_sequence'][int(x['start_juxta_before_%s' % TMD]):int(x['end_juxta_before_%s' % TMD])]

def slice_juxta_after_TMD_in_match(x,TMD):
    return x['match_alignment_sequence'][int(x['start_juxta_after_%s' % TMD]):int(x['end_juxta_after_%s' % TMD])]
    
def find_last_TMD():
    for n in range (1,24):
        if r_utils.isNaN(dfs['TM%.2d_start_in_SW_alignment'%n]):
            last_TMD = n