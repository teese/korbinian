
import pandas as pd

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
