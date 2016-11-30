import pandas as pd
import numpy as np
import ast
import korbinian
import matplotlib.pyplot as plt
#%matplotlib inline
plt.rcParams["savefig.dpi"] = 300


#########################
#                       #
#   load summary file   #
#                       #
#########################

# specify used list of proteins
proteinlist = 97
max_num_homologues = 4000
min_num_homologues = 10
# load summary files
# df1 = pd.read_csv("/Volumes/Musik/Databases/summaries/{a}/List{a}_summary.csv".format(a=proteinlist), index_col=0)
df_summary = pd.read_csv("/Volumes/Musik/Databases/summaries/{a}/List{a}_cr_summary.csv".format(a=proteinlist),
                         index_col=0)

# for every acc in df_summary add uniprot_entry_name from df1 to df_summary
# data could be already existent for new processed datasets (starting from Nov 28 2016)
# for acc in df_summary.loc[df_summary['list_of_TMDs'].notnull()].loc[df_summary['list_of_TMDs'] != 'nan'].index:
#    df_summary.loc[acc, 'uniprot_entry_name'] = df1.loc[acc, 'uniprot_entry_name']

# filter summary file for min and max number of homologues
print('Dropped:')
for acc in df_summary.loc[df_summary['list_of_TMDs'].notnull()].loc[df_summary['list_of_TMDs'] != 'nan'].index:
    if df_summary.loc[acc, 'TM01_AAIMON_n_homol'] > max_num_homologues:
        df_summary = df_summary.drop([acc])
        print(acc, sep=' ', end='  ')
for acc in df_summary.loc[df_summary['list_of_TMDs'].notnull()].loc[df_summary['list_of_TMDs'] != 'nan'].index:
    if df_summary.loc[acc, 'TM01_AAIMON_n_homol'] < min_num_homologues:
        df_summary = df_summary.drop([acc])
        print(acc, sep=' ', end='  ')

# save relevant parts for navigation through file system (database) in dictionaries
dict_TMDs = {}
for acc in df_summary.loc[df_summary['list_of_TMDs'].notnull()].loc[df_summary['list_of_TMDs'] != 'nan'].index:
    dict_TMDs[acc] = df_summary.loc[acc, 'list_of_TMDs']

dict_folder = {}
for key in dict_TMDs.keys():
    dict_folder[key] = key[:2]

dict_uniprot_entry = {}
for acc in df_summary.loc[df_summary['list_of_TMDs'].notnull()].loc[df_summary['list_of_TMDs'] != 'nan'].index:
    dict_uniprot_entry[acc] = df_summary.loc[acc, 'uniprot_entry_name']

#########################
#                       #
#       load data       #
#                       #
#########################

# initiate empty numpy array
data = np.empty([0,3])
# navigate through filesystem and open pickles from .zip
for key in dict_folder.keys():
    in_zipfile = "/Volumes/Musik/Databases/homol/{a}/{b}_{c}_cr_ratios.zip".format(a=dict_folder[key], b=key, c=dict_uniprot_entry[key])
    #print ('current directory: {a}' .format(a=in_zipfile))
    for TMD in ast.literal_eval(dict_TMDs[key]):
        filename = "{a}_{c}_{b}_cr_df.pickle".format(a=key, b=TMD, c=dict_uniprot_entry[key])
        filename_csv = "{a}_{c}_AAIMON_normalisation_data.csv".format(a=key, c=dict_uniprot_entry[key])
        #print ('current file: {a}' .format(a=filename))
        # generate column names necessary for current file
        columns = ['FASTA_gapped_identity', '{a}_AAIMON_ratio'.format(a=TMD), '{a}_AAIMON_ratio_n'.format(a=TMD)]
        # open dataframe  with function from korbinian, extract required columns, convert to np array
        df = korbinian.utils.open_df_from_pickle_zip(in_zipfile, filename)[columns].as_matrix()
        # join output data file with currently opened dataframe
        data = np.concatenate((data, df))
# drop every row with nan
data = data[~np.isnan(data).any(axis=1)]
# create real percentage values
data[:,0] = data[:,0]*100


#########################
#                       #
#      create bins      #
#                       #
#########################


linspace_binlist = np.linspace(1,100,100)
binned_data = np.empty([0,3])
for percentage in linspace_binlist:
    bin_for_mean = np.empty([0,3])
    for line in data:
        if line[0] < percentage and line [0] > percentage - 1:
            bin_for_mean = np.concatenate((bin_for_mean, line.reshape(1,3)))
    mean_data_in_bin = np.array ([percentage, bin_for_mean[:,1].mean(), bin_for_mean[:,2].mean()])
    binned_data = np.concatenate ((mean_data_in_bin.reshape(1,3), binned_data))
binned_data = binned_data[~np.isnan(binned_data).any(axis=1)]


#########################
#                       #
#       plot data       #
#                       #
#########################

#backgroundcolour = '0.5'
plt.style.use('ggplot')
# set default font size for plot
fontsize = 12
datapointsize = 0.5
alpha = 0.05
linewidth = 1
color_nonnorm = "#454545"  # grey
color_norm = "#0076B8"     # TUM-Blue
fig, ax = plt.subplots()

# set color of axis label to black
ax.tick_params(axis='x', colors='black')
ax.tick_params(axis='y', colors='black')
ax.yaxis.label.set_color('black')
ax.xaxis.label.set_color('black')

# pylab.rcParams['figure.figsize'] = (50.0, 40.0)
x = data[:,0] # FASTA_gapped_identity
y = data[:,1] # AAIMON for each TMD
ax.scatter(x=x, y=y, color=color_nonnorm, alpha=alpha, s=datapointsize) # color="#003366" is TUM-blue
plt.ylim(ymin = 0, ymax = 2.0)
plt.xlim(xmin = 60, xmax = 100)
# label the x-axis for each plot, based on the TMD
ax.set_xlabel('% identity', fontsize=fontsize)
# move the x-axis label closer to the x-axis
#ax.xaxis.set_label_coords(0.45, -0.085)
ax.set_ylabel('AAIMON ratio', fontsize=fontsize)
# change axis font size
ax.tick_params(labelsize=fontsize)
x_line = binned_data[:,0]
y_line = binned_data[:,1]
plt.plot(x_line, y_line, linewidth=linewidth, color=color_nonnorm) # plot linegraph

# plot normalised data
x = data[:,0] # FASTA_gapped_identity
y = data[:,2] # AAIMON_n for each TMD
ax.scatter(x=x, y=y, color=color_norm, alpha=alpha, s=datapointsize) # color="#FF6633" is TUM-orange
x_line = binned_data[:,0]
y_line = binned_data[:,2]
plt.plot(x_line, y_line, linewidth=linewidth, color=color_norm) # plot linegraph

fig.savefig('FigX18-Test_85000_TMDs_AAIMON_vs_AAIMON_n.png')