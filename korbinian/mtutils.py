#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Utilities file containing useful functions.
More recent functions are at the top.
Authors: Mark Teese, Rimma Jenske
Created on Fri Nov  8 15:45:06 2013
"""
from Bio import SeqIO
import csv
import ctypes
import logging
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import platform
import re
import subprocess, threading, time, sys
import tarfile

'''
************************************************************The uniprot functions start here.***************************************************************
'''

def aaa(df_or_series):
    """ Function for use in debugging.
    Saves pandas Series or Dataframes to a user-defined csv file.
    """
     # convert any series to dataframe
    if isinstance(df_or_series, pd.Series):
        df_or_series = df_or_series.to_frame()
    csv_out = r"D:\data\000_aaa_temp_df_out.csv"
    df_or_series.to_csv(csv_out, sep=",", quoting=csv.QUOTE_NONNUMERIC)


# Fibonacci numbers module. Use this to test that the utilities are working
def fib(n):    # write Fibonacci series up to n
    a, b = 0, 1
    while b < n:
        print(b),
        a, b = b, a+b

def import_amino_acid_substitution_matrices():
    """
    imports several aa sub matrices from Bio.SubsMat.MatrixInfo
    """
    from Bio.SubsMat.MatrixInfo import blosum30 as blosum30_matrix
    from Bio.SubsMat.MatrixInfo import blosum62 as blosum62_matrix
    from Bio.SubsMat.MatrixInfo import blosum95 as blosum95_matrix
    from Bio.SubsMat.MatrixInfo import ident as ident_matrix
    from Bio.SubsMat.MatrixInfo import gonnet as gonnet_matrix
    from Bio.SubsMat.MatrixInfo import pam250 as pam250_matrix
    from Bio.SubsMat.MatrixInfo import pam120 as pam120_matrix
    from Bio.SubsMat.MatrixInfo import pam30 as pam30_matrix
    from Bio.SubsMat.MatrixInfo import levin as levin_matrix

def all_df_in_list_contain_data(df_list_KW, title = '', KW = '', data_names_list = []):
    '''Function that takes a list of dataframes, and checks them all to make sure none are empty.
    Useful when slicing dataframes based on a list of keywords.
    Input: df_list_KW (list of dataframes)
    Input2: title of figure, KW (keyword to be examined), data_names_list
    Output: boolean object "both_df_contain_data"
    '''
    #simply create a list that numbers the dataframes, if they are not explicitly named
    if data_names_list == []:
        data_names_list = ['%i' % i for i in range(len(df_list_KW))]
    
    #first assume that they contain data
    both_df_contain_data = True
    #if any of them don't contain data, return false
    for n, dfK in enumerate(df_list_KW):
        if dfK.empty:
            logging.info('dataframe is empty for %s, %s, %s. Graph will not be shown' % (title,  KW, data_names_list[n]))
            both_df_contain_data = False
    return both_df_contain_data

def create_df_with_mean_AAIMON_each_TM(df):
    '''Takes a dataframe containing a list of proteins, and the average AAIMON ratio
    calculated for each TMD.
    Returns a dataframe with the average AAIMON for each TMD (TM01, TM02, etc) in the dataset.
    Returns the max_num_TMDs from all the proteins examined.
    Returns a list of the TMDs that can be used for iteration, etc.
    '''
    #find max number of TMDs in the whole dateset
    max_num_TMDs = df.number_of_TMDs.max()
    #create empty dictionaries to hold the mean values etc
    nested_dict_hist_data_AAIMON_each_TM = {}
    nested_dict_mean_AAIMON_each_TM = {}
    #create list to use as the legend in the later figures
    full_list_TMDs = []
    
    for i in range(1, max_num_TMDs + 1):
        TM = 'TM%02d' % i
        full_list_TMDs.append(TM)
        #create new series with the data (each datapoint is the mean for all homologues, 
        #for a single protein)
        hist_data_AAIMON_each_TM = df['TM%02d_AAIMON_ratio_mean' % i].dropna()
        #add data to nested dict
        nested_dict_hist_data_AAIMON_each_TM[TM] = hist_data_AAIMON_each_TM
        #calculate mean, std etc
        dict_mean_AAIMON_each_TM = {}
        dict_mean_AAIMON_each_TM['mean'] = hist_data_AAIMON_each_TM.mean()
        dict_mean_AAIMON_each_TM['std'] = hist_data_AAIMON_each_TM.std()
        dict_mean_AAIMON_each_TM['median'] = np.percentile(hist_data_AAIMON_each_TM, 50)
        dict_mean_AAIMON_each_TM['5_percentile'] = np.percentile(hist_data_AAIMON_each_TM, 5)
        dict_mean_AAIMON_each_TM['95_percentile'] = np.percentile(hist_data_AAIMON_each_TM, 50)
        #add mean etc to nested dict
        nested_dict_mean_AAIMON_each_TM[TM] = dict_mean_AAIMON_each_TM
    #convert nested dict to dataframe
    df_mean_AAIMON_each_TM = df.from_dict(nested_dict_hist_data_AAIMON_each_TM)
    
    #create new column to hold the last TMD
    df['last_TM_AAIMON_ratio_mean'] = np.nan
    #obtain data from last TMD for all proteins
    for acc in df.index:
        df.loc[acc,'last_TM_AAIMON_ratio_mean'] = df.loc[acc,'TM%02d_AAIMON_ratio_mean' % df.loc[acc,'number_of_TMDs']]
    AAIMON_last_TM = df['last_TM_AAIMON_ratio_mean'].dropna()
    #add the data for the last TMD to the dataframe
    df_mean_AAIMON_each_TM['last_TM_AAIMON_ratio_mean'] = AAIMON_last_TM
    
    return df_mean_AAIMON_each_TM, max_num_TMDs, full_list_TMDs


def improve_ggplot_for_4_plots(axarr,row_nr,col_nr,backgroundcolour,legend_obj):
    ''' Function designed to improve the appearance of matplotlib plots in the ggplot style
    Removes unnecssary ticks, changes background colour, etc.
    '''
    # Remove top axes and right axes ticks
    axarr[row_nr, col_nr].get_xaxis().tick_bottom()
    axarr[row_nr, col_nr].get_yaxis().tick_left()
    #change the position of the axis ticklabels so they are closer to the axis
    axarr[row_nr, col_nr].tick_params(direction = 'out', pad = 0.4)
    #change background colour of graph
    axarr[row_nr, col_nr].set_axis_bgcolor(backgroundcolour)
    #change background colour of legend box
    legend_obj.get_frame().set_facecolor(backgroundcolour)


    
def KW_list_contains_any_desired_KW(KW_list_to_search,list_desired_KW):
    ''' Determine if two lists contain any common values.
    Used to determine, for example, if a list of keywords contains any
    enzyme related words, from another list
    input:
    KW_list_to_search
    list_desired_KW
    note: in theory, this function could be updated to use set(), which should be slightly quicker
    '''
    is_enzyme = False
    for KW in list_desired_KW:
        if KW in KW_list_to_search:
            is_enzyme = True
            break
    return is_enzyme



def get_signif_symbol(number):
    '''
    Takes a number, and returns the approprate symbol for a graph. representing the statistical significance
    '''
    output_signif_string = ''
    signif_dict = {'ns' : (0.05,1.0), '*' : (0.01,0.05),'**' : (0.001, 0.01), '***' : (0.0, 0.001)}
    for key, val in signif_dict.items():
        if val[0] < number < val[1]:
            output_signif_string = key
    return output_signif_string


##define function to obtain regex output (start, stop, etc) as a tuple
#def get_start_and_end_of_TMD_in_query(x):
#    m = re.search(TMD_for_regular_expression_search, x)
#    if m:
#        #if the tmd is in the query, return True, start, stop
#        return (bool(m), m.start(), m.end())
#    else:
#        #if the tmd is not in the query, return False, NaN, NaN
#        return (bool(m), np.nan, np.nan)
def round_sig(x, sig=1):
    from math import floor
    from math import log10
    ''' Rounds a number roughly to a certain number of significant figures.
    Note that the result are floating point numbers that often don't exactly match the expected decimal.
    Note also that 4.5 will be rounded DOWN to 4, not up to 5. np.ceil could be incorporated somehow to fix this.
    http://stackoverflow.com/questions/3410976/how-to-round-a-number-to-significant-figures-in-python
    also see the more sophisticated to_precision script here:
    http://randlet.com/blog/python-significant-figures-format/
    '''
    return round(x, sig-int(floor(log10(x)))-1)

 
def create_hist_from_df_col(df,title,axarr,row_nr,col_nr,settings,data_column,color,alpha,col_width_value,fontsize,xlabel,ylabel,legend):
    #filter to remove sequences where no TMDs are found, 
    df = df.loc[df['list_of_TMDs'].notnull()]
    #filter to remove sequences where no TMDs are found (if string)
    df.loc[df['list_of_TMDs'] != 'nan']
    #filter to remove sequences where no TMDs are found (if string)
    df = df.loc[df['list_of_TMDs'] != 'nan']
    #iterate over the dataframe. Note that acc = uniprot accession here.    
    linspace_binlist = np.linspace(settings["hist_settings_mult_proteins"]["smallest_bin"],
                                   settings["hist_settings_mult_proteins"]["largest_bin"],
                                   settings["hist_settings_mult_proteins"]["number_of_bins"])

    #add 30 as the last bin, to make sure 100% of the data is added to the histogram, including major outliers
    binlist = np.append(linspace_binlist, settings["hist_settings_mult_proteins"]["final_highest_bin"])   
    #create numpy array of membranous over nonmembranous conservation ratios (identity)
    hist_data = np.array(df[data_column].dropna())
    #use numpy to create a histogram
    freq_counts_I, bin_array_I = np.histogram(hist_data, bins=binlist)
    #assuming all of the bins are exactly the same size, make the width of the column equal to 70% of each bin
    col_width = float('%0.3f' % (col_width_value * (bin_array_I[1] - bin_array_I[0])))
    #when align='center', the central point of the bar in the x-axis is simply the middle of the bins ((bin_0-bin_1)/2, etc)
    centre_of_bar_in_x_axis = (bin_array_I[:-2] + bin_array_I[1:-1]) / 2
    #add the final bin, which is physically located just after the last regular bin but represents all higher values
    bar_width = centre_of_bar_in_x_axis[3] - centre_of_bar_in_x_axis[2]
    centre_of_bar_in_x_axis = np.append(centre_of_bar_in_x_axis, centre_of_bar_in_x_axis[-1] + bar_width)
    barcontainer = axarr[row_nr, col_nr].bar(left=centre_of_bar_in_x_axis, height=freq_counts_I,
                                                         align='center', width=col_width, facecolor=color,
                                                         alpha=alpha, edgecolor='black', linewidth=0.1)  # edgecolor='black',
    #label the x-axis for each plot, based on the TMD
    axarr[row_nr, col_nr].set_xlabel(xlabel, fontsize=fontsize)
    #pylab.rcParams['figure.figsize'] = (50.0, 40.0)
    #pylab.rcParams['figure.figsize'] = (20.0, 16.0)
    #xlim_min = settings["hist_settings_mult_proteins"]["xlim_min01"]
    ##take x-axis max from settings
    #xlim_max = settings["hist_settings_mult_proteins"]["xlim_max01"]
    ##set x-axis min
    #axarr[row_nr, col_nr].set_xlim(xlim_min, xlim_max)
    #set x-axis ticks
    #use the slide selection to select every second item in the list as an xtick(axis label)
    axarr[row_nr, col_nr].set_xticks([float('%0.1f' % c) for c in centre_of_bar_in_x_axis[::3]])
    axarr[row_nr, col_nr].set_ylabel(ylabel, rotation='vertical', fontsize=fontsize)
    #change axis font size
    axarr[row_nr, col_nr].tick_params(labelsize=fontsize)
    #create legend?#http://stackoverflow.com/questions/9834452/how-do-i-make-a-single-legend-for-many-subplots-with-matplotlib
    axarr[row_nr, col_nr].legend(legend, loc='upper right',fontsize=fontsize)
    #add title
    #axarr[row_nr, col_nr].set_title(title,fontsize=fontsize)
    #add background grid
    #axarr[row_nr, col_nr].grid(True, color='0.75', alpha=0.5)
        

def create_line_hist_from_df_col(df,title,axarr,row_nr,col_nr,settings,data_column,color,alpha,col_width_value,fontsize,xlabel,ylabel,legend):
    #filter to remove sequences where no TMDs are found, 
    df = df.loc[df['list_of_TMDs'].notnull()]
    #filter to remove sequences where no TMDs are found (if string)
    df.loc[df['list_of_TMDs'] != 'nan']
    #filter to remove sequences where no TMDs are found (if string)
    df = df.loc[df['list_of_TMDs'] != 'nan']
    #iterate over the dataframe. Note that acc = uniprot accession here.    
    linspace_binlist = np.linspace(settings["hist_settings_mult_proteins"]["smallest_bin"],
                                   settings["hist_settings_mult_proteins"]["largest_bin"],
                                   settings["hist_settings_mult_proteins"]["number_of_bins"])

    #add 30 as the last bin, to make sure 100% of the data is added to the histogram, including major outliers
    binlist = np.append(linspace_binlist, settings["hist_settings_mult_proteins"]["final_highest_bin"])   
    #create numpy array of histogram data
    hist_data = np.array(df[data_column].dropna())
    #use numpy to create a histogram
    freq_counts_S, bin_array_S = np.histogram(hist_data, bins=binlist)
    #barcontainer_S = axarr[row_nr,col_nr].bar(left=centre_of_bar_in_x_axis, height=freq_counts_S, align='center', width=col_width, color="#0101DF", edgecolor="#0101DF", alpha = 0.5)
    centre_of_bar_in_x_axis = (bin_array_S[:-2] + bin_array_S[1:-1]) / 2
    #add the final bin, which is physically located just after the last regular bin but represents all higher values
    bar_width = centre_of_bar_in_x_axis[3] - centre_of_bar_in_x_axis[2]
    centre_of_bar_in_x_axis = np.append(centre_of_bar_in_x_axis, centre_of_bar_in_x_axis[-1] + bar_width)        
    #create a line graph rather than a bar graph for the AAISON (ident + similarity)
    linecontainer_AAISON_mean = axarr[row_nr, col_nr].plot(centre_of_bar_in_x_axis, freq_counts_S, color=color,
                                                           alpha=alpha)
    #other colours that are compatible with colourblind readers: #8A084B Dark red, #B45F04 deep orange, reddish purple #4B088A 
    #http://html-color-codes.info/
    #label the x-axis for each plot, based on the TMD
    axarr[row_nr, col_nr].set_xlabel(xlabel, fontsize=fontsize)
    #pylab.rcParams['figure.figsize'] = (50.0, 40.0)
    #pylab.rcParams['figure.figsize'] = (20.0, 16.0)
    #plt.show()
    xlim_min = settings["hist_settings_mult_proteins"]["xlim_min01"]
    #take x-axis max from settings
    xlim_max = settings["hist_settings_mult_proteins"]["xlim_max01"]
    #set x-axis min
    axarr[row_nr, col_nr].set_xlim(xlim_min, xlim_max)
    #set x-axis ticks
    #use the slide selection to select every second item in the list as an xtick(axis label)
    axarr[row_nr, col_nr].set_xticks([float('%0.1f' % c) for c in centre_of_bar_in_x_axis[::3]])
    axarr[row_nr, col_nr].set_ylabel(ylabel, rotation='vertical', fontsize=fontsize)
    #change axis font size
    axarr[row_nr, col_nr].tick_params(labelsize=fontsize)
    #create legend?#http://stackoverflow.com/questions/9834452/how-do-i-make-a-single-legend-for-many-subplots-with-matplotlib
    axarr[row_nr, col_nr].legend(legend, loc='upper right',fontsize=fontsize)
    #add background grid
    axarr[row_nr, col_nr].grid(True, color='0.75', alpha=0.5)

def create_hist_from_df_col_with_auto_binlist(df,title,axarr,num_bins,row_nr,col_nr,settings,data_column,color,alpha,col_width_value,fontsize,xlabel,ylabel,legend):
    #create numpy array of data
    hist_data = np.array(df[data_column].dropna())
    '''
    Calculated the bins for a histogram, even for highly non-normal data
    '''
    #calculate 5th percentile
    percentile_5 = np.percentile(hist_data, 5)#hist_data.min()
    #calculate 9th percentile
    percentile_95 = np.percentile(hist_data, 95)#hist_data.max()
    #calculate difference
    percentile_95_minus_5 = percentile_95 - percentile_5
    #create buffer for bins
    extra_xaxis_range = percentile_95_minus_5 / 4
    #lowest bin is the 5th percentile minus the buffer, except where that is below zero
    data_min = percentile_5 - extra_xaxis_range#hist_data.min()
    #ata_min = 0 if data_max < 0 else data_max
    #highest bin is the 95th percentile
    data_max = percentile_95 + extra_xaxis_range#hist_data.max()
    #create bins using the min and max
    binlist = np.linspace(data_min,data_max,num_bins)
    #print(binlist)
    #use numpy to create a histogram
    freq_counts_I, bin_array_I = np.histogram(hist_data, bins = binlist)
    #print(freq_counts_I)
    #print(bin_array_I)
    #assuming all of the bins are exactly the same size, make the width of the column equal to 70% of each bin
    col_width = float('%0.3f' % (col_width_value * (bin_array_I[1] - bin_array_I[0])))
    #when align='center', the central point of the bar in the x-axis is simply the middle of the bins ((bin_0-bin_1)/2, etc)
    centre_of_bar_in_x_axis = (bin_array_I[:-2] + bin_array_I[1:-1]) / 2
    #add the final bin, which is physically located just after the last regular bin but represents all higher values
    bar_width = centre_of_bar_in_x_axis[3] - centre_of_bar_in_x_axis[2]
    centre_of_bar_in_x_axis = np.append(centre_of_bar_in_x_axis, centre_of_bar_in_x_axis[-1] + bar_width)
    barcontainer = axarr[row_nr, col_nr].bar(left=centre_of_bar_in_x_axis, height=freq_counts_I,
                                                         align='center', width=col_width, facecolor=color,
                                                         alpha=alpha, edgecolor='black', linewidth=0.1)  # edgecolor='black',
    #label the x-axis for each plot, based on the TMD
    axarr[row_nr, col_nr].set_xlabel(xlabel, fontsize=fontsize)
    axarr[row_nr, col_nr].set_ylabel(ylabel, rotation='vertical', fontsize=fontsize)
    #change axis font size
    axarr[row_nr, col_nr].tick_params(labelsize=fontsize)
    #create legend?#http://stackoverflow.com/questions/9834452/how-do-i-make-a-single-legend-for-many-subplots-with-matplotlib
    axarr[row_nr, col_nr].legend(legend, loc='upper right',fontsize=fontsize)
    #add title
    #axarr[row_nr, col_nr].set_title(title,fontsize=fontsize)

def save_fig_with_subplots(fig, axarr, base_filename, fig_nr, fontsize):
    '''saves figure using the given base filename
    tightens using the subplots_adjust function, as an error often occurs with tighten_layout
    removes spines at top and right
    '''
    for ax in axarr.flat:                 
        #change axis font size
        ax.tick_params(labelsize = fontsize)
        #hide spines
        ax.spines["top"].set_visible(False)  
        ax.spines["right"].set_visible(False) 
        ax.tick_params(
            axis='x',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            bottom='off',      # ticks along the bottom edge are off
            right='off',
            top='off',         # ticks along the top edge are off
            labelbottom='on') # labels along the bottom edge are off    
    #automatically tighten the layout of plots in the figure
#    fig.subplots_adjust(bottom = 0)
#    fig.subplots_adjust(top = 1)
#    fig.subplots_adjust(right = 1)
#    fig.subplots_adjust(left = 0)
    fig.tight_layout()
    #save files
    fig.savefig(base_filename + '_%01d.png' % fig_nr, format='png', dpi=400)
    fig.savefig(base_filename + '_%01d.pdf' % fig_nr, format='pdf')

'''
This small function can be used to retrieve the TMD sequence from the full uniprot sequence by applying
the slice function to all rows on a pandas dataframe. In comparison to a loop, the pandas method applied to all rows
simultaneously should be at least 8x faster. Note that the notnull function should remove all np.nan values, but there seems to
be a bug, and it will cause an error when used as a function
For an unknown reason, this is only necessary when the .apply is used in a function.

'''
def slice_uniprot_TMD_seq(x, TMD):
   return x['uniprot_seq'][int(x['%s_start' % TMD] - 1):int(x['%s_end' % TMD])]

def slice_uniprot_TMD_plus_surr_seq(x, TMD):
    return x['uniprot_seq'][int(x['start_surrounding_seq_in_query_%s' % TMD] - 1):int(x['end_surrounding_seq_in_query_%s' % TMD])]

#create small throwaway functions to slice all sequences in dataframe simultaneously
#slice_SW_query_TMD_seq = lambda x: x['query_alignment_sequence'][int(x['%s_start_in_SW_alignment' % TMD]):int(x['%s_end_in_SW_alignment' % TMD])]
#slice_SW_markup_TMD = lambda x: x['alignment_markup'][int(x['%s_start_in_SW_alignment' % TMD]):int(x['%s_end_in_SW_alignment' % TMD])]
#slice_SW_match_TMD_seq = lambda x: x['match_alignment_sequence'][int(x['%s_start_in_SW_alignment' % TMD]):int(x['%s_end_in_SW_alignment' % TMD])]
#slice_SW_query_TMD_seq_plus_surr = lambda x: x['query_alignment_sequence'][int(x['%s_start_in_SW_alignment_plus_surr' % TMD]):int(x['%s_end_in_SW_alignment_plus_surr' % TMD])]
#slice_SW_markup_TMD_plus_surr = lambda x: x['alignment_markup'][int(x['%s_start_in_SW_alignment_plus_surr' % TMD]):int(x['%s_end_in_SW_alignment_plus_surr' % TMD])]
#slice_SW_match_TMD_seq_plus_surr = lambda x: x['match_alignment_sequence'][int(x['%s_start_in_SW_alignment_plus_surr' % TMD]):int(x['%s_end_in_SW_alignment_plus_surr' % TMD])]

#create small throwaway functions to slice all sequences in dataframe simultaneously
def slice_SW_query_TMD_seq(x, TMD):
    return x['query_alignment_sequence'][int(x['%s_start_in_SW_alignment' % TMD]):int(x['%s_end_in_SW_alignment' % TMD])]
def slice_SW_markup_TMD(x, TMD):
    return x['alignment_markup'][int(x['%s_start_in_SW_alignment' % TMD]):int(x['%s_end_in_SW_alignment' % TMD])]
def slice_SW_match_TMD_seq(x, TMD):
    return x['match_alignment_sequence'][int(x['%s_start_in_SW_alignment' % TMD]):int(x['%s_end_in_SW_alignment' % TMD])]
def slice_SW_query_TMD_seq_plus_surr(x, TMD):
    return x['query_alignment_sequence'][int(x['%s_start_in_SW_alignment_plus_surr' % TMD]):int(x['%s_end_in_SW_alignment_plus_surr' % TMD])]
def slice_SW_markup_TMD_plus_surr(x, TMD):
    return x['alignment_markup'][int(x['%s_start_in_SW_alignment_plus_surr' % TMD]):int(x['%s_end_in_SW_alignment_plus_surr' % TMD])]
def slice_SW_match_TMD_seq_plus_surr(x, TMD):
    return x['match_alignment_sequence'][int(x['%s_start_in_SW_alignment_plus_surr' % TMD]):int(x['%s_end_in_SW_alignment_plus_surr' % TMD])]
def find_indices_longer_than_prot_seq(df, TMD):
    return df['end_surrounding_seq_in_query_%s' % TMD] > df['uniprot_seqlen']

def create_indextuple_nonTMD_last(x):
    ''' Joins two columns into a tuple. Used to create the last tuple of the nonTMD region in the sequence.
    '''
    return (int(x['nonTMD_index_tuple_last0']), int(x['nonTMD_index_tuple_last1']))

def slice_with_tuple(string, tup):
    '''A function to slice a single string, taking the start and stop indices from a tuple
    '''
    #print(tup[0])
    return string[int(tup[0]):int(tup[1])]

def slice_with_nested_tuple(string, nested_tuple):
    '''A function to slice a sequence multiple times, using the indices from nested tuples
    '''
    #convert nested tuple from string to tuple 
    nested_tuple = eval(nested_tuple)
    #for each tuple, slice the input string. Make a list of all the sliced strings. Join list with no gaps
    return ''.join([slice_with_tuple(string, tup) for tup in nested_tuple])

def retrieve_selected_uniprot_records_from_flatfile(input_accession_list, large_input_uniprot_flatfile, output_flatfile):
    '''
    Function to select records from a large uniprot flatfile, and save them as a smaller flatfile of selected records.
    Input = list of uniprot accessions
    '''
    #accession_list = [line.strip() for line in open(input_accession_list, "r")]
    #open input flatfile
    uniprot_index_handle = SeqIO.index(large_input_uniprot_flatfile, "swiss")  
    #create an empty list to hold all the accessions that are not in the flatfile
    list_acc_not_in_flatfile = []
    with open(output_flatfile, "wb") as output:    
        for acc in input_accession_list:
            try:                  
                #get the raw uniprot record, and write to the output file
                output.write(uniprot_index_handle.get_raw(acc))              
            #if the record is not in the file
            except KeyError:
                list_acc_not_in_flatfile.append(acc)
    print("SwissProt records not found in %s:\n%s." % (large_input_uniprot_flatfile, list_acc_not_in_flatfile))


def get_start_and_end_of_TMD_in_query(x, regex_string):
    '''
    Returns a tuple containing (bool, start, stop) showing the location of a regex pattern
    in the target string.
    To be used with Pandas Dataframes. The aim is to conduct the computationally expensive regex
    search only once to obtain both the start and stop indices.
    Variables:
    x = target sequence (e.g. unicode string, DNA, Protein Seq)
    TMD_for_regular_expression_search = regex search string
    '''
    m = re.search(regex_string, x)
    if m:
        #if the tmd is in the query, return True, start, stop
        return (bool(m), m.start(), m.end())
    else:
        #if the tmd is not in the query, return False, NaN, NaN
        return (bool(m), np.nan, np.nan)

        
def find_disallowed_words(description, words_not_allowed_in_description):
    '''Finds disallowed words in the description (Patent, Synthetic, etc). 
    Returns the disallowed words, or an empty list. The lists are converted to strings.
    settings must be in globals
    '''
    #words_not_allowed_in_description = settings["simap_match_filters"]["words_not_allowed_in_description"]
    list_of_list_disallowed_words_in_descr = []
    for disallowed_word in words_not_allowed_in_description:
        if disallowed_word in description:
            #print(disallowed_word)
            list_of_list_disallowed_words_in_descr.append(disallowed_word)
    return str(list_of_list_disallowed_words_in_descr)

def create_dict_organising_subplots(n_plots_per_fig,n_rows,n_cols):
    '''
    Function to help organise the creation of figures that contain multiple plots.
    For example, 15 histograms printed in figures with 8 histograms per figure/page.
    Returns a dict that gives a tuple for each plot/graph. 
    newfig, savefig, fig_nr, plot_nr_in_fig, row_nr, col_nr
    row_nr and col_nr are used to index pyplot subplots as follows
    fig, axarr = plt.subplots(2,2)
    _ = axarr[row_nr,col_nr].plot(x, y)
    '''
    dict_organising_subplots = {}
    #figure number
    fig_nr = 0
    #plot number in figure
    plot_nr_in_fig = 0
    #row number in figure
    row_nr = 0
    #column number in figure
    col_nr = 0
    #whether the figure needs to be saved
    savefig = False
    #whether a new figure needs to be created
    newfig = True
    
    for plotnr in range(1, 500):
        #add current counters to dict
        dict_organising_subplots[plotnr] = (newfig, savefig, fig_nr, plot_nr_in_fig, row_nr, col_nr)
        plot_nr_in_fig += 1
        row_nr += 1
        newfig = False
        savefig = False
        #if plot_nr_in_fig is the last one before the new figure, then savefig = True
        if plot_nr_in_fig % (n_plots_per_fig - 1) == 0 and plot_nr_in_fig != 0:
            savefig = True
        #if plot_nr_in_fig is in a multiple of n_rows, then the plot goes to the second column
        if plot_nr_in_fig % n_rows == 0 and plot_nr_in_fig != 0:
            #print(plot_nr_in_fig)
            col_nr += 1
            row_nr = 0
        #if the plotnr is in a multple of n_plots_per_fig, then a new figure needs to created, and everything else reset
        if plotnr % n_plots_per_fig == 0 and plotnr != 0:
            #go to second figure
            fig_nr += 1
            #reset values
            plot_nr_in_fig = 0
            row_nr = 0
            col_nr = 0
            newfig = True
    return dict_organising_subplots

def create_new_fig_if_necessary(newfig, fig, axarr, nrows_in_each_fig, ncols_in_each_fig, dpi = 300):
    if newfig:
        #close any open figures
        plt.close('all')
        #create a new figure
        fig_new, axarr_new = plt.subplots(nrows=nrows_in_each_fig, ncols=ncols_in_each_fig, dpi=dpi)
        #if a new fig needs to be created, return new fig and axarr objects
        return fig_new, axarr_new
    else:
        #if a new fig does not need to be created, return the original fig and axarr objects
        return fig, axarr


def check_tarfile(df,i):
    '''
    Checks the tarball that contains the SIMAP output. 
    Looks to see if the tarball exists, if it is corrupted, if it contains the feature table and homologues from simap.
    '''
    if os.path.isfile(df.loc[i,'SIMAP_feature_table_XML_file_path']):
        feature_table_XML_exists = True
    else:
        feature_table_XML_exists = False
    if os.path.isfile(df.loc[i, 'SIMAP_homologues_XML_file_path']):
        homologues_XML_exists = True
    else:
        homologues_XML_exists = False
    if os.path.isfile(df.loc[i, 'SIMAP_tarfile']):
        SIMAP_tarfile_exists = True
    else:
        SIMAP_tarfile_exists = False                
    #at the moment, we'll only create the tarfile if both the feature table and homologue XML files downloaded successfully, but this might change depending on preference
    #feature_table_in_tarfile = False
    #homologues_XML_in_tarfile = False
    if SIMAP_tarfile_exists:
        try:
            with tarfile.open(df.loc[i, 'SIMAP_tarfile'], mode = 'r:gz') as tar:
                if df.loc[i,'SIMAP_feature_table_XML_file'] in [tarinfo.name for tarinfo in tar]:
                    feature_table_in_tarfile = True
                #else:
                #    feature_table_in_tarfile = False
                if df.loc[i,'SIMAP_homologues_XML_file'] in [tarinfo.name for tarinfo in tar]:
                    homologues_XML_in_tarfile = True
                #else:
                #    homologues_XML_in_tarfile = False 
        except EOFError:
            SIMAP_tarfile_exists = False 
    return feature_table_XML_exists, homologues_XML_exists, SIMAP_tarfile_exists

# VERKSOIK
def score_pairwise(seq1, seq2, matrix, gap_open_penalty, gap_extension_penalty, prev_site_contained_gap = True):
    '''
    Calculates a score between two aligned sequences, based on the gap penalties and matrix applied. 
    The zip seems to be a fast method of comparing individual letters in two strings of the same length
    A, B are each paired amino acid in the pairwise alignment
    yield is a generator function that returns a result. See http://stackoverflow.com/questions/231767/the-python-yield-keyword-explained
    yield should be faster than iterators, because the result does not need to be held in memory to access a second time, it woll only be read once     
    '''
    for A,B in zip(seq1, seq2):
        #ORIG SCRIPT: if either A or B is a gap, gap_exists = True. BUT this can't be used for self score ratio calculation! The presence of gaps in both the query and the subject shouldn't be penalised!        
        gap_exists = ('-'==A) or ('-'==B)
        #MT addition to script: determine if both sequences contain a gap at this position, and if they do, yield 0
        gap_in_both_query_and_match = True if ('-'==A) and ('-'==B) else False
        if gap_in_both_query_and_match:
            yield 0
        else:
            #easiest if read backwards: return the matrix value for A to B, unless A or B is a gap: return the gap open penalty, unless the previous aa pair was also a gap, return the gap extension penalty
            try:
                yield (gap_extension_penalty if prev_site_contained_gap else gap_open_penalty) if gap_exists else matrix[(A,B)]
            #in some cases, B is used as Asp or Asn. These should be very rare. Sequences with X are already removed. 
            except KeyError:
                yield 0
                logging.info('sequence pair contains non-IUPAC character: %s to %s' % (A,B))
            #the last amino acid pair contained a gap, so an extension penalty should be used instead of an opening penalty
        prev_site_contained_gap = gap_exists

def create_dict_of_data_from_uniprot_record(record):
    '''
    For each uniprot record, collect the desired data into a dictionary.
    '''
    global uniprot_TMD_start, uniprot_TMD_end, uniprot_TMD_description, uniprot_TMD_sequence     
    global TRANSMEM_missing_from_uniprot_features, output_dict, accession, uniprot_TMD_sequence, record_entry_name, record_gene_name, record_description, uniprot_TMD_start, uniprot_TMD_end, uniprot_TMD_description, record_sequence_length, comments_subcellular_location, record_keywords, record_features_all
    
    #convert the comments in uniprot to a dictionary. 
    #the comments in uniprot are unordered and non-hierarchical and need to be processed with string manipulation     
    comments_dict = {} #create a new empty dictionary
    create_dictionary_of_comments(record, comments_dict)            
    
    #create an empty output dictionary to holnd the uniprot data for each record
    output_dict = {}
    
    #print accession number
    logging.info(record.accessions[0])

    #by default this is zero
    output_dict['uniprot_TMD_start'] = 0
    output_dict['uniprot_TMD_end'] = 0
    output_dict['uniprot_TMD_description'] = 0
    output_dict['uniprot_record_included_in_csv'] = True
 
    #add data to dictionary
    output_dict['accession_uniprot'] = record.accessions[0]
    output_dict['full_list_of_accessions_uniprot'] = record.accessions
    output_dict['record_entry_name_uniprot'] = record.entry_name
    output_dict['record_gene_name_uniprot'] = record.gene_name
    output_dict['record_description_uniprot'] = record.description
    output_dict['sequence_uniprot'] = record.sequence
    output_dict['organism_classification'] = record.organism_classification
    output_dict['organism'] = record.organism
    output_dict['record_keywords_uniprot'] = record.keywords
    output_dict['record_features_all_uniprot'] = record.features
    output_dict['record_sequence_length_uniprot'] = record.sequence_length
    output_dict['comments_subcellular_location_uniprot'] = comments_dict['SUBCELLULAR LOCATION']

    #create a list of all the feature types (signal, transmem, etc)
    list_of_feature_types_in_uniprot_record = []
    for sublist in record.features:
        list_of_feature_types_in_uniprot_record.append(sublist[0])
        #logging.info(sublist)
    
    #list of the features that we want in the final csv
    desired_features_in_uniprot = ['TRANSMEM', 'VARIANT', 'CONFLICT', 'VAR_SEQ','VARSPLIC']
    desired_features_in_uniprot_dict = {}
    
    location_of_tmds_in_feature_list = []  
    
    for feature in desired_features_in_uniprot:
        if feature in list_of_feature_types_in_uniprot_record:
            #find the features in the feature list. For polytopic membrane protoins, there will be more than one tmd.
            location_of_features_in_feature_list = [i for i,x in enumerate(list_of_feature_types_in_uniprot_record) if x == feature]
            desired_features_in_uniprot_dict[feature] = location_of_features_in_feature_list        
            if feature == 'TRANSMEM':
                location_of_tmds_in_feature_list = location_of_features_in_feature_list
    
    #determine if the TMD is actually annotated in the uniprot record
    if 'TRANSMEM' not in desired_features_in_uniprot_dict.keys():
        TRANSMEM_missing_from_uniprot_features = True
    else:
        TRANSMEM_missing_from_uniprot_features = False
    output_dict['TRANSMEM_missing_from_uniprot_features'] = TRANSMEM_missing_from_uniprot_features
    output_dict['more_than_one_tmd_in_uniprot_annotation'] = False     
    
    if not TRANSMEM_missing_from_uniprot_features:  
        #add the TMD to the dictionary for single-pass membrane proteins
        if(len(location_of_tmds_in_feature_list)) == 1:
            location = location_of_tmds_in_feature_list[0]
            output_dict['number_of_tmds_in_seq'] = 1
            output_dict['uniprot_TMD_start'] = record.features[location][1]
            output_dict['uniprot_TMD_end'] = record.features[location][2]
            output_dict['uniprot_TMD_description'] = record.features[location][3]  
        
        #add the TMD to the dictionary for multi-pass membrane proteins
        if(len(location_of_tmds_in_feature_list)) > 1:
            logging.info('more than one "TRANSMEM" feature in uniprot for %s' % output_dict['accession_uniprot'])
            output_dict['number_of_tmds_in_seq'] = len(location_of_tmds_in_feature_list)
            output_dict['more_than_one_tmd_in_uniprot_annotation'] = True
            output_dict['uniprot_record_included_in_csv'] = False
#            for i in range(len(location_of_tmds_in_feature_list)):
#                output_dict['uniprot_TMD%s_start' % i] = record.features[location[i]][1]
#                output_dict['uniprot_TMD%s_end' % i] = record.features[location[i]][2]
#                output_dict['uniprot_TMD%s_description' % i] = record.features[location[i]][3]      
        
        if output_dict['number_of_tmds_in_seq'] == 1:
            #create a numpy array of any sequence variants are in the TMD region
            list_of_variant_types_in_uniprot = ['VARIANT', 'CONFLICT','VARSPLIC', 'VAR_SEQ']
            #array_of_all_variants_in_tmd = np.zeros(4)
            array_of_all_variants_in_tmd = ([])
            for variant_type in list_of_variant_types_in_uniprot:
                if variant_type in desired_features_in_uniprot_dict.keys():
                    list_of_variant_locations = list(desired_features_in_uniprot_dict[variant_type])       
                    for i in range(len(list_of_variant_locations)):
                        start_of_variant_in_seq = record.features[list_of_variant_locations[i]][1]
                        end_of_variant_in_seq = record.features[list_of_variant_locations[i]][2]
                        variant_description = record.features[list_of_variant_locations[i]][3]
                        variant_feature_identifier = record.features[list_of_variant_locations[i]][4]
                       #check if the variant is in the tmd   
                        start_of_variant_is_after_start_of_tmd = True if start_of_variant_in_seq > output_dict['uniprot_TMD_start'] else False
                        end_of_variant_is_before_end_of_tmd = True if end_of_variant_in_seq < output_dict['uniprot_TMD_end'] else False
                        variant_is_in_tmd = True if all(start_of_variant_is_after_start_of_tmd and end_of_variant_is_before_end_of_tmd) else False
                        
                        #add to numpy array that contains all the variants in the tmd region
                        if variant_is_in_tmd:
                            variant_array = np.array([variant_type, start_of_variant_in_seq, end_of_variant_in_seq, variant_description,variant_feature_identifier])
                            if array_of_all_variants_in_tmd != ([]):
                                array_of_all_variants_in_tmd = np.row_stack((array_of_all_variants_in_tmd, variant_array))
                            else:
                                array_of_all_variants_in_tmd = variant_array
            #if there were variants added, add them to the output dictionary
            output_dict['array_of_all_variants_in_tmd'] = array_of_all_variants_in_tmd
    

    #if the tmd region is annotated, get the sequence
    if not TRANSMEM_missing_from_uniprot_features:
        if output_dict['number_of_tmds_in_seq'] == 1:       
            output_dict['uniprot_TMD_sequence'] = record.sequence[output_dict['uniprot_TMD_start'] - 1:output_dict['uniprot_TMD_end']]           
        else:
            output_dict['uniprot_record_included_in_csv'] = False
            output_dict['uniprot_TMD_start'] = 0
            output_dict['uniprot_TMD_end'] = 0
            output_dict['uniprot_TMD_description'] = 0
    else:    
        output_dict['uniprot_record_included_in_csv'] = False
        logging.info('%s: "TRANSMEM" not found in uniprot features, therefore not incuded in csv file for further analysis' % output_dict['accession_uniprot'])
        output_dict['uniprot_TMD_start'] = 0
        output_dict['uniprot_TMD_end'] = 0
        output_dict['uniprot_TMD_description'] = 0
            
    #decide if the sequence goes in the csv file (more filters will probably be added later)
    #output_dict['uniprot_record_included_in_csv'] = True if not all(TRANSMEM_missing_from_uniprot_features and output_dict['more_than_one_tmd_in_uniprot_annotation']) else False
    return output_dict


#def create_list_of_files_from_csv_with_uniprot_data(input_file, list_of_keys):
#    '''
#    Generate the list of filenames, assuming SIMAP has already run
#    '''
#    #The list of keys are the headers in the csv file that should be inserter into the dictionary
#    #nested_dict_with_uniprot_seq = create_nested_dict_from_csv(input_file, list_of_keys)
#    #I want to save the uniprot entries based on their domain in the tree of life
#    list_of_files_with_feature_tables = []
#    list_of_files_with_homologues = []
#    list_of_protein_names = []
#    list_of_org_domains = []
#    for i in range(len(df_csv_file_with_uniprot_data)):
#        organism_classification_string = (df_csv_file_with_uniprot_data.loc[i,'organism_classification'])
#        organism_domain = convert_stringlist_to_list(organism_classification_string)[0]
#        list_of_org_domains.append(organism_domain)        
#        protein_name = '%s_%s' % (df_csv_file_with_uniprot_data.loc[i,'accession_uniprot'], df_csv_file_with_uniprot_data.loc[i,'record_entry_name_uniprot']) 
#        list_of_protein_names.append(protein_name)
#        SIMAP_feature_table_XML_file = r"E:\Databases\simap\%s\%s_feature_table.xml" % (organism_domain, protein_name)
#        list_of_files_with_feature_tables.append(SIMAP_feature_table_XML_file)
#        SIMAP_homologues_XML_file = r"E:\Databases\simap\%s\%s_homologues.xml" % (organism_domain, protein_name)
#        list_of_files_with_homologues.append(SIMAP_homologues_XML_file) 
#    return list_of_files_with_feature_tables, list_of_files_with_homologues, list_of_protein_names, list_of_org_domains


def retrieve_simap_feature_table(input_sequence, max_memory_allocation, output_file):
    '''
    Uses the java program to access the simap database and download the small file containing information on that protein, called a "feature table".
    '''
    #global number_of_files_not_found    
    #jar_file = r'E:\\Stephis\\Projects\\Programming\\Python\\programs\\eaSimap.jar'
    #jar_file = r'/nas/teeselab/programs/eaSimap.jar'
    jar_file = r"D:\Schweris\Projects\Programming\Python\programs\eaSimap.jar"
    #prepare input sequence and settings as a "java_run_string"
    #run command
    string_to_run_as_java_command = 'java -Xmx%im -jar %s -s %s -o %s -f' % (max_memory_allocation, jar_file, input_sequence, output_file)
    logging.info(string_to_run_as_java_command)
#    for output_line in run_command(string_to_run_as_java_command):
#        print(output_line)                          #note that nothing is actually printed
    command = Command(string_to_run_as_java_command)
    command.run(timeout=100)
    logging.info("Output file:     %s\n" % output_file), 
    
    #Indent the XML file in a more readable format. Note that this adds empty lines. 
    #I need to use more complex code, or switch to lxml at some stage    
#    try:    
#        xml = minidom.parse(output_file)
#        xml.normalize()
#        with open(output_file, 'w') as result:
#            result.write(xml.toprettyxml(indent = '  '))
#        #number_of_files_not_found = 0
#    except FileNotFoundError:
##        if 'number_of_files_not_found' not in globals():
##            number_of_files_not_found = 1
##        else:
##            number_of_files_not_found += 1        
#        logging.info('********************FileNotFoundError: for %s***************' % output_file)
    
    #If there are more than 10 errors, stop the attempts to download and sleep for 6 hours   
#    if number_of_files_not_found == 20:
#        logging.warning('number of files not found = %s, sleeping for 24 hours' % number_of_files_not_found)
#        sleep_24_hours()
#    if number_of_files_not_found == 3:       
#        logging.warning('number of files not found = %s, sleeping for 24 hours' % number_of_files_not_found)
#        sleep_6_hours()
    sleep_15_seconds()
    if not os.path.exists(output_file):
        logging.info('********************SIMAP download failed for : %s***************' % output_file)
#    #sleep for 30 seconds to not overload the server 
#    sys.stdout.write("sleeping .")
#    for i in range(5):
#        time.sleep(3)
#        sys.stdout.write(" .")
#        sys.stdout.flush()
#    print(' .')

def retrieve_simap_homologues(input_sequence, output_file, database, max_hits, max_memory_allocation, taxid):
    '''
    Uses the java program to access the simap database and download the large file containing all homologues of that protein.
    '''
    #global number_of_files_not_found    
    #locate jar executable file
    #jar_file = r'E:\\Stephis\\Projects\\Programming\\Python\\programs\\eaSimap.jar'
    #jar_file = r'/nas/teeselab/programs/eaSimap.jar'
    jar_file = r"D:\Schweris\Projects\Programming\Python\programs\eaSimap.jar"
    #set parameters
    #database = '' #leave blank('') for SIMAP all database. 'uniprot_swissprot' 'uniprot_trembl' 'refseq' 'Escherichia coli' 'Homo sapiens' 'Hot springs metagenome'
    database_dictionary = {313: 'uniprot_swissprot', 314: 'uniprot_trembl', 595: 'refseq', 721: 'Escherichia coli', 1296: 'Homo sapiens', 4250: 'Hot springs metagenome'}
    database_dictionary_reversed = {}     
    for v,k in database_dictionary.items():
        database_dictionary_reversed[k] = v
    print("in retr simap, database = %s" % database)
    database_search_string = '' if database == '""' else '-d %s' % database_dictionary_reversed[database]  # -d in java interface is currently not working, but could be functional at some stage
    taxid_search_string = '' if taxid == '""' else '-i %s' % taxid
    #note that windows has a character limit in the command prompt in theory of 8191 characters, but the command line java command seems to cause errors with sequences above 3000 amino acids.
    #the 3000 character limit is currently applied in the main_simap script, rather than here
    #run command
    string_to_run_as_java_command = 'java -Xmx%im -jar %s -s %s -m %s -o %s -x %s%s' % (max_memory_allocation,jar_file, input_sequence,
                                                                                           max_hits, output_file, 
                                                                                          database_search_string, taxid_search_string)
    logging.info(string_to_run_as_java_command)
#    for output_line in run_command(string_to_run_as_java_command):
#        print(output_line)                           #note that nothing is actually printed
    command = Command(string_to_run_as_java_command)
    timeout = max_hits/5 if max_hits > 500 else 100
    command.run(timeout=timeout) #give 1000 for 5000 hits to download?   
    logging.info("Output file:     %s\n'file saved'" % output_file)
    sleep_120_seconds()
    if not os.path.exists(output_file):
        logging.info('********************SIMAP download failed for : %s***************' % output_file)
    #print without new line
#    sys.stdout.write("sleeping .")
#    for i in range(5):
#        time.sleep(3)
#        sys.stdout.write(" .")
#        sys.stdout.flush()
#    print('.')
    '''There are many homologue XML files with nodes missing! Could this be due to the minidom parse?? 
    Maybe it's better to leave this out, and only parse to a readable format for some example proteins???
    '''
    
    #Indent the XML file in a more readable format. Note that this adds empty lines. 
    #I need to use more complex code, or switch to lxml at some stage.
#    try:    
#        xml = minidom.parse(output_file)
#        xml.normalize()
#        with open(output_file, 'w') as result:
#            result.write(xml.toprettyxml(indent = '  '))
#        number_of_files_not_found = 0
#    except FileNotFoundError:
#        logging.info('********************FileNotFoundError: for %s********************' % output_file)
#        if 'number_of_files_not_found' not in globals():
#            number_of_files_not_found = 1
#        else:
#            number_of_files_not_found += 1
#        print('number_of_files_not_found = %s' % number_of_files_not_found)
#    number_of_files_not_found = 0
#    if os.path.exists(output_file):
#        if 'number_of_files_not_found' not in globals():
#            pass
#        else:
#            number_of_files_not_found += 1        
        #If there are more than 10 errors, stop the attempts to download and sleep for 6 hours   
#    if number_of_files_not_found > 30:
#        sleep_24_hours()
#    if number_of_files_not_found == 20:
#        sleep_24_hours()
#    if number_of_files_not_found == 15:       
#        sleep_6_hours()
#    if number_of_files_not_found == 10:       
#        sleep_6_hours()
#def retrieve_simap_from_multiple_fasta(input_file):
#    records = SeqIO.parse(input_file, "fasta")
#    global list_of_files_with_feature_tables, list_of_files_with_homologues
#    list_of_files_with_feature_tables = []
#    list_of_files_with_homologues = []
#    recordcounter = 0
#    for record in records:
#        name = record.name.replace('|', '_')[:30]
#        accession = 'Acc'
#        label = '%s_%s' % (accession, name)
#        print(label)
#        SIMAP_feature_table_XML_file = r"E:\\Stephis\\Projects\\Programming\\Python\\files\\learning\\simap\\%s_simap_feature_table.xml" % label
#        list_of_files_with_feature_tables.append(SIMAP_feature_table_XML_file)
#        SIMAP_homologues_XML_file = r"E:\\Stephis\\Projects\\Programming\\Python\\files\\learning\\simap\\%s_simap_homologues.xml" % label
#        list_of_files_with_homologues.append(SIMAP_homologues_XML_file) 
#        retrieve_simap_feature_table(input_sequence=record.seq, output_file=SIMAP_feature_table_XML_file)
#        retrieve_simap_homologues(input_sequence=record.seq, output_file=SIMAP_homologues_XML_file, database='', max_hits='10', taxid='7227', extra_search_string='')
#        recordcounter += 1    
#    logging.info('Download complete, %s SIMAP records saved.' % recordcounter)                       
#input_seqs_mult_fasta = r'E:\Stephis\Projects\Programming\Python\files\learning\simap\multiple_protein_seqs_in_fasta_format.txt'
#retrieve_simap_from_multiple_fasta(input_seqs_mult_fasta)
#throwaway functions, currently kept in main
#slice_TMD_seq = lambda x: x['uniprot_seq'][int(x['%s_start' % TMD_name]-1):int(x['%s_end' % TMD_name])]
#slice_TMD_plus_surrounding_seq = lambda x: x['uniprot_seq'][int(x['start_surrounding_seq_in_query_%s' % TMD_name]-1):int(x['end_surrounding_seq_in_query_%s' % TMD_name])]

class Command(object):
    '''
    subprocess for running shell commands in win and linux 
    This will run commands from python as if it was a normal windows console or linux terminal.
    taken from http://stackoverflow.com/questions/17257694/running-jar-files-from-python)'
    '''
    def __init__(self, cmd):
        self.cmd = cmd
        self.process = None

    def run(self, timeout):
        def target():
            logging.info('Thread started')
            self.process = subprocess.Popen(self.cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            #self.process.communicate()
            stdout, stderr = self.process.communicate() # from http://stackoverflow.com/questions/14366352/how-to-capture-information-from-executable-jar-in-python
            logging.info('JAVA OUTPUT:%s' % stdout.decode("utf-8"))
            logging.warning('JAVA FAULTS: %s' %stderr.decode("utf-8"))
            logging.info('Thread finished')

        thread = threading.Thread(target=target)
        thread.start()

        thread.join(timeout)
        if thread.is_alive():
            logging.info('Terminating process')
            self.process.terminate()
            thread.join()
        logging.info(self.process.returncode)

def run_command(command):
    #this stopped working for some reason. Did I mess up a path variable?    
    p = subprocess.Popen(command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT)
    return iter(p.stdout.readline, b'')





def sleep_15_seconds():
    # sleep for 30 seconds to not overload the server
    sys.stdout.write("sleeping .")
    sys.stdout.flush()
    for i in range(15):
        time.sleep(1)
        sys.stdout.write(" .")
        sys.stdout.flush()
    print(' .')

def sleep_120_seconds():
    #sleep for 30 seconds to not overload the server 
    sys.stdout.write("sleeping .")
    sys.stdout.flush()
    for i in range(12):
        time.sleep(15)
        sys.stdout.write(" .")
        sys.stdout.flush()
    print(' .\n')

def sleep_6_hours():
    #sleep for 30 seconds to not overload the server 
    sys.stdout.write("sleeping .")
    sys.stdout.flush()
    for i in range(6):
        time.sleep(3600)
        sys.stdout.write(" .")
        sys.stdout.flush()
    print(' .\n')

def sleep_24_hours():
    #sleep for 30 seconds to not overload the server 
    sys.stdout.write("sleeping .")
    sys.stdout.flush()
    for i in range(24):
        time.sleep(3600)
        sys.stdout.write(" .")
        sys.stdout.flush()
    print(' .\n')

def create_dictionary_of_comments(uniprot_record_handle, output_dictionary):
    try:
        for comment in uniprot_record_handle.comments:
            # splits comments based on first ":" symbol, creates a list called split_comment        
            split_comment = comment.strip().split(': ', 1)               
             # several comments have the same name. need to check if it is already in the dictionary        
            if split_comment[0] in output_dictionary:              
                # list the different comments, one after another            
                output_dictionary[split_comment[0]] += ", %s" % split_comment[1] 
            else:
                output_dictionary[split_comment[0]] = split_comment[1]
    except AttributeError:
        #there are no comments in this uniprot file!
        logging.info('no comments in Uniprot file')
        output_dictionary = {}
        

#set up a function for showing object names, in order to write the csv header from a list of objects
def name_of_object_in_list_of_global_objects(object):
    for name_of_object,oid in globals().items():
        if oid is object:
            return name_of_object

def create_list_of_object_names_from_list_of_objects(input_list_of_objects):
    output_list = []
    for i in range(len(input_list_of_objects)):
        objectname = name_of_object_in_list_of_global_objects(input_list_of_objects[i])
        output_list.append(objectname) 
    return output_list

def save_list_as_row_in_csv(input_list, output_csv, open_method):
    #import csv
    open_method_with_quotation_marks =   "%s" % open_method
    with open(output_csv, open_method_with_quotation_marks) as f:
        writer = csv.writer(f, delimiter=',', quotechar='"', lineterminator='\n', quoting=csv.QUOTE_NONNUMERIC, doublequote = True)
        writer.writerow(input_list)

def create_regex_string(inputseq):
    ''' adds '-*' between each aa or nt/aa in a DNA or protein sequence, so that a particular
    aligned sequence can be identified via a regex search, even if it contains gaps
    inputseq : 'LQQLWNA'
    output   : 'L-*Q-*Q-*L-*W-*N-*A'
    '''
    search_string = ''
    for letter in inputseq:
        letter_with_underscore = letter + '-*'
        search_string += letter_with_underscore
    return search_string[:-2]


def count_non_protein_characters(inputseq):
    number_of_non_protein_characters = 0
    accepted_protein_alphabet = 'ACDEFGHIKLMNPQRSTVWY-'
    for character in inputseq:
        if character not in accepted_protein_alphabet:
            number_of_non_protein_characters += 1
    return number_of_non_protein_characters

def create_csv_header_fieldnames(input_dict):
    pass
    #creates a list that starts with the important fields, but also includes any new fields inserted into the dictionary
#    list_of_important_fields = ['hit_number', 'TMD_seq_in_hsp_match', 'expectation', 'organism', 'description', 'database', 
#    'ratio_percentage_identity_of_TMD_to_rest_of_hsp', 'md5', 'databaseId', 'percentage_identity_of_TMD', 
#    'identity', 'is_TMD_in_hsp', 'match_TMD_added_to_FastA_alignment', 'non_protein_characters', 'number_of_gaps_in_match_TMD', 
#    'number_of_gaps_in_query_TMD', 'percentage_identity_of_rest_of_alignment', 
#    'ratio_length_of_TMD_to_rest_of_hsp', 'ratio_length_of_query_TMD_to_rest_of_match_protein', 'taxonomy_node_id']
    #the list of keys from the dictionary   
#    keylist00 = list(input_dict.keys())
#    keylist = sorted(keylist00)
#    #remove any keys that are already in the list above
#    for field in list_of_important_fields:
#        if field in keylist:
#            keylist.remove(field)
##    join dictionaries
#    csv_header_fieldnames = list_of_important_fields + keylist
##    logging.info('\n\n\n\n')    
##    logging.info(csv_header_fieldnames)    
#    return csv_header_fieldnames

def save_dict_keys_as_header_in_csv(input_dict, header_fieldnames, output_csv, open_method):
    #import csv
    open_method_with_quotation_marks = "%s" % open_method
    with open(output_csv, open_method_with_quotation_marks) as f:
        writer = csv.DictWriter(f, fieldnames=header_fieldnames, extrasaction='ignore', delimiter=',', quotechar='"', lineterminator='\n', quoting=csv.QUOTE_NONNUMERIC, doublequote = True)
        writer.writeheader()

def save_dict_values_as_row_in_csv(input_dict, header_fieldnames, output_csv, open_method):
    #import csv
    open_method_with_quotation_marks =   "%s" % open_method
    with open(output_csv, open_method_with_quotation_marks) as f:
        #the extrasaction='ignore' should avoid the error that the dictionary contains fields that are not going to be written
        writer = csv.DictWriter(f, fieldnames=header_fieldnames, extrasaction='ignore', delimiter=',', quotechar='"', lineterminator='\n', quoting=csv.QUOTE_NONNUMERIC, doublequote = True)
        writer.writerow(input_dict)

def create_nested_dict_from_csv(csvfile, fieldlist):
    global selected_dict
    '''
    Choose a list of fields that you want to include in your final dictionary.
    
    fieldlist='all'
        with this option, all columns will be included in the final dict 
        don't use this option for large data files! Try Numpy and load data as an array.
    
    For the nested dictionary, the row number is used as the key.
    {1: {dictionary from all values in row 1}, 2: {dictionary from all values in row 2}, 
    
    if this data contains for example a uniprot number, this can be accessed from the nested dictionary as follows:
        nested_dict_with_uniprot_seq = create_nested_dict_from_csv(csv_file_with_uniprot_data, list_of_keys_to_keep)
        uniprot_number_for_seq_in_row_1 = nested_dict_with_uniprot_seq[1]['accession_uniprot']
        print(uniprot_number_for_seq_in_row_1)
    '''   
    global dict1, output_dict, reader
    #dict1 = {}
    output_dict = {}
    with open(csvfile, mode='r') as infile:
        reader = csv.reader(infile)
        rownumber = 0    
        for row in reader:
            dict1 = {}
            # if the input doesn't have a header, simply use the column mumber as the dictionary key
            if fieldlist == 'all':
                cellnumber = 0
                for cell in row:
                    dict1[cellnumber] = cell           
                    cellnumber += 1
                #for the nested dictionary, the row number is used as the key
                output_dict[rownumber] = dict1
            else:
                #create a dictionary from only the required fields                
                if rownumber == 0:        
                    header = row
                else:
                    cellnumber = 0
                    for cell in row:
                        key = header[cellnumber]                        
                        dict1[key] = cell           
                        cellnumber += 1
                    selected_dict = create_new_dict_with_only_selected_keys(dict1, fieldlist)           
                    output_dict[rownumber] = selected_dict                     

            rownumber += 1
    return output_dict

def convert_stringlist_to_list(input_string):
    '''
    Will convert the string "['Eukaryota', 'Metazoa', 'Chordata']" to a list.
    '''
    string1 = input_string.strip("'[]")
    list1 = string1.split("', '")
    return list1

def create_new_dict_with_only_selected_keys(inputdict, keylist):
    global output_dict3
    for key in keylist:
        output_dict3 = { key: inputdict[key] for key in keylist }
    #print(output_dict3)
    return output_dict3
#    for key in keylist:
#        try:        
#            output_dict = { key: inputdict[key] for key in keylist }
#        except KeyError:
#            pass
#    return output_dict

def convert_string_to_boolean_value(boolean_string):
    """
       Converts 'something' to boolean. Raises exception for invalid formats
           Possible True  values: 1, True, "1", "TRue", "yes", "y", "t"
           Possible False values: 0, False, None, [], {}, "", "0", "faLse", "no", "n", "f", 0.0, ...
    """
    if str(boolean_string).lower() in ("yes", "y", "true",  "t", "1", "ja"): return True
    if str(boolean_string).lower() in ("no",  "n", "false", "f", "0", "0.0", "", "none", "nein"): return False # can also add "[]", "{}" for empty lists if desired
    raise Exception('Invalid value for boolean conversion: ' + str(boolean_string))


def setup_error_logging(logfile):
    import os
    import logging
    import logging.config
    import sys
    import json
    #import platform
    
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
    with open(logfile, 'w') as f:
        pass
    
    #clear any previous logging handlers that might have been previously run in the console
    logging.getLogger('').handlers = [] 
    #load the logging settings from the modified json string
    logging.config.dictConfig(config)
    
    #write system settings to logfile
    logging.warning('LOGGING LEVEL: %s' % config["loggers"]["simpleExample"]["level"])
    #logging.critical('Example of critical-level error. Current logging settings are level %s. At level DEBUG this logfile should also contain examples of WARNING and INFO level reports.' % config['handlers']['console']['level'])
    #logging.warning('Example of warning-level error')
    #logging.info('Example of info-level error\n')
    logging.info('SYSTEM INFORMATION')
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
    logging.warning(system_settings_dict)
    #save the logging settings in the logfile    
    #logging.info('LOGGING SETTINGS FOR THIS RUN IN JSON FORMAT')
    #logging.info("%s\n" % config)
    
    #test error message reporting
    #logging.warning('LOGGING TEST:')
    #try:
    #    open('/path/to/does/not/exist', 'rb')
    #except (SystemExit, KeyboardInterrupt):
    #    raise
    #except Exception:
    #    logging.error('Failed to open file', exc_info=True)
    logging.warning('LOGGING SETUP IS SUCCESSFUL\n\nLOG INFORMATION STARTS HERE:\n')


def save_structured_array_to_csv(array1, file1):
    #save the column names in the structured array as a header
    header = [x[0] for x in array1.dtype.descr]
    save_list_as_row_in_csv(header, file1, 'w')
    
    #save the rest of the data in the csv file
    with open(file1, 'ab') as f:
        np.savetxt(f, array1, fmt='%s', delimiter=',', newline='\n', header='', footer='', comments='#')


def load_structured_array_from_csv(file2, dtype2):
    '''
    The data type(dtype) for each column can be is listed in this format:
    dtype_for_my_array = [('number', '<i4'), ('query_name', '<U30')]
    Note that if the list of data types is shorter than the number of columns 
    in the csv, numpy will simply ignore the rest of the data. Note also that
    the format needs to match the data closely, or you will have an empty array, 
    or 'nan' values. In python 3.3, you should use U for the unicode format.
    '''
    
    loaded_array = np.genfromtxt(file2, delimiter=',', 
                               dtype=dtype2,
                               comments='#',
                               skiprows = 1)
    return loaded_array


class HardDriveSpaceException(Exception):
    def __init__(self, value):
        self.parameter = value
    def __str__(self):
        return repr(self.parameter)


def get_free_space(folder, format="MB"):
    """ 
        Return folder/drive free space 
    """
    fConstants = {"GB": 1073741824,
                  "MB": 1048576,
                  "KB": 1024,
                  "B": 1
                  }
    if platform.system() == 'Windows':
        free_bytes = ctypes.c_ulonglong(0)
        ctypes.windll.kernel32.GetDiskFreeSpaceExW(ctypes.c_wchar_p(folder), None, None, ctypes.pointer(free_bytes))
        return (int(free_bytes.value/fConstants[format.upper()]), format)
    else:
        return (int(os.statvfs(folder).f_bfree*os.statvfs(folder).f_bsize/fConstants[format.upper()]), format)


'''****************************Small Bioinformatics Functions*****************************************'''

class SmallBioinformaticsFunctions(object):
    pass


def get_phobius_TMD_region(feature_table_root):
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

def get_TMD_seq_of_query_from_homologue_root(root, input_dict_with_query_data):
    for parameters in root[0][0][0][0].iter('parameters'):
        for sequences in parameters.iter('sequences'):        
            query_seq = sequences[0][0][0].text
            TMDstart = input_dict_with_query_data['phobius_TMD_start']
            TMDend = input_dict_with_query_data['phobius_TMD_end']
            TMD_seq = query_seq[TMDstart - 1:TMDend]   
            TMD_seq_uppercase = TMD_seq.upper()
    return TMD_seq_uppercase

def print_query_details_from_homologue_XML(root, query_details_dict):
    global counter_XML_to_CSV
    counter_XML_to_CSV = 1
    for parameters in root[0][0][0][0].iter('parameters'):
        SIMAP_input_seq_details_dict = parameters[0][0].attrib
        query_details_dict['SIMAP_input_seq_details_dict'] = SIMAP_input_seq_details_dict
        for filter in parameters.iter('filter'):
            SIMAP_filter_string = filter.text
        query_details_dict['SIMAP_filter_string'] = SIMAP_filter_string
        for resultSpecification in parameters.iter('resultSpecification'):
            SIMAP_resultSpecification_dict = resultSpecification.attrib
        SIMAP_input_seq_details_dict['SIMAP_resultSpecification_dict'] = SIMAP_resultSpecification_dict 
        for databases in parameters.iter('databases'):
            database_details_dict = databases[0].attrib
        SIMAP_input_seq_details_dict['database_details_dict'] = database_details_dict
        simap_version = root[0][0][0][0][0].attrib['version']
        SIMAP_input_seq_details_dict['simap_version'] = simap_version
        SIMAP_total_hits = int(root[0][0][0][1][0].attrib['total'])
        SIMAP_input_seq_details_dict['SIMAP_total_hits'] = SIMAP_total_hits
#    if counter_XML_to_CSV == 1:
#        logging.warning('This SIMAP parser was developed for SIMAP version 4.0. This XML file is SIMAP version %s.' % simap_version)
    if simap_version != "4.0":
        logging.warning('WARNING! Your XML file is simap version %s, however this SIMAP parser was developed for SIMAP version 4.0.' % simap_version)
    counter_XML_to_CSV += 1    
    logging.info('%s homologous sequences analysed' % SIMAP_total_hits)
    return query_details_dict
