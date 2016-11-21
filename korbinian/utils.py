#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Utilities file containing useful functions.
More recent functions are at the top.
Authors: Mark Teese, Rimma Jenske
Created on Fri Nov  8 15:45:06 2013
"""
import ast
import csv
import ctypes
import errno
import logging
import os
import pickle
import platform
import re as re
import subprocess
import sys
import tarfile
import threading
import time
import zipfile
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio.SeqUtils.ProtParam import ProteinAnalysis

def aaa(df_or_series):
    """ Function for use in debugging.
    Saves pandas Series or Dataframes to a user-defined csv file.
    """
     # convert any series to dataframe
    if isinstance(df_or_series, pd.Series):
        df_or_series = df_or_series.to_frame()
    csv_out = r"D:\data\000_aaa_temp_df_out.csv"
    df_or_series.to_csv(csv_out, sep=",", quoting=csv.QUOTE_NONNUMERIC)

def import_amino_acid_substitution_matrices():
    """
    imports several aa sub matrices from Bio.SubsMat.MatrixInfo
    """

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
    
    for i in range(1, max_num_TMDs.astype(np.int64) + 1):
        TM = 'TM%02d_AAIMON_ratio_mean' % i
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
#    m = re.search(TMD_regex_ss, x)
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
    #create a line graph rather than a bar graph for the AASMON (ident + similarity)
    linecontainer_AASMON_mean = axarr[row_nr, col_nr].plot(centre_of_bar_in_x_axis, freq_counts_S, color=color,
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
    #use numpy to create a histogram
    freq_counts_I, bin_array_I = np.histogram(hist_data, bins = binlist)
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
def slice_uniprot_SP_seg(x, SP):
    if x['SP01_end'] != '?':
        return x['full_seq'][int(x['SP01_start'] - 1):int(x['SP01_end'])]

def slice_uniprot_TMD_seq(x, TMD):
   return x['full_seq'][int(x['%s_start'%TMD] - 1):int(x['%s_end'%TMD])]

def slice_uniprot_TMD_plus_surr_seq(x, TMD):
    return x['full_seq'][int(x['%s_start_plus_surr'%TMD] - 1):int(x['%s_end_plus_surr'%TMD])]

#create small throwaway functions to slice all sequences in dataframe simultaneously
#slice_SW_query_TMD_seq = lambda x: x['query_align_seq'][int(x['%s_start_in_SW_alignment'%TMD]):int(x['%s_end_in_SW_alignment'%TMD])]
#slice_SW_markup_TMD = lambda x: x['align_markup_seq'][int(x['%s_start_in_SW_alignment'%TMD]):int(x['%s_end_in_SW_alignment'%TMD])]
#slice_SW_match_TMD_seq = lambda x: x['match_align_seq'][int(x['%s_start_in_SW_alignment'%TMD]):int(x['%s_end_in_SW_alignment'%TMD])]
#slice_SW_query_TMD_seq_plus_surr = lambda x: x['query_align_seq'][int(x['%s_start_in_SW_alignment_plus_surr'%TMD]):int(x['%s_end_in_SW_alignment_plus_surr'%TMD])]
#slice_SW_markup_TMD_plus_surr = lambda x: x['align_markup_seq'][int(x['%s_start_in_SW_alignment_plus_surr'%TMD]):int(x['%s_end_in_SW_alignment_plus_surr'%TMD])]
#slice_SW_match_TMD_seq_plus_surr = lambda x: x['match_align_seq'][int(x['%s_start_in_SW_alignment_plus_surr'%TMD]):int(x['%s_end_in_SW_alignment_plus_surr'%TMD])]

#create small throwaway functions to slice all sequences in dataframe simultaneously
# def slice_SW_query_TMD_seq(x, TMD):
#     return x['query_align_seq'][int(x['%s_start_in_SW_alignment'%TMD]):int(x['%s_end_in_SW_alignment'%TMD])]
# def slice_SW_markup_TMD(x, TMD):
#     return x['align_markup_seq'][int(x['%s_start_in_SW_alignment'%TMD]):int(x['%s_end_in_SW_alignment'%TMD])]
# def slice_SW_match_TMD_seq(x, TMD):
#     return x['match_align_seq'][int(x['%s_start_in_SW_alignment'%TMD]):int(x['%s_end_in_SW_alignment'%TMD])]
# def slice_SW_query_TMD_seq_plus_surr(x, TMD):
#     return x['query_align_seq'][int(x['%s_start_in_SW_alignment_plus_surr'%TMD]):int(x['%s_end_in_SW_alignment_plus_surr'%TMD])]
# def slice_SW_markup_TMD_plus_surr(x, TMD):
#     return x['align_markup_seq'][int(x['%s_start_in_SW_alignment_plus_surr'%TMD]):int(x['%s_end_in_SW_alignment_plus_surr'%TMD])]
# def slice_SW_match_TMD_seq_plus_surr(x, TMD):
#     return x['match_align_seq'][int(x['%s_start_in_SW_alignment_plus_surr'%TMD]):int(x['%s_end_in_SW_alignment_plus_surr'%TMD])]
# def find_indices_longer_than_prot_seq(df, TMD):
#     return df['%s_end_plus_surr'%TMD] > df['seqlen']


def slice_SW_query_TMD_seq(x, TMD):
    return x['query_align_seq'][int(x['%s_start_in_SW_alignment'%TMD]):int(x['%s_end_in_SW_alignment'%TMD])]
def slice_SW_markup_TMD(x, TMD):
    return x['align_markup_seq'][int(x['%s_start_in_SW_alignment'%TMD]):int(x['%s_end_in_SW_alignment'%TMD])]
def slice_SW_match_TMD_seq(x, TMD):
    return x['match_align_seq'][int(x['%s_start_in_SW_alignment'%TMD]):int(x['%s_end_in_SW_alignment'%TMD])]
def slice_SW_query_TMD_seq_plus_surr(x, TMD):
    return x['query_align_seq'][int(x['%s_start_in_SW_alignment_plus_surr'%TMD]):int(x['%s_end_in_SW_alignment_plus_surr'%TMD])]
def slice_SW_markup_TMD_plus_surr(x, TMD):
    return x['align_markup_seq'][int(x['%s_start_in_SW_alignment_plus_surr'%TMD]):int(x['%s_end_in_SW_alignment_plus_surr'%TMD])]
def slice_SW_match_TMD_seq_plus_surr(x, TMD):
    return x['match_align_seq'][int(x['%s_start_in_SW_alignment_plus_surr'%TMD]):int(x['%s_end_in_SW_alignment_plus_surr'%TMD])]


def create_indextuple_nonTMD_last(x):
    ''' Joins two columns into a tuple. Used to create the last tuple of the nonTMD region in the sequence.
    '''
    return (int(x['nonTMD_index_tuple_last0']), int(x['nonTMD_index_tuple_last1']))

def slice_with_listlike(string, tup, start=0, end=1):
    '''A function to slice a single string, taking the start and stop indices from a tuple
    '''
    return string[int(tup[start]):int(tup[end])]

def slice_with_nested_tuple(string, nested_tuple):
    '''A function to slice a sequence multiple times, using the indices from nested tuples
    '''
    #convert nested tuple from string to tuple
    nested_tuple = ast.literal_eval(nested_tuple)
    #for each tuple, slice the input string. Make a list of all the sliced strings. Join list with no gaps
    return ''.join([slice_with_listlike(string, tup) for tup in nested_tuple])


def get_start_and_end_of_TMD_in_query(x, regex_string):
    '''
    Returns a tuple containing (bool, start, stop) showing the location of a regex pattern
    in the target string.
    To be used with Pandas Dataframes. The aim is to conduct the computationally expensive regex
    search only once to obtain both the start and stop indices.
    Variables:
    x = target sequence (e.g. unicode string, DNA, Protein Seq)
    TMD_regex_ss = regex search string
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


def check_SIMAP_tarfile(SIMAP_tar, ft_xml_path, homol_xml_path, acc, logging, delete_corrupt=False):
    ''' Checs
    Checks the tarball that contains the SIMAP output. 
    Looks to see if the tarball exists, if it is corrupted, if it contains the feature table and homologues from simap.
    '''

    ft_xml_filename = os.path.basename(ft_xml_path)
    homol_xml_filename = os.path.basename(homol_xml_path)

    if os.path.isfile(ft_xml_path):
        ft_XML_exists = True
    else:
        ft_XML_exists = False
    if os.path.isfile(homol_xml_path):
        homol_XML_exists = True
    else:
        homol_XML_exists = False
    if os.path.isfile(SIMAP_tar):
        SIMAP_tar_exists = True
    else:
        SIMAP_tar_exists = False
    # check if feature table and homologues XML files are in the simap tarball
    ft_in_tar = False
    homol_in_tar = False
    if SIMAP_tar_exists:
        try:
            with tarfile.open(SIMAP_tar, mode = 'r:gz') as tar:
                if ft_xml_filename in [tarinfo.name for tarinfo in tar]:
                    ft_in_tar = True
                if homol_xml_filename in [tarinfo.name for tarinfo in tar]:
                    homol_in_tar = True
        except EOFError:
            if delete_corrupt == True:
                logging.info("{} SIMAP_tar seems corrupt, will be deleted.".format(acc))
                os.remove(SIMAP_tar)
            else:
                SIMAP_tar_exists = False
                logging.info("{} SIMAP_tar seems corrupt.".format(acc))
    return ft_XML_exists, homol_XML_exists, SIMAP_tar_exists, ft_in_tar, homol_in_tar

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

def score_pairwise_gapless(seq1, seq2, matrix):
    '''
    Calculates a score between two aligned sequences without gaps, based on matrix applied.

    A, B are each paired amino acid in the pairwise alignment

    Usage:
    # import various matrices from biopython (many available)
    from Bio.SubsMat.MatrixInfo import ident, blosum62, pam120, levin
    a = "ACGEGGGFFFCCC"
    b = "ACFGGGTFFTCCC"
    c = score_pairwise_gapless(a, b, blosum62_matrix)
    score = sum(c)
    '''
    for A, B in zip(seq1, seq2):
        pair = (A, B)
        if pair not in matrix:
            yield matrix[(tuple(reversed(pair)))]
        else:
            yield matrix[pair]

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
            #logging.info('Thread started')
            self.process = subprocess.Popen(self.cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            #self.process.communicate()
            stdout, stderr = self.process.communicate() # from http://stackoverflow.com/questions/14366352/how-to-capture-information-from-executable-jar-in-python
            # Thus far, SIMAP has only ever given java faults, never java output. Don't bother showing.
            # if the console prints anything longer than 5 characters, log it
            if len(stderr.decode("utf-8")) > 5:
                logging.warning('FAULTS: %s' % stderr.decode("utf-8"))
            #logging.info('Thread finished')

        thread = threading.Thread(target=target)
        thread.start()

        thread.join(timeout)
        if thread.is_alive():
            logging.info('Terminating process')
            self.process.terminate()
            thread.join()
        # simply returns 0 every time it works. Waste of logging space! :)
        #logging.info(self.process.returncode)

def run_command(command):
    #this stopped working for some reason. Did I mess up a path variable?    
    p = subprocess.Popen(command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT)
    return iter(p.stdout.readline, b'')

def sleep_x_seconds(x, print_stuff=True):
    # sleep for several seconds to not overload a server, for example
    if print_stuff == True:
        sys.stdout.write("sleeping ")
    for i in range(x):
        time.sleep(1)
        if print_stuff == True:
            sys.stdout.write(" .")
            sys.stdout.flush()
    if print_stuff == True:
        sys.stdout.write(' .')

def sleep_x_hours(x):
    """Sleeps for a certain number of hours. Prints a dot each hour.

    Parameters
    ----------
    x : int
        Number of hours to sleep

    """
    #sleep for 30 seconds to not overload the server
    sys.stdout.write("sleeping .")
    sys.stdout.flush()
    for i in range(x):
        time.sleep(3600)
        sys.stdout.write(" .")
        sys.stdout.flush()
    sys.stdout.write(' .\n')

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

# Function to store Dataframes in an Excelfile; Converting lists etc. into strings
def df_to_excel(dataframe, path):
    temp_df = pd.DataFrame(dataframe, dtype=str)
    temp_df.to_excel(path)

# Creates subfolders, based on first to letters of files; Distributes files to these subfolders
def distribute_files_to_subfolders(path):
    list_of_filenames = os.listdir(path)
    for n in list(set(n[0:2] for n in list_of_filenames)):
        os.makedirs(path + "/%s" % n)
    for n in list_of_filenames:
        directory = path + "/" + n
        new_directory = path + "/" + n[0:2] + "/" + n
        os.replace(directory, new_directory)

# Moves files from a subfolder to the root folder;
# Length of subfoldername = amount of letters; necessary for path-identification
def move_files_from_subfolder_to_folder(path_of_subfolder, length_of_subfoldername):
    for n in os.listdir(path_of_subfolder):
        directory = path_of_subfolder + "/" + n
        new_directory = str(path_of_subfolder)[:-length_of_subfoldername] + "/" + n
        os.replace(directory, new_directory)

# this functions works exclusivly with dataframes; query, start, stop, and new_name refer to columns
# and as well, it does not work yet, still working on it
def slicing(df, columname_of_sequence, start, stop, columnname_for_spliced_sequence):
    # df was missing from this function!
    for n in df["%s" % columname_of_sequence]:
        df["%s" % columnname_for_spliced_sequence] = n[start, stop]

def getting_list_of_indices_of_M_in_a_string(string):
    ind_list = [i for i, element in enumerate(string) if element == "M"]  # find(Topo_data)
    return ind_list

def getting_list_of_gapindices_of_string(string):
    gap_list = [i for i, element in enumerate(string) if element == "-"]  # find(Topo_data)
    return gap_list

def getting_list_of_gapindices_of_string(string):
    gap_list = [i for i, element in enumerate(string) if element == "-"]  # find(Topo_data)
    return gap_list

def create_border_list(index_list):
    m_borders = []
    m_borders.append(index_list[0])
    for n in range(0, len(index_list) - 1):
        if index_list[n] + 1 != index_list[n + 1]:
            m_borders.append(index_list[n] + 1)
            m_borders.append(index_list[n + 1])
    m_borders.append(index_list[-1] + 1)
    return m_borders

def create_list_of_TMDs(amount_of_TMDs):
    list_of_TMDs = []
    for n in range(1, int(amount_of_TMDs) + 1):
        list_of_TMDs.append("TM%.2d" % n)
    return list_of_TMDs

def isEven(number):
    return number % 2 == 0

def isOdd(number):
    return number % 2 != 0

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
            frequency = frequency + 1
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

def get_end_juxta_before_TMD(x, input_TMD):
    TM_int = int(input_TMD[2:])
    if input_TMD == "TM01":
        x['end_juxta_before_%s_in_query' % input_TMD] = np.where(x['%s_start_in_SW_alignment' % input_TMD] == 0, 0,
                                                                 x['%s_start_in_SW_alignment' % input_TMD] - 1)
    else:
        x["end_juxta_before_%s_in_query" % input_TMD] = x["%s_start_in_SW_alignment" % input_TMD] - 1

def get_end_juxta_after_TMD(x, input_TMD, list_of_tmds):
    # list_of_tmds was missing from this function! added by MT 20.07.2016
    # this function contained dfs instead of x! added by MT 20.07.2016
    TM_int = int(input_TMD[2:])
    last_TMD = list_of_tmds[-1]
    if input_TMD == last_TMD:
        x["end_juxta_after_%s" % input_TMD] = np.where(
            (x["%s_end_in_SW_alignment"] + 30) < x["len_query_aligment_sequence"], x["%s_end_in_SW_alignment"] + 30,
            x["len_query_aligment_sequence"])
    else:
        x["end_juxta_after_%s" % input_TMD] = x["%s_end_in_SW_alignment" % input_TMD] + (
        (x["TM%.2d_start_in_SW_alignment" % (TM_int + 1)] - x["%s_end_in_SW_alignment" % input_TMD]) / 2).apply(
            lambda x: int(x) if not np.isnan(x) else np.nan)

        # else:
        #     x["end_juxta_after_%s" % input_TMD] = dfs["%s_end_in_SW_alignment" % input_TMD] + ((dfs["TM%.2d_start_in_SW_alignment" % (TM_int + 1)] - dfs["%s_end_in_SW_alignment" % input_TMD]) / 2).apply(    lambda x: int(x) if not np.isnan(x) else np.nan)


def get_start_and_end_of_TMD_in_query(x, TMD_regex_ss):
    '''
    define function to obtain regex output (start, stop, etc) as a tuple
    the functions were originally in MAIN, as I am not sure how to apply **args and **kwargs in functions that apply to pandas
    note that TMD_regex_ss is therefore a global variable
    '''
    m = re.search(TMD_regex_ss, x)
    if m:
        # if the tmd is in the query, return True, start, stop
        return [bool(m), m.start(), m.end()]
    else:
        # if the tmd is not in the query, return False, 0, 0
        return np.nan

def slice_juxta_before_TMD_in_query(x, TMD):
    return x['query_align_seq'][int(x['start_juxta_before_%s'%TMD]):int(x['end_juxta_before_%s'%TMD])]

def slice_juxta_after_TMD_in_query(x, TMD):
    return x['query_align_seq'][int(x['start_juxta_after_%s'%TMD]):int(x['end_juxta_after_%s'%TMD])]

def slice_juxta_before_TMD_in_match(x, TMD):
    return x['match_align_seq'][int(x['start_juxta_before_%s'%TMD]):int(x['end_juxta_before_%s'%TMD])]

def slice_juxta_after_TMD_in_match(x, TMD):
    return x['match_align_seq'][int(x['start_juxta_after_%s'%TMD]):int(x['end_juxta_after_%s'%TMD])]

def find_last_TMD(dfs):
    # dfs was missing from input, added by MT 20.07.2016
    for n in range(1, 24):
        if isNaN(dfs['TM%.2d_start_in_SW_alignment' % n]):
            last_TMD = n

def convert_truelike_to_bool(input_item, convert_int=False, convert_float=False, convert_nontrue=False):
    """Converts true-like values ("true", 1, True", "WAHR", etc) to python boolean True.

    Parameters
    ----------
    input_item : string or int
        Item to be converted to bool (e.g. "true", 1, "WAHR" or the equivalent in several languagues)
    convert_float: bool
        Convert floats to bool.
        If True, "1.0" will be converted to True
    convert_nontrue : bool
        If True, the output for input_item not recognised as "True" will be False.
        If True, the output for input_item not recognised as "True" will be the original input_item.

    Returns
    -------
    return_value : True, or input_item
        If input_item is True-like, returns python bool True. Otherwise, returns the input_item.

    Usage
    -----
    # convert a single value or string
    convert_truelike_to_bool("true")
    # convert a column in a pandas DataFrame
    df["column_name"] = df["column_name"].apply(convert_truelike_to_bool)
    """
    list_True_items = [True, 'True', "true","TRUE","T","t",'wahr', 'WAHR', 'prawdziwy', 'verdadeiro', 'sann', 'istinit',
                       'veritable', 'Pravda', 'sandt', 'vrai', 'igaz', 'veru', 'verdadero', 'sant', 'gwir', 'PRAWDZIWY',
                       'VERDADEIRO', 'SANN', 'ISTINIT', 'VERITABLE', 'PRAVDA', 'SANDT', 'VRAI', 'IGAZ', 'VERU',
                       'VERDADERO', 'SANT', 'GWIR', 'bloody oath', 'BLOODY OATH', 'nu', 'NU','damn right','DAMN RIGHT']

    # if you want to accept 1 or 1.0 as a true value, add it to the list
    if convert_int:
        list_True_items += ["1"]
    if convert_float:
        list_True_items += [1.0, "1.0"]
    # check if the user input string is in the list_True_items
    input_item_is_true = input_item in list_True_items
    # if you want to convert non-True values to "False", then nontrue_return_value = False
    if convert_nontrue:
        nontrue_return_value = False
    else:
        # otherwise, for strings not in the True list, the original string will be returned
        nontrue_return_value = input_item
    # return True if the input item is in the list. If not, return either False, or the original input_item
    return_value = input_item_is_true if input_item_is_true == True else nontrue_return_value
    # special case: decide if 1 as an integer is True or 1
    if input_item == 1:
        if convert_int == True:
            return_value = True
        else:
            return_value = 1
    return return_value

def convert_falselike_to_bool(input_item, convert_int=False, convert_float=False):
    """Converts false-like values ("false", 0, FALSE", "FALSCH", etc) to python boolean False.

    Parameters
    ----------
    input_item : string or int
        Item to be converted to bool (e.g. "FALSE", 0, "FALSCH" or the equivalent in several languagues)
    convert_float: bool
        Convert floats to bool.
        If True, "0.0" will be converted to True

    Returns
    -------
    return_value : False, or input_item
        If input_item is False-like, returns python bool False. Otherwise, returns the input_item.

    Usage
    -----
    # convert a single value or string
    convert_falselike_to_bool("false")
    # convert a column in a pandas DataFrame
    df["column_name"] = df["column_name"].apply(convert_falselike_to_bool)
    """
    list_False_items = [False, "False", "false", "FALSE", "F", "f", "falsch", "FALSCH", "valse", "lažna", "fals",
                        "NEPRAVDA", "falsk", "vals", "faux", "pa vre", "tsis tseeb", "hamis", "palsu", "uongo", "ngeb",
                        "viltus", "klaidinga", "falz", "falso", "USANN", "wartosc false", "falošné", "falskt", "yanlis",
                        "sai", "ffug", "VALSE", "LAŽNA", "FALS", "FALSK", "VALS", "FAUX", "PA VRE", "TSIS TSEEB",
                        "HAMIS", "PALSU", "UONGO", "NGEB", "VILTUS", "KLAIDINGA", "FALZ", "FALSO", "WARTOSC FALSE",
                        "FALOŠNÉ", "FALSKT", "YANLIS", "SAI", "FFUG"]

    # if you want to accept 0 or 0.0 as a false value, add it to the list
    if convert_int:
        list_False_items += [0, "0"]
    if convert_float:
        list_False_items += [0.0,"0.0"]
    # return boolean False if the input item is in the list. If not, return the original input_item
    return_value = False if input_item in list_False_items else input_item

    return return_value

def calc_hydrophob(seq):
    """ Calculates the average hydrophobicity of a sequence according to the Hessa biological scale.

    Hessa T, Kim H, Bihlmaier K, Lundin C, Boekel J, Andersson H, Nilsson I, White SH, von Heijne G. Nature. 2005 Jan 27;433(7024):377-81

    The Hessa scale has been calculated empirically, using the glycosylation assay of TMD insertion.
    Negative values indicate hydrophobic amino acids with favourable membrane insertion.

    Other hydrophobicity scales are in the settings folder. They can be generated as follows.
    hydrophob_scale_path = r"D:\korbinian\korbinian\settings\hydrophobicity_scales.xlsx"
    df_hs = pd.read_excel(hydrophob_scale_path, skiprows=2)
    df_hs.set_index("1aa", inplace=True)
    dict_hs = df_hs.Hessa.to_dict()
    hessa_scale = np.array([value for (key, value) in sorted(dict_hs.items())])
    ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K',
     'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V',
     'W', 'Y']

    Parameters:
    -----------
    seq : string
        Sequence to be analysed. Gaps (-) and unknown amino acids (x) should be ignored.

    Returns:
    --------
    mean hydrophobicity value for the sequence entered

    Usage:
    ------
    from korbinian.utils import calc_hydrophob
    # for a single sequence
    s = "SAESVGEVYIKSTETGQYLAG"
    calc_hydrophob(s)
    # for a series of sequences
    TMD_ser = df2.TM01_SW_match_seq.dropna()
    hydro = TMD_ser.apply(lambda x : calc_hydrophob(x))

    Notes:
    ------
    %timeit results:
    for a 20aa seq: 136 µs per loop
    for a pandas series with 852 tmds: 118 ms per loop
    """
    # hydrophobicity scale
    hessa_scale = np.array([0.11, -0.13, 3.49, 2.68, -0.32, 0.74, 2.06, -0.6, 2.71,
                            -0.55, -0.1, 2.05, 2.23, 2.36, 2.58, 0.84, 0.52, -0.31,
                            0.3, 0.68])
    # convert to biopython analysis object
    analysed_seq = ProteinAnalysis(seq)
    # biopython count_amino_acids returns a dictionary.
    aa_counts_dict = analysed_seq.count_amino_acids()
    # convert dictionary to array, sorted by aa
    aa_counts_arr = np.array([value for (key, value) in sorted(aa_counts_dict.items())])
    multiplied = aa_counts_arr * hessa_scale
    return multiplied.sum()

def make_sure_path_exists(path, isfile=False):
    """ If path to directory or folder doesn't exist, creates the necessary folders.

    Parameters
    ----------
    path : str
        Path to desired directory or file.
    isfile :
        If True, the path is to a file, and the subfolder will be created if necessary
    """
    if isfile:
        path = os.path.dirname(path)
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

def save_df_to_csv_zip(df,out_zipfile,open_method="w"):
    """ Save a pandas dataframe to a zipped csv (.csv.zip)

    Parameters
    ----------
    df : pd.DataFrame
        pandas dataframe to be saved
    out_zipfile : filepath
        Path of zipfile to be created or added to.
    open_method : str
        Method to open file, e.g. "w" for write mode

    Saved Files and Figures
    -----------------------
    out_zipfile : zipfile

    Note
    -------
    Much faster than saving to excel.
    """
    # create a temporary csv file path, equivalent to .csv.zip minus the .zip
    temp_csv = out_zipfile[:-4]
    # extract filename
    filename = os.path.basename(temp_csv)
    #save
    df.to_csv(temp_csv, quoting=csv.QUOTE_NONNUMERIC)
    # either create new zip and add ("w"), or open existing zip and add "a"
    with zipfile.ZipFile(out_zipfile,open_method, zipfile.ZIP_DEFLATED) as zipout:
        zipout.write(temp_csv, arcname=filename)
    # delete temporary csv file
    os.remove(temp_csv)

def open_df_from_csv_zip(in_zipfile, filename=None):
    """ Opens a pandas dataframe that is saved as a zipped csv (.csv.zip)

    Parameters
    ----------
    in_zipfile : str
        Path to zip file
    filename : str
        Filename. Default is "None", which will result in the opening of the first file in the zipfile.

    Returns
    -------
    df : pd.DataFrame
        pandas Dataframe

    Note
    -------
    Much faster than reading from excel.
    """
    if os.path.isfile(in_zipfile):
        with zipfile.ZipFile(in_zipfile, "r", zipfile.ZIP_DEFLATED) as openzip:
            if filename == None:
                # if a filename is not given, open the first file in the list
                filename = openzip.namelist()[0]
            # open the file
            csv_file_handle = openzip.open(filename)
            # read as pandas dataframe
            df = pd.read_csv(csv_file_handle, sep=",", quoting=csv.QUOTE_NONNUMERIC, index_col=0)
    else:
        raise FileNotFoundError("{} not found".format(in_zipfile))
    return df

def open_df_from_pickle_zip(in_zipfile, filename=None, delete_corrupt=False):
    """ Opens a pandas dataframe that is saved as a zipped pickle file (.pickle.zip)

    Parameters
    ----------
    in_zipfile : str
        Path to zip file
    filename : str
        Filename inside zipfile to open. Default is "None", which will result in the opening of the first .pickle file in the zipfile.

    Returns
    -------
    df_loaded : pd.DataFrame
        Output pandas Dataframe.

    Note
    -------
    Much faster than reading from excel.
    """
    # create bool deciding whether zip file will be deleted
    deletezip = False
    if os.path.isfile(in_zipfile):
        try:
            with zipfile.ZipFile(in_zipfile, "r", zipfile.ZIP_DEFLATED) as openzip:
                filenamelist = openzip.namelist()
                if filename == None:
                    # if a filename is not given, open the first file in the list
                    for file_in_zip in filenamelist:
                        if file_in_zip[-7:] == ".pickle":
                            filename = file_in_zip
                            # pickle is found, stop searching
                            break
                if filename is not None:
                    # if a filename is available, check if the file is in the zip
                    if  filename in filenamelist:
                        csv_file_handle = openzip.open(filename)
                        # read as pandas dataframe
                        df_loaded = pickle.load(csv_file_handle)
                        # make sure that the pickled object was REALLY a pandas object, and not some other python datatype that was pickled.
                        assert isinstance(df_loaded, (pd.Series, pd.DataFrame))
                    else:
                        # the desired file is not in the zip. Either delete the zip, or return an empty dataframe.
                        if delete_corrupt == True:
                            deletezip = True
                        else:
                            df_loaded = pd.DataFrame()
                else:
                    # if the zipfile doesn't contain ANY pickle files, something is seriously wrong. Either delete or raise Error.
                    if delete_corrupt == True:
                        deletezip = True
                    else:
                        raise FileNotFoundError("{} does not contain a valid pickle file".format(in_zipfile))
        except zipfile.BadZipFile:
            # the desired file is not in the zip. Either delete the zip, or return an empty dataframe.
            if delete_corrupt == True:
                deletezip = True
            else:
                df_loaded = pd.DataFrame()
    else:
        raise FileNotFoundError("{} not found".format(in_zipfile))
    if deletezip:
        logging.info("{} does not contain expected pickle file {}. File is old or damaged, and has been deleted".format(in_zipfile, filename))
        os.remove(in_zipfile)
        df_loaded = pd.DataFrame()
    return df_loaded

def create_colour_lists():
    '''
    Converts several lists of rgb colours to the python format (normalized to between 0 and 1)
    Returns a dictionary that contains dictionaries of palettes with named colours (eg. TUM blues)
    and also lists of unnamed colours (e.g. tableau20)
    (copied from tlabtools 2016.08.08)
    '''
    import numpy as np
    import matplotlib.colors as colors
    output_dict = {}

    matplotlib_150 = list(colors.cnames.values())
    output_dict['matplotlib_150'] = matplotlib_150

    #define colour dictionaries. TUM colours are based on the style guide.
    colour_dicts = {
                    'TUM_colours' : {
                                    'TUMBlue':(34,99,169),
                                    'TUM1':(100,160,200),
                                    'TUM2':(1,51,89),
                                    'TUM3':(42,110,177),
                                    'TUM4':(153,198,231),
                                    'TUM5':(0,82,147)
                                    },
                    'TUM_accents' : {
                                    'green':(162,183,0),
                                    'orange':(227,114,34),
                                    'ivory':(218,215,203),
                                    }
                    }

    #convert the nested dicts to python 0 to 1 format
    for c_dict in colour_dicts:
        for c in colour_dicts[c_dict]:
            #define r, g, b as ints
            r, g, b = colour_dicts[c_dict][c]
            #normalise r, g, b and add to dict
            colour_dicts[c_dict][c] = (r / 255., g / 255., b / 255.)
        #add normalised colours to output dictionary
        output_dict[c_dict] = colour_dicts[c_dict]

    #define colour lists
    colour_lists = {
                    'tableau20' : [
                                 (31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),
                                 (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),
                                 (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),
                                 (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
                                 (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)
                                    ],
                    'tableau20blind' : [
                                         (0, 107, 164), (255, 128, 14), (171, 171, 171), (89, 89, 89),
                                         (95, 158, 209), (200, 82, 0), (137, 137, 137), (163, 200, 236),
                                         (255, 188, 121), (207, 207, 207)
                                          ]
                    }
    #normalise the colours for the colour lists
    for rgb_list in colour_lists:
        colour_array = np.array(colour_lists[rgb_list])/255.
        colour_array_tup = tuple(map(tuple,colour_array))
        colour_lists[rgb_list] = colour_array_tup
        #add normalised colours to output dictionary
        output_dict[rgb_list] = colour_lists[rgb_list]
    #create a mixed blue/grey colour list, with greys in decreasing darkness
    TUM_colours_list_with_greys = []
    grey = 0.7
    for c in colour_dicts['TUM_colours'].values():
        TUM_colours_list_with_greys.append('%0.2f' % grey)
        TUM_colours_list_with_greys.append(c)
        grey -= 0.1
    output_dict['TUM_colours_list_with_greys'] = TUM_colours_list_with_greys
    return output_dict

def savefig_if_necessary(savefig, fig, fig_nr, base_filepath, tight_layout = False, formats = ['png','pdf'], dpi = 400):
    '''
    Function to save figure with multiple subplots. (i.e., a canvas containing multiple figures)
    Designed to work with the function create_dict_organising_subplots(), which creates a bool object "savefig".
    Automatically names the figure based on the figure number (fig_nr), using a previously defined file path as a base.
    '''
    if savefig:
        if 'png' in formats:
            fig.savefig(base_filepath + '_%01d.png' % fig_nr, format='png', dpi=dpi)
        if 'pdf' in formats:
            fig.savefig(base_filepath + '_%01d.pdf' % fig_nr, format='pdf')
        #close any open figures
        plt.close('all')

def save_figure(s, fig, Fig_Nr, base_filepath, dpi = 400):
    if not os.path.exists(base_filepath):
        os.makedirs(base_filepath)
    if s['save_fig_to_png']:
        fig.savefig(os.path.join(base_filepath, '_figs') + '_%01d.png'  % Fig_Nr, format='png', dpi=dpi)
    if s['save_fig_to_pdf']:
        fig.savefig(os.path.join(base_filepath, '_figs') + '_%01d.pdf' % Fig_Nr, format='pdf')
    # close any open figures
    plt.close('all')

class Log_Only_To_Console(object):
    def __init__(self):
        pass
    def info(self, message):
        print(message)
    def warning(self, message):
        print(message)
    def critical(self, message):
        print(message)

def convert_summary_csv_to_input_list(s, pathdict, logging, list_excluded_acc=None):
    # open dataframe with list of proteins
    df = pd.read_csv(pathdict["list_summary_csv"], sep=",", quoting=csv.QUOTE_NONNUMERIC, index_col=0)
    # exclude any proteins where there is no list_of_TMDs
    df = df.loc[df['list_of_TMDs'].notnull()].loc[df['list_of_TMDs'] != 'nan']
    # add the accession
    df["acc"] = df.index

    if list_excluded_acc != None:
        # remove any excluded acc from index
        non_excluded = set(df.index) - set(list_excluded_acc)
        # redefine df, skipping excluded acc
        df = df.loc[non_excluded, :]

    # convert to dict
    df_as_dict = df.to_dict(orient="index")
    # convert values to list
    list_p = list(df_as_dict.values())

    for p in list_p:
        # print("in for loop", p["acc"])
        p["s"] = s
        p["pathdict"] = pathdict
        p["logging"] = logging

    return list_p


def get_list_not_in_homol_db(pathdict):
    not_in_homol_db = []
    if os.path.isfile(pathdict["acc_not_in_homol_db_txt"]):
        # Extracts accession numbers out of file
        with open(pathdict["acc_not_in_homol_db_txt"], "r") as source:
            for line in source:
                line = line.strip()
                not_in_homol_db.append(line)
    return not_in_homol_db

def get_list_failed_downloads(pathdict):
    acc_list_failed_downloads = []
    if os.path.isfile(pathdict["failed_downloads_txt"]):
        # Extracts accession numbers out of file
        with open(pathdict["failed_downloads_txt"], "r") as source:
            for line in source:
                line = line.strip()
                acc_list_failed_downloads.append(line)
    return acc_list_failed_downloads