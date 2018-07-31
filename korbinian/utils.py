#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Utilities file containing useful functions.
"""
#Switch matplotlib library for possible commandline execution
import matplotlib

# commented out as it causes errors in dependencies.
#matplotlib.use('Agg')

import ast
import csv
import ctypes
import errno
import glob
import inspect
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
import matplotlib.colors as colors
import numpy as np
import pandas as pd
from shutil import copyfile
from time import strftime
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from inspect import currentframe, getframeinfo, stack
import io
import urllib.request
from tqdm import tqdm

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
        TM = 'TM%02d_AAIMON_mean' % i
        full_list_TMDs.append(TM)
        #create new series with the data (each datapoint is the mean for all homologues,
        #for a single protein)
        hist_data_AAIMON_each_TM = df['TM%02d_AAIMON_mean' % i].dropna()
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
    df['last_TM_AAIMON_mean'] = np.nan
    #obtain data from last TMD for all proteins
    for acc in df.index:
        df.loc[acc,'last_TM_AAIMON_mean'] = df.loc[acc,'TM%02d_AAIMON_mean' % df.loc[acc,'number_of_TMDs']]
    AAIMON_last_TM = df['last_TM_AAIMON_mean'].dropna()
    #add the data for the last TMD to the dataframe
    df_mean_AAIMON_each_TM['last_TM_AAIMON_mean'] = AAIMON_last_TM

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
    ''' Checks the tarball that contains the SIMAP output.
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
        sys.stdout.write(uniprot_number_for_seq_in_row_1)
    '''
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
        # name changed to allow p-r( to be unique to the print function
        canonical_string_representation = repr
        return canonical_string_representation(self.parameter)

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

def calc_lipophilicity(seq, method = "mean"):
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
    method : string
        Method to be used to average the hydrophobicity values over the whole sequence.
        The hydrophobicity score is positive for polar/charged aa, negative for hydrophobic aa.
            "sum" will return the sum of the hydrophobicity scores over the sequence
            "mean" will return the mean of the hydrophobicity scores over the sequence

    Returns:
    --------
    mean hydrophobicity value for the sequence entered

    Usage:
    ------
    from korbinian.utils import calc_lipophilicity
    # for a single sequence
    s = "SAESVGEVYIKSTETGQYLAG"
    calc_lipophilicity(s)
    # for a series of sequences
    TMD_ser = df2.TM01_SW_match_seq.dropna()
    hydro = TMD_ser.apply(lambda x : calc_lipophilicity(x))

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
    # get the number of AA residues used to calculated the hydrophobicity
    # this is not simply the sequence length, as the sequence could include gaps or non-natural AA
    aa_counts_excluding_gaps = np.array(list(aa_counts_dict.values()))
    number_of_residues = aa_counts_excluding_gaps.sum()
    # if there are no residues, don't attempt to calculate a mean. Return np.nan.
    if number_of_residues == 0:
        return np.nan
    # convert dictionary to array, sorted by aa
    aa_counts_arr = np.array([value for (key, value) in sorted(aa_counts_dict.items())])
    multiplied = aa_counts_arr * hessa_scale
    sum_of_multiplied = multiplied.sum()
    if method == "mean":
        return sum_of_multiplied / number_of_residues
    if method == "sum":
        return sum_of_multiplied

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

def open_df_from_csv_zip(in_zipfile, filename=None, delete_corrupt=False):
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
    # in case someone accidentally applies the csv function to a pickle file, etc, print a warning
    if filename is not None:
        if filename[-4:] != ".csv":
            sys.stdout.write("Warning: extension is not .csv. ({})".format(filename)), sys.stdout.flush()
    # create bool deciding whether zip file will be deleted
    deletezip = False
    if os.path.isfile(in_zipfile):
        try:
            with zipfile.ZipFile(in_zipfile, "r", zipfile.ZIP_DEFLATED) as openzip:
                filenamelist = openzip.namelist()
                if filename is None:
                    # if a filename is not given, open the first file in the list
                    for file_in_zip in filenamelist:
                        if file_in_zip[-4:] == ".csv":
                            filename = file_in_zip
                            # stop searching after finding the first csv file
                            break
                    if filename is None:
                        # code should only continue here if no csv files found
                        # if the zipfile doesn't contain ANY csv files, something is seriously wrong. Either delete or raise Error.
                        if delete_corrupt == True:
                            deletezip = True
                            df_loaded = pd.DataFrame()
                        else:
                            raise FileNotFoundError("{} does not contain a valid csv file".format(in_zipfile))
                if filename is not None:
                    # if a filename is available, check if the file is in the zip
                    if filename in filenamelist:
                        csv_file_handle = openzip.open(filename)
                        # read as pandas dataframe
                        df_loaded = pd.read_csv(csv_file_handle, sep=",", quoting=csv.QUOTE_NONNUMERIC, index_col=0)
                    else:
                        # the desired file is not in the zip. Either delete the zip, or return an empty dataframe.
                        if delete_corrupt == True:
                            deletezip = True
                        else:
                            df_loaded = pd.DataFrame()
        except zipfile.BadZipFile:
            # the desired file is not in the zip. Either delete the zip, or return an empty dataframe.
            if delete_corrupt == True:
                deletezip = True
            else:
                df_loaded = pd.DataFrame()
    else:
        raise FileNotFoundError("{} not found".format(in_zipfile))
    if deletezip:
        logging.info("{} does not contain expected csv file {}. File is old or damaged, and has been deleted".format(in_zipfile, filename))
        os.remove(in_zipfile)
        df_loaded = pd.DataFrame()
    return df_loaded


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
        If the zip was corrupted or file could not be found, an empty dataframe will be returned.

    Note
    -------
    Much faster than reading from excel.
    To check that the file was successfully opened:
        df = open_df_from_pickle_zip(in_zipfile, filename)
        if df.empty:
            raise ValueError("corrupt zip or file not found")
    """
    # in case someone accidentally applies the csv function to a pickle file, etc, print a warning
    if filename is not None:
        if filename[-7:] != ".pickle":
            sys.stdout.write("Warning: extension is not .pickle. ({})".format(filename)), sys.stdout.flush()
    # create bool deciding whether zip file will be deleted
    deletezip = False
    if os.path.isfile(in_zipfile):
        try:
            with zipfile.ZipFile(in_zipfile, "r", zipfile.ZIP_DEFLATED) as openzip:
                filenamelist = openzip.namelist()
                if filename is None:
                    # if a filename is not given, open the first file in the list
                    for file_in_zip in filenamelist:
                        if file_in_zip[-7:] == ".pickle":
                            filename = file_in_zip
                            # pickle is found, stop searching
                            break
                    if filename is None:
                        # code should only continue here if no pickle files found
                        # if the zipfile doesn't contain ANY pickle files, something is seriously wrong. Either delete or raise Error.
                        if delete_corrupt == True:
                            deletezip = True
                            df_loaded = pd.DataFrame()
                        else:
                            raise FileNotFoundError("{} does not contain a valid pickle file".format(in_zipfile))

                if filename is not None:
                    # if a filename is available, check if the file is in the zip
                    if  filename in filenamelist:
                        with openzip.open(filename) as pickle_file_handle:
                            # read as pandas dataframe
                            # TEMP try/except function to understand the compatibility issues
                            try:
                                df_loaded = pd.read_pickle(pickle_file_handle)
                                # DEPRECATED METHOD
                                # df_loaded = pickle.load(pickle_file_handle)
                            except (EOFError, ModuleNotFoundError, io.UnsupportedOperation) as e:
                                # pickle files created by an older version of python have compatibility issues
                                if delete_corrupt:
                                    deletezip = True
                                    df_loaded = pd.DataFrame()
                                else:
                                    # assume they can be recreated. Return empty dataframe, usually a sign that the file should be deleted.
                                    df_loaded = pd.DataFrame()
                                    #raise IOError("{} in {} was created by an older version of pandas and cannot be opened.".format(filename, in_zipfile))
                            # make sure that the pickled object was REALLY a pandas object, and not some other python datatype that was pickled.
                            assert isinstance(df_loaded, (pd.Series, pd.DataFrame))
                    else:
                        # the desired file is not in the zip. Either delete the zip, or return an empty dataframe.
                        if delete_corrupt == True:
                            deletezip = True
                        else:
                            df_loaded = pd.DataFrame()
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
                    'TUM_oranges': {
                        'TUM0': (202, 101, 10),
                        'TUM1': (213, 148, 96),
                        'TUM2': (102, 49, 5),
                        'TUM3': (220, 108, 11),
                        'TUM4': (247, 194, 148),
                        'TUM5': (160, 78, 8)
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

    # custom colour lists
    output_dict['HTML_list01'] = ['#808080', '#D59460', '#005293', '#A1B11A', '#9ECEEC', '#0076B8', '#454545', "#7b3294", "#c2a5cf", "#008837", "#a6dba0"]
    output_dict['BGO_rgb'] = [(0,93,151), (83,167,214),(158,206,236),(94,152,49),(105,190,158),(166,213,196),(238,119,4),(214,168,82)] + colour_lists["tableau20"]
    output_dict['BGO_arr'] = [np.array(x)/255 for x in output_dict['BGO_rgb']]
    output_dict['BGO_dark_arr'] = [darken_or_lighten(c, -0.33) for c in output_dict['BGO_arr']]

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

def darken_or_lighten(col_array, fraction):
    col_array = col_array + fraction
    col_array[col_array > 1] = 1
    col_array[col_array < 0] = 0
    return col_array

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

def save_figure(fig, Fig_name, base_filepath, save_png, save_pdf, dpi = 400, close=True, transparent=True):
    """

    Parameters
    ----------
    fig
    Fig_name
    base_filepath
    save_png
    save_pdf
    dpi
    close
    transparent

    Returns
    -------

    """
    if not os.path.exists(base_filepath):
        os.makedirs(base_filepath)
    if save_png:
        fig.savefig(os.path.join(base_filepath, '{a}.png'.format(a=Fig_name)), format='png', dpi=dpi)
    if save_pdf:
        base_filepath_pdf = os.path.join(base_filepath, 'pdf')
        if not os.path.exists(base_filepath_pdf):
            os.makedirs(base_filepath_pdf)
        fig.savefig(os.path.join(base_filepath_pdf, '{a}.pdf'.format(a=Fig_name)), format='pdf', transparent=transparent)
    # close any open figures
    if close == True:
        plt.close('all')
    sys.stdout.write('Figure processed: {}\n'.format(Fig_name)), sys.stdout.flush()

class Log_Only_To_Console(object):
    def __init__(self):
        pass
    def info(self, message):
        sys.stdout.write("\n{}".format(message))
    def warning(self, message):
        sys.stdout.write("\n{}".format(message))
    def critical(self, message):
        sys.stdout.write("\n{}".format(message))

def convert_summary_csv_to_input_list(s, pathdict, logging, list_excluded_acc=None):
    # open dataframe with list of proteins
    df = pd.read_csv(pathdict["list_csv"], sep=",", quoting=csv.QUOTE_NONNUMERIC, index_col=0)
    # exclude any proteins where there is no list_of_TMDs
    df = df.loc[df['list_of_TMDs'].notnull()]
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

    # # extract values from settings file based on entered list number
    # rand_TM = s["rand_TM"]
    # rand_nonTM = s["rand_nonTM"]

    for p in list_p:
        p["s"] = s
        p["pathdict"] = pathdict
        p["logging"] = logging
        # p["rand_TM"] = rand_TM
        # p["rand_nonTM"] = rand_nonTM

    return list_p


def get_acc_list_from_txt(txtfile):
    """Opens a list of uniprot accessions and converts to a python list.

    Parameters
    ----------
    txtfile : str
        text file with list of acc, each with a new line

    Returns
    -------
    acc_list : list
        list of each uniprot accession.
        empty lines will be ignored
    """
    acc_list = []
    if os.path.isfile(txtfile):
        # Extracts accession numbers out of file
        with open(txtfile, "r") as source:
            for line in source:
                line = line.strip()
                if line is not "":
                    acc_list.append(line)
    # remove any redundant accessions
    acc_list = list(set(acc_list))
    return acc_list

def send_email_when_finished(s, pathdict):
    """ Sends an email to specified address when job is finished

    Parameters
    ----------
    s : dict
        Settings dictionary extracted from excel settings file.

    Returns
    -------
    nothing but sends an email

    """
    import smtplib
    from email.mime.multipart import MIMEMultipart
    from email.mime.text import MIMEText
    from email.mime.base import MIMEBase
    from email import encoders

    fromaddr = s["email_address"]
    toaddr = s['send_email_to']

    msg = MIMEMultipart()

    msg['From'] = fromaddr
    msg['To'] = toaddr
    msg['Subject'] = "korbinian run is finished"

    # if s['multiple_lists_to_analyse']:
    #     body = '{a}\n\nmultiple lists activated! lists to analyse: {c}\nprocessed list: {b}'.format(a=s['email_message'], b=list_number, c=s['multiple_lists_to_analyse'])
    #
    # else:
    #     body = '{a}\n\n processed list: {b}'.format(a=s['email_message'], b=list_number)

    body = '{a}\n\n processed list: {b}'.format(a=s['email_message'], b=s["list_number"])
    msg.attach(MIMEText(body, 'plain'))

    email_fig_list = []
    for settings_parameter in s.keys():
        if settings_parameter[-5:] == "email":
            if s[settings_parameter] == True:
                email_fig_list.append(settings_parameter)

    if email_fig_list != []:
        for email_fig in email_fig_list:
            Fig_name = email_fig[:-6]
            filepath = os.path.join(pathdict["single_list_fig_path"], Fig_name + ".png")
            sys.stdout.write("filepath : {}".format(filepath))
            if os.path.isfile(filepath):
                attachment = open(filepath, "rb")
                part = MIMEBase('application', 'octet-stream')
                part.set_payload((attachment).read())
                encoders.encode_base64(part)
                part.add_header('Content-Disposition', "attachment; filename= %s" % Fig_name + ".png")
                msg.attach(part)

    server = smtplib.SMTP('smtp.gmail.com', 587)
    server.starttls()
    server.login(s["email_address"], s["email_p_w0rd"])
    text = msg.as_string()
    server.sendmail(fromaddr, toaddr, text)
    server.quit()
    sys.stdout.write('\nEmail sent to {}\n'.format(toaddr))

def filter_for_truncated_sequences(min_perc_nonTMD, df_cr):
    if min_perc_nonTMD != 1:
        list_of_hits_to_keep = []
        list_of_hits_to_drop = []
        number_of_hits = len(df_cr.index)
        for hit in df_cr.index:
            if df_cr.loc[hit, 'perc_nonTMD_coverage'] >= min_perc_nonTMD:
                list_of_hits_to_keep.append(hit)
            else:
                list_of_hits_to_drop.append(hit)
        # keep only hits that were not excluded due to truncation
        df_cr = df_cr.loc[list_of_hits_to_keep, :]
        sys.stdout.write('Truncated alignments; homologues dropped: -- {}/{} --\n'.format(len(list_of_hits_to_drop), number_of_hits))
    else:
        sys.stdout.write('min_perc_nonTMD = 1 ; no filtering for truncated sequences \n')
    return df_cr

def calc_alpha_from_datapoints(data):
    if len(data) > 500:
        alpha = 500 / len(data) - 0.1
        if alpha < 0.1:
            alpha = 0.1
    else:
        alpha = 0.9
    return float(alpha)

def get_publication_colors():
    color_list = ['#808080', '#D59460', '#005293']
    # add other HTML colours
    color_list = color_list + ['#A1B11A', '#9ECEEC', '#0076B8', '#454545']
    return color_list

def varname(p):
    """Returns the variable name, e.g. "df" from the globals().

    This is somewhat of a hack, and may not be future proof. But debuggers are slow, and python objects don't hold their own names.
    see Stackoverflow page : http://stackoverflow.com/questions/592746/how-can-you-print-a-variable-name-in-python

    Parameters
    ----------
    p : python object

    Returns
    -------
    variable_name : str
        Variable name for the object entered.
    """
    for line in inspect.getframeinfo(inspect.currentframe().f_back)[3]:
        m = re.search(r'\bvarname\s*\(\s*([A-Za-z_][A-Za-z0-9_]*)\s*\)', line)
        if m:
            return m.group(1)

def pc(p):
    """Prints the variable name, followed by the value, separated by a comma.

    Use only when debugging.

    This is somewhat of a hack, and may not be future proof. But debuggers are slow, and python objects don't hold their own names.
    see Stackoverflow page : http://stackoverflow.com/questions/592746/how-can-you-print-a-variable-name-in-python

    Parameters
    ----------
    p : python object
    """
    variable_name = None
    for line in inspect.getframeinfo(inspect.currentframe().f_back)[3]:
        m = re.search(r'\bpc\s*\(\s*([A-Za-z_][A-Za-z0-9_]*)\s*\)', line)
        if m:
            variable_name = m.group(1)
            break
    sep = ", "
    frameinfo = getframeinfo(currentframe())
    sys.stdout.write("\nline {}, {}\n".format(frameinfo.lineno, os.path.basename(frameinfo.filename)))
    sys.stdout.write("{}{}{}\n".format(variable_name, sep, p))
    sys.stdout.flush()

def pn(p):
    """Prints the variable name, followed by the value, separated by a newline.

    Use only when debugging.

    This is somewhat of a hack, and may not be future proof. But debuggers are slow, and python objects don't hold their own names.
    see Stackoverflow page : http://stackoverflow.com/questions/592746/how-can-you-print-a-variable-name-in-python

    Parameters
    ----------
    p : python object
    """
    variable_name = None
    for line in inspect.getframeinfo(inspect.currentframe().f_back)[3]:
        m = re.search(r'\bpn\s*\(\s*([A-Za-z_][A-Za-z0-9_]*)\s*\)', line)
        if m:
            variable_name = m.group(1)
            break
    sep = "\n"
    frameinfo = getframeinfo(currentframe())
    caller = getframeinfo(stack()[1][0])
    line = caller.lineno
    scriptname = os.path.basename(caller.filename)
    sys.stdout.write("\nline {}, {}\n".format(line, scriptname))
    sys.stdout.write("{}{}{} ({})\n".format(variable_name, sep, p, type(p)))
    sys.stdout.flush()

def pr(p):
    """Shortened version of the print function

    Part of the debugging set of functions (pr, pc, pn)

    Note: Uses sys.stdout.write rather than print, so that debugging print functions
    in the code can be found and removed quickly through find-and-replace.


    """
    sys.stdout.write("{}\n".format(p))
    sys.stdout.flush()

# simple flatten function
flatten = lambda x: [item for sublist in x for item in sublist]

def read_signalp_output(filepath):
    """read signalp output file

    Note: SignalP uses a non-standard tab-formatted output file. Copy everything INCLUDING THE FIRST LINE that shows the SignalP version.

    # SignalP-4.1 gram- predictions
    # name                     Cmax  pos  Ymax  pos  Smax  pos  Smean   D     ?  Dmaxcut    Networks-used
    O83335                     0.280  25  0.284  25  0.467   2  0.351   0.309 Y  0.300      SignalP-TM
    A0A0J3ZJ69                 0.494  43  0.473  43  0.544  38  0.253   0.392 Y  0.300      SignalP-TM
    and so on

    This is better than the gff.txt, which only contains the proteins deemed to be signal peptides, and not the full list of proteins including those without signals.

    Parameters
    ----------
    filepath : str
        Full path to .txt or .csv file with signalp output


    Returns
    -------
    dfsp : pd.DataFrame
        Pandas dataframe with signalp output data.
    """
    cols = ['Cmax', 'pos', 'Ymax', 'pos', 'Smax', 'pos', 'Smean', 'D', 'Q', 'Dmaxcut', 'Networks-used']
    dfsp = pd.read_csv(filepath, skiprows=2, sep=r"\s*", header=None, index_col = 0, engine="python")
    dfsp.columns = cols
    return dfsp

def read_signalp_gff(filepath):
    """Returns the signalP gff file as a pandas dataframe

    Parameters
    ----------
    filepath : str
        Full path to file.

    Returns
    -------
    df_gff : pd.DataFrame
        Dataframe with gff data. Note that this only contains the proteins expected to be signal peptides.
    """
    df_gff = pd.read_csv(filepath, skiprows=3, sep="\t", index_col=0, header=None)
    cols = ["version", "signal", "start", "end", "score", ".", "..", "is_signal"]
    df_gff.columns = cols
    return df_gff


def HTMLColorToRGB(colorstring):
    """ convert #RRGGBB to an (R, G, B) tuple
    from http://code.activestate.com/recipes/266466-html-colors-tofrom-rgb-tuples/
    """
    colorstring = colorstring.strip()
    if colorstring[0] == '#': colorstring = colorstring[1:]
    if len(colorstring) != 6:
        raise ValueError ("input #%s is not in #RRGGBB format" % colorstring)
    r, g, b = colorstring[:2], colorstring[2:4], colorstring[4:]
    r, g, b = [int(n, 16) for n in (r, g, b)]
    return (r, g, b)

def file_is_old(filepath, oldest_acceptable_file_date_str):
    file_mtime = os.path.getmtime(filepath)
    oldest_acceptable_file_date_timeobject = time.strptime(oldest_acceptable_file_date_str, '%d_%m_%Y')
    oldest_acceptable_file_date_epochtime = time.mktime(oldest_acceptable_file_date_timeobject)
    if file_mtime < oldest_acceptable_file_date_epochtime:
        return True
    else:
        return False

def file_is_old_and_will_be_removed(filepath, oldest_acceptable_file_date_str, logging, acc):
    if file_is_old(filepath, oldest_acceptable_file_date_str):
        os.remove(filepath)
        message = "{} skipped, filepath is old and has been deleted".format(acc)
        logging.info(message)
    else:
        message = ""
    return message


concat_file = r"D:\databases\predictions\topology_fasta_smaller.txt"
predictions_dir = r"D:\databases\predictions"


def gen_fasta_records(concat_file):
    """Generator yields individual fasta records from file with concatenated records.

    The "record" is a string that contains any lines until the next ">" symbol.

    The first record will be an empty string.

    Parameters
    ----------
    concat_file : str
        Path to TMSEG file with concatenated predictions in fasta format

    Returns
    -------
    record : str
        Fasta string, e.g. ">P62258|1433E\nMDDREDLVYQAKLAEQAERYDEMVE"

    Usage
    -----
    concat_file = r"D:\data\concatfasta.fas"
    records = gen_fasta_records(concat_file)
    # skip first
    next(records)
    for n, record in enumerate(records):
        do_stuff_with_record(record)
    """
    with open(concat_file) as f:
        record = ""
        for line in f.readlines():
            if line[0] == ">":
                yield record
                record = ""
            record += line

def split_TMSEG_fasta_into_separate_files(concat_file, predictions_dir):
    """Splits concatenated TMSEG fasta-like predictions into separate files.

    Saves is subdirectory A1, A2, Q0 etc based on first two letters of uniprot accession.

    Parameters
    ----------
    concat_file : str
        Path to file with concat TMSEG fasta records.
    predictions_dir : str
        Parent path in which files are saved
    """
    records = gen_fasta_records(concat_file)
    next(records)
    for n, record in enumerate(records):
        acc = record[1:].split("|")[0]
        outfile = os.path.join(predictions_dir, acc[0:2], "{}_TMSEG_fastastyle.txt".format(acc))
        make_sure_path_exists(outfile, isfile=True)
        with open(outfile, "w") as f:
            f.write(record)
        sys.stdout.write("."), sys.stdout.flush()
    sys.stdout.write("\n{} records processed".format(n))

def downloaderInterface(url, outputPath, logging):
    """Downloads a file from the web (url) and saves it in the provided outputPath.

    Parameters
    ----------
    url : string
        String containing a valid url to a file which should be downloaded.
    outputPath : string
        String containing a valid path where the downloaded file should be saved.
    logging : boolean
        Logging variable (boolean) decides if a progressbar should be printed to console or not.

    Saved Files and Figures
    -----------------------
    Downloaded file saved at the provided output path
    """
    # response = requests.get(url, stream=True)
    # total_size = int(response.headers.get('content-length', 0));
    # block_size = 1024
    # with open(outputPath, 'wb') as f:
    #     if logging:
    #         for data in tqdm(response.iter_content(block_size), unit='KB',\
    #                     unit_scale=True, total=math.ceil(total_size//block_size)):
    #             f.write(data)
    #     else:
    #         for data in response.iter_content(block_size):
    #             f.write(data)

    #code obtained from "https://pypi.org/project/tqdm/"
    class TqdmUpTo(tqdm):
        """Provides `update_to(n)` which uses `tqdm.update(delta_n)`."""
        def update_to(self, b=1, bsize=1, tsize=None):
            """
            b  : int, optional
                Number of blocks transferred so far [default: 1].
            bsize  : int, optional
                Size of each block (in tqdm units) [default: 1].
            tsize  : int, optional
                Total size (in tqdm units). If [default: None] remains unchanged.
            """
            if tsize is not None:
                self.total = tsize
            self.update(b * bsize - self.n)  # will also set self.n = b * bsize

    if logging:
        with TqdmUpTo(unit='B', unit_scale=True, miniters=1,
                      desc=url.split('/')[-1]) as t:  # all optional kwargs
            urllib.request.urlretrieve(url, filename=outputPath,
                               reporthook=t.update_to, data=None)
    else:
        urllib.request.urlretrieve(url, filename=outputPath)
