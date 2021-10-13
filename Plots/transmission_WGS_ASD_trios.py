# -*- coding: utf-8 -*-
"""
Created on Mon Oct 11 19:03:36 2021

@author: Martin
"""

#Import modules

import sys, os, gzip
import pandas as pd
import numpy as np
import math
import scipy as sp
import seaborn as sns
import matplotlib.pyplot as plt
import pingouin as pg
from astropy import modeling

from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle, Circle, Polygon, Patch
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D as Line

plt.rcParams['svg.fonttype'] = 'none'
plt.ioff()

#Import modules done
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Functions

def remove_wrong_ids(df):
    '''remove the ID string that is wrong for Blood; and replace - with _'''
    
    df['ID'] = df.apply(lambda row: '_'.join(row.ID[4:].split('-')), axis=1)


def simplify_MOSAIC_CLASS(df):
    '''simplify MOSAIC_CLASS_SSB and drop the more complex current one.'''
    
    dictionary = {'Sperm' : '1_sperm', 'Shared' : '2_shared',
                  'Blood-A' : '3_soma'}
    
    df['MOSAIC_CLASS_SSS'] = df.apply(lambda row:
                                          dictionary[row.MOSAIC_CLASS_SSByBa],
                                      axis=1)
    df.drop('MOSAIC_CLASS_SSByBa', axis=1, inplace=True)

#------------------------------------------------------------------------------
#add_id_columns

def add_info_columns(df):
    '''add info columns like unique identifier CHR_POS_REF_ALT; change things
    to be easily sortable; add sort levels; sort'''
    
    df['CHR_POS_REF_ALT'] = df.apply(cpra_ID, axis =1)
    df['TYPE'] = df.apply(lambda row: row.ID.split('_')[-1], axis=1)
    df['SD_SPERM_AF'] = df.apply(lambda row: int(row.TYPE == 'Sperm') *
                                             (row.SET_SPERM_MOSAIC) *
                                             (row.MAF),
                                 axis=1)
    df['SPERM_MOSAICISM_AF'] = df.groupby('CHR_POS_REF_ALT')['SD_SPERM_AF']\
                                 .transform('max')
    df.sort_values(by=['CHR_POS_REF_ALT', 'TYPE'], inplace=True)
    df.reset_index(inplace=True, drop=True)

def cpra_ID(row):
    '''make the CHR_POS_REF_ALT ID.'''
    
    if row.CHROM in ['X', 'Y']:
        chrom = row.CHROM
    else:
        if int(row.CHROM) > 9:
            chrom = str(row.CHROM)
        else:
            chrom = '0' + str(row.CHROM)
    
    return '_'.join([chrom, str(row.POS), row.REF, row.ALT])

#------------------------------------------------------------------------------    '''make the CHR_POS_REF_ALT ID.'''

def child_genotype(df):
    '''assign genotype based on the MAF; use 0.3 as cutoff based on the
    overview.'''
    
    df['CHILD_STATUS'] = df.apply(child_status, axis=1)
    
def child_status(row):
    
    dictionary = {0 : 'neg', 1 : 'het'}
    
    if 'Child' in row.TYPE:
        return dictionary[int(row.MAF > 0.3)]
    else:
        return 'father'



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Plotting

#---
#Transmission plot

def transmission_plot(df, first=True):
    '''use the Blast column and annotate it with rgb colors to then use imshow
    to get a tile plot. do it for all three families independently.
    first argument splits the df into half'''
    
    f, axs = plt.subplots(3, 1, gridspec_kw = {'height_ratios':[4,1,24]},
                          sharex=True)
    
    df = df[(df.TYPE.str.contains('Child')) & (df.SET_SPERM_MOSAIC)]
    
    if first == True:
        lst = ['F01', 'F02', 'F03', 'F04']
    else:
        lst = ['F05', 'F06', 'F07', 'F08']
    
    df = df[df.FAMILY.str.contains('|'.join(lst))]
    
    df.sort_values(['FAMILY', 'SPERM_MOSAICISM_AF'],
                   ascending=[True, False], inplace=True)
    
    sp_AF = df[~df.CHR_POS_REF_ALT.duplicated()]
    sp_AF = sp_AF.reset_index(drop=True).reset_index()
    
    IDs = list(set(sp_AF.FAMILY.to_list()))
    IDs.sort()
    dictionary_col = {ID : str(int(ID[-1]) * 1/8) for ID in IDs}
    
    colors = [dictionary_col[ID] for ID in sp_AF.FAMILY]
    colors_ss = ['g' if sss == '1_sperm' else 'xkcd:brown'
                   for sss in sp_AF.MOSAIC_CLASS_SSS]
    
    axs[0].bar((sp_AF.iloc[:,0]), sp_AF.SPERM_MOSAICISM_AF**0.5,
               color=colors,edgecolor='k', width=0.75)
    
    axs[1].bar((sp_AF.iloc[:,0]), [1 for i in sp_AF.iloc[:,0]],
               color=colors_ss, edgecolor='k', width=0.75)
    
    df.sort_values(['TYPE', 'FAMILY', 'SPERM_MOSAICISM_AF'],
                    ascending=[True, True, False], inplace=True)
    
    dictionary_status = {'het' : (255,0,0), 'neg' : (255,255,220)}
    
    if first == True:
        sp_n = 62
    else:
        sp_n = 69
    
    image = np.full((3, sp_n, 3), 255)
    offset = 0
    
    for ID in df.FAMILY.unique():
        
        sub_df = df[df.FAMILY == ID]
        
        child_n = len(sub_df.TYPE.unique())
        samp_n = len(sub_df.CHR_POS_REF_ALT.unique())
        
        sub_img_lst = []
        
        for samp in sub_df.TYPE.unique():
            
            df_sub_samp = sub_df[sub_df.TYPE == samp]
            l = [dictionary_status[stat] for stat in df_sub_samp.CHILD_STATUS]
            
            sub_img_lst.append(l)
        
        sub_img = np.array(sub_img_lst)
        
        image[0:child_n, offset:(offset + samp_n)] = sub_img
        
        offset += samp_n
    
    axs[2].imshow(image)
    
    if first == True:
        X = [18.5, 29.5, 43.5]
    else:
        X = [10.5, 29.5, 43.5]
    
    for x in X:
        axs[2].axvline(x=x, ymin=0, ymax=1, linestyle='--', color='k')
    
    axs[0].set_xticks([])
    axs[0].set_yticks([0.0, 0.05**0.5, 0.25**0.5, 0.5**0.5])
    axs[0].set_yticklabels(['0.00', '0.05', '0.25', '0.50'])
    axs[0].set_ylabel('Sperm AF (sqrt-t)')
    
    sns.despine(bottom=True, trim=True, offset=5, ax=axs[0])
    sns.despine(bottom=True, left=True, offset=5, ax=axs[1])
    sns.despine(bottom=True, left=True, offset=5, ax=axs[2])
    
    axs[1].set_xticks([])
    axs[1].set_yticks([])
    
    axs[2].set_xticks([], minor=True)
    axs[2].set_xlabel('Mosaic Variants')
    axs[2].set_yticks([i for i in range(4)])
    axs[2].set_yticklabels(['#' + str(i) for i in range(1,4)])
    axs[2].set_ylabel('Children')
    
    axs[2].set_xticks(np.arange(sp_n) - 0.5, minor=True)
    axs[2].set_yticks(np.arange(4) - 0.5, minor=True)
    axs[2].grid(which="minor", color="w", linestyle='-', linewidth=3)
    axs[2].tick_params(which='minor', bottom=False, left=False)
    axs[0].tick_params(which='minor', bottom=False, left=False)
    axs[1].tick_params(which='minor', bottom=False, left=False)
    
    plt.show()

#---

#---
#Expected/observed

def plot_expected_observed(df, x=60, y=0.10):
    '''to check results; use overview as an easier way to represent data.'''
    
    X = len(df[df.CHILD_STATUS == 'het'])
    df = get_final_prob_list(df)
    
    f, axs = plt.subplots()
    xticks = df.NUMBER.tolist()[::10]
    #xticks.append(618)
    #xticks.insert(1, 19)
    
    mean_init = df[df.PROBABILITY == df.PROBABILITY.max()]['NUMBER']\
                .to_list()[0]
    amp_init = df.PROBABILITY.max()
    fitter = modeling.fitting.LevMarLSQFitter()
    model = modeling.models.Gaussian1D(mean=mean_init, amplitude=amp_init)
    fitted_model = fitter(model, df.NUMBER, df.PROBABILITY)
    mean = fitted_model.mean.value
    std = fitted_model.stddev.value
    
    sns.barplot(x='NUMBER', y='PROBABILITY', data=df, color='0.8', ax=axs)
    plt.plot(np.linspace(0, len(df), 10000),
             fitted_model(np.linspace(0, len(df), 10000)),
             color='g', linestyle='--')
    axs.axvline(x=X, color='r')
    axs.axvline(x=(mean - 1.96 * std), color='g', linestyle='--')
    axs.axvline(x=(mean + 1.96 * std), color='g', linestyle='--')
    axs.set_xticks(xticks)
    axs.set_xticklabels(xticks)
    plt.xlim(0, x)
    plt.ylim(-0.001, y)
    plt.xlabel('Number of Transmitted Variants')
    plt.ylabel('Probability of Transmission')
    
    sns.despine(offset=5, trim=True)
    plt.show()
    
    
def plot_expected_observed_all(df, s=0):
    
    lst = []
    Xs = []
    
    IDs = ['F0' + str(i) for i in range(1, 9)]
    IDs.append('ALL')
    
    for ID in IDs:
        
        if ID == 'ALL':
            X = len(df[df.CHILD_STATUS == 'het'])
            df_ = get_final_prob_list(df)
        else:
            X = len(df[(df.CHILD_STATUS == 'het') &
                       (df.FAMILY == ID)])
            df_ = get_final_prob_list(df[df.FAMILY == ID])
        
        mean_init = df_[df_.PROBABILITY == df_.PROBABILITY.max()]['NUMBER']\
                    .to_list()[0]
        amp_init = df_.PROBABILITY.max()
        fitter = modeling.fitting.LevMarLSQFitter()
        model = modeling.models.Gaussian1D(mean=mean_init, amplitude=amp_init)
        fitted_model = fitter(model, df_.NUMBER, df_.PROBABILITY)
        mean = fitted_model.mean.value
        std = fitted_model.stddev.value
        
        lst.append([mean, std])
        Xs.append(X)
    
    data = np.array(lst)
    
    plt.errorbar(y=[i - s for i in range(9)], x=data[:, 0].astype(float),
                 xerr=data[:, 1].astype(float) * 1.96, marker='o', ms=5,
                 mfc='k', mec ='k', ls='', capsize=5, ecolor='k')
    
    plt.errorbar(y=[i + s for i in range(9)], x=Xs,
                 marker='d', ms=10, mfc='r', mec ='k', ls='', capsize=5,
                 ecolor='k')
    
    plt.xlabel('Transmitted Variants')
    plt.xlim(-0.4, 25)
    
    plt.yticks([0, 1, 2, 3, 4, 5, 6, 7, 8], IDs)
    
    sns.despine(offset=5, trim=True)
    #plt.yticks(ticks=[0, 2, 1, 3], labels=['F01', 'F02', 'F03', 'All'])
    plt.ylim(8.5, -0.5)

    plt.yticks()
    
    plt.show()


def get_final_prob_list(df):
    
    df_list = combine_all_children(df)
    prob_df = calculate_final_prob_df(df_list)
    
    return prob_df

def calculate_final_prob_df(df_list):
    '''using the df_list from combine_all_children, calculate the number of
    transmitted events blast by child. after each combination of two children,
    summarize to obtain a single probability df.'''
    
    combined = np.array([(0,1)])
    
    for c_trans in df_list:
        
        c_array = zip_df_to_array(c_trans)
        
        tmp_combine = np.array([(0,0)])
        
        for line in combined:
            arr_1 = line[0] + c_array[:,0]
            arr_2 = line[1] * c_array[:,1]
            arr_combined = np.stack((arr_1, arr_2), axis=1)
            tmp_combine = np.concatenate((tmp_combine, arr_combined))
        
        tmp_df = pd.DataFrame(tmp_combine, columns=['NUMBER', 'PROBABILITY'])
        
        tmp_df = tmp_df.groupby(by='NUMBER').sum().reset_index()
        
        combined = zip_df_to_array(tmp_df)
        
    df = pd.DataFrame(combined, columns=['NUMBER', 'PROBABILITY'])
    df['NUMBER'] = df['NUMBER'].astype('int64')
    
    return df

def zip_df_to_array(df):
    '''utility to obtain a numpy array from a 2 column data frame.'''
    
    zipped = list(zip(df.iloc[:,0], df.iloc[:,1]))
    array = np.array(zipped)
    return array

def combine_all_children(df):
    '''Number based on the children.'''
    
    IDs = df.FAMILY.unique()
    
    child_lst = []
    
    for ID in IDs:
        
        c_nr = len(df[(df.FAMILY == ID) &
                      (df.TYPE.str.contains('Child'))].TYPE.unique())
            
        child_lst.extend([ID for i in range(c_nr)])
    
    
    df_list = []
    
    for ID in child_lst:
        
        prob_series = df[(df.SET_SPERM_MOSAIC == True) & (df.FAMILY == ID) &
                         ~(df.CHR_POS_REF_ALT.duplicated())]\
                         ['SPERM_MOSAICISM_AF']
        trans_df = number_of_transmissions_per_child(prob_series)
        df_list.append(trans_df)
    
    return df_list

def number_of_transmissions_per_child(prob_series):
    '''uses dynamic programming to obtain the matrix of probability for n
    transmissions using the prob_series passed as an argument. returns the
    last column and number of events as a pd.DataFrame. prob_series is the
    sperm mosaic variants for the individual'''
    
    #get a tuple list with number, prob, and anti-prob
    probs = prob_series.tolist()
    antiprobs = [(1-p) for p in probs]
    all_p = list(zip(range(len(probs)), probs, antiprobs))
    
    #initialize the matrix with all 0, add the 1 at the position[1,0]
    mat = np.zeros(((len(probs) +2),(len(probs) + 1)))
    mat[1,0] = 1
    
    #dynamically calculate each field by adding up the left * anti-prob and
    #the upper left * prob fields
    for i in range(len(probs) + 1):
        i += 1
        for j, p, ap in all_p:
            j += 1
            mat[i,j] = mat[i, (j-1)] * ap + mat[(i-1), (j-1)] * p
    
    dist = pd.Series(mat[1:, -1])
    df = dist.reset_index()
    df.columns = ['NUMBER', 'PROBABILITY']
    
    return df


























