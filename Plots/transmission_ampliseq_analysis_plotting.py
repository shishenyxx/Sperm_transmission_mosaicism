# -*- coding: utf-8 -*-
"""
Created on Wed Nov 18 14:11:30 2020

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

#Annotate and prepare

#------------------------------------------------------------------------------
#add_id_columns

def add_id_columns(df):
    '''add id columns like unique identifier CHR_POS_REF_ALT, the data types 
    from ID; change things to be easily sortable; add sort levels; sort'''
    
    df['CHR_POS_REF_ALT'] = df.apply(cpra_ID, axis =1)
    df['FAMILY_ID'] = df.apply(lambda row: return_info(row)[0], axis=1)
    df['SAMPLE'] = df.apply(lambda row: return_info(row)[1], axis=1)
    df['TYPE'] = df.apply(lambda row: return_info(row)[2], axis=1)
    df['SORT_LVL_0'] = df.apply(lambda row: return_info(row)[3], axis=1)
    
    df.sort_values(by=['CHR_POS_REF_ALT', 'FAMILY_ID','SORT_LVL_0', 'SAMPLE'],
                   inplace=True)
    df.reset_index(inplace=True, drop=True)

def cpra_ID(row):
    '''make the CHR_POS_REF_ALT ID.'''
    
    if row.CHROM == 'X':
        chrom = row.CHROM
    else:
        if int(row.CHROM) > 9:
            chrom = str(row.CHROM)
        else:
            chrom = '0' + str(row.CHROM)
    
    return '_'.join([chrom, str(row.POS), row.REF, row.ALT])

def return_info(row):
    '''splits the ID field to get the family ID, samp_ID, and samp_type.'''
    
    ids = row.ID.split('_')
    
    f_id = ids[0]
    
    if len(ids) == 2 and ids[1] == 'Blood': #Control 1 and 2
        f_id = 'Control'
        samp_ID = ids[0]
        samp_type = ids[1]
        level = 3
    
    elif len(ids) == 2 and ids[1] != 'Blood':
        samp_ID = add_0(ids[1])
        samp_type = 'Blastocyst'
        level = 2
    
        
    elif len(ids) ==3:
        samp_type = ids[-1]
        if ids[1] == 'ED':
            samp_ID = 'EggDonor_' + ids[-1]
            level = 1
        elif ids[1] == 'Fa':
            samp_ID = 'SpermDonor_' + ids[-1]
            level = 0
        
    else:
        samp_ID, samp_type, level = 'zonk', 'zonk', 'zonk'
    
    return f_id, samp_ID, samp_type, level

def add_0(samp_ID):
    '''adds 0 to single digits for sorting.'''
    
    ID = samp_ID
    
    if len(str(ID)) == 2:
        ID = ID[0] + '0' + ID[1]
    
    return ID

#------------------------------------------------------------------------------

def add_indel_info(df):
    '''add info whether a variant is an indel or a long indel (more than 1bp
    difference).'''
    
    df['INDEL'] = (df.REF.str.len() != df.ALT.str.len())
    df['LONG_INDEL'] = ((df.REF.str.len() > 2) | (df.ALT.str.len() > 2))


def add_depth(df):
    '''add depth'''
    
    df['DEPTH'] = df.ALT_COUNT + df.REF_COUNT
    df['DEPTH'] = df.DEPTH.fillna(0)


def add_genotype(df, het_hom):
    '''use the het_hom_annotation_only table to annotate the table with the
    genotypes. make sure that the FAMILY_ID is not int64, but str object.'''
    
    mrg = pd.merge(df, het_hom, how='left', on=['CHR_POS_REF_ALT', 'FAMILY_ID',
                                                'SAMPLE'])
    mrg['GENOTYPE'] = mrg.GENOTYPE.fillna('no_WGS_genotype')
    
    return mrg

#------------------------------------------------------------------------------
#annotate hets and flag noise

def annotate_noise_het(df):
    '''annotate noise with threshold of {s.b.} for lower to not be considered
    noise, in addition to being higher than the control upper. het is assigned
    for Blasts if not noise, for tissues if upper is higher than 0.40. These
    values are based on empirical analysis from het variants also uses a 100x
    depth filter for both noise and het.'''
    
    df['FLAG_DEPTH_20'] = (df.DEPTH < 20)
    df['FLAG_DEPTH_100'] = (df.DEPTH < 100)
    df['FLAG_NOISE'] = df.apply(flag_noise, axis=1)
    df['FLAG_HET'] = df.apply(flag_het, axis=1)

def flag_noise(row):
    '''for the threshold analysis the following command was used to determine
    the 95% CI across 444 ref hom events that were depth > 100
    
    data[(data.GENOTYPE == 'ref_hom') & (data.DEPTH > 100)].LOWER_CI.\
    describe(percentiles=[0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975])
    
    count    4.440000e+02
    mean     3.773686e-04
    std      6.221819e-04
    min     -2.710000e-20
    2.5%     0.000000e+00
    5%       0.000000e+00
    25%      5.955000e-05
    50%      2.171115e-04
    75%      4.386203e-04
    95%      1.368370e-03
    97.5%    1.638222e-03
    max      8.011659e-03
    Name: LOWER_CI, dtype: float64'''
    
    above_ctrl1 = row.LOWER_CI > row.NORMAL1_UPPER_CI
    above_ctrl2 = row.LOWER_CI > row.NORMAL2_UPPER_CI
    norm_thresh = ((row.NORMAL1_LOWER_CI < 0.001368370) &
                  (row.NORMAL2_LOWER_CI < 0.001368370))
    
    if row.TYPE == 'Blastocyst':
        depth = (row.DEPTH >= 20)
        above_thresh = row.LOWER_CI > 0.05
    else:
        depth = (row.DEPTH >= 100)
        above_thresh = row.LOWER_CI > 0.001368370
    
    alt = (row.ALT_COUNT >= 3)
    not_ctrl = (row.FAMILY_ID != 'Control')

    #note that this will not work for some of the SNPs where the ctrls have
    #them, therefore this is not used for SNPs.
    return not bool(above_thresh and above_ctrl1 and above_ctrl2 and
                    norm_thresh and depth and alt and not_ctrl)

def flag_het(row):
    '''define het as described above.'''

    if row.TYPE == 'Blastocyst':
        return not (row.FLAG_NOISE)
    
    else:
        return bool((row.UPPER_CI >= 0.4) & (row.DEPTH >= 20))
#---


def add_seed_info(df):
    '''add a recurrence flag, so analysis on the variant level that do not
    need organ information can be done.'''
    
    df['SEED'] = ~df.duplicated('CHR_POS_REF_ALT', keep='first')


def define_helper_columns(df):
    
    '''define additional columns that will help with determination of mosaicism
    status. for now, this includes whether a sample was used to call mosaicism,
    whether it matches the sample origin. Note that data.INDIVIDUAL detected
    has to be cast as string for MATCHED_ID to work.'''
    
    df['WGS_SAMPLE'] = df.SAMPLE.str.contains('SpermDonor')
    df['MATCHED_ID'] = (df.INDIVIDUAL_DETECTED == df.FAMILY_ID)
    df['MOSAIC_SD'] = (df.WGS_SAMPLE & df.MATCHED_ID & ~(df.FLAG_NOISE) &
                       ~(df.FLAG_DEPTH_100) & ~(df.FLAG_HET)
                       )
    df['MOSAIC_SD_SPERM'] = (df.MOSAIC_SD &
                             df.apply(lambda row: row.TYPE == 'Sperm', axis=1))
    df['SNP_IN_SAMPLE'] = (df.GENOTYPE != 'no_WGS_genotype')

#------------------------------------------------------------------------------
#annotate variant level

def define_mosaic_class(df):
    '''use the helper columns to get mosaic variants and to define the mosaic
    class.'''
    
    grp = df.groupby('CHR_POS_REF_ALT')
    
    df['SET_MOSAIC'] = grp['MOSAIC_SD'].transform(flag_signal)
    df['SET_SPERM_MOSAIC'] = grp['MOSAIC_SD_SPERM'].transform(flag_signal)
    df['SET_BOTH_MOSAIC'] = grp['MOSAIC_SD'].transform(flag_signal_2)
    df['AMPLI_MOSAIC_CLASS'] = df.apply(mosaic_class_ampli, axis=1)
    df['SET_SNP'] = grp['SNP_IN_SAMPLE'].transform(flag_signal)

def flag_signal(col):
    
    return sum(col) > 0

def flag_signal_2(col):
    
    return sum(col) > 1

def mosaic_class_ampli(row):
    '''actuall function extracting the class.'''
    
    sss = 'zonk'
    
    if row.SET_SPERM_MOSAIC == True:
        
        if row.SET_BOTH_MOSAIC == True:
            sss = '2_shared'
            
        else:
            sss = '1_sperm'
            
    elif row.SET_MOSAIC == True:
        sss = '3_soma'
    
    return sss

#---

#work on mosaic table
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def make_mosaic_snp_table(df):
    '''take only variants that are mosaic or were annotated SNPs. Also rename
    columns that have to renamed.'''
    
    mos = df[(df.SET_MOSAIC == True) | (df.SET_SNP == True)]
    
    mos = mos[['ID', 'CHR_POS_REF_ALT', 'CHROM', 'POS', 'REF', 'ALT', 'INDEL',
               'LONG_INDEL', 'ANNO', 'GENE', 'REF_SEQ', 'FAMILY_ID', 'SAMPLE',
               'TYPE', 'SORT_LVL_0', 'SEED', 'DEPTH', 'REF_COUNT', 'ALT_COUNT',
               'MAF', 'LOWER_CI', 'UPPER_CI', 'FLAG_DEPTH_20',
               'FLAG_DEPTH_100', 'FLAG_NOISE', 'FLAG_HET', 'WGS_SAMPLE',
               'GENOTYPE', 'SET_MOSAIC', 'SET_SPERM_MOSAIC', 'SET_BOTH_MOSAIC',
               'AMPLI_MOSAIC_CLASS', 'INDIVIDUAL_DETECTED', 'MATCHED_ID',
               'SET_SNP']]
    
    mos.rename({'AMPLI_MOSAIC_CLASS' : 'MOSAIC_CLASS_SSS',
                'FLAG_HET' : 'FLAG_HET_HOM', 'GENOTYPE' : 'WGS_GENOTYPE'},
               axis=1, inplace=True)
    
    return mos

#------------------------------------------------------------------------------
#SNP genotyper
def define_SNP_genotypes(df):
    '''define het/hom for parental genotypes. first step is to remove those
    SNPs that were absent across all data.'''
    
    df['TOTAL_DEPTH'] = df.groupby('CHR_POS_REF_ALT')['DEPTH'].transform('sum')
    
    df.drop(df.loc[df.TOTAL_DEPTH == 0].index, inplace=True)
    df.drop('TOTAL_DEPTH', axis=1, inplace=True)
    
    df['SNP_GENOTYPE'] = df.apply(lambda row: genotyper(row)[0], axis=1)
    df['SNP_GENOTYPE_N'] = df.apply(lambda row: genotyper(row)[1], axis=1)
    
    grp = df.groupby(['CHR_POS_REF_ALT', 'FAMILY_ID'])
    
    df['SNP_PARENTS'] = grp['SNP_GENOTYPE_N'].transform(parental_gts)
    df['SNP_PARENTS'] = df.apply(fix_8090, axis=1)
    
    
    
def genotyper(row):
    '''use logic to determine SNP_GENOTYPE.'''
    
    gt = 'zonk'
    n = -1
    
    if row.SET_SNP == True:
    
        if row.TYPE == 'Blastocyst':
            if row.LOWER_CI > 0.05:
                gt = 'Blastocyst_pos'
                n = 0
            elif row.FLAG_HET_HOM == False and row.FLAG_DEPTH_20 == False:
                gt = 'Blastocyst_neg'
                n = 0
            else:
                gt = 'LOW_DEPTH'
                n = 0
    
        elif 'Sperm' in row.SAMPLE or 'Egg' in row.SAMPLE:
            #this is based on visual inspection;
            #good separation at these levels
            if row.MAF <0.2:
                gt = 'ref'
                n = 0
            elif row.MAF >= 0.2 and row.MAF < 0.8:
                if row.TYPE != 'Sperm':
                    gt = 'het'
                    n = 1
                else:
                    gt = 'het'
                    n = 0
            elif row.MAF >= 0.8:
                if row.TYPE != 'Sperm':
                    gt = 'hom'
                    n = 10
                else:
                    gt = 'hom'
                    n = 0
        
        elif 'Control' in row.SAMPLE:
            if row.MAF <0.2:
                gt = 'ref_ctrl'
                n = 0
            elif row.MAF >= 0.2 and row.MAF < 0.8:
                gt = 'het_ctrl'
                n = 0
            elif row.MAF >= 0.8:
                gt = 'hom_ctrl'
                n = 0
    
    return (gt, n)

def parental_gts(col):
    '''use the groupby transform logic to determine SNP status in both parents.
    Uses negative values to identify non-SNPs and the combination of genotypes
    to determine the rest.'''
    
    if sum(col) == 0:
        return 'refxref'
    elif sum(col) == 1:
        return 'refxhet'
    elif sum(col) == 2:
        return 'hetxhet'
    elif sum(col) == 10:
        return 'refxhom'
    elif sum(col) == 11:
        return 'hetxhom'
    elif sum(col) == 20:
        return 'homxhom'
    elif sum(col) < 0:
        return 'not_SNP'
    else:
        return 'zonk'

def fix_8090(row):
    '''genotype of mother is unknown for 8090.'''
    
    dictionary = {'refxref' : 'refxunknown', 'refxhet' : 'hetxunknown',
                  'refxhom' : 'homxunknown'}
    
    if row.FAMILY_ID not in [8090, '8090']:
        return row.SNP_PARENTS
    else:
        return dictionary.get(row.SNP_PARENTS, row.SNP_PARENTS)
#---

#------------------------------------------------------------------------------
#transmission genotyper

#prep functions

def add_sperm_AF(df):
    '''make a column showing the related MAF in the appropriate father.'''
    
    #will make sure that only the father in which mosaicism is detected is
    #populated with the actual number. enables a max function.
    df['SD_SPERM_AF'] = df.apply(lambda row: int(row.TYPE == 'Sperm') *
                                             (row.MATCHED_ID) *
                                             (row.SET_SPERM_MOSAIC) *
                                             (row.SET_SNP == False) *
                                             (row.MAF),
                                 axis=1)
    
    df['SPERM_MOSAICISM_AF'] = df.groupby('CHR_POS_REF_ALT')['SD_SPERM_AF']\
                                 .transform('max')

#---
def sperm_mosaic_detected_in_blastocyst(df):
    '''use logic columns to determine if sperm mosaic variant is present in the
    blastocyst. also includes a control.'''
    
    df['BLAST_STATUS'] = df.apply(blast_status, axis=1)


def blast_status(row):
    
    stat = 'zonk'
    
    #is blast
    if row.TYPE == 'Blastocyst' and row.SET_SPERM_MOSAIC:
        
        #is in the potentially transmitting family
        if row.MATCHED_ID == False:
            
            if row.FLAG_DEPTH_20 == True:
                stat = 'cross_low_depth'
            elif row.FLAG_HET_HOM == True:
                stat = 'cross_pos'
            else:
                stat = 'cross_neg'
        
        else:
            
            if row.FLAG_DEPTH_20 == True:
                stat = 'low_depth'
            elif row.FLAG_HET_HOM == True:
                stat = 'pos'
            else:
                stat = 'neg'
    
    return stat
#---

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Plotting

#------------------------------------------------------------------------------
#parental mosaicism analysis; similar to ASD and CC papers

def pivot_spermdonor(df):
    '''restrict analysis table to SpermDonor and combine into one line.'''
    
    df = df[(df.SAMPLE.str.contains('Sperm')) & ~(df.SET_SNP) &
            (df.MATCHED_ID)]
    
    dictionary = {'Sperm' : 'SPERM', 'Blood' : 'SOMA', 'Saliva' : 'SOMA'}
    
    df['TYPE_SS'] = df.apply(lambda row: dictionary[row.TYPE], axis=1)
    
    piv = pd.pivot_table(df, values=['MAF', 'LOWER_CI', 'UPPER_CI',
                                     'REF_COUNT', 'ALT_COUNT'],
                             index=['CHR_POS_REF_ALT', 'CHROM', 'POS', 'REF',
                                    'ALT', 'FAMILY_ID', 'MOSAIC_CLASS_SSS',
                                    'INDEL'],
                             columns=['TYPE_SS'])
    
    piv.columns = ['_'.join(col) for col in piv.columns.values]
    
    piv.reset_index(inplace=True)
    
    return piv


def plot_with_custom_ci_sp_sh_bl(df, x_s=35, x_sh=25, x_b=140,
                                 y_s=0.21, y_sh=0.51, y_b=0.21):
    
    df = pivot_spermdonor(df)
    
    for organ in ['_SOMA', '_SPERM']:
    
        df['LOWER_ERROR' + organ] = df['MAF' + organ] - df['LOWER_CI' + organ]
        df['UPPER_ERROR' + organ] = df['UPPER_CI' + organ] - df['MAF' + organ]
    
    dict_div = {str(n) : 1 for n in range(1,24)}
    dict_div['X'] = 2
    dict_div['Y'] = 2
    
    for cat in ['MAF_SPERM', 'MAF_SOMA', 'LOWER_ERROR_SOMA',
                'UPPER_ERROR_SOMA', 'UPPER_CI_SPERM', 'LOWER_CI_SPERM',
                'UPPER_CI_SOMA', 'LOWER_CI_SOMA']:
        
        df[cat] = df.apply(lambda row: row[cat]/dict_div[str(row['CHROM'])],
                           axis=1)
    
    df_s = df[df.MOSAIC_CLASS_SSS == '1_sperm']\
             .sort_values(by=['MAF_SPERM'], ascending=False)
    df_s = df_s.reset_index(drop=True).reset_index()
    
    df_sh = df[df.MOSAIC_CLASS_SSS == '2_shared']\
             .sort_values(by=['MAF_SPERM'], ascending=False)
    df_sh = df_sh.reset_index(drop=True).reset_index()
    
    df_b = df[df.MOSAIC_CLASS_SSS == '3_soma']\
             .sort_values(by=['MAF_SOMA'], ascending=False)
    df_b = df_b.reset_index(drop=True).reset_index()
    
    f, axs = plt.subplots(nrows=3)
    
    
    axs[0].fill_between((df_s.iloc[:,0] + 1), df_s['UPPER_CI_SPERM'],
                                           df_s['LOWER_CI_SPERM'], color='0.8')
    axs[0].plot((df_s.iloc[:,0] + 1), df_s['MAF_SPERM'], color='g')
    
    
    axs[1].fill_between((df_sh.iloc[:,0] + 1), df_sh['UPPER_CI_SPERM'],
                                          df_sh['LOWER_CI_SPERM'], color='0.8')
    axs[1].plot((df_sh.iloc[:,0] + 1), df_sh['MAF_SPERM'],
                color='xkcd:khaki')
    axs[1].errorbar(x=(df_sh.iloc[:,0] + 1), y=df_sh['MAF_SOMA'],
               yerr=df_sh[['LOWER_ERROR_SOMA', 'UPPER_ERROR_SOMA']].T.values,
               ecolor='xkcd:burnt orange',
               capsize=0, elinewidth=2, marker='o', alpha=0.25,
               linestyle='None', mec='None', mfc='xkcd:black', markersize=4)
    
    
    axs[2].fill_between((df_b.iloc[:,0] + 1), df_b['UPPER_CI_SOMA'],
                                            df_b['LOWER_CI_SOMA'], color='0.8')
    axs[2].plot((df_b.iloc[:,0] + 1), df_b['MAF_SOMA'],
                color='xkcd:bright orange')
    
    axs[0].set_xlim(0,x_s)
    axs[0].set_ylim(0,y_s)
    sns.despine(offset=5, trim=True, ax=axs[0])
    axs[0].set_xlabel('Ranked Mosaic Variants')
    axs[0].set_ylabel('Norm. Sperm AF')
    
    axs[1].set_xlim(0,x_sh)
    axs[1].set_ylim(0,y_sh)
    axs[1].set_yticks([0, 0.25, 0.5])
    sns.despine(offset=5, trim=True, ax=axs[1])
    axs[1].set_xlabel('Ranked Mosaic Variants')
    axs[1].set_ylabel('Norm. Sperm/Soma AF')
    
    axs[2].set_xlim(0,x_b)
    axs[2].set_ylim(0,y_b)
    sns.despine(offset=5, trim=True, ax=axs[2])
    axs[2].set_xlabel('Ranked Mosaic Variants')
    axs[2].set_ylabel('Norm. Soma AF')
    
    plt.show()


def ratio_overlay_combine(df, x=60, n=20):
    
    f, axs = plt.subplots(2, 1, gridspec_kw = {'height_ratios':[3,1]},
                          sharex=True)
    
    df = pivot_spermdonor(df)
    
    for organ in ['_SOMA', '_SPERM']:
    
        df['LOWER_ERROR' + organ] = df['MAF' + organ] - df['LOWER_CI' + organ]
        df['UPPER_ERROR' + organ] = df['UPPER_CI' + organ] - df['MAF' + organ]
    
    dict_div = {str(n) : 1 for n in range(1,24)}
    dict_div['X'] = 2
    dict_div['Y'] = 2
    
    for cat in ['MAF_SPERM', 'MAF_SOMA', 'LOWER_ERROR_SOMA',
                'UPPER_ERROR_SOMA', 'UPPER_CI_SPERM', 'LOWER_CI_SPERM']:
        
        df[cat] = df.apply(lambda row: row[cat]/dict_div[str(row['CHROM'])],
                           axis=1)
    
    df = df[df.MOSAIC_CLASS_SSS != '3_soma']
    
    df['RATIO'] = df.apply(lambda row:
                                math.log((row['MAF_SPERM'] + 10**(-8))/
                                         (row['MAF_SOMA'] + 10**(-8))),
                           axis=1)
    df = df.sort_values(by=['MAF_SPERM'], ascending=False)
    df = df.reset_index(drop=True).reset_index()
    
    axs[0].errorbar(x=(df.iloc[:,0] + 1), y=df['MAF_SOMA'],
                 yerr=df[['LOWER_ERROR_SOMA', 'UPPER_ERROR_SOMA']].T.values,
                 ecolor='xkcd:bright orange',
                 capsize=0, elinewidth=2, marker='o', alpha=0.25,
                 linestyle='None', mec='None', mfc='xkcd:black', markersize=4)
    
    axs[0].fill_between((df.iloc[:,0] + 1), df['UPPER_CI_SPERM'],
                                           df['LOWER_CI_SPERM'], color='0.8')
    axs[0].plot((df.iloc[:,0] + 1), df['MAF_SPERM'], color='g')
    #axs[0].set_ylim(-0.01,0.6)
    axs[0].set_xlim(0,x)
    sns.despine(offset=5, trim=True, ax=axs[0])
    axs[0].set_ylabel('Norm. Sperm/Blood AF')
    
    axs[1].plot((df.iloc[:,0] + 1),df['RATIO'], color='k', linestyle='None',
             marker='o', markersize=3, mfc='None')
    roav = np.convolve(df['RATIO'], np.ones((n,))/n, mode='same')
    axs[1].plot((df.iloc[:,0] + 1), roav, color='0.5')
    axs[1].set_ylim(-5, 21)
    sns.despine(offset=5, trim=True, ax=axs[1])
    axs[1].set_ylabel('log(AF Ratio)')
    axs[1].set_xlabel('Ranked Mosaic Variants')
    
    plt.show()


def mosaic_number_per_individual(df, y=105):
    
    df = pivot_spermdonor(df)
    
    colors=['g', 'xkcd:brown', 'xkcd:bright orange']
    df = df.sort_values(by=['FAMILY_ID', 'MOSAIC_CLASS_SSS'])
    
    sns.countplot(x='FAMILY_ID', hue='MOSAIC_CLASS_SSS', data=df,
                  palette=colors)
    plt.ylim(0,y)
    sns.despine(bottom=True, offset=5, trim=True)
    plt.ylabel('Number of Mosaic Variants')

    plt.show()


def mosaic_number_indel_entirecohort(df, y=140):
    '''Breaking down SNV/INDELs'''
    
    df = pivot_spermdonor(df)
    
    colors=['g', 'xkcd:brown', 'xkcd:bright orange']
    df = df.sort_values(by=['MOSAIC_CLASS_SSS', 'INDEL'])
    
    sns.countplot(x='INDEL', hue='MOSAIC_CLASS_SSS', data=df,
                  palette=colors)
    
    plt.ylim(0,y)
    sns.despine(bottom=True, offset=5, trim=True)
    plt.ylabel('Number of Mosaic Variants')

    plt.show()


def mini_violins(df):
    '''make af violins for all individuals for sperm, shared-sperm,
    shared-soma, soma.'''
    
    df = pivot_spermdonor(df)
    
    #make squre root-transformed column for use in violin
    #not, that to reduce changes in function, old monikers are used
    df['MAF_SPERM_A_'] = df.MAF_SPERM**0.5
    df['MAF_BLOOD_'] = df.MAF_SOMA**0.5
    
    ids = df.FAMILY_ID.sort_values().unique()
    n = len(ids)
    
    f, axs = plt.subplots(nrows=4, ncols=len(ids))
    plt.subplots_adjust(wspace=0.)
    
    for i, sid in enumerate(ids):
        
        #sperm-specific
        sns.violinplot(x='FAMILY_ID', y='MAF_SPERM_A_',
                       data=df[(df.FAMILY_ID == sid) &
                               (df.MOSAIC_CLASS_SSS == '1_sperm')],
                       palette=['g'], ax=axs[0,i], cut=0., inner=None)
        
        #shared sp-af
        sns.violinplot(x='FAMILY_ID', y='MAF_SPERM_A_',
                       data=df[(df.FAMILY_ID == sid) &
                               (df.MOSAIC_CLASS_SSS == '2_shared')],
                       palette=['xkcd:khaki'], ax=axs[1,i], cut=0., inner=None)
    
        #shared bl-af
        sns.violinplot(x='FAMILY_ID', y='MAF_BLOOD_',
                       data=df[(df.FAMILY_ID == sid) &
                               (df.MOSAIC_CLASS_SSS == '2_shared')],
                       palette=['xkcd:burnt orange'], ax=axs[2,i], cut=0.,
                       inner=None)
        
        #blood-specific
        sns.violinplot(x='FAMILY_ID', y='MAF_BLOOD_',
                       data=df[(df.FAMILY_ID == sid) &
                               (df.MOSAIC_CLASS_SSS == '3_soma')],
                       palette=['xkcd:bright orange'], ax=axs[3,i], cut=0.,
                       inner=None)
    
    for i in range(4):
        for j in range(n):
        
            axs[i,j].set_ylim(0.,1.)
            axs[i,j].set_xlabel('')
            axs[i,j].set_ylabel('')
            sns.despine(bottom=True, trim=True, offset=5, ax=axs[i,j])
    
    for i in range(3):
        for j in range(n):
        
            axs[i,j].set_xticks([])
    
    for i in range(4):
        for j in range(1, n):
        
            sns.despine(left=True, trim=True, ax=axs[i,j])
            axs[i,j].set_yticks([])
    
    for i in range(4):
        
        axs[i,0].set_yticks([0., 0.2**0.5, 1.])
        axs[i,0].set_yticklabels(['0.0', '0.2', '1.0'])
    
    axs[2,0].set_ylabel('Sperm/Blood AF (square root-transformed)')
    
    plt.show()


def loliplots_sperm_mosaic(df):
    '''show just the sperm mosaic variants per individual.'''
    
    df = pivot_spermdonor(df)
    df = df[df.MOSAIC_CLASS_SSS != '3_soma'].groupby('FAMILY_ID').count()\
                                            .reset_index()
    
    plt.hlines(xmin=[0, 0, 0], xmax=df.POS, y=[3, 2, 1])
    plt.scatter(x=df.POS, y=[3, 2, 1], color='0.5', s=100, edgecolors='k',
                zorder=100)
    
    plt.xlim(0, 30)
    plt.ylim(0.75, 3.5)
    
    plt.xticks(ticks=[0, 15, 30])
    plt.yticks(ticks=[1, 2, 3], labels=['8090', '7670', '7465'])
    
    plt.xlabel('')
    plt.ylabel('FAMILY_ID')
    
    sns.despine(trim=True, offset=5)
    plt.show()


def stacked_bar_sperm_mosaic(df):
    '''show total number of mosaic variants per person.'''
    
    df = pivot_spermdonor(df)
    
    df.sort_values('FAMILY_ID', inplace=True)
    
    df_sp = df[df.MOSAIC_CLASS_SSS == '1_sperm'].groupby('FAMILY_ID').count()\
                                            .reset_index()
    df_sh = df[df.MOSAIC_CLASS_SSS == '2_shared'].groupby('FAMILY_ID').count()\
                                            .reset_index()
    
    plt.bar(x=[0, 1, 2], height=df_sh.POS, color='xkcd:brown', ec='k')
    plt.bar(x=[0, 1, 2], height=df_sp.POS, bottom=df_sh.POS, color='g',
            ec='k')
    
    plt.ylabel('Number of Sperm Mosaic Variants')
    plt.xticks(ticks=[0, 1, 2], labels=['7465', '7670', '8090'])
    sns.despine(trim=True, offset=5)
    
    plt.show()
    
    
    


#---
#Transmission plot

def transmission_plot(df):
    '''use the Blast column and annotate it with rgb colors to then use imshow
    to get a tile plot. do it for all three families independently.'''
    
    f, axs = plt.subplots(3, 1, gridspec_kw = {'height_ratios':[4,1,24]},
                          sharex=True)
    
    df = df[(df.TYPE == 'Blastocyst') & (df.SET_SPERM_MOSAIC) &
            (df.MATCHED_ID)]
    df.sort_values(['FAMILY_ID', 'SPERM_MOSAICISM_AF'],
                   ascending=[True, False], inplace=True)
    
    sp_AF = df[~df.CHR_POS_REF_ALT.duplicated()]
    sp_AF = sp_AF.reset_index(drop=True).reset_index()
    
    dictionary_col = {'7465' : '0.25', '7670' : '0.5', '8090' : '0.75'}
    colors = [dictionary_col[ID] for ID in sp_AF.FAMILY_ID]
    colors_ss = ['g' if sss == '1_sperm' else 'xkcd:brown'
                   for sss in sp_AF.MOSAIC_CLASS_SSS]
    
    axs[0].bar((sp_AF.iloc[:,0]), sp_AF.SPERM_MOSAICISM_AF**0.5,
               color=colors,edgecolor='k', width=0.75)
    
    axs[1].bar((sp_AF.iloc[:,0]), [1 for i in sp_AF.iloc[:,0]],
               color=colors_ss, edgecolor='k', width=0.75)
    
    #visually not nice; maybe elongate line from below in AI
    #for x in [11.5, 24.5]:
    #    axs[0].axvline(x=x, ymin=0, ymax=1, linestyle='--', color='k')
    
    df.sort_values(['SAMPLE', 'FAMILY_ID', 'SPERM_MOSAICISM_AF'],
                    ascending=[True, True, False], inplace=True)
    
    dictionary_status = {'pos' : (255,0,0), 'neg' : (255,255,220),
                         'low_depth' : (224,224,224)}
    
    image = np.full((14, 55, 3), 255)
    offset = 0
    
    for ID in df.FAMILY_ID.unique():
        
        sub_df = df[df.FAMILY_ID == ID]
        
        bl_n = len(sub_df.SAMPLE.unique())
        samp_n = len(sub_df.CHR_POS_REF_ALT.unique())
        
        sub_img_lst = []
        
        for samp in sub_df.SAMPLE.unique():
            
            df_sub_samp = sub_df[sub_df.SAMPLE == samp]
            l = [dictionary_status[stat] for stat in df_sub_samp.BLAST_STATUS]
            
            sub_img_lst.append(l)
        
        sub_img = np.array(sub_img_lst)
        
        image[0:bl_n, offset:(offset + samp_n)] = sub_img
        
        offset += samp_n
    
    axs[2].imshow(image)
    
    for x in [11.5, 24.5]:
        axs[2].axvline(x=x, ymin=0, ymax=1, linestyle='--', color='k')
    
    sns.despine(bottom=True, trim=True, offset=5, ax=axs[0])
    sns.despine(bottom=True, left=True, offset=5, ax=axs[1])
    sns.despine(bottom=True, left=True, offset=5, ax=axs[2])

    
    axs[0].set_xticks([])
    axs[0].set_yticks([0.0, 0.02**0.5, 0.1**0.5, 0.25**0.5])
    axs[0].set_yticklabels(['0.00', '0.02', '0.10', '0.25'])
    axs[0].set_ylabel('Sperm AF (sqrt-t)')
    
    axs[1].set_xticks([])
    axs[1].set_yticks([])
    
    axs[2].set_xticks([], minor=True)
    axs[2].set_xlabel('Mosaic Variants')
    axs[2].set_yticks([i for i in range(15)])
    axs[2].set_yticklabels(['#' + str(i) for i in range(1,15)])
    axs[2].set_ylabel('Blastocysts')
    
    axs[2].set_xticks(np.arange(56) - 0.5, minor=True)
    axs[2].set_yticks(np.arange(15) - 0.5, minor=True)
    axs[2].grid(which="minor", color="w", linestyle='-', linewidth=3)
    axs[2].tick_params(which='minor', bottom=False, left=False)
    axs[0].tick_params(which='minor', bottom=False, left=False)
    axs[1].tick_params(which='minor', bottom=False, left=False)
    
    plt.show()


#---
#SNP plots

def plot_refxhet_hetxhet(df):
    
    lst = []
    
    for gt in [('refxhet', 0.5), ('hetxhet', 0.75)]:
        
        gt_df = df[(df.SNP_PARENTS == gt[0]) & (df.TYPE == 'Blastocyst') &
                   (df.SNP_GENOTYPE != 'LOW_DEPTH')]
        
        counts = gt_df.SNP_GENOTYPE.value_counts().sort_values(ascending=False)
        
        x = counts[0]
        n = counts.sum()
        p = x/n
        ci = 1.96*(p*(1-p)/n)**0.5
        P = sp.stats.binom_test(x, n, p=gt[1])
    
        lst_temp = [p, ci, x, n, P, gt[0]]
        
        lst.append(lst_temp)
    
    data = np.array(lst)
    
    plt.errorbar(x=[0, 1], y=data[:, 0].astype(float),
                 yerr=data[:, 1].astype(float), marker='d', ms=10, mfc='k',
                 mec ='k', ls='', capsize=5, ecolor='k')
    
    plt.axhline(xmin=0, xmax=1, y=0.5, linestyle='--', color='r')
    plt.axhline(xmin=0, xmax=1, y=0.75, linestyle='--', color='b')
    
    plt.ylim(0, 1)
    plt.ylabel('Relative SNP Transmission')
    
    sns.despine(offset=5, trim=True)
    plt.xticks(ticks=[0, 1], labels=data[:, -1])
    plt.yticks(ticks=[0, 0.25, 0.5, 0.75, 1.0])
    
    plt.xticks()
    
    plt.show()
    
    return data


def plot_false_neg_SNP_corrected(df):
    
    fpr_fnr = []
    
    mos_neg_df = df[df.BLAST_STATUS.isin(['cross_pos', 'cross_neg'])]
    
    counts =mos_neg_df.BLAST_STATUS.value_counts()
    
    rate = min(counts)/counts.sum()
    
    fpr_fnr.append(rate)
    
    gt_df = df[(df.SNP_PARENTS == 'refxhom') & (df.TYPE == 'Blastocyst') &
               (df.SNP_GENOTYPE != 'LOW_DEPTH')]
    
    counts = gt_df.SNP_GENOTYPE.value_counts()
    
    rate = min(counts)/counts.sum()
    
    fpr_fnr.append(rate)
    
    
    plt.bar(x=[0, 1], height=fpr_fnr, color=['0.4', '0.6'], ec='k')
    
    plt.ylim(0, 1)
    plt.ylabel('Fraction of Blastocysts with Indicated Genotype')
    
    sns.despine(offset=5, trim=True)
    plt.xticks(ticks=[0, 1], labels=['FPR', 'FNR'])
    plt.yticks(ticks=[0, 0.2, 0.4, 0.6, 1.0])
    
    plt.xticks()
    
    plt.show()
    
    return fpr_fnr

#------------------------------------------------------------------------------
#Expected/Observed

def plot_expected_observed(df, x=60, y=0.10):
    '''7465, 7670: x=10, y=0.251
    8090: x=50, y=0.10'''
    
    X = len(df[df.BLAST_STATUS == 'pos'])
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
    
    for ID in [7465, 7670, 8090, 'all']:
        
        if ID == 'all':
            X = len(df[df.BLAST_STATUS == 'pos'])
            df_ = get_final_prob_list(df)
        else:
            X = len(df[(df.BLAST_STATUS == 'pos') &
                       (df.INDIVIDUAL_DETECTED == ID)])
            df_ = get_final_prob_list(df[df.INDIVIDUAL_DETECTED == ID])
        
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
    
    plt.errorbar(y=[i - s for i in range(4)], x=data[:, 0].astype(float),
                 xerr=data[:, 1].astype(float) * 1.96, marker='o', ms=5,
                 mfc='k', mec ='k', ls='', capsize=5, ecolor='k')
    
    plt.errorbar(y=[i + s for i in range(4)], x=Xs,
                 marker='d', ms=10, mfc='r', mec ='k', ls='', capsize=5,
                 ecolor='k')
    
    plt.xlabel('Transmitted Variants')
    plt.xlim(-2, 50)
    
    plt.yticks([0, 1, 2, 3], ['F01', 'F02', 'F03', 'All'])
    
    sns.despine(offset=5, trim=True)
    #plt.yticks(ticks=[0, 2, 1, 3], labels=['F01', 'F02', 'F03', 'All'])
    plt.ylim(3.5, -0.5)

    plt.yticks()
    
    plt.show()
    
    
    


def get_final_prob_list(df):
    
    df_list = combine_all_blasts(df)
    prob_df = calculate_final_prob_df(df_list)
    
    return prob_df

def calculate_final_prob_df(df_list):
    '''using the df_list from combine_all_blasts, calculate the number of
    transmitted events blast by blast. after each combination of two blasts,
    summarize to obtain a single probability df.'''
    
    combined = np.array([(0,1)])
    
    for b_trans in df_list:
        
        b_array = zip_df_to_array(b_trans)
        
        tmp_combine = np.array([(0,0)])
        
        for line in combined:
            arr_1 = line[0] + b_array[:,0]
            arr_2 = line[1] * b_array[:,1]
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

def combine_all_blasts(df):
    '''Number based on the blasts that were of sufficient quality.'''
    
    IDs = df[df.FAMILY_ID != 'Control'].FAMILY_ID.unique()
    
    blast_lst = []
    
    for ID in IDs:
        
        if ID == '7670':
            #manually remove bad Blastocysts in 7670
            b_nr = len(df[(df.FAMILY_ID == ID) &
                           ~(df.SAMPLE.isin(['B04', 'B07'])) &
                          (df.TYPE == 'Blastocyst')].SAMPLE.unique())
        else:
        
            b_nr = len(df[(df.FAMILY_ID == ID) &
                          (df.TYPE == 'Blastocyst')].SAMPLE.unique())
            
        blast_lst.extend([ID for i in range(b_nr)])
    
    
    df_list = []
    
    for ID in blast_lst:
        
        prob_series = df[(df.SET_SPERM_MOSAIC == True) & (df.SEED) &
                         (df.INDIVIDUAL_DETECTED.astype('str') == ID)]\
                         ['SPERM_MOSAICISM_AF']
        trans_df = number_of_transmissions_per_blast(prob_series)
        df_list.append(trans_df)
    
    return df_list

def number_of_transmissions_per_blast(prob_series):
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



















