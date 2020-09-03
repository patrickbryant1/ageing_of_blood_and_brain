#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import sys
import os
import glob
import pandas as pd
import numpy as np

import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt

from sklearn.decomposition import PCA
from scipy.stats import pearsonr
from statsmodels.stats.multitest import fdrcorrection
from scipy.stats import ttest_ind
import pdb



#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''Analyze all statistically significant methylation markers btw age groups''')
parser.add_argument('--joined_betas', nargs=1, type= str, default=sys.stdin, help = 'Path to joined_betas.')
parser.add_argument('--sample_sheet', nargs=1, type= str, default=sys.stdin, help = 'Path to sample sheet.')
parser.add_argument('--gene_annotations', nargs=1, type= str, default=sys.stdin, help = 'Path to gene annotations.')
parser.add_argument('--agelabel', nargs=1, type= str, default=sys.stdin, help = 'Agelabel.')
parser.add_argument('--outdir', nargs=1, type= str, default=sys.stdin, help = 'Path to outdir.')


#######FUNCTIONS#######
def get_ages(sample_sheet, sample_names, agelabel):
    '''Get the age for each participant
    '''
    sample_ages = []


    for name in sample_names:
        try:
            sample_ages.append(float(sample_sheet[sample_sheet['Source Name']==name+' 1'][agelabel].values[0]))
        except:
            pdb.set_trace()
            if sample_sheet[sample_sheet['Source Name']==name+' 1'][agelabel].values[0] == '>90':
                sample_ages.append(90.0)
            else:
                pdb.set_trace()

    return np.array(sample_ages)

def group_by_age(ages):

    age_indices = []
    bins = [19,30,40,50,60,70,80,max(ages)+1]
    for i in range(len(bins)-1):
        sel = np.where((ages>=bins[i])&(ages<bins[i+1]))
        age_indices.append(sel[0])

    return age_indices

def compare_probes(joined_betas, sample_sheet, gene_annotations, outdir):
    '''Analyze correaltion of probes to age
    '''


    #Merge on probe id
    merged = pd.merge(joined_betas, gene_annotations, left_on='Reporter Identifier', right_on='Name')
    #Check how many are zeros (unquantified, beta can't be 0)
    zeros = (merged[merged.columns[2:-34]] == 0).astype(int).sum(axis=1)
    zero_indices = np.where(zeros<10)
    #Get ages
    ages = get_ages(sample_sheet, joined_betas.columns[2:], agelabel)
    #Save ages
    age_df = pd.DataFrame()
    age_df['Sample'] = joined_betas.columns[2:]
    age_df['Age'] = ages
    age_df.to_csv(outdir+'ages.csv')

    #Get age groups
    age_indices = group_by_age(ages)
    #Save age groups
    np.save(outdir+'age_groups.npy',np.array(age_indices))

    #Looking at single marker comparisons
    markers = np.array(merged['Reporter Identifier'])
    markers = markers[zero_indices]


    #Methylation values
    X = np.array(merged[merged.columns[2:-34]])
    #Check how many are zeros (unquantified, beta can't be 0)
    print(np.round(100*X[X==0].shape[0]/(X.shape[0]*X.shape[1]),2), '% zeros')
    #Take all samples with less than 10 zeros
    X = X[zero_indices,:][0]
    print('Removed ', len(merged)-len(X), 'markers that had over 10 missing values (Beta=0)')
    for ai in range(len(age_indices)-1):
        for bi in range(ai+1,len(age_indices)): #Compare all combinations
            print(ai,bi)
            #R = np.zeros(X.shape[0])
            #p = np.zeros(X.shape[0])

            i1 = age_indices[ai]
            i2 = age_indices[bi]
            X1 = X[:,i1]
            X2 = X[:,i2]
            stats, pvals = ttest_ind(X1,X2,axis=1)
            #Volcano plot
            fold_change = np.average(X2,axis=1)/np.average(X1,axis=1)
            volcano_plot(fold_change, pvals, ai, bi, outdir)
            #Plot pvals
            plot_pvals(pvals, ai, bi, outdir)
            #agesel = np.append(age_indices[ai],age_indices[ai+1])
            #Xsel = X[:,agesel]

            # for xi in range(X.shape[0]):
            #     if xi%1000==0:
            #         print(xi)
            #     R[xi], p[xi] = pearsonr(Xsel[xi,:], agesel)
            #Save
            df = pd.DataFrame()
            df['Reporter Identifier']=markers
            df['stat']=stats
            df['p']=pvals
            df['fold_change']=fold_change
            #Save df
            df.to_csv(outdir+str(ai)+'_'+str(bi)+'_corr_results.csv')


    return None

def plot_pvals(pvals, ai, bi, outdir):
    #Plot pvalue distribution
    agebins = ['19-30','30-40','40-50','50-60','60-70','70-80','80+']
    fig,ax = plt.subplots(figsize=(4.5/2.54, 4.5/2.54))
    sns.distplot(pvals)
    plt.title(agebins[ai]+' vs '+agebins[bi]+'\n'+str(len(pvals[pvals<0.05/len(pvals)]))+ ' out of '+str(len(pvals))+' sig on 0.05')
    plt.xlabel('p-value')
    plt.tight_layout()
    plt.savefig(outdir+'pval_'+str(ai)+'_'+str(bi)+'.png', format='png', dpi=300)
    plt.close()

def volcano_plot(fold_change, pvals, ai, bi, outdir):
    '''Do a volcano plot
    '''
    agebins = ['19-30','30-40','40-50','50-60','60-70','70-80','80+']
    fig,ax = plt.subplots(figsize=(4.5/2.54, 4.5/2.54))
    log2fc = np.log2(fold_change)
    neglog10pval = -np.log10(pvals)
    plt.scatter(log2fc,neglog10pval, s=0.2, color='lightsteelblue')
    bonferroni_t  =-np.log10(0.05/fold_change.shape[0])
    plt.axhline(bonferroni_t, min(log2fc),max(log2fc), linewidth=0.3, linestyle ="--", color = 'firebrick', label='bonferroni threshold')
    #Plot those with fold change less than 2 (log2=1) (or half = -1)
    fc_t = np.log2(1.5)
    high_fc_i = np.where((log2fc<-fc_t)|(log2fc>fc_t))
    plt.scatter(log2fc[high_fc_i[0]],neglog10pval[high_fc_i[0]], s=0.2, color='midnightblue', label='FC>1.5')
    plt.legend()
    plt.title(agebins[ai]+' vs '+agebins[bi])
    plt.xlabel('log2 fold change')
    plt.ylabel('-log10 p-value')
    plt.tight_layout()
    plt.savefig(outdir+'volcano_'+str(ai)+'_'+str(bi)+'.png', format='png', dpi=300)
    plt.close()
###########MAIN###########
#Plt
plt.rcParams.update({'font.size': 7})
#Args
args = parser.parse_args()
agelabels = {'blood':"Characteristics [age y]"}
agelabel=agelabels[args.agelabel[0]]
gene_annotations = pd.read_csv(args.gene_annotations[0],low_memory=False)
joined_betas = pd.read_csv(args.joined_betas[0], low_memory=False)
print('Read betas')
sample_sheet = pd.read_csv(args.sample_sheet[0], sep = '\t')
outdir = args.outdir[0]

#Compare probes btw age stratified samples
compare_probes(joined_betas, sample_sheet, gene_annotations, outdir)

#Plot white box
fig,ax = plt.subplots(figsize=(4.5/2.54, 4.5/2.54))
plt.axis('off')
plt.savefig(outdir+'whitebox.png',format='png',dpi=300)
