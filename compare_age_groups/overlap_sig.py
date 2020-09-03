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
from statsmodels.stats.multitest import multipletests
from scipy.stats import ttest_ind
import pdb



#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''Analyze all statistically significant methylation markers btw age groups''')
parser.add_argument('--joined_betas', nargs=1, type= str, default=sys.stdin, help = 'Path to joined_betas.')
parser.add_argument('--sample_sheet', nargs=1, type= str, default=sys.stdin, help = 'Path to sample sheet.')
parser.add_argument('--gene_annotations', nargs=1, type= str, default=sys.stdin, help = 'Path to gene annotations.')
parser.add_argument('--agelabel', nargs=1, type= str, default=sys.stdin, help = 'Agelabel.')
parser.add_argument('--indir', nargs=1, type= str, default=sys.stdin, help = 'indir.')
parser.add_argument('--outdir', nargs=1, type= str, default=sys.stdin, help = 'Path to outdir.')


#######FUNCTIONS#######
def adjust_pvals(corr_df):
    '''Adjust the pvalues
    '''
    res = multipletests(corr_df['p'], 0.05)
    rej, cor_pval = res[0], res[1]
    corr_df['Rejection on 0.05']=rej
    corr_df['qval']=cor_pval
    num_sig = rej[rej==True].shape[0]
    print(num_sig, 'out of', len(corr_df))


    return corr_df, num_sig

#corr_df[corr_df['Rejection on 0.05']==True]

def plot_qvals(qvals, ai, bi, outdir):
    #Plot pvalue distribution
    agebins = ['19-30','30-40','40-50','50-60','60-70','70-80','80+']
    fig,ax = plt.subplots(figsize=(4.5/2.54, 4.5/2.54))
    sns.distplot(qvals)
    plt.title(agebins[ai]+' vs '+agebins[bi])
    plt.xlabel('q-value')
    plt.tight_layout()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.savefig(outdir+'qval_'+str(ai)+'_'+str(bi)+'.png', format='png', dpi=300)
    plt.close()

def volcano_plot(fold_change, qvals, ai, bi, outdir):
    '''Do a volcano plot
    '''
    agebins = ['19-30','30-40','40-50','50-60','60-70','70-80','80+']
    fig,ax = plt.subplots(figsize=(4.5/2.54, 4.5/2.54))
    log2fc = np.log2(fold_change)
    neglog10qval = -np.log10(qvals)
    plt.scatter(log2fc,neglog10qval, s=0.2, color='lightsteelblue')
    t  =-np.log10(0.05)
    plt.axhline(t, min(log2fc),max(log2fc), linewidth=0.3, linestyle ="--", color = 'firebrick', label='bonferroni threshold')
    #Plot those with fold change less than 2 (log2=1) (or half = -1)
    fc_t = np.log2(1.5)
    high_fc_i = np.where((log2fc<-fc_t)|(log2fc>fc_t))
    plt.scatter(log2fc[high_fc_i[0]],neglog10qval[high_fc_i[0]], s=0.2, color='midnightblue', label='FC>1.5')
    plt.title(agebins[ai]+' vs '+agebins[bi])
    plt.xlabel('log2 fold change')
    plt.ylabel('-log10 p-value')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig(outdir+'qval_volcano_'+str(ai)+'_'+str(bi)+'.png', format='png', dpi=300)
    plt.close()

###########MAIN###########
#Plt
plt.rcParams.update({'font.size': 7})
#Args
args = parser.parse_args()
agelabels = {'blood':"Characteristics [age y]"}
agelabel=agelabels[args.agelabel[0]]
gene_annotations = args.gene_annotations[0]
joined_betas = args.joined_betas[0]
sample_sheet = args.sample_sheet[0]
indir = args.indir[0]
outdir = args.outdir[0]


#Get age group correlations

joined_dfs = pd.DataFrame()
ids1 = []
ids2 = []
num_sig = []
for i in range(7):
    for j in range(i+1,7):
        corr_df = pd.read_csv(indir+str(i)+'_'+str(j)+'_corr_results.csv')
        #Adjust pvals
        adjusted_corr_df,sig = adjust_pvals(corr_df)
        adjusted_corr_df['id1']=i
        adjusted_corr_df['id2']=j
        #Volcano plot
        volcano_plot(np.array(adjusted_corr_df['fold_change']), np.array(adjusted_corr_df['qval']), i, j, outdir)
        #qval plot
        plot_qvals(np.array(adjusted_corr_df['p']), i, j, outdir)
        num_sig.append(sig)
        #adjusted_corr_df['id']=id
        joined_dfs = joined_dfs.append(adjusted_corr_df)

pdb.set_trace()
