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
    markers = merged['Reporter Identifier']
    #Save
    ages_compared = []
    pvals = []
    genes = []

    #Methylation values
    X = np.array(merged[merged.columns[2:-34]])

    for ai in range(len(age_indices)-1):
        print(ai)
        R = np.zeros(X.shape[0])
        p = np.zeros(X.shape[0])
        i1 = age_indices[ai]
        i2 = age_indices[ai+1]
        X1 = X[:,i1]
        X2 = X[:,i2]
        stats, pvals = ttest_ind(X1,X2,axis=1)

        fig,ax = plt.subplots(figsize=(6/2.54, 6/2.54))
        sns.distplot(pvals)
        plt.title(str(len(pvals[pvals<0.05/len(pvals)]))+ ' out of '+str(len(pvals))+' sig on 0.05')
        plt.xlabel('p-value')
        plt.tight_layout()
        plt.savefig(outdir+'pval'+str(ai)+'.png', format='png', dpi=300)
        plt.close()
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
        #Save df
        df.to_csv(outdir+str(ai)+'_corr_results.csv')


    return None

def adjust_pvals(results, outdir):
    '''Adjust the pvalues
    '''
    rej, cor_pval = fdrcorrection(results['pval'], 0.001)
    results['Rejection on 0.001']=rej
    results['qval']=cor_pval
    print(rej[rej==True].shape[0], 'out of', len(results))
    fig,ax = plt.subplots(figsize=(6/2.54, 6/2.54))
    sns.distplot(results['qval'])
    plt.title(str(rej[rej==True].shape[0])+ ' out of '+str(len(results))+' sig on 0.001')
    plt.xlabel('qval')
    plt.tight_layout()
    plt.savefig(outdir+'qval.png', format='png', dpi=300)
    plt.close()
    return results



###########MAIN###########
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
