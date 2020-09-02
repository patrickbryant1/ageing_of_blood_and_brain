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


def corr_probe_to_age(joined_betas, sample_sheet, gene_annotations, outdir):
    '''Analyze correaltion of probes to age
    '''

    #Merge on probe id
    merged = pd.merge(joined_betas, gene_annotations, left_on='Reporter Identifier', right_on='Name')

    #Get ages
    y = get_ages(sample_sheet, joined_betas.columns[2:], agelabel)

    #Looking at single marker correlations
    markers = merged['Reporter Identifier']
    #Save
    pvals = []
    genes = []
    pdb.set_trace()
    for index, row in merged.iterrows():
        #Get CpG vals
        X = np.array(row[merged.columns[2:-34]])
        pdb.set_trace()
        R,p = pearsonr(X,y)
        #Save
        pvals.append(p)
        genes.append(row['UCSC_RefGene_Name'])
        print(index, row['UCSC_RefGene_Name'])

    #Create df with results
    results = pd.DataFrame()
    results['Reporter Identifier'] = markers
    results['pval'] = pvals
    results['gene'] = genes
    #Save df
    results.to_csv(outdir+'single_corr_results.csv')


    return results

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

try:
    single_results = pd.read_csv(outdir+'single_corr_results.csv')
except:
    single_results = corr_probe_to_age(joined_betas, sample_sheet, gene_annotations, outdir)

#Set fontsize
plt.rcParams.update({'font.size': 7})
#Adjust pvals
single_results = adjust_pvals(single_results, outdir)
#Look at the significant markers
plot_sig(joined_betas, sample_sheet, single_results, agelabel, outdir)
