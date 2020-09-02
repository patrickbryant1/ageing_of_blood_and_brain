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
parser.add_argument('--indirl', nargs=1, type= str, default=sys.stdin, help = 'indir.')
parser.add_argument('--outdir', nargs=1, type= str, default=sys.stdin, help = 'Path to outdir.')


#######FUNCTIONS#######
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
joined_betas = args.joined_betas[0]
sample_sheet = pd.read_csv(args.sample_sheet[0], sep = '\t')
indir = args.indir[0]
outdir = args.outdir[0]

plt.rcParams.update({'font.size': 7})
#Adjust pvals
single_results = adjust_pvals(single_results, outdir)
#Look at the significant markers
plot_sig(joined_betas, sample_sheet, single_results, agelabel, outdir)
