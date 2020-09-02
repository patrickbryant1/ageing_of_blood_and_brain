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
def adjust_pvals(corr_df, id):
    '''Adjust the pvalues
    '''
    res = multipletests(corr_df['p'], 0.05)
    rej, cor_pval = res[0], res[1]
    corr_df['Rejection on 0.05']=rej
    corr_df['qval']=cor_pval
    print(rej[rej==True].shape[0], 'out of', len(corr_df))


    return corr_df[corr_df['Rejection on 0.05']==True]



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
ids = []
num_sig = []
for i in range(6):
    corr_df = pd.read_csv(indir+str(i)+'_corr_results.csv')

    id = str(i)
    #Adjust pvals
    adjusted_corr_df = adjust_pvals(corr_df, id)
    ids.append(int(id))
    num_sig.append(len(adjusted_corr_df))
    #adjusted_corr_df['id']=id
    joined_dfs = joined_dfs.append(adjusted_corr_df)


#Plot num significant against grouping
age_groups = np.load(indir+'age_groups.npy', allow_pickle=True)
agebins = ['19-30 vs 30-40','30-40 vs 40-50','40-50 vs 50-60','50-60 vs 60-70','60-70 vs 70-80','70-80 vs 80+']
fig,ax = plt.subplots(figsize=(12/2.54, 12/2.54))
plt.bar(ids,num_sig)
plt.xticks(ids,agebins, rotation='vertical')
plt.ylabel('Number significant markers')
plt.title('Significant markers (FDR 0.05)')
plt.tight_layout()
plt.savefig(outdir+'num_sig_agegroups.png', format = 'png', dpi=300)
pdb.set_trace()
