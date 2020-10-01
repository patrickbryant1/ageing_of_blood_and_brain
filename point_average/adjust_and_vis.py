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
from scipy.stats import ttest_ind
import matplotlib.pylab as pl
import pdb



#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''Adjust the pvals and visualize the change in selected markers vs age''')
parser.add_argument('--gene_annotations', nargs=1, type= str, default=sys.stdin, help = 'Path to gene annotations.')
parser.add_argument('--running_averages', nargs=1, type= str, default=sys.stdin, help = 'Path to marker running averages.')
parser.add_argument('--max_fold_change_df', nargs=1, type= str, default=sys.stdin, help = 'Path to marker max fold changes and pvals.')
parser.add_argument('--outdir', nargs=1, type= str, default=sys.stdin, help = 'Path to outdir.')


#######FUNCTIONS#######


def plot_probes(X,markers,ages,age_indices,overlapping_probes):
    '''Plot the change in probe vs age.
    '''


    u_probes = overlapping_probes['Reporter Identifier'].unique()
    for u_probe in u_probes:
        sel = overlapping_probes[overlapping_probes['Reporter Identifier']==u_probes[0]]
        sel = sel.reset_index()
        i = np.where(markers==u_probe)[0]
        vals = X[i,:][0,:]
        #Plot ages vs vals
        fig,ax = plt.subplots(figsize=(6/2.54, 6/2.54))
        #Get ra
        x_av,y_av = running_average(ages,age_indices,vals)
        plt.plot(x_av,np.log10(y_av), color = 'k', linewidth=1)
        plt.scatter(ages, np.log10(vals), color = 'midnightblue', s=0.1)

        #Format plot
        plt.title(u_probe)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        plt.ylabel('log Beta value')
        plt.xlabel('Age')
        plt.tight_layout()
        plt.savefig(outdir+'fold_changes/markers/'+u_probe+'_vs_age.png', format='png', dpi=300)
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
overlapping_probes = pd.read_csv(args.overlapping_probes[0]) #Probes with up/down regulation in different age group comparisons
overlapping_genes = pd.read_csv(args.overlapping_genes[0]) #Genes with up/down regulation in different age group comparisons
diff_probes = np.load(args.diff_probes[0], allow_pickle=True)#Probes not overlapping btw correlation analysis and age group comprison
top10_corr = pd.read_csv(args.top10_corr[0]) #Top 10 marker correlations with age
overlap_10_year_gaps = pd.read_csv(args.overlap_10_year_gaps[0])

#pdb.set_trace()
#Get data
X, markers, ages, age_indices = format_probes(joined_betas, sample_sheet, gene_annotations, outdir)
#plot_probes(X,markers,ages,age_indices,overlapping_probes)
#plot_genes(X,markers,ages,age_indices,overlapping_genes)
#plot_diff_probes(X,markers,ages,age_indices,diff_probes)
#plot_top10_probes(X,markers,ages,age_indices,top10_corr['Reporter Identifier'], top10_corr['R'])

plot_overlap_10_year_gaps(X,markers,ages,age_indices,overlap_10_year_gaps['Reporter Identifier'])
