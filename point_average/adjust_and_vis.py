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
def adjust_pvals(comparison_df):
    '''Adjust the pvalues
    '''
    res = multipletests(comparison_df['p'], 0.05, method='fdr_bh')
    rej, cor_pval = res[0], res[1]
    comparison_df['Rejection on 0.05']=rej
    comparison_df['qval']=cor_pval
    return comparison_df

def calc_derivatives(max_fold_change_df, running_averages):
    '''Calculate the derivatives for all significant probes
    with FC >2 (or less than 1/2)
    '''
    
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
gene_annotations = pd.read_csv(args.gene_annotations[0],low_memory=False)
running_averages = np.load(args.running_averages[0], allow_pickle=True)
max_fold_change_df = pd.read_csv(args.max_fold_change_df[0])
outdir = args.outdir[0]

#Adjust pvals
max_fold_change_df = adjust_pvals(max_fold_change_df)
