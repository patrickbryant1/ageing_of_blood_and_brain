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

from statsmodels.stats.multitest import multipletests
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

def vis_pvals(comparison_df):
    '''Visualize the pvals
    '''

    fig,ax = plt.subplots(figsize=(6/2.54, 6/2.54))
    sns.distplot(comparison_df['p'])
    #Format plot
    plt.title('p-value distribution')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.ylabel('Density')
    plt.xlabel('p-value')
    plt.tight_layout()
    plt.savefig(outdir+'pval_distribution.png', format='png', dpi=300)
    plt.close()


def calc_derivatives(sel, running_averages):
    '''Calculate the derivatives for all significant probes
    with FC >2 (or less than 1/2)
    '''
    gradients = np.zeros((len(sel),running_averages.shape[1])) #Save gradients
    max_grad_diff = np.zeros(len(sel))
    sel_indices = np.array(sel.index) #Indices
    #Plot the selected ra as well
    fig1,ax1 = plt.subplots(figsize=(6/2.54, 6/2.54))
    fig2,ax2 = plt.subplots(figsize=(6/2.54, 6/2.54))
    for i in range(len(sel)):
        si = sel_indices[i] #Get index
        gradients[i,:]=np.gradient(running_averages[si,:]) #Calc gradient
        max_grad_diff[i] = (max(gradients[i,:])-min(gradients[i,:]))

        ax1.plot(np.arange(19,102),running_averages[si,:]/max(running_averages[si,:]),color='b', linewidth=0.1,alpha=0.1)
        ax2.plot(np.arange(19,102),gradients[i,:],color='b', linewidth=0.1,alpha=0.1)

    #Format plot
    ax1.set_title('Running averages')
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.set_ylabel('Normalized beta value')
    ax1.set_xlabel('Age')
    fig1.tight_layout()
    fig1.savefig(outdir+'ra.png', format='png', dpi=300)

    ax2.set_title('Gradients of running averages')
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.set_ylabel('Gradient')
    ax2.set_xlabel('Age')
    fig2.tight_layout()
    fig2.savefig(outdir+'gradients.png', format='png', dpi=300)
    plt.close()

    fig,ax = plt.subplots(figsize=(6/2.54, 6/2.54))
    #Plot gradients with diff >0.1
    sel_grad = gradients[np.where(max_grad_diff>0.1)]
    for i in range(len(sel_grad)):
        plt.plot(np.arange(19,102),sel_grad[i,:],color='b', linewidth=0.1)
    plt.show()
    #Plot distribution of the max grad diff
    fig,ax = plt.subplots(figsize=(6/2.54, 6/2.54))
    sns.distplot(max_grad_diff)
    #Format plot
    plt.title('Max gradient difference')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.ylabel('Density')
    plt.xlabel('Max gradient difference')
    plt.tight_layout()
    plt.savefig(outdir+'grad_diff_distribution.png', format='png', dpi=300)
    plt.close()


    #Plot max grad diff vs fold change
    fig,ax = plt.subplots(figsize=(6/2.54, 6/2.54))
    matplotlib.rc('lines', linewidth=0.5, linestyle='--')
    plt.scatter(max_grad_diff, np.log10(sel['fold_change']),s=0.1,color='cornflowerblue')
    sns.kdeplot(max_grad_diff, np.log10(sel['fold_change']),shade_lowest =False,color="w", ax=ax)
    #Format plot
    plt.title('Max grad. diff. vs FC')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.ylabel('log10(FC)')
    plt.xlabel('Max gradient difference')
    plt.tight_layout()
    plt.savefig(outdir+'grad_diff_vs_FC.png', format='png', dpi=300)
    plt.close()
    pdb.set_trace()

    #Plot the top 10 gradient changes


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
#Visualize pvals
vis_pvals(max_fold_change_df)
#Adjust pvals
max_fold_change_df = adjust_pvals(max_fold_change_df)
#Select significant probes (FDR<0.05) with FC >2 (or less than 1/2)
sel = max_fold_change_df[max_fold_change_df['Rejection on 0.05']==True]
sel = sel[np.absolute(sel['fold_change'])>2]
#Print the number selected
print(len(sel),'selected markers out of', len(max_fold_change_df))
#Calculate derivatives
calc_derivatives(sel, running_averages)
