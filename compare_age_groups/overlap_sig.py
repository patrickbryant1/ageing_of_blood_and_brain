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
    return corr_df

def plot_pvals(adjusted_corr_df, ai, bi, outdir):
    #Plot pvalue distribution
    agebins = ['19-30','30-40','40-50','50-60','60-70','70-80','80+']
    fig,ax = plt.subplots(figsize=(4.5/2.54, 4.5/2.54))
    pvals = np.array(adjusted_corr_df['p'])
    num_sig = len(adjusted_corr_df[adjusted_corr_df['Rejection on 0.05']==True])
    sns.distplot(pvals)
    plt.title(agebins[ai]+' vs '+agebins[bi]+'\n'+str(num_sig)+' sig on FDR 0.05')
    plt.xlabel('p-value')
    plt.tight_layout()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.savefig(outdir+'qval_'+str(ai)+'_'+str(bi)+'.png', format='png', dpi=300)
    plt.close()

def volcano_plot(adjusted_corr_df, ai, bi, outdir):
    '''Do a volcano plot
    '''
    agebins = ['19-30','30-40','40-50','50-60','60-70','70-80','80+']
    fold_change = np.array(adjusted_corr_df['fold_change'])
    pvals = np.array(adjusted_corr_df['p'])
    sig_pvals_index = adjusted_corr_df[adjusted_corr_df['Rejection on 0.05']==True].index

    #Log transform
    log2fc = np.log2(fold_change)
    neglog10pval = -np.log10(pvals)
    neglog10sigpval = neglog10pval[sig_pvals_index]
    log2fcsig = log2fc[sig_pvals_index]
    #Plot
    fig,ax = plt.subplots(figsize=(4.5/2.54, 4.5/2.54))
    plt.scatter(log2fc,neglog10pval, s=0.2, color='lightsteelblue')
    t  =-np.log10(0.05/len(log2fc))
    #plt.axhline(t, min(log2fc),max(log2fc), linewidth=0.3, linestyle ="--", color = 'firebrick', label='bonferroni threshold')
    #Plot those with fold change less than 1.5 (or 1/1.5)
    fc_t = np.log2(1.5)
    #high_fc_i = np.where((log2fc<-fc_t)|(log2fc>fc_t))
    #plt.scatter(log2fc[high_fc_i[0]],neglog10pval[high_fc_i[0]], s=0.2, color='midnightblue', label='FC>1.5')
    #Plot the ones sig on FDR 0.05 with FC >1.5
    high_fc_i = np.where((log2fcsig<-fc_t)|(log2fcsig>fc_t))
    plt.scatter(log2fcsig[high_fc_i[0]],neglog10sigpval[high_fc_i[0]], s=0.2, color='midnightblue', label='FC>1.5')
    plt.title(agebins[ai]+' vs '+agebins[bi])
    plt.xlabel('log2 fold change')
    plt.ylabel('-log10 p-value')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig(outdir+'qval_volcano_'+str(ai)+'_'+str(bi)+'.png', format='png', dpi=300)
    plt.close()

def group_genes(unique_genes):
    '''Go through all unique genes and group them according to overlaps
    '''

    single_genes = {}
    for gene in unique_genes:
        gene = gene.split(';')
        if len(gene)>1: #If more than one gene associated to the marker
            fetched_singles = [] #Many genes are repeated
            for g in gene: #Loop through all individual genes
                if g in fetched_singles:
                    continue #Make sure each gene in the split is only recorded once
                if g in single_genes.keys(): #If already found
                    single_genes[g].append(gene)
                else:
                    single_genes[g] = [gene]
                fetched_singles.append(g)
        else:
            if gene[0] in single_genes.keys(): #If already found
                single_genes[gene[0]].append(gene[0])
            else:
                single_genes[gene[0]] = [gene[0]]

    unique_genes_grouped=single_genes
    return unique_genes_grouped

def find_overlap(joined_dfs, agebins):
    '''Find markers that overlap significantly and have FC>1.5 across comparisons
    '''

    unique_probes = joined_dfs['Reporter Identifier'].unique()
    # for u_probe in unique_probes:
    #     overlap = joined_dfs[joined_dfs['Reporter Identifier']==u_probe]
    #
    #     log2fc = np.log2(overlap['fold_change'])
    #     if min(log2fc)<0 and max(log2fc)>0:
    #         print('FOUND!!!!')
    #         pdb.set_trace()
        #fig,ax = plt.subplots(figsize=(2.5/2.54, 2.5/2.54))
        # plt.bar(np.arange(len(overlap)), log2fc)
        # plt.title(u_probe)
        # ax.spines['top'].set_visible(False)
        # ax.spines['right'].set_visible(False)
        # plt.tight_layout()
        # plt.savefig(outdir+'fold_changes/'+u_probe+'.png', format='png', dpi=300)
        # plt.close()

    gene_annotations = pd.read_csv(args.gene_annotations[0])
    #Join on probe id
    joined_dfs = pd.merge(joined_dfs, gene_annotations, left_on='Reporter Identifier', right_on='Name', how='left')
    #Group genes
    unique_genes = joined_dfs['UCSC_RefGene_Name'].dropna()
    unique_genes = unique_genes.unique()
    unique_genes_grouped = group_genes(unique_genes)

    #Loop through all gene groups
    for gene_group in unique_genes_grouped:
        gene_df = pd.DataFrame()
        for gene in unique_genes_grouped[gene_group]:
            if type(gene) == list and len(gene)>1:
                gene = ';'.join(gene)
            gene_df = pd.concat([gene_df, joined_dfs[joined_dfs['UCSC_RefGene_Name']==gene]])
        #Check fold change
        log2fc = np.log2(gene_df['fold_change'])
        if min(log2fc)<0 and max(log2fc)>0:

                 print('FOUND!!!!', gene_group)
                 #Get age bins
                 age_comparisons = []
                 for index, row in gene_df.iterrows():
                     age_comparisons.append(agebins[row['id1']]+' vs '+agebins[row['id2']])

                 fig,ax = plt.subplots(figsize=(6/2.54, 6/2.54))
                 plt.bar(np.arange(len(log2fc)), log2fc, color = 'midnightblue')
                 plt.title(gene_group)
                 plt.xticks(np.arange(len(log2fc)), gene_df['Reporter Identifier'].values+'\n'+age_comparisons, rotation='vertical')
                 ax.spines['top'].set_visible(False)
                 ax.spines['right'].set_visible(False)
                 plt.ylabel('log2 fold change')
                 plt.tight_layout()
                 plt.savefig(outdir+'fold_changes/'+gene_group+'.png', format='png', dpi=300)
                 plt.close()

                 pdb.set_trace()
###########MAIN###########
#Plt
plt.rcParams.update({'font.size': 7})
#Args
args = parser.parse_args()
agelabels = {'blood':"Characteristics [age y]"}
#agelabel=agelabels[args.agelabel[0]]
#gene_annotations = args.gene_annotations[0]
#joined_betas = args.joined_betas[0]
#sample_sheet = args.sample_sheet[0]
indir = args.indir[0]
outdir = args.outdir[0]


#Get age group correlations
agebins = ['19-30','30-40','40-50','50-60','60-70','70-80','80+']
try:
    joined_dfs = pd.read_csv(outdir+'sig_markers.csv')
except:
    joined_dfs = pd.DataFrame()
    age1 = []
    age2 = []
    num_sig = []

    for i in range(7):
        for j in range(i+1,7):
            corr_df = pd.read_csv(indir+str(i)+'_'+str(j)+'_corr_results.csv')
            #Adjust pvals
            adjusted_corr_df = adjust_pvals(corr_df)
            adjusted_corr_df['id1']=i
            adjusted_corr_df['id2']=j
            #Volcano plot
            volcano_plot(adjusted_corr_df, i, j, outdir)
            #qval plot
            plot_pvals(adjusted_corr_df, i, j, outdir)
            #Select on qval
            sel = adjusted_corr_df[adjusted_corr_df['Rejection on 0.05']==True]
            sel = sel.reset_index()
            log2fc = np.log2(sel['fold_change'])
            t = np.log2(1.5)
            #Select on FC above or below log2(1.5)/-log2(1.5)
            sel_i = np.where((log2fc>t)|(log2fc<-t))[0]
            sel = sel.loc[sel_i]

            #Save
            age1.append(agebins[i])
            age2.append(agebins[j])
            num_sig.append(len(sel))
            print(len(sel),'significant markers with fold change >1.5 out of',len(corr_df),'on FDR<0.05')
            joined_dfs = joined_dfs.append(sel)


    #Save
    results = pd.DataFrame()
    results['Age group 1']=age1
    results['Age group 2']=age2
    results['Significant markers with fold change >1.5']=num_sig
    results.to_csv(outdir+'num_sig_fc1_5.csv')
    #Save joined dfs
    joined_dfs.to_csv(outdir+'sig_markers.csv')


#Find overlapping markers
find_overlap(joined_dfs, agebins)
pdb.set_trace()
