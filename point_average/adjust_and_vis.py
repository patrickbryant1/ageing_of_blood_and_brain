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
parser.add_argument('--marker_values', nargs=1, type= str, default=sys.stdin, help = 'Path to marker values.')
parser.add_argument('--ages', nargs=1, type= str, default=sys.stdin, help = 'Path to sample ages.')
parser.add_argument('--age_points', nargs=1, type= str, default=sys.stdin, help = 'Path to age points.')
parser.add_argument('--sample_sheet', nargs=1, type= str, default=sys.stdin, help = 'Path to sample sheet.')
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

    fig,ax = plt.subplots(figsize=(9/2.54, 9/2.54))
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

def vis_age_distr(ages, age_points, sample_sheet):
    '''Visualize the distribution of ages and the age points
    '''
    #Merge dfs
    ages['Sample'] = ages['Sample']+ ' 1'
    merged = pd.merge(ages, sample_sheet, left_on='Sample', right_on = 'Source Name', how='left')
    ages = np.array(ages['Age'])
    fig,ax = plt.subplots(figsize=(9/2.54, 6/2.54))
    sns.distplot(ages,color='grey', label='All')
    sns.distplot(merged[merged['Characteristics [sex]']=='female']['Age'],color='darkred', hist=False,label='Female')
    sns.distplot(merged[merged['Characteristics [sex]']=='male']['Age'], color='royalblue', hist=False,label='Male')
    plt.title('Age distribution')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.ylabel('Density')
    plt.xlabel('Age')
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir+'ag_distribution.png', format='png', dpi=300)
    plt.close()

    #Plot cutoffs
    fig,ax = plt.subplots(figsize=(9/2.54, 6/2.54))
    y=0.03/len(age_points)
    sns.distplot(ages,color='grey')
    age=19
    for i in range(len(age_points)):
        agesel = ages[np.array(age_points[i,:],dtype='int32')]
        plt.plot([min(agesel),max(agesel)],[y,y],alpha=0.5, color='royalblue',linewidth=1)
        plt.scatter(age,y,color='k',marker='|',s=1)
        y+=0.03/len(age_points)
        age+=1


    plt.title('Running average groups')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.ylabel('Density')
    plt.xlabel('Age')
    plt.tight_layout()
    plt.savefig(outdir+'ag_point_cutoffs.png', format='png', dpi=300)
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


def calc_derivatives(sel, ages, running_averages, marker_values):
    '''Calculate the derivatives for all significant probes
    with FC >2 (or less than 1/2)
    '''
    gradients = np.zeros((len(sel),running_averages.shape[1])) #Save gradients
    max_grad_diff = np.zeros(len(sel))
    sel_indices = np.array(sel.index) #Indices
    #Save the positive and neg corr in different lists
    pos_sel_ra = []
    pos_sel_marker_values = []
    neg_sel_ra = []
    neg_sel_marker_values = []

    #Plot the selected gradients as well
    fig,ax = plt.subplots(figsize=(6/2.54, 6/2.54))
    for i in range(len(sel)):
        si = sel_indices[i] #Get index
        gradients[i,:]=np.gradient(running_averages[si,:]) #Calc gradient
        #Save normalized selected ra
        if np.sum(gradients[i,:]) >0:
            pos_sel_ra.append(running_averages[si,:]/max(marker_values[si,:]))
            pos_sel_marker_values.append(marker_values[si,:]/max(marker_values[si,:]))
        else:
            neg_sel_ra.append(running_averages[si,:]/max(marker_values[si,:]))
            neg_sel_marker_values.append(marker_values[si,:]/max(marker_values[si,:]))
        #Calculate the maximal gradient difference
        max_grad_diff[i] = (max(gradients[i,:])-min(gradients[i,:]))
        ax.plot(np.arange(19,102),gradients[i,:],color='royalblue', linewidth=0.1,alpha=0.2)

    #Format plot
    #Plot total gradients
    ax.plot(np.arange(19,102),np.average(gradients,axis=0),color='k', linewidth=1)
    ax.set_title('Gradients of running averages')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_ylabel('Gradient')
    ax.set_xlabel('Age')
    fig.tight_layout()
    fig.savefig(outdir+'gradients.png', format='png', dpi=300)
    plt.close()

    #Plot running averages with pos and neg gradients
    #Positive
    fig,ax = plt.subplots(figsize=(9/2.54, 9/2.54))
    for pi in range(len(pos_sel_ra)):
        ax.plot(np.arange(19,102),pos_sel_ra[pi],color='royalblue', linewidth=0.5,alpha=0.2)
    #Plot total ra
    pos_sel_ra = np.array(pos_sel_ra)
    print('Positively correlated markers:', len(pos_sel_ra))
    pos_sel_marker_values = np.array(pos_sel_marker_values)
    ax.plot(np.arange(19,102),np.average(pos_sel_ra,axis=0),color='k', linewidth=1)
    ax.scatter(ages,np.average(pos_sel_marker_values,axis=0),color='k',s=0.2)
    ax.set_title('Positive running averages')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_ylabel('Normalized beta value')
    ax.set_xlabel('Age')
    fig.tight_layout()
    fig.savefig(outdir+'pos_ra.png', format='png', dpi=300)

    #Negative
    fig,ax = plt.subplots(figsize=(9/2.54, 9/2.54))
    for pi in range(len(neg_sel_ra)):
        ax.plot(np.arange(19,102),neg_sel_ra[pi],color='royalblue', linewidth=0.5,alpha=0.5)
    #Plot total ra
    neg_sel_ra = np.array(neg_sel_ra)
    print('Negatively correlated markers:', len(neg_sel_ra))
    neg_sel_marker_values = np.array(neg_sel_marker_values)
    ax.plot(np.arange(19,102),np.average(neg_sel_ra,axis=0),color='k', linewidth=1)
    ax.scatter(ages,np.average(neg_sel_marker_values,axis=0),color='k',s=0.2)
    ax.set_title('Negative running averages')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_ylabel('Normalized beta value')
    ax.set_xlabel('Age')
    fig.tight_layout()
    fig.savefig(outdir+'neg_ra.png', format='png', dpi=300)
    # #Plot selected markers
    # fig,ax = plt.subplots(figsize=(6/2.54, 6/2.54))
    # #Plot ra with diff >0.1
    # sel_ra = running_averages[sel_indices[np.where((max_grad_diff>0.04)&(max_grad_diff<0.05))[0]]]
    # sel_markers = marker_values[sel_indices[np.where((max_grad_diff>0.04)&(max_grad_diff<0.05))[0]]]
    # for i in range(len(sel_ra)):
    #     plt.scatter(ages,sel_markers[i,:])
    #     plt.plot(np.arange(19,102),sel_ra[i,:],color='b', linewidth=1)
    #     plt.show()

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


def group_markers_by_gene(sel, unique_genes_grouped):
    '''Group the selected markers per gene and return the
    genes that are regulated by at least 2 markers
    '''
    multi_marker_gene_df = pd.DataFrame()
    for gene_group in unique_genes_grouped:
        gene_df = pd.DataFrame()
        for gene in unique_genes_grouped[gene_group]:
            if type(gene) == list and len(gene)>1:
                gene = ';'.join(gene)
            gene_df = pd.concat([gene_df, sel[sel['UCSC_RefGene_Name']==gene]])
            gene_df['gene_group']=gene_group

        if len(gene_df)>1:
            multi_marker_gene_df = pd.concat([multi_marker_gene_df,gene_df])
            print(gene_group)

    return multi_marker_gene_df


def plot_multi_markers(multi_marker_gene_df,running_averages,marker_values,ages):
    '''Plot the running averages of the genes regulated by multiple markers
    '''

    u_genes = multi_marker_gene_df['gene_group'].unique()
    for gene in u_genes:
        multi_markers = multi_marker_gene_df[multi_marker_gene_df['gene_group']==gene]
        multi_markers = multi_markers.reset_index()
        #Plot
        fig,ax = plt.subplots(figsize=(6/2.54, 6/2.54))
        #Colors
        colors = pl.cm.viridis(np.linspace(0,1,len(multi_markers)))
        for i in range(len(multi_markers)):
            row = multi_markers.loc[i]
            index = row['Unnamed: 0_x']
            plt.plot(np.arange(19,102),running_averages[index,:],color = colors[i], linewidth=1,label=row['Reporter Identifier'])
            plt.scatter(ages, marker_values[index,:], color=colors[i],s=0.1)

        #Format plot
        plt.title(gene)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        plt.ylabel('Beta value')
        plt.xlabel('Age')
        plt.tight_layout()
        plt.savefig(outdir+'genes/'+gene+'.png', format='png', dpi=300)
        plt.close()


###########MAIN###########
#Plt
plt.rcParams.update({'font.size': 7})
#Args
args = parser.parse_args()
gene_annotations = pd.read_csv(args.gene_annotations[0],low_memory=False)
running_averages = np.load(args.running_averages[0], allow_pickle=True)
max_fold_change_df = pd.read_csv(args.max_fold_change_df[0])
marker_values = np.load(args.marker_values[0], allow_pickle=True)
ages = pd.read_csv(args.ages[0])
age_points = np.load(args.age_points[0],allow_pickle=True)
sample_sheet = pd.read_csv(args.sample_sheet[0],sep='\t')
outdir = args.outdir[0]

#Visualize pvals
#vis_pvals(max_fold_change_df)
#Visualize age distribution and cutoffs
vis_age_distr(ages, age_points, sample_sheet)
pdb.set_trace()
#Adjust pvals
max_fold_change_df = adjust_pvals(max_fold_change_df)
#Select significant probes (FDR<0.05) with FC >2 (or less than 1/2)
sel = max_fold_change_df[max_fold_change_df['Rejection on 0.05']==True]
sel = sel[np.absolute(sel['fold_change'])>2]
#Print the number selected
print(len(sel),'selected markers out of', len(max_fold_change_df))
#Get the gene annotations for the selected markers
sel = pd.merge(sel,gene_annotations,left_on='Reporter Identifier',right_on='Unnamed: 0', how='left')
unique_genes_grouped = group_genes(sel['UCSC_RefGene_Name'].unique()[1:]) #The first is nan

#Calculate derivatives
#calc_derivatives(sel, ages['Age'], running_averages, marker_values)
multi_marker_gene_df = group_markers_by_gene(sel, unique_genes_grouped)
#Plot
plot_multi_markers(multi_marker_gene_df,running_averages,marker_values,ages['Age'])
pdb.set_trace()
