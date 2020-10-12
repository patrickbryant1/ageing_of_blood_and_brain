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
from scipy.signal import savgol_filter

from collections import Counter
from sklearn.cluster import KMeans
from sklearn.manifold import TSNE
import pdb



#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''Adjust the pvals and visualize the change in selected markers vs age''')
parser.add_argument('--gene_annotations', nargs=1, type= str, default=sys.stdin, help = 'Path to gene annotations.')
parser.add_argument('--running_averages', nargs=1, type= str, default=sys.stdin, help = 'Path to marker running averages.')
parser.add_argument('--max_fold_change_df', nargs=1, type= str, default=sys.stdin, help = 'Path to marker max fold changes and pvals.')
parser.add_argument('--marker_values', nargs=1, type= str, default=sys.stdin, help = 'Path to marker values.')
parser.add_argument('--ages', nargs=1, type= str, default=sys.stdin, help = 'Path to sample ages.')
parser.add_argument('--point_indices', nargs=1, type= str, default=sys.stdin, help = 'Path to age points.')
parser.add_argument('--sample_sheet36194', nargs=1, type= str, default=sys.stdin, help = 'Path to sample sheet with accession 36194.')
parser.add_argument('--sample_sheet1575', nargs=1, type= str, default=sys.stdin, help = 'Path to sample sheet with accession 1575.')
parser.add_argument('--horvath_markers', nargs=1, type= str, default=sys.stdin, help = 'Path to horvath markers (71).')
parser.add_argument('--correlation_results', nargs=1, type= str, default=sys.stdin, help = 'Path to correlation results.')
parser.add_argument('--n_clusters', nargs=1, type= int, default=sys.stdin, help = 'n_clusters for kmeans clustering.')
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

def vis_age_distr(age_df, point_indices):
    '''Visualize the distribution of ages and the age points
    '''
    #Merge dfs
    ages = np.array(age_df['Age'])
    fig,ax = plt.subplots(figsize=(9/2.54, 6/2.54))
    sns.distplot(ages,color='grey', label='All')
    sns.distplot(age_df[age_df['Sex']=='female']['Age'],color='darkred', hist=False,label='Female')
    sns.distplot(age_df[age_df['Sex']=='male']['Age'], color='royalblue', hist=False,label='Male')
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
    y=0.02/len(point_indices)
    sns.distplot(ages,color='grey')
    age=0
    for i in range(len(point_indices)):
        agesel = ages[np.array(point_indices[i,:],dtype='int32')]
        plt.plot([min(agesel),max(agesel)],[y,y],alpha=0.5, color='royalblue',linewidth=1)
        plt.scatter(age,y,color='k',marker='|',s=1)
        y+=0.02/len(point_indices)
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


def calc_derivatives(sel, ages, running_averages, marker_values, point_indices,n_clusters):
    '''Calculate the derivatives for all significant probes
    with FC >2 (or less than 1/2)
    '''
    gradients = np.zeros((len(sel),running_averages.shape[1])) #Save gradients
    max_grad_diff = np.zeros(len(sel))
    sel_indices = np.array(sel['Unnamed: 0_x']) #Indices

    #Save the corr in different lists
    sel_ra = []
    sel_marker_values = []
    #Loop through the significant markers
    keep_indices = [] #keep the markers with sufficiently small std compared to the median
    norm=True
    for i in range(len(sel)):
        si = sel_indices[i] #Get index
        if norm == True:
            divider = max(marker_values[si,:])
        else:
            divider=1

        #Calculate the standard deviation for the running medians
        #Go through all point indices
        running_std = []
        for pi in range(len(point_indices)):
            age_points = point_indices[pi]
            running_std.append(np.std(marker_values[si,:][np.array(age_points,dtype='int32')]))

        #Check how big the std is compared to the mean on average
        running_std = np.array(running_std)
        rel_std_size = running_std/running_averages[si,:]

        if np.average(rel_std_size) >0.5: #If the average rel std is above 0.5
            continue

        keep_indices.append(i)
        #Calculate the maximal gradient difference
        gradients[i,:]=np.gradient(running_averages[si,:]/divider) #Calc gradient with norm
        max_grad_diff[i] = (max(gradients[i,:])-min(gradients[i,:]))

        #Save normalized selected ra
        sel_ra.append(running_averages[si,:]/divider)
        sel_marker_values.append(marker_values[si,:]/divider)

    #Select the gradients associated with the markers that had sufficiently low spread
    gradients = gradients[keep_indices]

    #Cluster the gradients
    k=n_clusters #Number of clusters
    kmeans = KMeans(n_clusters=k, random_state=0).fit(gradients)
    #In practice, the k-means algorithm is very fast (one of the fastest clustering algorithms available),
    #but it falls in local minima.
    #That’s why it can be useful to restart it several times.
    cluster_labels = kmeans.labels_
    #Visualize the gradient clustering
    colors = pl.cm.viridis(np.linspace(0,1,k))
    fig2,ax2 = plt.subplots(figsize=(6/2.54, 6/2.54))
    for cl in range(k):
        fig,ax = plt.subplots(figsize=(6/2.54, 6/2.54))
        cluster_indices = np.where(cluster_labels==cl)[0]
        for i in cluster_indices:
            ax.plot(np.arange(19,102),sel_ra[i],color=colors[cl],linewidth=1,alpha=0.5)

        ax.plot(np.arange(19,102), np.median(np.array(sel_ra)[cluster_indices],axis=0),color='k', linewidth=3)
        ax.scatter(ages,np.median(np.array(sel_marker_values)[cluster_indices],axis=0),color='k',s=1)
        #Plot gradients
        ax2.plot(np.arange(19,102), np.median(np.array(gradients)[cluster_indices],axis=0),color=colors[cl], linewidth=1)
        #Format plot
        ax.set_title('Cluster '+str(cl+1)+'|'+str(len(cluster_indices))+' markers')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_ylabel('Normalized beta value')
        ax.set_xlabel('Age')
        fig.tight_layout()
        fig.savefig(outdir+'/clustering/'+str(cl+1)+'.png', format='png', dpi=300)


    #Format plot
    ax2.set_title('Cluster gradients')
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.set_ylabel('Normalized beta value per year')
    ax2.set_xlabel('Age')
    fig2.tight_layout()
    fig2.savefig(outdir+'/clustering/gradients.png', format='png', dpi=300)
    plt.close()
    #Look at the clusters cin tsne
    #Tsne
    gradients_embedded = TSNE(n_components=2,random_state=0).fit_transform(gradients)
    fig,ax = plt.subplots(figsize=(6/2.54, 6/2.54))
    for cl in range(k):
        cluster_indices = np.where(cluster_labels==cl)[0]
        sel_emb_grads = gradients_embedded[cluster_indices]
        plt.scatter(sel_emb_grads[:,0],sel_emb_grads[:,1],color=colors[cl],label=cl+1,s=1)
    #Format plot
    plt.title('T-sne of gradient clusters')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.ylabel('Component 2')
    plt.xlabel('Component 1')
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir+'/clustering/tsne.png', format='png', dpi=300)
    plt.close()
    pdb.set_trace()
    #Select the keep indices
    sel = sel.loc[keep_indices]


    return sel

def group_markers_by_gene(sel, unique_genes_grouped):
    '''Group the selected markers per gene and return the
    genes that are regulated by at least 2 markers
    '''
    multi_marker_gene_df = pd.DataFrame()
    unique_gene_file = open(outdir+'genes/unique_genes.txt','w')
    for gene_group in unique_genes_grouped:
        unique_gene_file.write(gene_group+',\n')
        gene_df = pd.DataFrame()
        for gene in unique_genes_grouped[gene_group]:
            if type(gene) == list and len(gene)>1:
                gene = ';'.join(gene)
            gene_df = pd.concat([gene_df, sel[sel['UCSC_RefGene_Name']==gene]])
            gene_df['gene_group']=gene_group

        if len(gene_df)>1:
            multi_marker_gene_df = pd.concat([multi_marker_gene_df,gene_df])
            print(gene_group)
    unique_gene_file.close()
    return multi_marker_gene_df


def plot_multi_markers(multi_marker_gene_df,running_averages,marker_values,ages):
    '''Plot the running averages of the genes regulated by multiple markers
    '''

    u_genes = multi_marker_gene_df['gene_group'].unique()
    for gene in u_genes:
        multi_markers = multi_marker_gene_df[multi_marker_gene_df['gene_group']==gene]
        multi_markers = multi_markers.reset_index()
        #Plot
        fig,ax = plt.subplots(figsize=(6/2.54, 5.5/2.54))
        #Colors
        colors = pl.cm.viridis(np.linspace(0,1,len(multi_markers)))
        for i in range(len(multi_markers)):
            row = multi_markers.loc[i]
            index = row['Unnamed: 0_x']
            plt.plot(np.arange(0,103),running_averages[index,:],color = colors[i], linewidth=1,label=row['Reporter Identifier'])
            plt.scatter(ages, marker_values[index,:], color=colors[i],s=0.1)

        #Format plot
        plt.title(gene)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        plt.legend()
        plt.ylabel('Beta value')
        plt.xlabel('Age')
        plt.tight_layout()
        plt.savefig(outdir+'genes/'+gene+'.png', format='png', dpi=300)
        plt.close()

def analyze_horvath(horvath_markers,sel):
    '''Analyze which of the significant markers are used in the horvath clock
    '''
    overlap = pd.merge(horvath_markers,sel,left_on='Marker',right_on='Reporter Identifier', how='inner')
    print(len(overlap), 'of the selected markers overlap with the Horvath markers.')
    fig,ax = plt.subplots(figsize=(6/2.54, 6/2.54))
    sns.distplot(horvath_markers['Coefficient'],bins=30,label='Horvath', color='cornflowerblue')
    sns.distplot(overlap['Coefficient'],bins=30,label='Running median',color='darkgreen')
    plt.legend()
    plt.xlabel('Horvath coeffecient')
    plt.ylabel('Density')
    plt.title('Horvath markers and coefficients')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig(outdir+'horvath_markers.png', format='png', dpi=300)
    plt.close()

def correlation_overlap(correlation_results, sel):
    '''Check the overlap with the markers that have sig correlations
    with age on FDR=0.05
    '''

    #Adjust pvals for correlation results
    correlation_results = adjust_pvals(correlation_results)
    #Get only the sig
    sig_correlation_results =  correlation_results[correlation_results['Rejection on 0.05']==True]
    print(len(sig_correlation_results), 'are significant on FDR 0.05.')
    #Get  correlations for sel
    sel = pd.merge(sel,sig_correlation_results,on='Reporter Identifier', how='left')

    #See how many markers overlap
    print(len(sel[sel['Rejection on 0.05_y']==True]),'markers out of',len(sel),'were found in the correlation analysis')

    #Plot
    fig,ax = plt.subplots(figsize=(6/2.54, 6/2.54))

    sns.distplot(sig_correlation_results['R'],color='cornflowerblue', label='Significant correlations')
    sns.distplot(sel['R'], color='darkgreen', label='Running median')
    plt.xlabel('Pearson R')
    plt.ylabel('Density')
    plt.title('Marker correlations')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir+'correlations.png', format='png', dpi=300)
    plt.close()


def reg_feature_groups(sel):
    '''Check Regulatory_Feature_Group for the curves with pos and neg gradients respectively
    '''

    pos = sel[sel['pos_neg_grad']=='pos']
    neg = sel[sel['pos_neg_grad']=='neg']
    #Plot
    fig,ax = plt.subplots(figsize=(6/2.54, 6/2.54))
    pos_counts = Counter(pos['Regulatory_Feature_Group'])

    plt.xlabel('Pearson R')
    plt.ylabel('Density')
    plt.title('Marker correlations')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir+'correlations.png', format='png', dpi=300)
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
age_df = pd.read_csv(args.ages[0])
point_indices = np.load(args.point_indices[0],allow_pickle=True)
sample_sheet36194 = pd.read_csv(args.sample_sheet36194[0], sep = '\t')
sample_sheet1575 = pd.read_csv(args.sample_sheet1575[0], sep = '\t')
horvath_markers = pd.read_csv(args.horvath_markers[0])
correlation_results = pd.read_csv(args.correlation_results[0])
n_clusters = args.n_clusters[0]
outdir = args.outdir[0]

#Visualize pvals
vis_pvals(max_fold_change_df)
#Visualize age distribution and cutoffs
vis_age_distr(age_df, point_indices)

#Adjust pvals
max_fold_change_df = adjust_pvals(max_fold_change_df)
#Select significant probes (FDR<0.05) with FC >2 (or less than 1/2)
sel = max_fold_change_df[max_fold_change_df['Rejection on 0.05']==True]
sel = sel[np.absolute(sel['fold_change'])>2]

#Get the gene annotations for the selected markers
###NOTE!!! This resets the index!!!

sel = pd.merge(sel,gene_annotations,left_on='Reporter Identifier',right_on='Name', how='left')
#Calculate derivatives
sel = calc_derivatives(sel, age_df['Age'], running_averages, marker_values, point_indices,n_clusters)
#Print the number selected
print(len(sel),'selected markers out of', len(max_fold_change_df))
#Group genes
unique_genes_grouped = group_genes(sel['UCSC_RefGene_Name'].dropna().unique())
print(len(unique_genes_grouped.keys()),'unique genes')
#Get the genes regulated by multiple markers
multi_marker_gene_df = group_markers_by_gene(sel, unique_genes_grouped)
#Plot
plot_multi_markers(multi_marker_gene_df,running_averages,marker_values,age_df['Age'])

#Analyze horvath markers
#Get gene annotations
horvath_markers = pd.merge(horvath_markers, gene_annotations,left_on='Marker', right_on='Name', how='left')
#Save Horvath gene annotations
horvath_markers.to_csv(outdir+'genes/horvath/horvath_genes.csv')
#Group horvath markers
unique_genes_grouped = group_genes(horvath_markers['UCSC_RefGene_Name'].dropna().unique())
analyze_horvath(horvath_markers,sel)

#Analyze overlap with correlations
correlation_overlap(correlation_results, sel)

#Analyze 'Regulatory_Feature_Group' in relation to pos/neg gradients, sel['pos_neg_grad']
#reg_feature_groups(sel)
#Save sel
sel.to_csv(outdir+'ra_sig_markers.csv')
