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
import pdb



#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''Visualize the change in selected markers vs age''')
parser.add_argument('--joined_betas', nargs=1, type= str, default=sys.stdin, help = 'Path to joined_betas.')
parser.add_argument('--sample_sheet', nargs=1, type= str, default=sys.stdin, help = 'Path to sample sheet.')
parser.add_argument('--gene_annotations', nargs=1, type= str, default=sys.stdin, help = 'Path to gene annotations.')
parser.add_argument('--agelabel', nargs=1, type= str, default=sys.stdin, help = 'Agelabel.')
parser.add_argument('--overlapping_probes', nargs=1, type= str, default=sys.stdin, help = 'Path to overlapping probes.')
parser.add_argument('--overlapping_genes', nargs=1, type= str, default=sys.stdin, help = 'Path to overlapping probes.')
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


def format_probes(joined_betas, sample_sheet, gene_annotations, outdir):
    '''Format the probes, their values and the sample ages
    '''


    #Merge on probe id
    merged = pd.merge(joined_betas, gene_annotations, left_on='Reporter Identifier', right_on='Name')
    #Check how many are zeros (unquantified, beta can't be 0)
    zeros = (merged[merged.columns[2:-34]] == 0).astype(int).sum(axis=1)
    zero_indices = np.where(zeros<10)
    #Get ages
    ages = get_ages(sample_sheet, joined_betas.columns[2:], agelabel)
    #Get age groups
    age_indices = group_by_age(ages)

    #Save ages
    age_df = pd.DataFrame()
    age_df['Sample'] = joined_betas.columns[2:]
    age_df['Age'] = ages
    age_df.to_csv(outdir+'ages.csv')

    #Looking at single marker comparisons
    markers = np.array(merged['Reporter Identifier'])
    markers = markers[zero_indices]


    #Methylation values
    X = np.array(merged[merged.columns[2:-34]])
    #Check how many are zeros (unquantified, beta can't be 0)
    print(np.round(100*X[X==0].shape[0]/(X.shape[0]*X.shape[1]),2), '% zeros')
    #Take all samples with less than 10 zeros
    X = X[zero_indices,:][0]


    return X, markers, ages, age_indices


def plot_probes(X,markers,ages,age_indices,overlapping_probes):
    '''Plot the change in probe vs age.
    '''

    def running_average(ages,age_indices,vals):
        x = []
        y = []
        for i in range(len(age_indices)):
            x.append(np.average(ages[age_indices[i]]))
            y.append(np.average(vals[age_indices[i]]))

        return x,y
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
    pdb.set_trace()
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
overlapping_probes = pd.read_csv(args.overlapping_probes[0])
overlapping_genes = pd.read_csv(args.overlapping_genes[0])
#Get data
X, markers, ages, age_indices = format_probes(joined_betas, sample_sheet, gene_annotations, outdir)
plot_probes(X,markers,ages,age_indices,overlapping_probes)
