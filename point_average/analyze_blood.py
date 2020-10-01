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

def point_indices(ages):

    #Get unique ages
    unique_ages = np.unique(ages)
    unique_ages = np.sort(unique_ages)

    #Go through each unique age and see how many there are of each
    num_per_age = np.zeros(len(unique_ages))
    indices_per_age = [] #Save age indices per age
    for i in range(len(unique_ages)):
        #Check where ages = ages[i]
        num_per_age[i]=len(np.where(ages==ages[i])[0])
        indices_per_age.append(np.where(ages==ages[i])[0])


    #Get 5% of points for each year
    mina = min(unique_ages)
    maxa = max(unique_ages)
    point_indices = [] #save point indices
    target = int(len(ages)*0.05) #Number of points to fetch for each age
    for age in np.arange(mina,maxa+1):
        offset=0
        num_fetched = 0
        fetched_indices = [] #Save the indices fetched per age
        while num_fetched<target:
            #Start with the exact age match
            pos_match_i = np.where(unique_ages==age+offset)[0]
            neg_match_i = np.where(unique_ages==age-offset)[0]

            #Count the number fetched
            if pos_match_i and offset!=0: #don't want to sample twice on offset=0
                num_fetched+=num_per_age[pos_match_i]
                pos_indices = indices_per_age[pos_match_i]
            else:
                pos_indices = np.array([]) #Add empty array

            if neg_match_i:
                num_fetched+=num_per_age[neg_match_i])
                neg_indices = indices_per_age[neg_match_i]
            #All indices
            all_indices = np.concatenate([pos_indices,neg_indices])
            pdb.set_trace()
            #Want the same number of points per age - adjustments are therefore needed
            if num_fetched>target:
                diff = num_fetched-target

                #Take away points so the ages are equally represented, meaning removing
                #more points from the bigger group



            #Increase offset
            offset+=1


            pdb.set_trace()




    return age_indices

def compare_probes(joined_betas, sample_sheet, gene_annotations, outdir):
    '''Analyze the relationship between probe beta values and age
    '''


    #Merge on probe id
    merged = pd.merge(joined_betas, gene_annotations, left_on='Reporter Identifier', right_on='Name')
    #Check how many are zeros (unquantified, beta can't be 0)
    zeros = (merged[merged.columns[2:-34]] == 0).astype(int).sum(axis=1)
    zero_indices = np.where(zeros<10)
    #Get ages
    ages = get_ages(sample_sheet, joined_betas.columns[2:], agelabel)
    #Save ages
    age_df = pd.DataFrame()
    age_df['Sample'] = joined_betas.columns[2:]
    age_df['Age'] = ages
    age_df.to_csv(outdir+'ages.csv')

    #Get point indices
    age_indices = point_indices(ages)
    #Save age groups
    np.save(outdir+'age_points.npy',np.array(age_indices))

    #Looking at single marker comparisons
    markers = np.array(merged['Reporter Identifier'])
    markers = markers[zero_indices]


    #Methylation values
    X = np.array(merged[merged.columns[2:-34]])
    #Check how many are zeros (unquantified, beta can't be 0)
    print(np.round(100*X[X==0].shape[0]/(X.shape[0]*X.shape[1]),2), '% zeros')
    #Take all samples with less than 10 zeros
    X = X[zero_indices,:][0]

    print('Removed ', len(merged)-len(X), 'markers that had over 10 missing values (Beta=0)')

    #Compare age groups
    for ai in range(len(age_indices)-1):
        for bi in range(ai+1,len(age_indices)): #Compare all combinations
            print(ai,bi)


            i1 = age_indices[ai]
            i2 = age_indices[bi]
            X1 = X[:,i1]
            X2 = X[:,i2]
            stats, pvals = ttest_ind(X1,X2,axis=1)
            #agesel = np.append(age_indices[ai],age_indices[ai+1])
            #Xsel = X[:,agesel]


            fold_change = np.average(X2,axis=1)/np.average(X1,axis=1)
            #Save
            df = pd.DataFrame()
            df['Reporter Identifier']=markers
            df['stat']=stats
            df['p']=pvals
            df['fold_change']=fold_change
            df['beta1'] = np.average(X1,axis=1)
            df['beta2'] = np.average(X2,axis=1)
            #Save df
            df.to_csv(outdir+str(ai)+'_'+str(bi)+'_age_comparison_results.csv')
    #Correlate probe values with age
    R = np.zeros(X.shape[0])
    p = np.zeros(X.shape[0])
    for xi in range(X.shape[0]):
        if xi%1000==0:
            print(xi)
        R[xi], p[xi] = pearsonr(X[xi,:], ages)

    #Save marker-age correlations
    corr_df = pd.DataFrame()
    corr_df['Reporter Identifier']=markers
    corr_df['R']=R
    corr_df['p']=p
    corr_df.to_csv(outdir+'correlation_results.csv')





    return None

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

#Compare probes btw age stratified samples
compare_probes(joined_betas, sample_sheet, gene_annotations, outdir)
