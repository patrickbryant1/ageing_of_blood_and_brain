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
from scipy.interpolate import Rbf #Radian basis function

from collections import Counter
import pdb



#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''Analyze all statistically significant methylation markers btw age groups''')
parser.add_argument('--joined_betas', nargs=1, type= str, default=sys.stdin, help = 'Path to joined_betas.')
parser.add_argument('--sample_sheet36194', nargs=1, type= str, default=sys.stdin, help = 'Path to sample sheet with accession 36194.')
parser.add_argument('--sample_sheet1575', nargs=1, type= str, default=sys.stdin, help = 'Path to sample sheet with accession 1575.')
parser.add_argument('--gene_annotations', nargs=1, type= str, default=sys.stdin, help = 'Path to gene annotations.')
parser.add_argument('--outdir', nargs=1, type= str, default=sys.stdin, help = 'Path to outdir.')


#######FUNCTIONS#######
def get_ages(joined_betas, sample_sheet1575, sample_sheet36194):
    '''Get the age for each participant
    '''
    sample_names = joined_betas.columns[2:]
    sample_ages = []
    sample_sexes = []
    sample_tissues = []
    for name in sample_names:
        sheet, id = name.split('_')
        if sheet=='E-GEOD-36194':
            sample_sheet = sample_sheet36194
            agelabel='Characteristics[age (y)]'
            sexlabel='Characteristics[sex]'
            tissue_label='Characteristics[organism part]'
        else:
            sample_sheet = sample_sheet1575
            agelabel='Characteristics[age]'
            sexlabel='Characteristics[gender]'
            tissue_label='Characteristics[tissue]'

        #Save sex
        sample_sexes.append(sample_sheet[sample_sheet['Source Name']==id+' 1'][sexlabel].values[0])
        #Save tissue type
        sample_tissues.append(sample_sheet[sample_sheet['Source Name']==id+' 1'][tissue_label].values[0])
        #Save age
        try:
            sample_ages.append(float(sample_sheet[sample_sheet['Source Name']==id+' 1'][agelabel].values[0]))

        except:
            if sample_sheet[sample_sheet['Source Name']==id+' 1'][agelabel].values[0] == '>90':
                sample_ages.append(90.0)
            else:
                pdb.set_trace()
    return sample_ages, sample_sexes, sample_tissues

def get_point_indices(ages):


    #Get unique ages
    unique_ages = np.unique(ages)
    unique_ages = np.sort(unique_ages)

    #Go through each unique age and see how many there are of each
    num_per_age = np.zeros(len(unique_ages))
    indices_per_age = [] #Save age indices per age
    for i in range(len(unique_ages)):
        #Check where ages = ages[i]
        num_per_age[i]=len(np.where(ages==unique_ages[i])[0])
        indices_per_age.append(np.where(ages==unique_ages[i])[0])


    #Get 10% of points for each year
    mina = min(unique_ages)
    maxa = max(unique_ages)
    point_indices = [] #save point indices
    target = int(len(ages)*0.1) #Number of points to fetch for each age
    print('Target:', target)
    for age in np.arange(mina,maxa+1):
        offset=0
        num_fetched = 0
        fetched_indices = [] #Save the indices fetched per age
        while num_fetched<target:
            #Start with the exact age match
            pos_match_i = np.where(unique_ages==age+offset)[0]
            neg_match_i = np.where(unique_ages==age-offset)[0]

            #Count the number fetched
            if len(pos_match_i)>0 and offset!=0: #don't want to sample twice on offset=0
                num_fetched+=num_per_age[pos_match_i][0]
                pos_indices = indices_per_age[pos_match_i[0]]
            else:
                pos_indices = np.array([]) #Add empty array

            if len(neg_match_i)>0:
                num_fetched+=num_per_age[neg_match_i][0]
                neg_indices = indices_per_age[neg_match_i[0]]
            else:
                neg_indices = np.array([]) #Add empty array
            #All indices
            all_indices = np.concatenate([pos_indices,neg_indices])

            #Want the same number of points per age - adjustments are therefore needed
            if num_fetched>target:
                diff = int(num_fetched-target)
                #Take away points so the target is not exceeded.
                #By making a random selection on all indices
                remove_indices = np.random.choice(all_indices,diff,replace=False)
                keep_indices = np.setdiff1d(all_indices,remove_indices)
                all_indices = keep_indices
            #Save indices
            fetched_indices.extend(all_indices)

            #Increase offset
            offset+=1

        point_indices.append(fetched_indices)

    return point_indices

def compare_probes(joined_betas, sample_sheet1575, sample_sheet36194, gene_annotations, outdir):
    '''Analyze the relationship between probe beta values and age
    '''


    #Merge on probe id
    merged = pd.merge(joined_betas, gene_annotations, left_on='Reporter Identifier', right_on='Name',how='left')
    #Check how many are zeros (unquantified, beta can't be 0)
    zeros = (merged[merged.columns[2:-3]] == 0).astype(int).sum(axis=1)
    zero_indices = np.where(zeros<10)
    #Get ages
    ages, sexes, sample_tissues = get_ages(joined_betas, sample_sheet1575, sample_sheet36194)

    #Round ages
    ages = np.round(ages)

    #Save ages
    age_df = pd.DataFrame()
    age_df['Sample'] = joined_betas.columns[2:]
    age_df['Age'] = ages
    age_df['Sex'] = sexes
    age_df['Tissue'] = sample_tissues
    age_df['Tissue']=age_df['Tissue'].replace('FCTX','frontal cortex') #Rename FCTX to frontal cortex
    age_df.to_csv(outdir+'ages.csv')

    #Looking at single marker comparisons
    markers = np.array(merged['Reporter Identifier'])
    markers = markers[zero_indices]


    #Methylation values
    X = np.array(merged[merged.columns[2:-3]])
    #Check how many are zeros (unquantified, beta can't be 0)
    print(np.round(100*X[X==0].shape[0]/(X.shape[0]*X.shape[1]),3), '% zeros')
    #Take all samples with less than 10 zeros
    X = X[zero_indices,:][0]

    print('Removed ', len(merged)-len(X), 'markers that had over 10 missing values (Beta=0)')
    #Go through the tissues
    for tissue in [ 'cerebellum','frontal cortex']:
        #Get frontal cortex ages
        tissue_indices = age_df[age_df['Tissue']==tissue].index
        tissue_ages = ages[tissue_indices]
        print(tissue, len(tissue_indices),'samples')

        #Save the tissue ages
        age_df.loc[tissue_indices].to_csv(outdir+tissue+'_ages.csv')

        #Get point indices
        point_indices = get_point_indices(tissue_ages)
        #Save point_indices
        np.save(outdir+tissue+'_age_points.npy',np.array(point_indices))

        #Get tissue marker values
        X_tissue = X[:,tissue_indices]
        #Save X
        np.save(outdir+tissue+'_marker_values.npy',X_tissue)
        #Min and max age
        minage = min(tissue_ages)
        maxage = max(tissue_ages)
        point_ages = np.arange(minage,maxage+1)
        #Save running averages
        running_averages = np.zeros((X_tissue.shape[0],len(point_ages)))
        max_fold_changes = np.zeros(X_tissue.shape[0])
        max_fold_change_pvals = np.zeros(X_tissue.shape[0])
        #Calculate running point average
        for xi in range(len(X_tissue)):
            if xi%1000==0: #Print if congruent with 1000
                print(xi)
            Xsel = X_tissue[xi,:]
            #Go through all point indices
            for pi in range(len(point_indices)):
                age_points = point_indices[pi]
                running_averages[xi,pi]=np.average(Xsel[np.array(age_points,dtype='int32')])

            maxi = np.where(running_averages[xi,:]==max(running_averages[xi,:]))[0][0]
            mini = np.where(running_averages[xi,:]==min(running_averages[xi,:]))[0][0]
            max_fold_changes[xi] = running_averages[xi,maxi]/running_averages[xi,mini]

            #Calculate p-value between samples belonging to max/min fold change
            xmax = Xsel[np.array(point_indices[maxi],dtype='int32')]
            xmin = Xsel[np.array(point_indices[mini],dtype='int32')]
            stat, pval = ttest_ind(xmax,xmin)
            max_fold_change_pvals[xi]=pval

        #Save running averages and fold changes
        np.save(outdir+tissue+'_running_averages.npy', running_averages)
        np.save(outdir+tissue+'_max_fold_changes.npy', max_fold_changes)
        np.save(outdir+tissue+'_max_fold_change_pvals.npy', max_fold_change_pvals)
        df = pd.DataFrame()

        df['Reporter Identifier']=markers
        df['p']=max_fold_change_pvals
        df['fold_change']=max_fold_changes
        #Save df
        df.to_csv(outdir+tissue+'_marker_max_FC_pval.csv')

        #Correlate probe values with age
        R = np.zeros(X_tissue.shape[0])
        p = np.zeros(X_tissue.shape[0])
        for xi in range(X_tissue.shape[0]):
            if xi%1000==0:
                print(xi)
            R[xi], p[xi] = pearsonr(X_tissue[xi,:], tissue_ages)

        #Save marker-age correlations
        corr_df = pd.DataFrame()
        corr_df['Reporter Identifier']=markers
        corr_df['R']=R
        corr_df['p']=p
        corr_df.to_csv(outdir+tissue+'_correlation_results.csv')





    return None

###########MAIN###########
#Plt
plt.rcParams.update({'font.size': 7})
#Seed
np.random.seed(42) #Answer to everything
#Args
args = parser.parse_args()
gene_annotations = pd.read_csv(args.gene_annotations[0],low_memory=False)
joined_betas = pd.read_csv(args.joined_betas[0], low_memory=False)
joined_betas = joined_betas.fillna(0) #Fill nans with 0
print('Read betas')
sample_sheet36194 = pd.read_csv(args.sample_sheet36194[0], sep = '\t')
#The sex annotation is screqed up, this fixes that
sample_sheet36194['Characteristics[sex]'][470:]=np.array(sample_sheet36194['Characteristics[sex].1'][470:])
sample_sheet1575 = pd.read_csv(args.sample_sheet1575[0], sep = '\t')
outdir = args.outdir[0]

#Compare probes btw age stratified samples
compare_probes(joined_betas,  sample_sheet1575, sample_sheet36194, gene_annotations, outdir)
