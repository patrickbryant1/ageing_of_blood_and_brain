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
from sklearn.metrics import mutual_info_score
import pdb



#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''Analyze all statistically significant methylation markers btw age groups''')
parser.add_argument('--joined_betas', nargs=1, type= str, default=sys.stdin, help = 'Path to joined_betas.')
parser.add_argument('--sample_sheet', nargs=1, type= str, default=sys.stdin, help = 'Path to sample sheet.')
parser.add_argument('--gene_annotations', nargs=1, type= str, default=sys.stdin, help = 'Path to gene annotations.')
parser.add_argument('--agelabel', nargs=1, type= str, default=sys.stdin, help = 'Agelabel.')
parser.add_argument('--median_range', nargs=1, type=str, default=sys.stdin, help = 'Range for which group median is the current age.')
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
            if sample_sheet[sample_sheet['Source Name']==name+' 1'][agelabel].values[0] == '>90':
                sample_ages.append(90.0)
            else:
                pdb.set_trace()

    return np.array(sample_ages)

def clean_outliers(X, outdir, tissue):
    '''Remove the outlier samples by investigating the entropy btw the mean beta value distribution
    and each sample's beta value distribution
    '''

    #Calculate KL-divergence
    entropies = []
    mean_beta_distr = np.average(X, axis = 1)
    for i in range(X.shape[1]):
        entropies.append(mutual_info_score(mean_beta_distr, X[:,i]))
    entropies = np.array(entropies)
    #Plot
    sns.distplot(entropies)
    plt.savefig(outdir+tissue+'_entropies.png', format='png')
    plt.close()
    #Plot betas, color in deviations
    dev_samples = np.where(entropies<(np.mean(entropies)-(np.std(entropies)*3)))[0] #select all 3 stds away
    for i in range(X.shape[1]):
        if i in dev_samples:
            color = 'r'
        else:
            color = 'b'
        sns.distplot(X[:,i], color = color, hist = False, kde_kws={"lw": 1},)
    plt.title('Beta values colored by mutual information to mean.\nblue=kept|red=removed')
    plt.savefig(outdir+tissue+'_betaplot.png', format='png', dpi=300)

    print(len(dev_samples), 'samples had MI 3 stds aeay from average.')
    return np.where(entropies>=(np.mean(entropies)-(np.std(entropies)*3)))[0]

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

def compare_probes(joined_betas, sample_sheet, gene_annotations, median_range, outdir):
    '''Analyze the relationship between probe beta values and age
    '''


    #Merge on probe id
    merged = pd.merge(joined_betas, gene_annotations, left_on='Reporter Identifier', right_on='Name')
    #Check how many are zeros (unquantified, beta can't be 0)
    zeros = (merged[merged.columns[2:-34]] == 0).astype(int).sum(axis=1)
    zero_indices = np.where(zeros<10)
    #Get ages
    ages = get_ages(sample_sheet, joined_betas.columns[2:], agelabel)

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

    #Clean outliers
    remain_indices = clean_outliers(X, outdir, 'blood')
    X = X[:,remain_indices]
    ages = ages[remain_indices]
    #Get point indices
    point_indices = get_point_indices(ages)
    #Save point_indices
    np.save(outdir+'age_points.npy',np.array(point_indices))

    #Save ages
    age_df = pd.DataFrame()
    age_df['Sample'] = joined_betas.columns[2:][remain_indices]
    age_df['Age'] = ages
    age_df.to_csv(outdir+'ages.csv')

    #Save X
    np.save(outdir+'marker_values.npy',X)
    #Min and max age
    minage = min(ages)
    maxage = max(ages)
    point_ages = np.arange(minage,maxage+1)
    #Save running averages
    running_averages = np.zeros((X.shape[0],len(point_ages)))
    max_fold_changes = np.zeros(X.shape[0])
    max_fold_change_pvals = np.zeros(X.shape[0])
    max_abs_changes = np.zeros(X.shape[0])
    #Calculate running point average
    for xi in range(len(X)):
        if xi%1000==0: #Print if congruent with 1000
            print(xi)
        Xsel = X[xi,:]
        #Go through all point indices
        for pi in range(len(point_indices)):
            age_points = point_indices[pi]
            running_averages[xi,pi]=np.median(Xsel[np.array(age_points,dtype='int32')])

        maxi = np.where(running_averages[xi,:]==max(running_averages[xi,:][median_range[0]-20:median_range[1]-19]))[0][0] #Starts at year 19 --> first index is age-20
        mini = np.where(running_averages[xi,:]==min(running_averages[xi,:][median_range[0]-20:median_range[1]-19]))[0][0]
        max_fold_changes[xi] = running_averages[xi,maxi]/running_averages[xi,mini]
        max_abs_changes[xi] = running_averages[xi,maxi]-running_averages[xi,mini]
        #Calculate p-value between samples belonging to max/min fold change
        xmax = Xsel[np.array(point_indices[maxi],dtype='int32')]
        xmin = Xsel[np.array(point_indices[mini],dtype='int32')]
        stat, pval = ttest_ind(xmax,xmin)
        max_fold_change_pvals[xi]=pval

    #Save running averages and fold changes
    np.save(outdir+'running_averages.npy', running_averages)
    np.save(outdir+'max_fold_changes.npy', max_fold_changes)
    np.save(outdir+'max_fold_change_pvals.npy', max_fold_change_pvals)
    df = pd.DataFrame()

    df['Reporter Identifier']=markers
    df['p']=max_fold_change_pvals
    df['fold_change']=max_fold_changes
    df['abs_change']=max_abs_changes
    #Save df
    df.to_csv(outdir+'marker_max_FC_pval.csv')


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
#Seed
np.random.seed(42) #Answer to everything
#Args
args = parser.parse_args()
agelabels = {'blood':"Characteristics [age y]"}
agelabel=agelabels[args.agelabel[0]]
gene_annotations = pd.read_csv(args.gene_annotations[0],low_memory=False)
joined_betas = pd.read_csv(args.joined_betas[0], low_memory=False)
print('Read betas')
sample_sheet = pd.read_csv(args.sample_sheet[0], sep = '\t')
median_range = np.array(args.median_range[0].split(','),dtype='int32')
outdir = args.outdir[0]

#Compare probes btw age stratified samples
compare_probes(joined_betas, sample_sheet, gene_annotations, median_range, outdir)
