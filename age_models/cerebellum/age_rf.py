#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import sys
import os
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import KFold
from scipy.stats import pearsonr

import numpy as np



import pdb
#Arguments for argparse module
parser = argparse.ArgumentParser(description = '''Create a rf predictor for CA using marker values from selected markers''')
parser.add_argument('--selected_markers', nargs=1, type= str, default=sys.stdin, help = 'Path to df with selected markers.')
parser.add_argument('--marker_values', nargs=1, type= str, default=sys.stdin, help = 'Path to marker values.')
parser.add_argument('--ages', nargs=1, type= str, default=sys.stdin, help = 'Path to sample ages.')
parser.add_argument('--horvath_markers', nargs=1, type= str, default=sys.stdin, help = 'Path to Horvath marker values.')
parser.add_argument('--max_fold_change_df', nargs=1, type= str, default=sys.stdin, help = 'Path to marker max fold changes and pvals. This df contains the marker ids in order.')
parser.add_argument('--outdir', nargs=1, type= str, default=sys.stdin, help = 'Path to outdir.')

def create_horvath_df(age_df, all_marker_values,all_marker_ids, horvath_markers, outdir):
    '''Get Horvath marker values
    '''

    samples = age_df['Sample'].values
    df_to_horvath_clock = pd.DataFrame()
    df_to_horvath_clock['Name']=all_marker_ids
    for i in range(len(samples)):
        df_to_horvath_clock[samples[i]]=all_marker_values[:,i]


    #Merge
    df_to_horvath_clock = pd.merge(horvath_markers,df_to_horvath_clock,on='Name',how='left')
    #Drop unwanted cols
    df_to_horvath_clock = df_to_horvath_clock.drop(columns={'Gene_ID', 'GenomeBuild', 'Chr', 'Accession', 'overallMeanByCpGacross50data', 'CoefficientHannum'})

    #Save df
    df_to_horvath_clock.to_csv(outdir+'df_to_horvath_clock.csv',index=False)

    return None


def rf_fit(sel_marker_values, ages, horvath_preds, outdir):
    '''5 fold CV
    '''

    fig,ax = plt.subplots(figsize=(6/2.54, 6/2.54))
    #Plot Horvath
    horvath_error = np.average(np.absolute(ages-horvath_preds))
    horvath_corr = pearsonr(ages,horvath_preds)[0]
    plt.scatter(ages,horvath_preds,color='cornflowerblue',s=1,label='Horvath')
    #5-fold CV
    kf = KFold(n_splits=5, random_state=42, shuffle=True)
    errors = []
    corrs = []
    fold = 0
    for ti, vi in kf.split(sel_marker_values):
        fold+=1 #Increase fold
        X_train = sel_marker_values[ti]
        y_train = ages[ti]
        X_valid = sel_marker_values[vi]
        y_valid = ages[vi]
        #Define model
        #All default params are used - no optimization is performed
        regr = RandomForestRegressor(n_jobs=-1, random_state=42)
        regr.fit(X_train, y_train)
        pred = regr.predict(X_valid)
        errors.append(np.average(np.absolute(pred-y_valid)))
        corrs.append(pearsonr(pred,y_valid)[0])
        if fold ==5:
            plt.scatter(y_valid,pred,s=1,color='darkgreen',label='Random forest',alpha=0.5)
        else:
            plt.scatter(y_valid,pred,s=1,color='darkgreen',alpha=0.5)

    #Plot diagonal line
    plt.plot([min(ages),max(ages)],[min(ages),max(ages)],color='k',linewidth=0.5)
    plt.xlabel('True age')
    plt.ylabel('Predicted age')
    plt.title('Cerebellum')
    plt.legend(frameon = False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig(outdir+'cv_results.png', format='png', dpi=300)
    plt.close()

    print('Horvath',horvath_error,horvath_corr)
    print('Random forest',np.average(errors),np.std(errors),np.average(corrs),np.std(corrs))

###########MAIN###########
#Plt
plt.rcParams.update({'font.size': 7})
#Args
args = parser.parse_args()
selected_markers = pd.read_csv(args.selected_markers[0],low_memory=False)
marker_values = np.load(args.marker_values[0], allow_pickle=True)
age_df = pd.read_csv(args.ages[0])
horvath_markers = pd.read_csv(args.horvath_markers[0])
max_fold_change_df = pd.read_csv(args.max_fold_change_df[0])
outdir = args.outdir[0]

#Select marker values
#the column Unnamed: 0_x contains the indx of the selected marker
sel_marker_values = marker_values[selected_markers['Unnamed: 0_x'].values]

try:
    horvath_preds = pd.read_csv(outdir+'norm_df_to_horvath_clock.output.csv')
except:
    #Format Horvath markers to run clock
    create_horvath_df(age_df, marker_values,max_fold_change_df['Reporter Identifier'].values, horvath_markers,outdir)
    raise IOError('No Horvath preds')

#Select ages
ages = age_df['Age'].values
#Fit a rf model
rf_fit(sel_marker_values.T, ages, horvath_preds['DNAmAge'].values, outdir)
