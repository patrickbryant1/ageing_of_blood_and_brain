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
parser.add_argument('--hannum_markers', nargs=1, type= str, default=sys.stdin, help = 'Path to hannum marker values.')
parser.add_argument('--max_fold_change_df', nargs=1, type= str, default=sys.stdin, help = 'Path to marker max fold changes and pvals. This df contains the marker ids in order.')
parser.add_argument('--mode', nargs=1, type= str, default=sys.stdin, help = 'Mode.')
parser.add_argument('--outdir', nargs=1, type= str, default=sys.stdin, help = 'Path to outdir.')

def rf_fit(sel_marker_values, ages, hannum_preds, horvath_preds, mode):
    '''5 fold CV
    '''


    fig,ax = plt.subplots(figsize=(6/2.54, 6/2.54))
    #Plot hannum
    hannum_error = np.average(np.absolute(ages-hannum_preds))
    hannum_corr = pearsonr(ages,hannum_preds)[0]
    plt.scatter(ages,hannum_preds,color='magenta',s=1,label='Hannum')
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
    titles = {'FC':'FC','abs':'Absolute value','overlap':'Overlap'}
    plt.plot([min(ages),max(ages)],[min(ages),max(ages)],color='k',linewidth=0.5)
    plt.xlabel('True age')
    plt.ylabel('Predicted age')
    plt.title(titles[mode])
    plt.legend(frameon = False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig(outdir+mode+'_cv_results.png', format='png', dpi=300)
    plt.close()

    print('Hannum',hannum_error,hannum_corr)
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
hannum_markers = pd.read_csv(args.hannum_markers[0])
max_fold_change_df = pd.read_csv(args.max_fold_change_df[0])
mode = args.mode[0]
outdir = args.outdir[0]

#Get Horvath ages
horvath_preds = pd.read_csv(outdir+'horvath_ages.csv')
horvath_preds=pd.merge(age_df,horvath_preds,on='Sample',how='left')
#Select marker values
#the column Unnamed: 0_x contains the indx of the selected marker
sel_marker_values = marker_values[selected_markers['Unnamed: 0_x'].values]
#Select hannum marker values
hannum_markers =  pd.merge(hannum_markers,max_fold_change_df,left_on='Marker',right_on='Reporter Identifier', how='left').dropna()
hannum_marker_values = marker_values[np.array(hannum_markers['Unnamed: 0'].values,dtype='int')]
hannum_coefs = hannum_markers['Coefficient'].values
hannum_preds = np.dot(hannum_marker_values.T, hannum_coefs)

#Select ages
ages = age_df['Age'].values
#Fit a rf model

rf_fit(sel_marker_values.T, ages, hannum_preds,horvath_preds.DNAmAge.values, mode)
