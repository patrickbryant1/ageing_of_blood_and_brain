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

def rf_fit(sel_marker_values, ages, horvath_preds, outdir):
    '''5 fold CV
    '''

    kf = KFold(n_splits=5, random_state=42, shuffle=True)
    fig,ax = plt.subplots(figsize=(6/2.54, 6/2.54))
    errors = []
    for ti, vi in kf.split(sel_marker_values):
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
        plt.scatter(y_valid,pred,s=1,color='cornflowerblue')
    #Plot Horvath
    plt.scatter(ages,horvath_preds,color='r',s=1,label='Horvath')
    plt.xlabel('True age')
    plt.ylabel('Predicted age')
    plt.title('Average error '+str(np.round(np.average(errors),2))+' +/- '+str(np.round(np.std(errors),2)))
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir+'cv_results.png', format='png', dpi=300)
    plt.close()


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
#Select Horvath marker values
horvath_markers =  pd.merge(horvath_markers,max_fold_change_df,left_on='Marker',right_on='Reporter Identifier', how='left')
horvath_marker_values = marker_values[horvath_markers['Unnamed: 0'].values]
horvath_coefs = horvath_markers['Coefficient'].values
horvath_preds = np.dot(horvath_marker_values.T, horvath_coefs)

#Select ages
ages = age_df['Age'].values
#Fit a rf model
rf_fit(sel_marker_values.T, ages, horvath_preds, outdir)
