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

parser.add_argument('--outdir', nargs=1, type= str, default=sys.stdin, help = 'Path to outdir.')

def rf_fit(sel_marker_values, ages, outdir):
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
        plt.scatter(pred,y_valid,s=1,color='cornflowerblue')
    plt.xlabel('Predicted age')
    plt.ylabel('True age')
    plt.title('Average error '+str(np.round(np.average(errors),2))+' +/- '+str(np.round(np.std(errors),2)))
    plt.tight_layout()
    plt.savefig(outdir+'cv_results.png', format='png', dpi=300)
    plt.close()


###########MAIN###########
#Plt
plt.rcParams.update({'font.size': 7})
#Args
args = parser.parse_args()
selected_markers = pd.read_csv(args.selected_markers[0],low_memory=False) #the column Unnamed: 0_x contains the indx of the selected marker
marker_values = np.load(args.marker_values[0], allow_pickle=True)
age_df = pd.read_csv(args.ages[0])
outdir = args.outdir[0]

#Select marker values
sel_marker_values = marker_values[selected_markers['Unnamed: 0_x'].values]
#Select ages
ages = age_df['Age'].values
#Fit a rf model
rf_fit(sel_marker_values.T, ages, outdir)
