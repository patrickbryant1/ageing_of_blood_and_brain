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

import pdb



#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''Visualize the gene ontlogy''')
parser.add_argument('--ra_GO', nargs=1, type= str, default=sys.stdin, help = 'Path to running average GO.')
parser.add_argument('--hannum_GO', nargs=1, type= str, default=sys.stdin, help = 'Path to Hannum GO.')
parser.add_argument('--outdir', nargs=1, type= str, default=sys.stdin, help = 'Path to outdir.')


#######FUNCTIONS#######

def plot_GO(go_df):
    '''Plot the GO as bar charts
    '''

    #Plot
    fig,ax = plt.subplots(figsize=(6/2.54, 6/2.54))

    sns.distplot(sig_correlation_results['R'],color='cornflowerblue', label='Significant correlations')

    plt.ylabel('Count')
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
ra_GO = pd.read_csv(args.ra_GO[0],low_memory=False)
hannum_GO = pd.read_csv(args.hannum_GO[0],low_memory=False)
outdir = args.outdir[0]

pdb.set_trace()
#Plot ra GO
plot_GO(ra_go)
