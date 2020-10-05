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
import matplotlib.pylab as pl

import pdb



#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''Visualize the gene ontlogy''')
parser.add_argument('--ra_GO', nargs=1, type= str, default=sys.stdin, help = 'Path to running average GO.')
parser.add_argument('--hannum_GO', nargs=1, type= str, default=sys.stdin, help = 'Path to Hannum GO.')
parser.add_argument('--outdir', nargs=1, type= str, default=sys.stdin, help = 'Path to outdir.')


#######FUNCTIONS#######

def plot_GO(go_df,colors,title,outname):
    '''Plot the GO as bar charts
    '''

    #Plot
    fig,ax = plt.subplots(figsize=(18/2.54, 12/2.54))

    plt.bar(np.arange(len(go_df)),go_df[2],color=colors)
    ax.set_xticks(np.arange(len(go_df)))
    ax.set_xticklabels(go_df[1], rotation='vertical')
    plt.ylabel('Count')
    plt.title(title)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig(outname, format='png', dpi=300)
    plt.close()

###########MAIN###########
#Plt
plt.rcParams.update({'font.size': 7})
#Args
args = parser.parse_args()
ra_GO = pd.read_csv(args.ra_GO[0],sep='\t', header=None)
hannum_GO = pd.read_csv(args.hannum_GO[0],sep='\t', header=None)
outdir = args.outdir[0]


#Plot ra GO
colors = pl.cm.tab20b(np.linspace(0,1,len(ra_GO)))
plot_GO(ra_GO,colors[:len(ra_GO)], 'Running average gene ontology enrichment', outdir+'genes/go_ra.png')
colors = pl.cm.tab20b(np.linspace(0,1,len(hannum_GO)))
plot_GO(hannum_GO,colors,'Hannum gene ontology enrichment', outdir+'genes/Hannum/go_hannum.png')
