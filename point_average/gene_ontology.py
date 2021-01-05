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
parser.add_argument('--blood1', nargs=1, type= str, default=sys.stdin, help = 'Path to blood GO1.')
parser.add_argument('--blood2', nargs=1, type= str, default=sys.stdin, help = 'Path to blood GO2.')
parser.add_argument('--frontal_cortex', nargs=1, type= str, default=sys.stdin, help = 'Path to FCTX GO.')
parser.add_argument('--cerebellum', nargs=1, type= str, default=sys.stdin, help = 'Path to cerebellum GO.')
parser.add_argument('--horvath', nargs=1, type= str, default=sys.stdin, help = 'Path to Horvath GO.')
parser.add_argument('--hannum', nargs=1, type= str, default=sys.stdin, help = 'Path to Hannum GO.')
parser.add_argument('--outdir', nargs=1, type= str, default=sys.stdin, help = 'Path to outdir.')


#######FUNCTIONS#######

def plot_GO(go_df,colors,all_categories,title,outname):
    '''Plot the GO as bar charts
    '''
    #Get colors
    sel_colors = []
    for key in go_df[1]:
        sel_colors.append(colors[np.where(all_categories==key)[0]][0])

    #Plot
    fig,ax = plt.subplots(figsize=(18/2.54, 12/2.54))

    plt.bar(np.arange(len(go_df)),go_df[2],color=np.array(sel_colors))
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
blood1 = pd.read_csv(args.blood1[0],sep='\t', header=None)
blood2 = pd.read_csv(args.blood2[0],sep='\t', header=None)
frontal_cortex = pd.read_csv(args.frontal_cortex[0],sep='\t', header=None)
cerebellum = pd.read_csv(args.cerebellum[0],sep='\t', header=None)
horvath = pd.read_csv(args.horvath[0],sep='\t', header=None)
hannum = pd.read_csv(args.hannum[0],sep='\t', header=None)
outdir = args.outdir[0]

all_categories = np.concatenate([np.array(blood1[1].unique()),np.array(blood2[1].unique()),
                np.array(frontal_cortex[1].unique()),
                np.array(cerebellum[1].unique()),np.array(horvath[1].unique()),np.array(hannum[1].unique())])
all_categories = np.unique(all_categories)

#Plot ra GO
colors = pl.cm.tab20b(np.linspace(0,1,len(all_categories)))


plot_GO(blood1,colors,all_categories,'Blood GO enrichment 1',outdir+'blood/genes/go1.png')
plot_GO(blood2,colors,all_categories,'Blood GO enrichment 2',outdir+'blood/genes/go2.png')
plot_GO(frontal_cortex,colors,all_categories,'Frontal cortex GO enrichment',outdir+'brain/frontal_cortex/genes/go.png')
plot_GO(cerebellum,colors,all_categories,'Cerebellum GO enrichment',outdir+'brain/cerebellum/genes/go.png')

plot_GO(hannum,colors,all_categories,'Hannum GO enrichment', outdir+'blood/genes/Hannum/go.png')
plot_GO(horvath,colors,all_categories,'Horvath GO enrichment', outdir+'brain/go.png')
