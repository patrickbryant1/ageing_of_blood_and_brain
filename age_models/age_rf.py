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


###########MAIN###########
#Plt
plt.rcParams.update({'font.size': 7})
#Args
args = parser.parse_args()
selected_markers = pd.read_csv(args.selected_markers[0],low_memory=False)

marker_values = np.load(args.marker_values[0], allow_pickle=True)
age_df = pd.read_csv(args.ages[0])

outdir = args.outdir[0]
