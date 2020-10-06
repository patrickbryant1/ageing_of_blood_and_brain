#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import sys
import os
import glob
import numpy as np
import math

import pdb

#Arguments for argparse module:
parser = argparse.ArgumentParser(description = '''Analyze all statistically significant methylation markers btw age groups''')
parser.add_argument('--x', nargs=1, type= int, default=sys.stdin, help = 'Number of overlapping markers.')
parser.add_argument('--n', nargs=1, type= int, default=sys.stdin, help = 'Number of markers in group 1.')
parser.add_argument('--D', nargs=1, type= int, default=sys.stdin, help = 'Number of markers in group 2.')
parser.add_argument('--N', nargs=1, type= int, default=sys.stdin, help = 'Number of total markers')

#####FUNCTIONS#####
def expected_n(n,D,N):
    '''Calculate the expected number of genes
     p(A&B)=p(A)*p(B|A)=n/N*D/N=(n*D)/N
    '''
    return (n*D)/N

def representation_factor(exp_n,x):
    '''Calculate the representation factor,
    the overlap compared to the expected overlap
    x / expected # of genes
    '''

    return x/exp_n

def striling_aprroximation(n):
    '''Calculate the approzimation of n! according to the
    striling formula
    n! = sqrt(2*pi*n)*(n/e)^n
    '''
    return np.sqrt(2*math.pi*n)*(n/math.e)**n

def comb_a_b(a,b):
    '''Calculate the possible number of ways to choose b from a
    C(a,b) = a! / ((a - b)! * b!)
    '''
    a_fac = striling_aprroximation(a)
    a_b_fac = striling_aprroximation(a-b)
    b_fac = striling_aprroximation(b)
    return a_fac/(a_b_fac*b_fac)

def hypergeometric_pmf(i,n,D,N,C_N_n):
    '''Calculate the probability of i markers from n and D overlapping
    when n and D are selected from N
    p(i) = C(D, i) * C(N-D, n-i) / C(N,n)
    ,where C(a,b) = a! / ((a - b)! * b!)
    '''

    C_Di = comb_a_b(D,i)
    C_ND_ni = comb_a_b(N-D,n-i)

    return C_Di*C_ND_ni/C_N_n


#############MAIN#############
#Args
args = parser.parse_args()
x = args.x[0]
n = args.n[0]
D = args.D[0]
N = args.N[0]
#Calculate the expected overlap and the representation factor
exp_n = expected_n(n,D,N)
rf = representation_factor(exp_n,x)
print('Expexted number of overlapping markers:',exp_n)
print('Representation factor', rf)

#Calculate the probability of x overlapping btw n and D chosen from N
C_N_n = comb_a_b(N,n)
for i in range(x+1):
    pi = hypergeometric_pmf(i,n,D,N,C_N_n)
    pdb.set_trace()
