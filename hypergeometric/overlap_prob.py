#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import sys
import os
import glob
import numpy as np
import math

import pdb

def striling_aprroximation(n):
    '''Calculate the approzimation of n! according to the
    striling formula
    n! = sqrt(2*pi*n)*(n/e)^n
    '''
    return np.sqrt(2*math.pi*n)*(n/math.e)**n

    
