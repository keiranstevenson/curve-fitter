# -*- coding: utf-8 -*-
"""
Created on Mon Nov  7 19:05:42 2016

@author: S1002426
"""

import numpy as np
import pandas as pd
from fitderiv import fitderiv
from scipy import signal
import matplotlib.pyplot as plt


indata = pd.read_csv('Datafiles/160705_2308.CSV',header= None,skiprows= 6)
deletewells = 0
labelcols = 3

cols = indata.iloc[:,0]
colindex = (cols == 'A') | (cols == 'H')
cols = indata.iloc[:,1]
colindex = (cols == 1) | (cols == 12) | colindex
if deletewells == 1:
    colindex = colindex == False
    data = indata.loc[colindex]
    data = data.copy()
else:        
    data = indata.copy()
    for i in range(0,len(colindex)):
        if 1 == colindex[i]:
            data.iloc[[i],3:] = np.zeros(data.shape[1]-3)
