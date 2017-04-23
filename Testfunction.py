# -*- coding: utf-8 -*-
"""
Created on Fri Apr 21 17:12:28 2017

@author: s1002426
"""

from glob import glob
import pandas as pd
import numpy as np
import os
import platereaderimport as pltim


filename = 'Datafiles/truncatedtest.CSV'
skiprows= 6
labelcols= 3
header= None

infile = pd.read_csv(filename, header=header, skiprows=skiprows)
normvalue= 0.05
alignvalue= 0.1

pltim.normalisetraces(infile)
dataset = np.array(infile.iloc[2:11,labelcols:])
alignpoint= normvalue+alignvalue
startindexes = np.int_(dataset[:,0])

for i in range(0,dataset.shape[0]):
    for ii in range(0,dataset.shape[1]):
        if (dataset[i,ii]> alignpoint) & (dataset[i,ii+1]> alignpoint) & (dataset[i,ii+2]> alignpoint):
            x = ii
            break
    startindexes[i] = np.int_(x)
    
minindex = startindexes.min()
maxindex = startindexes.max()
maxrowlength = dataset.shape[1]-maxindex-1
for i in range(0,dataset.shape[0]):
    rowindex = startindexes[i]
    newdata = dataset[i,rowindex-minindex:rowindex+maxrowlength-1]
    if i == 0:
        stack= newdata
    else:
        stack= np.vstack((stack,newdata))


