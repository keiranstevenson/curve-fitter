# -*- coding: utf-8 -*-
"""
Created on Wed Jan 20 18:17:34 2016

@author: s1002426
"""

import numpy 
from fitderiv import fitderiv
import os
from scipy import signal
import pandas

## Location setup and file retrieval
#os.chdir("C:\Users\s1002426\Local Documents\Google Drive\LAB\Python\Peters program 5")
dirfiles = os.listdir()
files = [ num for num in dirfiles if num.endswith(".csv") ]
iterations = len(files)

## ouput setup
gr = numpy.zeros((iterations,2))
notes = [None] * iterations

## access files in sequence
for i in range(0,iterations-1):
    
    file = files[i]
    print(file)
        
    ld= numpy.genfromtxt(file, delimiter=',')
    t= ld[0,:]
    od= numpy.transpose(ld[1:,:])
    od = numpy.ravel(od)
    od = signal.medfilt(od,5)
    
    # Check for growth
    maxi = numpy.amax(od)
    mini = numpy.amin(od)
    diff = maxi-mini
    if diff > 0.1: 
        
        fitty = fitderiv(t,od,bd= {0: [-5,8],1: [-6,-1], 2: [-5,2]},exitearly= False) # Peters program
        gr[i,0] = fitty.ds['max df']
        gr[i,1] = fitty.ds['max df var']
        notes[i] = ["Normal"]
        
    else:
        gr[i,0] = 0
        gr[i,1] = 0
        notes[i] = ["No Growth"]

numpy.savetxt("output.csv", gr, delimiter=",")