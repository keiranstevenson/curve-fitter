# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 16:09:32 2016

@author: s1002426
"""

import numpy as np
import pandas as pd
from fitderiv import fitderiv
from scipy import signal

infile = pd.read_csv('MV26-low-1.csv')

growthrates = infile.ix[1:,0:7]

time = infile.iloc[0]
time = time[3:]
time = np.float64(time)

for i in range(1,61):
    print(i)    
    
    od = infile.iloc[i]
    
    od = od[3:]
    od = od+0.1
    
#    od2 = infile.iloc[10+i]
#    od2 = od2[3:]
#    od2 = od2+0.1
    
#    od3 = infile.iloc[20+i]
#    od3 = od3[3:]
#    od3 = od3+0.1
    
#    odboth = pd.concat([od,od2,od3],axis=1)
    odboth2 = np.float64(od)  
    
    # Removes nans from data that matlab put in
    odboth3 =  [ num for num in odboth2 if (False ==(np.isnan(num)).any())] 
    t = time[:len(odboth3)]
   

    odboth3 = signal.medfilt(odboth3,5)
    
    # Check for growth
    maxi = np.amax(od)
    mini = np.amin(od)
    diff = maxi-mini
    if diff > 0.05: 
        # Runs fitderiv only if growth is over 0.05
        fitty = fitderiv(t,odboth2,bd= {0: [-5,8],1: [-6,-1], 2: [-5,2]},exitearly= False,nosamples= 20) # Peters program
        Gr = fitty.ds['max df']
        Err = fitty.ds['max df var']
        Lag = fitty.ds['lag time']
        TimeofmaxGr = fitty.ds['time of max df']
        note = ["Normal"]
        
    else:
        # Returns no growth if none detected
        Gr = 0
        Err = 0
        Lag = 0
        note = ["No Growth"]
        TimeofmaxGr = 0
    
    
    growthrates.ix[i,2] = note    
    growthrates.ix[i,3] = Gr
    growthrates.ix[i,4] = Err
    growthrates.ix[i,5] = Lag
    growthrates.ix[i,6] = TimeofmaxGr