# -*- coding: utf-8 -*-
"""
Created on Wed Jan 20 18:17:34 2016

@author: s1002426
"""

import numpy as np
import genutils as gu
import pandas as pd
import string as string
from fitderiv import fitderiv

namelist = ['LB','MV36','MV26','osmoLB','osmoMV36','osmoMV26'];

Gr = np.zeros(shape=(60,1))
Err = np.zeros(shape=(60,1))

for i in range(1,7):
    for ii in range(1,11):
#filename = "LB_0.csv"
        firstname = namelist[i-1]
        iii = ii-1
        filename = str(firstname)+'_'+ str(iii)+'.csv'
        ld= np.genfromtxt(filename, delimiter=',')
        t= ld[0,:]
        d= np.transpose(ld[1:,:])
        fitty = fitderiv(t,d,bd= {0: [-5,8],1: [-6,-1], 2: [-5,2]},exitearly= True)

        ind = 10*(i-1) + ii
        Gr[ind-1] =fitty.ds['max df']
        Err[ind-1] = fitty.ds['max df var']
        print([ind])






