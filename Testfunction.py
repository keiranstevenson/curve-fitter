# -*- coding: utf-8 -*-
"""
Created on Fri Apr 21 17:12:28 2017

@author: s1002426
"""

import platereaderimport as pi
from glob import glob
import pandas as pd
header = None
skiprows = 6
filedirectory = 'C:\Data\Gitstuff\Curve-fitter\Datafiles'


files = glob(filedirectory + '/*.csv')
#files = [file.split('/')[-1] for file in files]
for i in range(0,(len(files)-1)):
    if i == 0:
        stackfile = pd.read_csv(files[i],header= header,skiprows= skiprows)
    else:
        newfile = pd.read_csv(files[i],header= header,skiprows= skiprows+1)
        stack = [stackfile,newfile]
        stackfile = pd.concat(stack)
