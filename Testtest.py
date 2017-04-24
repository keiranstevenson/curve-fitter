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

a = np.array(cols)
for i in range(cols.size): 
    if re.match('Sample*',cols[i]):
        print x
        a[i]= False
    else:
        a[i]= True
        
