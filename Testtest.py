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

infile = indata.iloc[:,0:]
data = indata.iloc[:,3:]

time = data.iloc[0]
od = data.iloc[16]

plt.figure()
plt.subplot(2,1,1)
plt.plot(time,od,'r.')
plt.ylabel('OD')
plt.xlabel('Time [h]')
plt.subplot(2,1,2)
plt.plot(time,od,'b.')
plt.ylabel('GR [Hr$^{-1}$]')
plt.xlabel('Time [h]')
