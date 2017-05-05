# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 16:09:32 2016

@author: Keiran Stevenson
"""
 
from curvefitter import curvefitter as cfit

## This test takes a while due to high number of replicates
#filename = 'Datafiles/Multireptest/'
#cfit(filename, predefinedinput= 'BMG', replicates=True, waterwells=True)

filename = 'Datafiles/Shorttest.CSV'
cfit(filename, labelcols = 3, replicols = 3, skiprows=6, waterwells=True, normalise = 0.1, noruns=1, nosamples=4)

filename = 'Datafiles/Singlereplicate ex.CSV'
cfit(filename, predefinedinput= 'BMG', waterwells=True, normalise = 0.1, noruns=1, nosamples=4, showplots= False)

filename = 'Datafiles/Difformat/Tecanstyle.xlsx'
cfit(filename, showplots=False, skiprows=63, labelcols=1, growthmin = 0.101, fitparams={0:[-5,8], 1:[-6,2], 2:[-5,2]})

filename = 'Datafiles/'
cfit(filename, predefinedinput= 'BMG', waterwells=True)

print('YAY')