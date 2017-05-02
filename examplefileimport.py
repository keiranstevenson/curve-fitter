# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 16:09:32 2016

@author: s1002426
"""
 
from curvefitter import curvefitter as cfit

filename = 'Datafiles/Alextest/'
cfit(filename, predefinedinput= 'BMG', noruns=1, nosamples=5, replicates=True, waterwells=True)

#filename = 'Datafiles/truncatedtest.CSV'
#cfit(filename, predefinedinput= 'BMG', noruns=1, nosamples=5, waterwells=True)

#filename = 'Datafiles/EtOHtestdata.csv'
#cfit(filename, predefinedinput= 'BMG', noruns=1, nosamples=5, waterwells=True)

filename = 'Datafiles/'
cfit(filename, predefinedinput= 'BMG', noruns=1, nosamples=5, waterwells=True)

print('YAY')