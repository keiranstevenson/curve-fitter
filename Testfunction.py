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


filename = 'Datafiles/Alextest/Alextest.csv'
skiprows= 6
labelcols= 3
header= None

infile = pd.read_csv(filename, header=header, skiprows=skiprows)
normvalue= 0.05
alignvalue= 0.1
