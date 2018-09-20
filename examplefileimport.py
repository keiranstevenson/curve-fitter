# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 16:09:32 2016

@author: Keiran Stevenson
"""

from curvefitter import CurveFitter as cfit

filename = 'Datafiles/Shorttest.csv'

fitter = cfit(filename, skip_rows=6, label_columns=3, replicate_column=3, no_runs=4,no_samples=10, show_plots=False)
# fitter.fit_data()

## This test takes a while due to high number of replicates
fitter.replicates_exist = True
fitter.replicate_ignore = 'Sample X*'
fitter.file_import('Datafiles/multireptest/')
# fitter.file_import('Datafiles/Singlereplicate ex.csv')
fitter.fit_data()

fitter.file_import('Datafiles/Shorttest.CSV')
fitter.replicates_exist=False
fitter.replicate_ignore= None
fitter.fit_data()

filename = 'Datafiles/Singlereplicate ex.csv'
fitter.replicates_exist=True
fitter.predefined_input='BMG'
fitter.file_import(filename)
fitter.fit_data()

# filename = 'Datafiles/Difformat/Tecanstyle.xlsx'
# cfit(filename, show_plots=False, skip_rows=63, label_columns=1, growth_minimum=0.101, fiting_parameters={0: [-3, 8], 1: [-6, 2], 2: [-5, 2]})
#
# filename = 'Datafiles/'
# cfit(filename, predefined_input='BMG', show_plots=False)

print('YAY! No errors!')


