# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 16:09:32 2016

@author: Keiran Stevenson
"""

from curvefitter import CurveFitter as cfit

filename = 'Datafiles/Shorttest.csv'

fitter = cfit(filename, skip_rows=6, label_columns=3, replicate_column=3)
fitter.fit_data()

## This test takes a while due to high number of replicates
# filename = 'Datafiles/multireptest/'

# cfit(filename, predefined_input= 'BMG', replicates_exist=True, no_runs=1, no_samples=4)

# filename = 'Datafiles/Shorttest.CSV'
# cfit(filename, label_columns=3, replicate_column=3, skip_rows=6, no_runs=1, no_samples=4)

# filename = 'Datafiles/Singlereplicate ex.csv'
# cfit(filename, predefined_input='BMG', no_runs=1, no_samples=4, show_plots=False,
#      replicates_exist=True)

# filename = 'Datafiles/Difformat/Tecanstyle.xlsx'
# cfit(filename, sho_wplots=False, skip_rows=63, label_columns=1, growth_min=0.101, fiting_parameters={0: [-3, 8], 1: [-6, 2], 2: [-5, 2]})

#filename = 'Datafiles/'
# cfit(filename, predefined_input='BMG', show_plots=False)

print('YAY! No errors!')


