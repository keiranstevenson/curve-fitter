from curvefitter import CurveFitter as cfit
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os


# logistic
def funclog(x, a, u, l):
    y = a/(1+np.exp(l-u*x))
    return y

# bio logistic
def funcbiolog(x, a, u, l):
    sub = (((4*u)/a)*(l-x))+2
    y = a/(1+np.exp(sub))
    # y = a/(1+exp((())))
    return y

# Gompertz
def funcgompertz(x,)

# Bio Gompertz

#exponential
def funcexp(x, a, u, c):
    l = 20 # flattens first 20 time units
    y = a*np.exp(u*x)+c
    y[:l] = y[l-1]
    y[l:l*5] = 0
    return y

#linear
def funclinear(x, a, u, c):
    y = x.copy
    y = u*x+c
    y[:20] = y[20+1]
    return y


if os.path.exists('testfunctions') is False:
    os.mkdir('testfunctions')

time = np.arange(0,21,0.1)
a = 1
u = 2
c = 10
y = funcexp(time, a, u, c)
# plt.plot(time,y)
# plt.show()


frame = pd.DataFrame([time,y])
frame.index = ['time','exponential u{}'.format(u)]
print(frame)
outputname = 'testfunctions\\exp.xlsx'
writer = pd.ExcelWriter(outputname, engine='xlsxwriter')
frame.to_excel(writer)
writer.close()
fitter = cfit(outputname, growth_minimum=0.001, skip_rows=1, label_columns=1, replicate_column=1, no_runs=10,no_samples=50,
              show_plots=True,replicates_exist=False, logdata=True)
fitter.fiting_parameters = {0: [-3, 8], 1: [-6, 8], 2: [-5, 2]}
fitter.logdata=True
fitter.fit_data()

a = 1
y = funcbiolog(time, a, u, c)
# plt.plot(time,y)
# plt.show()

frame = pd.DataFrame([time,y])
frame.index = ['time','biologistic u{}'.format(u)]
print(frame)
outputname = 'testfunctions\\biologistic.xlsx'
writer = pd.ExcelWriter(outputname, engine='xlsxwriter')
frame.to_excel(writer)
writer.close()
fitter.file_import(outputname)
fitter.fit_data()

y = funclog(time, a, u, c + 10)
# plt.plot(time,y)
# plt.show()

frame = pd.DataFrame([time,y])
frame.index = ['time','logistic u{}'.format(u)]
print(frame)
outputname = 'testfunctions\\logistic.xlsx'
writer = pd.ExcelWriter(outputname, engine='xlsxwriter')
frame.to_excel(writer)
writer.close()
fitter.file_import(outputname)
fitter.fit_data()

y = funclinear(time, a, u, c + 10)

frame = pd.DataFrame([time,y])
frame.index = ['time','linear u{}'.format(u)]
print(frame)
outputname = 'testfunctions\\linear.xlsx'
writer = pd.ExcelWriter(outputname, engine='xlsxwriter')
frame.to_excel(writer)
writer.close()
fitter.file_import(outputname)
fitter.fit_data()
