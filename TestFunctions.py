from curvefitter import CurveFitter as cfit
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


# logistic
def fitlog(x,a,u,l):
    y = a/(1+np.exp(l-u*x))
    return y

# bio logistic
def fitbiolog(x,a,u,l):
    sub = (((4*u)/a)*(l-x))+2
    y = a/(1+np.exp(sub))
    # y = a/(1+exp((())))
    return y

#exponential
def fitexp(x,a,u,c):
    l = 20 # flattens first 20 time units
    y = a*np.exp(u*x)+c
    y[:l] = y[l-1]
    y[l:l*5] = 0
    return y



time = np.arange(0,21,0.1)
a = 1
u = 0.5
c = 10
y = fitexp(time,a,u,c)
# plt.plot(time,y)
# plt.show()


frame = pd.DataFrame([time,y])
frame.index = ['time','exponential u{}'.format(u)]
print(frame)
outputname = 'testfunctions\\exp.xlsx'
writer = pd.ExcelWriter(outputname, engine='xlsxwriter')
frame.to_excel(writer)
writer.close()
fitter = cfit(outputname, growth_minimum=0.001, skip_rows=1, label_columns=1, replicate_column=1, no_runs=4,no_samples=10, show_plots=True,replicates_exist=False, logdata=True)
fitter.logdata=True
fitter.fit_data()

a = 1
y = fitbiolog(time,a,u,c)
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

u=u*4
y = fitlog(time,a,u,c+10)
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
