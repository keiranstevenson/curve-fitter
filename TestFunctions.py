from curvefitter import CurveFitter as cfit
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os


# logistic
def funclog(x, a, b, c):
    y = a/(1 + np.exp(c - b * x))
    return y

# bio logistic
def funcbiolog(x, a, u, l):
    sub = (((4*u)/a)*(l-x))+2
    y = a/(1+np.exp(sub))
    # y = a/(1+exp((())))
    return y

# Gompertz
def funcgompertz(x,a,c,b):
    y = a*np.exp(-np.exp(b-(c*x)))
    return y

# Bio Gompertz
def funcbiogompertz(x,a,u,l):
    y = a*np.exp(-np.exp((((u*np.e)/a)*(l-x))+1))
    return y

#exponential
def funcexp(x, a, u, c):
    l = 20 # flattens first 20 time units
    y = a*np.exp(u*x)+c
    y[:l] = y[l-1]
    y[l:l*5] = 0
    y[-20:]=y[-21]
    return y

#linear
def funclinear(x, a, b, c):
    y = x.copy
    y = b * x + c
    y[:20] = y[20+1]
    return y


if os.path.exists('testfunctions') is False:
    os.mkdir('testfunctions')

time = np.arange(0,40,0.1)
a = 1
u = 0.8
c = 10
y = funcexp(time, a, u, c)
plt.plot(time,y)
plt.title('exponential u{}'.format(u))
plt.show()


alloutputs=[]

frame = pd.DataFrame([time,y])
frame.index = ['time','exponential u{}'.format(u)]
outputname = 'testfunctions\\exp.xlsx'
alloutputs.append(outputname)
writer = pd.ExcelWriter(outputname, engine='xlsxwriter')
frame.to_excel(writer)
writer.close()

##
a = 1
y = funcbiolog(time, a, u, c)
plt.plot(time,y)
plt.title('biologistic u{}'.format(u))
plt.show()

frame = pd.DataFrame([time,y])
frame.index = ['time','biologistic u{}'.format(u)]
outputname = 'testfunctions\\biologistic.xlsx'
alloutputs.append(outputname)
writer = pd.ExcelWriter(outputname, engine='xlsxwriter')
frame.to_excel(writer)
writer.close()


# ##
# y = funclog(time, a, u, c + 10)
# plt.plot(time,y)
# plt.title('logistic u{}'.format(u))
# plt.show()
#
# frame = pd.DataFrame([time,y])
# frame.index = ['time','logistic u{}'.format(u)]
# print(frame)
# outputname = 'testfunctions\\logistic.xlsx'
# alloutputs.append(outputname)
# writer = pd.ExcelWriter(outputname, engine='xlsxwriter')
# frame.to_excel(writer)
# writer.close()

##
y = funclinear(time, a, u, c + 10)
plt.plot(time,y)
plt.title('linear u{}'.format(u))
plt.show()
frame = pd.DataFrame([time,y])
frame.index = ['time','linear u{}'.format(u)]
outputname = 'testfunctions\\linear.xlsx'
alloutputs.append(outputname)
writer = pd.ExcelWriter(outputname, engine='xlsxwriter')
frame.to_excel(writer)
writer.close()

# ##
# y = funcgompertz(time, a, u, c+10)
# plt.plot(time,y)
# plt.title('gompertz u{}'.format(u))
# plt.show()
# frame = pd.DataFrame([time,y])
# frame.index = ['time','gompertz u{}'.format(u)]
# print(frame)
# outputname = 'testfunctions\\gompertz.xlsx'
# alloutputs.append(outputname)
# writer = pd.ExcelWriter(outputname, engine='xlsxwriter')
# frame.to_excel(writer)
# writer.close()

##
y = funcbiogompertz(time, a, u, c)
plt.plot(time,y)
plt.title('biogompertz u{}'.format(u))
plt.show()
frame = pd.DataFrame([time,y])
frame.index = ['time','biogompertz u{}'.format(u)]
outputname = 'testfunctions\\biogompertz.xlsx'
alloutputs.append(outputname)
writer = pd.ExcelWriter(outputname, engine='xlsxwriter')
frame.to_excel(writer)
writer.close()






fitter = cfit(outputname, growth_minimum=0.001, skip_rows=1, label_columns=1, replicate_column=1, no_runs=10,no_samples=50,
              show_plots=True,replicates_exist=False, logdata=True)
fitter.fiting_parameters = {0: [-3, 8], 1: [-6, 8], 2: [-5, 2]}
fitter.logdata=False
for file in alloutputs:
    fitter.file_import(file)
    fitter.fit_data()