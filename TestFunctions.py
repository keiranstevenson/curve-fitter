from curvefitter import CurveFitter as cfit
import pandas as pd
import numpy as np
# import matplotlib.pyplot as plt


# logistic
def fitlog(x,a,b,c):
    y = a/(1+np.exp(b-c*x))
    return y

# bio logistic
def fitbiolog(x,a,u,l):
    sub = ((4*u)/a)*(l-x)+2
    y = a/(1+np.exp(sub))
    # y = a/(1+exp((())))
    return y

#exponential
def fitexp(x,a,u,c):
    y = a*np.exp(u*x)+c
    return y



time = np.arange(0,21,0.1)
a = 1
u = 0.5
c = 2
# y = fitbiolog(time,a,u,c)
y = fitexp(time,a,u,c)
# plt.plot(time,y)
# plt.show()


frame = pd.DataFrame([time,y])
print(frame)
outputname = 'temp.xlsx'
writer = pd.ExcelWriter(outputname, engine='xlsxwriter')
frame.to_excel(writer)
writer.close()
fitter = cfit(writer, growth_minimum=0.001, skip_rows=1, label_columns=1, replicate_column=1, no_runs=4,no_samples=10, show_plots=True,replicates_exist=False, logdata=True)
fitter.fit_data()
