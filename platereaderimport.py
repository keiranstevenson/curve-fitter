# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 16:09:32 2016

@author: Keiran Stevenson
"""

import numpy as np
import pandas as pd
from fitderiv import fitderiv
import matplotlib.pyplot as plt
import os
from glob import glob


def platereaderimport(filename,header= None, manufacturer= 'BMG', skiprows= 0, labelcols = 1, 
                      waterwells = 1, replicates = 0, growthmin = 0.05,
                      fitparams= {0: [-5,8],1: [-6,-1], 2: [-5,2]}, noruns= 5, nosamples= 20,
                      makeplots = 1):
    # filename or folder location, if latter then all files imported and processed
    # waterwells; ignores wells on the outside
    # replicates; indicates presense of replicates, denoted by same sample no in 3rd column BMG files
    # Manufacturer; Software that the data was exported for, sets skiprows and labelcol based on that
    # skiprows; lines of input file to skip before retrieving data. First row assumed as time with data following immediately below
    # labelcol; assumes first n rows are labels/text and are used to populate the output
    # fitparams; fitparameters used by the deODoriser
    # noruns; number of fitting attempts made
    # nosamples; number of samples used to calculate error
    
    if (manufacturer == 'BMG'):
        skiprows= 6
        labelcols = 3
    elif (manufacturer == 'Tecan'):
        skiprows= 63
        labelcols = 1
        
    if makeplots == 1: # makes directory for plots to be put in
        if not os.path.isdir(filename.split('/')[-2] + '/plots/'):
            os.makedirs(filename.split('/')[-2] + '/plots/') 
        
    # Process files before inputting
    if replicates == 1 & os.path.isdir(filename)==True:
        infile = multifilerepimport(filename, header, skiprows, labelcols, waterwells)
    elif os.path.isdir(filename):
        replicates == 1
    else:
        infile = pd.read_csv(filename, header=header, skiprows=skiprows)
    
        if waterwells == 1:
            infile = removewaterwells(infile,labelcols,0)    
        
    
    # Prepare output variables with names
    growthrates = infile.iloc[1:,0:7]
    growthcurves = infile.copy()
    growthcurveserr = infile.copy()
    growthcurvesder = infile.copy()
    growthcurvesdererr = infile.copy()
    varnames = ['Row','Column','Note','GR','GR Err', 'Lag', 'Time of max GR']
    growthrates.columns = varnames
    # Separate time variable
    time = infile.iloc[0]
    time = time[labelcols:]
    time = np.float64(time)
    # Gather info about raw numerical data
    dataheight = infile.shape[0]-1
    datalength = infile.shape[1]    
    
    for i in range(1,dataheight):
        location = '++++++++++ Processing row ' + str(i) + ' of ' + str(dataheight) + ' ++++++++++'
        print(location)     
        
        od = infile.iloc[i]
        
        od = od[labelcols:]
        od = od
        
        # Converts columns to float format for fitderiv    
        odfloat = np.float64(od)  
        
        # Removes nans from data that matlab put in
        odfloat =  odfloat[np.isfinite(odfloat)] #  [ num for num in odfloat if (False ==(np.isnan(num)).any())] 
        datalength = len(odfloat)
        t = time[:datalength]
        
        
        # Check for growth
        diff = np.amax(od)-np.amin(od)        
        
        if diff > growthmin: 
            # Runs fitderiv only if growth is over 0.05
            fitty = fitderiv(t,odfloat,bd= fitparams,exitearly= False,nosamples= nosamples, noruns= noruns) # Peters program
            Gr = fitty.ds['max df']
            Err = fitty.ds['max df var']
            Lag = fitty.ds['lag time']
            TimeofmaxGr = fitty.ds['time of max df']
            note = ["Normal"]
            
            fitcurve = fitty.f
            fitcurveerr = fitty.fvar
            fitdercurve = fitty.df
            fitdercurveerr = fitty.dfvar
            
            if makeplots == 1:
                plt.figure()
                plt.subplot(2,1,1)
                plt.plot(t,fitty.d,'r.')
                plt.plot(t,fitcurve,'b')
                plt.ylabel('log OD')
                plt.xlabel('Time [h]')
                plt.subplot(2,1,2)
                plt.plot(t,fitdercurve,'b')
                plt.fill_between(t, fitdercurve-np.sqrt(fitdercurveerr),fitdercurve+np.sqrt(fitdercurveerr), facecolor= 'blue', alpha=0.2)
                plt.ylabel('GR [Hr$^{-1}$]')
                plt.xlabel('Time [h]')

                if replicates == 1:
                    picname = str(infile.iloc[i,2])
                        
                elif manufacturer == 'BMG':
                    picname = str(infile.iloc[i,0]) + str(int(infile.iloc[i,1]))
                    
                else:
                    picname = str(infile.iloc[i,0])
                
                picname = filename.split('/')[0] + '/plots/' + picname + '.PNG'
                plt.savefig(picname)
                plt.show()
                        
                
        else:
            # Returns no growth if none detected
            Gr = 0
            Err = 0
            Lag = 0
            note = ["No Growth"]
            TimeofmaxGr = 0
            
            fitcurve = np.zeros(datalength)
            fitcurveerr = np.zeros(datalength)
            fitdercurve = np.zeros(datalength)
            fitdercurveerr = np.zeros(datalength)
            if diff == 0:
                print("Empty well, skipping analysis.")
            else:
                print("No growth found! Less than " + str(growthmin) + ' change in OD detected')
        
        # Sticks into the output variable (allows individual debugging)
        growthrates.ix[i,2] = note    
        growthrates.ix[i,3] = Gr
        growthrates.ix[i,4] = Err
        growthrates.ix[i,5] = Lag
        growthrates.ix[i,6] = TimeofmaxGr
        
        growthcurves.ix[i,labelcols:datalength+labelcols-1] = fitcurve
        growthcurveserr.ix[i,labelcols:datalength+labelcols-1] = fitcurveerr
        growthcurvesder.ix[i,labelcols:datalength+labelcols-1] = fitdercurve
        growthcurvesdererr.ix[i,labelcols:datalength+labelcols-1] = fitdercurveerr
        
    Outputname = filename.split('.')[-2] + '_GRoutput.xlsx'
    writer = pd.ExcelWriter(Outputname, engine='xlsxwriter')
    growthrates.to_excel(writer, sheet_name='Stats')
    growthcurves.to_excel(writer, sheet_name='fit')
    growthcurveserr.to_excel(writer, sheet_name='fit err')
    growthcurvesder.to_excel(writer, sheet_name='Derivative')
    growthcurvesdererr.to_excel(writer, sheet_name='Derivative err')
    writer.save()
    
    
    
def removewaterwells(indata,labelcols,deletewells):
       
    cols = indata.iloc[:,0]
    colindex = (cols == 'A') | (cols == 'H')
    cols = indata.iloc[:,1]
    colindex = (cols == 1) | (cols == 12) | colindex
    if deletewells == 1:
        colindex = colindex == False
        data = indata.loc[colindex]
        data = data.copy()
    else:        
        data = indata.copy()
        for i in range(0,len(colindex)):
            if 1 == colindex[i]:
                data.iloc[[i],3:] = np.zeros(data.shape[1]-3)
    return data
            
def multifilerepimport(filedirectory, header, skiprows, labelcols, waterwells):
    files = glob(filedirectory + '/*.csv')
    #files = [file.split('/')[-1] for file in files]
    for i in range(0,(len(files)-1)):
        if i == 0:
            stackfile = pd.read_csv(files[i],header= header,skiprows= skiprows)
        else:
            newfile = pd.read_csv(files[i],header= header,skiprows= skiprows+1)
            stack = [stackfile,newfile]
            stackfile = pd.concat(stack)
    if waterwells == 1:
        removewaterwells(stackfile,labelcols,1)
   
    return stackfile
    
