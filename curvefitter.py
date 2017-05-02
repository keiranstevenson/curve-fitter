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
import re
import sys


def curvefitter(filename,header= None, predefinedinput= None, skiprows= 0, labelcols= 3, replicols= 3, 
                      waterwells= False, replicates= False, repignore= None, normalise= 0.05, growthmin= 0.05, alignvalue= 0.1,
                      fitparams= {0:[-5,8], 1:[-6,-1], 2:[-5,2]}, noruns= 5, nosamples= 20,
                      makeplots = True, showplots= True):
    '''
    filename: filename or folder location (only if replicates==1)
    predefinedinput: 'BMG' sets skiprows, labelcols, replicols and repignore based on standard format for BMG platereader files
    skiprows; lines of input file to skip before retrieving data. First row assumed as time with data following immediately below
    labelcols; first n columns are labels/text and are used to populate the output
    replicols; column containing the strings used to match replicates
    waterwells; ignores wells on the outside of 96 well plate
    replicates; indicates presense of replicates to be used for data sorting automatically runs normalise and align on replicates to ensure most accurate GR
    repignore; regex string that defines replicates to be ignored ie 'Sample *' for BMG files
    normalise; value that data is normalised to at the start DO NOT USE 0 or log function will fail
    growthmin; minimum value required for growth to be counted and fitted, fitting purely flat functions consumes time for fitting and produces unreliable results
    alignvalue; aligns replicates so that this value is reached at the same time for all reps
    fitparams; fitparameters used by the deODoriser
    noruns; number of fitting attempts made
    nosamples; number of samples used to calculate error
    makeplots; determines if program makes plots and saves to output folder
    showplots; displays plots during processing
    '''
    
    if (predefinedinput == 'BMG'):
        skiprows= 6
        labelcols = 3
        replicols = 3
        repignore = 'Sample X*'

        
    replicols = replicols-1 # convers to index from number
    # Process files before inputting
    filename = os.path.realpath(filename)
    if replicates & os.path.isdir(filename)==True:
        infile = multifilerepimport(filename, header, skiprows, labelcols, waterwells)
        filepath = os.path.join(filename, os.path.split(filename)[-1] + ' outputdata')
        filename = os.path.split(filename)[-1]
        infile = cleannonreps(infile, replicols, repignore)
        
        reps = infile.iloc[1:,replicols]
        uniquereps = np.unique(reps)
        
        # Provide info about datashape
        dataheight = uniquereps.shape[0]
        datalength = infile.shape[1]
        firstline = infile.iloc[0]
        
        print('++++++++++ Found '+str(dataheight)+' replicates ++++++++++')
        for x in uniquereps: print(x)
        sys.stdout.flush()
    elif os.path.isdir(filename):
        files = glob(os.path.join(filename, '*.[cC][sS][vV]')) + glob(os.path.join(filename, '*.[xX][lL][sS][xX]'))
        print('++++++++++Detected folder. Processing ' + str(len(files)) + ' files++++++++++')
        for i in range(0,len(files)):
            filename = files[i]
            print('++++++++++ Processing file ++++++++++')
            print(filename)
            curvefitter(filename, header, predefinedinput, skiprows, labelcols, replicols, 
                      waterwells, replicates, repignore, normalise, growthmin, alignvalue,
                      fitparams, noruns, nosamples, makeplots, showplots)
        return
    else:
        try:
            infile = pd.read_csv(filename, header=header, skiprows=skiprows)
        except pd.parser.CParserError:
            infile = pd.read_excel(filename, header=header, skiprows=skiprows)
        
        infile = infile.iloc[:,0:-1]
        filepath = os.path.split(filename)[0]
        filename = os.path.split(filename)[1]
        filename = filename.split('.')[-2]
        filepath = os.path.join(filepath, filename + ' outputdata')
        if waterwells:
            infile = removewaterwells(infile,labelcols)    
        
        # Gather info about raw numerical data
        dataheight = infile.shape[0]-1 #ignore time row
        datalength = infile.shape[1]
        firstline = infile.iloc[0]

    infile = normalisetraces(infile, normalise, labelcols)
    # Checks and makes output directories
    if not os.path.isdir(filepath):
        os.makedirs(filepath)       
    if makeplots:
        if not os.path.isdir(os.path.join(filepath, 'plots', )):
            os.makedirs(os.path.join(filepath, 'plots', )) 
    
    # Separate time variable
    time = infile.iloc[0]
    time = time.iloc[labelcols:]
    time = np.float64(time)
    
    # sanity check time
    timecheck = np.gradient(time)
    if any(timecheck<0):
        print('Time does not always increase along its length!!!')
        print('Estimating missing values')
        timegradient = np.diff(time)
        meanstep = np.mean(timegradient[0:10])
        
        for i in range(len(time)-1):
            if abs(time[i]-time[i+1]) > meanstep*1.5:
                time[i+1] = time[i] + meanstep
    
    try:
        for i in range(1,dataheight+1):
            if replicates:
                location = '++++++++++ Processing replicate set ' + str(i) + ' of ' + str(dataheight) + ' ++++++++++'
                print(location)
            else:
                location = '++++++++++ Processing row ' + str(i) + ' of ' + str(dataheight) + ' ++++++++++'
                print(location) 
            if replicates:
                repset = uniquereps[i-1]
                repselection = repset == infile.iloc[:,replicols]
                od = infile.loc[repselection]
                labels = od.iloc[0,0:labelcols]
                labels = labels.copy()
                od = (od.iloc[:,labelcols:]).copy()
                # Converts columns to float format for fitderiv    
                odfloat = np.array(od,dtype='float64')  
                odfloat = alignreplicates(odfloat, normalise, alignvalue)
                repinfo = 'Fitting ' + str(odfloat.shape[0]) + ' replicates'
                print(repinfo)
                
            else:
                od = infile.iloc[i]
                od = od[labelcols+1:]
                
                labels = infile.iloc[i,0:labelcols]
                labels = labels.copy()
                # Converts columns to float format for fitderiv    
                odfloat = np.array(od,dtype='float64')  
            sys.stdout.flush()
            # Removes nans from data that matlab put in
            datalength = odfloat.shape[-1]
            t = time[:datalength]
            
            # Check for growth
            #diff = np.amax(np.ndarray.flatten(odfloat))-np.amin(np.ndarray.flatten(odfloat))
            odfloat2 = np.ndarray.flatten(odfloat)
            growthminnorm = np.amin(odfloat2) + growthmin
                
            for ii in range(0,odfloat2.shape[0]-3):
                if (odfloat2[ii]> growthminnorm) & (odfloat2[ii+1]> growthminnorm) & (odfloat2[ii+2]> growthminnorm):
                    diff = True
                    break
                else:
                    diff = False
            
            if diff: 
                # Runs fitderiv only if growth is over 0.05
                for attemptno in range(5):
                    try:
                        fitty = fitderiv(t,np.transpose(odfloat),bd= fitparams,exitearly= False,nosamples= nosamples, noruns= noruns) # Peters program
                        break
                    except KeyboardInterrupt:
                        raise(KeyboardInterrupt)
                    except:
                        if attemptno == 4:
                            raise
                Gr = fitty.ds['max df']
                Err = fitty.ds['max df var']
                Lag = fitty.ds['lag time']
                TimeofmaxGr = fitty.ds['time of max df']
                note = ["Normal"]
                
                fitcurve = fitty.f
                fitcurveerr = fitty.fvar
                fitdercurve = fitty.df
                fitdercurveerr = fitty.dfvar
                functime = fitty.t
                
                if makeplots:
                    plt.figure()
                    plt.subplot(2,1,1)
                    plt.plot(functime,fitty.d,'r.')
                    plt.plot(functime,fitcurve,'b')
                    plt.ylabel('log OD')
                    plt.xlabel('Time [h]')
                    plt.subplot(2,1,2)
                    plt.plot(functime,fitdercurve,'b')
                    plt.fill_between(functime, fitdercurve-np.sqrt(fitdercurveerr),fitdercurve+np.sqrt(fitdercurveerr), facecolor= 'blue', alpha=0.2)
                    plt.ylabel('GR [Hr$^{-1}$]')
                    plt.xlabel('Time [h]')
    
                    if replicates:
                        picname = labels.iloc[-1]
                    elif predefinedinput == 'BMG':
                        picname = labels.iloc[0] + str(labels.iloc[1])
                    else:
                        picname = str(labels.iloc[i,0])
                    picname = picname + '.PNG'
                    picname = os.path.join(filepath, 'plots',picname)
                    plt.savefig(picname)
                    if showplots:
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
                print("No growth found! Less than " + str(growthmin) + ' change in OD detected')
            
            # Sticks into the output variable (allows individual debugging)
            
            if i == 1:
                growthrates         = (pd.concat([labels,pd.DataFrame([Gr,Err,Lag,TimeofmaxGr])],ignore_index = True)).transpose()
                growthcurves        = (pd.DataFrame(firstline)).transpose()
                growthcurveserr     = (pd.DataFrame(firstline)).transpose()
                growthcurvesder     = (pd.DataFrame(firstline)).transpose()
                growthcurvesdererr  = (pd.DataFrame(firstline)).transpose()       
            else:    
                growthratesin          = (pd.concat([labels,pd.DataFrame([Gr,Err,Lag,TimeofmaxGr])],ignore_index = True)).transpose()
                growthrates         = pd.concat([growthrates,growthratesin],ignore_index = True)
    
            growthcurvesin          = (pd.concat([labels,pd.DataFrame(fitcurve)],ignore_index = True)).transpose()
            growthcurveserrin       = (pd.concat([labels,pd.DataFrame(fitcurveerr)],ignore_index = True)).transpose()
            growthcurvesderin       = (pd.concat([labels,pd.DataFrame(fitdercurve)],ignore_index = True)).transpose()
            growthcurvesdererrin    = (pd.concat([labels,pd.DataFrame(fitdercurveerr)],ignore_index = True)).transpose()
            
            growthcurves        = pd.concat([growthcurves,growthcurvesin],ignore_index = True)
            growthcurveserr     = pd.concat([growthcurveserr,growthcurveserrin],ignore_index = True)
            growthcurvesder     = pd.concat([growthcurvesder,growthcurvesderin],ignore_index = True)
            growthcurvesdererr  = pd.concat([growthcurvesdererr,growthcurvesdererrin],ignore_index = True)   
    except:
        if i>1:
            Outputname = os.path.join(filepath, filename + ' Analysed.xlsx')
            writer = pd.ExcelWriter(Outputname, engine='xlsxwriter')
            growthrates.to_excel(writer, sheet_name='Stats')
            growthcurves.to_excel(writer, sheet_name='fit')
            growthcurveserr.to_excel(writer, sheet_name='fit err')
            growthcurvesder.to_excel(writer, sheet_name='Derivative')
            growthcurvesdererr.to_excel(writer, sheet_name='Derivative err')
            writer.save()
        raise
            
    varnames = ['Row','Column','Note','GR','GR Err', 'Lag', 'Time of max GR']
    growthrates.columns = varnames        
        
    Outputname = os.path.join(filepath, filename + ' Analysed.xlsx')
    writer = pd.ExcelWriter(Outputname, engine='xlsxwriter')
    growthrates.to_excel(writer, sheet_name='Stats')
    growthcurves.to_excel(writer, sheet_name='fit')
    growthcurveserr.to_excel(writer, sheet_name='fit err')
    growthcurvesder.to_excel(writer, sheet_name='Derivative')
    growthcurvesdererr.to_excel(writer, sheet_name='Derivative err')
    writer.save()
    return
    
    
def removewaterwells(indata,labelcols,deletewells=1):
    # deletewells removes from tables else forces points to zero
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
    
    
def cleannonreps(indata, replicol, repignore):
    if repignore == None:
        return indata
    if isinstance(repignore,str):
        cols = indata.iloc[:,replicol]
        a = np.array(cols)
        for i in range(cols.size): 
            try:
                if re.match(repignore,cols[i]):
                    a[i]= False
                else:
                    a[i]= True
            except:
                print('WARNING: unparsable replicate label, skipping')
                a[i]= False
                pass
        indata = (indata.loc[a]).copy()
        
      
    return indata
            
            
def multifilerepimport(filedirectory, header, skiprows, labelcols, waterwells):
    files = glob(os.path.join(filedirectory, '*.[cC][sS][vV]'))
    print('++++++++++ Found '+ str(len(files)) +' files ++++++++++')
    # First need to determine minimum file length to use as first input
    lengths = np.array(np.zeros([np.shape(files)[0],1]),dtype='int')
    for i in range(0,(len(files))):
        testfile = pd.read_csv(files[i], header= header, skiprows= skiprows)
        lengths[i] = testfile.shape[1]
    
    minlength = np.amin(lengths)
    
    # Assembles all files into a single large dataset
    for i in range(0,(len(files))):
        if i == 0:
            stackfile = pd.read_csv(files[i], header= header, skiprows= skiprows)
            if stackfile.shape[1] > minlength:
                stackfile = stackfile.iloc[:,0:minlength-1]
                
        else:
            newfile = pd.read_csv(files[i], header= header, skiprows= skiprows+1)
            if newfile.shape[1] > minlength:
                newfile = newfile.iloc[:,0:minlength-1]
            stack = [stackfile,newfile]
            stackfile = pd.concat(stack,ignore_index = True)
    if waterwells:
        removewaterwells(stackfile,labelcols)
    return stackfile
    
    
def normalisetraces(indata, normvalue= 0.05, labelcols= 3, timerow= 1):

    data = indata.iloc[timerow:,labelcols:]
    # Normalises line by line on points 5:15
    for i in range(0,data.shape[0]):
        datasnip = data.iloc[i,:]
        datasnip = np.array(datasnip)
        zeroingvalue= np.mean(datasnip[4:14])
        zeroingvalue= normvalue - zeroingvalue
        newnumline = datasnip + zeroingvalue
        
        if i == 0:
            stack= newnumline
        else:
            stack= np.vstack((stack,newnumline))
    indata.iloc[1:,labelcols:] = stack 
    return indata    


def alignreplicates(dataset, normvalue= 0.05, alignvalue= 0.1):
    alignpoint= normvalue+alignvalue
    startindexes = np.int_(dataset[:,0])
    
    # Finds where data > alignpoint for 3 consecutive points
    for i in range(0,dataset.shape[0]):
        for ii in range(0,dataset.shape[1]):
            if (dataset[i,ii]> alignpoint) & (dataset[i,ii+1]> alignpoint) & (dataset[i,ii+2]> alignpoint):
                x = ii
                break
        startindexes[i] = np.int_(x)
        
    # Aligns data into array by dropping cells in everything but lowest lag set
    minindex = startindexes.min()
    maxindex = startindexes.max()
    maxrowlength = dataset.shape[1]-maxindex-1
    for i in range(0,dataset.shape[0]):
        rowindex = startindexes[i]
        newdata = dataset[i,rowindex-minindex:rowindex+maxrowlength-1]
        if i == 0:
            stack= newdata
        else:
            stack= np.vstack((stack,newdata))
    return stack

