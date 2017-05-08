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
import platform

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

    if (predefinedinput == 'BMG'): # For BMG output from TP lab reader
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

        # Provide info about datashape for interation to use
        dataheight = uniquereps.shape[0]
        datalength = infile.shape[1]
        firstline = infile.iloc[0]

        print('++++++++++ Found '+str(dataheight)+' replicates ++++++++++')
        for x in uniquereps: print(x)
        sys.stdout.flush()

    elif replicates & os.path.isfile(filename):
        try:
            infile = pd.read_csv(filename, header=header, skiprows=skiprows)
        except pd.parser.CParserError:
            infile = pd.read_excel(filename, header=header, skiprows=skiprows)

        infile = cleannonreps(infile, replicols, repignore)
        reps = infile.iloc[1:,replicols]
        uniquereps = np.unique(reps)

        dataheight = uniquereps.shape[0]
        datalength = infile.shape[1]
        firstline = infile.iloc[0]

        print('++++++++++ Found '+str(dataheight)+' replicates ++++++++++')
        for x in uniquereps: print(x)
        sys.stdout.flush()

        filepath = os.path.split(filename)[0]
        filename = os.path.split(filename)[1]
        filename = filename.split('.')[-2]
        filepath = os.path.join(filepath, filename + ' outputdata')

    elif os.path.isdir(filename):
        files = glob(os.path.join(filename, '*.[cC][sS][vV]')) + glob(os.path.join(filename, '*.[xX][lL][sS][xX]'))
        print('++++++++++Detected folder. Processing ' + str(len(files)) + ' files++++++++++')
        for i in range(0,len(files)):
            filename = files[i]
            print('++++++++++ Processing file ++++++++++')
            print(filename)
            # Yay recursion
            curvefitter(filename, header, predefinedinput, skiprows, labelcols, replicols,
                      waterwells, replicates, repignore, normalise, growthmin, alignvalue,
                      fitparams, noruns, nosamples, makeplots, showplots)
        return

    elif os.path.isfile(filename):
        try:
            infile = pd.read_csv(filename, header=header, skiprows=skiprows)
        except pd.parser.CParserError:
            infile = pd.read_excel(filename, header=header, skiprows=skiprows)

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
    else:
        raise ImportError('File or directory not found')

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
        print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! \nTime does not always increase along its length \n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
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
                repset = uniquereps[i-1]
                repselection = repset == infile.iloc[:,replicols]
                od = infile.loc[repselection]
                # Defines names for results
                labels = pd.DataFrame([repset])
                repinfo1 = 'Fitting ' + repset
                print(repinfo1)

                od = (od.iloc[:,labelcols:]).copy()
                # Converts columns to float format for fitderiv
                odfloat = np.array(od,dtype='float64')
                odfloat = normalisetraces(odfloat, normalise)
                odfloat = alignreplicates(odfloat, normalise, alignvalue)
                noofreps = odfloat.shape[0]
                repinfo2 = 'Found ' + str(noofreps) + ' replicates'
                print(repinfo2)

                nantest = np.isnan(odfloat)
                for ii in range(0,nantest.shape[1]):
                    if any(nantest[:,ii]):
                        x = ii-1
                        break
                odfloat = odfloat[:,:ii]
                t = time[:ii]
            else:
                location = '++++++++++ Processing row ' + str(i) + ' of ' + str(dataheight) + ' ++++++++++'
                print(location)
                od = infile.iloc[i]
                od = od[labelcols+1:]
                noofreps = 1
                # Defines names for results
                labels = infile.iloc[i,0:labelcols]
                labels = labels.copy()

                # Converts columns to float format for fitderiv
                odfloat = np.array(od,dtype='float64')
                odfloat = normalisetraces(odfloat, normalise)
                nantest = np.isnan(odfloat)

                for ii in range(0,len(nantest)):
                    if nantest[ii]:
                        x = ii-1
                        break
                odfloat = odfloat[:ii]
                t = time[:ii]
            sys.stdout.flush() # Forces prints to display immediately



            # Check for growth
            diff = checkforgrowth(odfloat,growthmin)
            if (type(diff) == np.ndarray) | (type(diff) == list):
                invdiff = [not x for x in diff]
                if np.any(invdiff) & replicates:
                    print('Error, no growth detected in replicate, dropping from analysis')
                    odfloat = odfloat[np.array(diff,bool),:]
                diff=True

            if diff:
                # Runs fitderiv only if growth is over growthmin
                for attemptno in range(5):
                    try:
                        fitty = fitderiv(t,np.transpose(odfloat),bd= fitparams,exitearly= False,nosamples= nosamples, noruns= noruns) # Peters program
                        break
                    except KeyboardInterrupt:
                        raise KeyboardInterrupt('User aborted run')
                    except MemoryError:
                        version = platform.architecture()[0]
                        if version == '32bit':
                            raise MemoryError('Out of Memory while fitting. Try installing 64-bit python or using fewer replicates')
                        elif version == '64bit':
                            raise MemoryError('Out of memory while fitting. Try using fewer replicates')
                        else:
                            raise MemoryError('Out of memory while fitting. Unable to determine version, try making more memory available or using fewer replicates')

                    except:
                        print('Fitting failure, retrying')
                        if attemptno == 4:
                            raise
                # Pulls stats and data from Peters routine
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
                        picname = repset
                    elif predefinedinput == 'BMG':
                        picname = labels.iloc[0] + str(labels.iloc[1])
                    elif len(labels) == 1:
                         picname = labels.astype('str')
                         picname = picname.str.cat()
                    else:
                         picname = labels.astype('str')
                         picname = picname.str.cat()
                    picname = picname + '.PNG'
                    picname = os.path.join(filepath, 'plots',picname)
                    plt.savefig(picname)
                    if showplots:
                        plt.ion()
                        plt.show()
                    else:
                        plt.close()

            else:
                # Returns no growth if none detected
                Gr = 0
                Err = 0
                Lag = 0
                note = ["No Growth"]
                TimeofmaxGr = 0
                noofreps = 0

                fitcurve = np.zeros(datalength)
                fitcurveerr = np.zeros(datalength)
                fitdercurve = np.zeros(datalength)
                fitdercurveerr = np.zeros(datalength)
                print("No growth found! Less than " + str(growthmin) + ' change in OD detected')

            # Sticks into the output variable (allows individual debugging)
            if i == 1:
                growthrates         = (pd.concat([labels,pd.DataFrame([Gr,Err,Lag,TimeofmaxGr,noofreps])],ignore_index = True)).transpose()
                growthcurves        = (pd.DataFrame(firstline)).transpose()
                growthcurveserr     = (pd.DataFrame(firstline)).transpose()
                growthcurvesder     = (pd.DataFrame(firstline)).transpose()
                growthcurvesdererr  = (pd.DataFrame(firstline)).transpose()
            else:
                growthratesin          = (pd.concat([labels,pd.DataFrame([Gr,Err,Lag,TimeofmaxGr,noofreps])],ignore_index = True)).transpose()
                growthrates         = pd.concat([growthrates,growthratesin],ignore_index = True)

            growthcurvesin          = (pd.concat([labels,pd.DataFrame(fitcurve)],ignore_index = True)).transpose()
            growthcurveserrin       = (pd.concat([labels,pd.DataFrame(fitcurveerr)],ignore_index = True)).transpose()
            growthcurvesderin       = (pd.concat([labels,pd.DataFrame(fitdercurve)],ignore_index = True)).transpose()
            growthcurvesdererrin    = (pd.concat([labels,pd.DataFrame(fitdercurveerr)],ignore_index = True)).transpose()

            growthcurves        = pd.concat([growthcurves,growthcurvesin],ignore_index = True)
            growthcurveserr     = pd.concat([growthcurveserr,growthcurveserrin],ignore_index = True)
            growthcurvesder     = pd.concat([growthcurvesder,growthcurvesderin],ignore_index = True)
            growthcurvesdererr  = pd.concat([growthcurvesdererr,growthcurvesdererrin],ignore_index = True)
    except KeyboardInterrupt:
        raise
    except:
        print('ERROR DURING FIT: DUMPING DATA TO OUTPUT FILE')
        raise
    finally:
        # Always runs data saving in case of error
        if 'growthrates' in locals():
            if replicates:
                varnames = ['Replicate Name','GR','GR Err', 'Lag', 'Time of max GR','no. of replicates']
            else:
                x = list(range(labels.shape[0]))
                for labelindex in range(labels.shape[0]):
                    x[labelindex] = 'Label'
                varnames = x + ['GR','GR Err', 'Lag', 'Time of max GR','no. of replicates']
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
    # Removes wells matching a particular regex
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
    files = glob(os.path.join(filedirectory, '*.[cC][sS][vV]')) + glob(os.path.join(filedirectory, '*.[xX][lL][sS][xX]'))
    print('++++++++++ Found '+ str(len(files)) +' files ++++++++++')

    # First need to determine max time length to use as first input
    for i in range(0,(len(files))):
        if i == 0:
            testfile = pd.read_csv(files[i], header= header, skiprows= skiprows)
            time = testfile.iloc[0,:]
            time = pd.DataFrame(time)
        else:
            testfile = pd.read_csv(files[i], header= header, skiprows= skiprows)
            if testfile.shape[1] > len(time):
                time = testfile.iloc[0,:]
                time = pd.DataFrame(time)

    stackfile = time.transpose()

    for i in range(0,len(files)):
        newfile = pd.read_csv(files[i], header= header, skiprows= skiprows+1)
        stackfile = stackfile.append(newfile, ignore_index=True)

    if waterwells:
        removewaterwells(stackfile,labelcols)
    return stackfile


def normalisetraces(dataset, normvalue= 0.05):
    # Normalises line by line on points 5:15
    try:
        x = dataset.shape[1]
        for i in range(0,dataset.shape[0]):
            zeroingvalue= np.mean(dataset[i,4:14])
            zeroingvalue= normvalue - zeroingvalue
            dataset[i,:] = dataset[i,:] + zeroingvalue

        return dataset
    except IndexError:
        zeroingvalue= np.mean(dataset[4:14])
        zeroingvalue= normvalue - zeroingvalue
        dataset = dataset + zeroingvalue
        return dataset
    except:
        raise


def alignreplicates(dataset, normvalue= 0.05, alignvalue= 0.1):
    diff = checkforgrowth(dataset,alignvalue)
    invdiff = [not x for x in diff]
    if np.any(invdiff):
        print('Error, replicate does not reach alignvalue, dropping from alignment')
        dataset = dataset[np.array(diff,bool),:]

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

def checkforgrowth(data, growthmin):
    try:
        test = np.array(range(0,data.shape[0]))
        for i in range(0,data.shape[0]):
            norm = np.mean(data[i,1:10])
            growthminnorm = growthmin + norm
            for ii in range(0,data.shape[1]-3):
                if (data[i,ii]> growthminnorm) & (data[i,ii+1]> growthminnorm) & (data[i,ii+2]> growthminnorm):
                    growth = True
                    break
                else:
                    growth = False
            test[i] = growth
        return test
    except IndexError:
        norm = np.mean(data[1:10])
        growthminnorm = growthmin + norm
        for i in range(0,data.shape[0]-3):
            if (data[i]> growthminnorm) & (data[i+1]> growthminnorm) & (data[i+2]> growthminnorm):
                growth = True
                break
            else:
                growth = False
        test = growth
        return test
    except:
        raise
