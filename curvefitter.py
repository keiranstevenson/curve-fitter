# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 16:09:32 2016

@author: Keiran Stevenson
"""

import os
import platform
import re
import sys
from glob import glob

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from fitderiv import fitderiv


def curvefitter(filename, header=None, predefinedinput=None, skiprows=0, labelcolumns=3, replicatecolumn=3,
                replicatesexist=False, replicateignore=None, 
                normalise=0.05, normby=(4,14), growthmin=0.05, alignvalue=0.1,
                fitparams={0: [-5, 8], 1: [-6, -1], 2: [-5, 2]}, noruns=5, nosamples=20, logdata=True,
                makeplots=True, showplots=False, startnormalise=1):
    '''
    filename: filename or folder location (only if replicatesexist==1)
    predefinedinput: 'BMG' sets skiprows, labelcolumns, replicatecolumn and replicateignore based on standard format for BMG platereader files
    skiprows; lines of input file to skip before retrieving data. First row assumed as time with data following immediately below
    labelcolumns; first n columns are labels/text and are used to populate the output
    replicatecolumn; column containing the strings used to match replicatesexist

    waterwells; ignores wells on the outside of 96 well plate
    replicatesexist; indicates presence of replicatesexist to be used for data sorting automatically runs normalise and align on replicatesexist to ensure most accurate GR
    replicateignore; regex string that defines replicatesexist to be ignored ie 'Sample *' for BMG files
    normalise; value that data is normalised to at the start DO NOT USE 0 or log function will fail
    growthmin; minimum value required for growth to be counted and fitted, fitting purely flat functions consumes time for fitting and produces unreliable results
    alignvalue; aligns replicatesexist so that this value is reached at the same time for all reps
    fitparams; fitparameters used by the deODoriser
    noruns; number of fitting attempts made
    nosamples; number of samples used to calculate error
    makeplots; determines if program makes plots and saves to output folder
    showplots; displays plots during processing
    '''

    if predefinedinput == 'BMG':  # For BMG output from TP lab reader
        skiprows = 6
        labelcolumns = 3
        replicatecolumn = 3
        replicateignore = 'Sample X*'
    elif predefinedinput == 'Tecan':
        skiprows = 63
        labelcolumns = 1

    waterwells = False  # no longer necessary but left
    replicatecolumn = replicatecolumn - 1  # converts to index from number
    startnormalise = startnormalise - 1

    # Process files before inputting
    filename = os.path.realpath(filename)
    if replicatesexist & os.path.isdir(filename) is True:
        infile = multifilerepimport(filename, header, skiprows, labelcolumns)
        filepath = os.path.join(filename, 'curvefitter' + ' outputdata')
        filename = os.path.split(filename)[-1]
        infile = cleannonreps(infile, replicatecolumn, replicateignore)

        reps = infile.iloc[1:, replicatecolumn]
        unique_replicates = np.unique(reps)


        # Provide info about datashape for interation to use
        dataheight = unique_replicates.shape[0]
        datalength = infile.shape[1]
        firstline = infile.iloc[0, labelcolumns - 1:].copy()
        firstline = firstline.reset_index(drop=True)


        print('++++++++++ Found ' + str(dataheight) + ' replicates ++++++++++')
        for x in unique_replicates: print(x)
        sys.stdout.flush()

    elif os.path.isdir(filename):
        files = glob(os.path.join(filename, '*.[cC][sS][vV]')) + glob(os.path.join(filename, '*.[xX][lL][sS][xX]'))
        print('++++++++++Detected folder. Processing ' + str(len(files)) + ' files++++++++++')
        for i in range(0, len(files)):
            filename = files[i]
            print('++++++++++ Processing file ' + filename + '++++++++++')
            print(filename)
            # Yay recursion
            curvefitter(filename=filename, header=header, predefinedinput=predefinedinput, skiprows=skiprows,
                        labelcolumns=labelcolumns, replicatecolumn=replicatecolumn + 1, replicatesexist=replicatesexist,
                        replicateignore=replicateignore, normalise=normalise, growthmin=growthmin,
                        alignvalue=alignvalue, fitparams=fitparams, noruns=noruns, nosamples=nosamples, logdata=logdata,
                        makeplots=makeplots, showplots=showplots)
        return

    elif replicatesexist & os.path.isfile(filename):
        try:
            infile = pd.read_csv(filename, header=header, skiprows=skiprows)
        except pd.parser.CParserError:
            infile = pd.read_excel(filename, header=header, skiprows=skiprows)
        except pd.parser.ParserError:
            infile = pd.read_excel(filename, header=header, skiprows=skiprows)


        infile = cleannonreps(infile, replicatecolumn, replicateignore)
        reps = infile.iloc[1:, replicatecolumn]
        unique_replicates = np.unique(reps)

        dataheight = unique_replicates.shape[0]
        datalength = infile.shape[1]

        firstline = infile.iloc[0, labelcolumns - 1:].copy()
        firstline = firstline.reset_index(drop=True)
        print('++++++++++ Processing file ' + filename + '++++++++++')
        print('++++++++++ Found ' + str(dataheight) + ' replicates ++++++++++')
        for x in unique_replicates: print(x)

        sys.stdout.flush()

        filepath = os.path.split(filename)[0]
        filename = os.path.split(filename)[1]
        filename = filename.split('.')[-2]
        filepath = os.path.join(filepath, 'curvefitter' + ' outputdata')

    elif os.path.isfile(filename) is True:
        try:
            infile = pd.read_csv(filename, header=header, skiprows=skiprows)
        except pd.parser.CParserError:
            infile = pd.read_excel(filename, header=header, skiprows=skiprows)
        except pd.parser.ParserError:
            infile = pd.read_excel(filename, header=header, skiprows=skiprows)

        filepath = os.path.split(filename)[0]
        filename = os.path.split(filename)[1]
        filename = filename.split('.')[-2]
        filepath = os.path.join(filepath, 'curvefitter' + ' outputdata')

        # Gather info about raw numerical data
        dataheight = infile.shape[0] - 1  # ignore time row
        datalength = infile.shape[1]
        firstline = infile.iloc[0]

    else:
        raise ImportError('File or directory not found')

    # Checks and makes output directories
    if not os.path.isdir(filepath):
        os.makedirs(filepath)
    if makeplots:
        if not os.path.isdir(os.path.join(filepath, filename + ' plots', )):
            os.makedirs(os.path.join(filepath, filename + ' plots', ))

    # Separate time variable
    time = infile.iloc[0]
    time = time.iloc[labelcolumns:]
    time = np.float64(time)


    sanitychecks(infile, labelcolumns, normalise, normby, time)

    try:
        for i in range(1, dataheight + 1):
            if replicatesexist:
                location = '++++++++++ Processing replicate set ' + str(i) + ' of ' + str(dataheight) + ' ++++++++++'
                print(location)
                replicate_chosen = unique_replicates[i - 1]
                replicate_index = replicate_chosen == infile.iloc[:, replicatecolumn]
                od = infile.loc[replicate_index]

                # Defines names for results
                labels = pd.DataFrame([replicate_chosen])
                replicate_info1 = 'Fitting ' + replicate_chosen
                print(replicate_info1)
                od = (od.iloc[:, labelcolumns:]).copy()

                # Converts columns to float format for fitderiv
                od_float = np.array(od, dtype='float64')
                od_float = normalise_traces(od_float, normalise, normby)
                od_float = align_replicates(od_float, normalise, alignvalue)
                noofreps = od_float.shape[0]
                replicate_info2 = 'Found ' + str(noofreps) + ' replicates'
                print(replicate_info2)

                # Removes NaNs from analysis
                nantest = np.isnan(od_float)
                for ii in range(0, nantest.shape[1]):
                    if any(nantest[:, ii]):
                        x = ii - 1
                        break
                od_float = od_float[:, :ii]
                t = time[:ii]
                if noofreps == 0:
                    growthfound = False
                else:
                    growthfound = True
            else:
                location = '++++++++++ Processing row {:03} of {:03} ++++++++++'.format(i, dataheight)
                print(location)
                od = infile.iloc[i]
                od = od[labelcolumns + 1:]
                noofreps = 1

                # Defines names for results
                labels = infile.iloc[i, 0:labelcolumns]
                labels = labels.copy()

                # Converts columns to float format for fitderiv
                od_float = np.array(od, dtype='float64')
                od_float = normalise_traces(od_float, normalise, normby)

                # Removes NaNs from analysis
                nantest = np.isnan(od_float)
                for ii in range(0, len(nantest)):
                    if nantest[ii]:
                        x = ii - 1
                        break
                od_float = od_float[:ii]
                t = time[:ii]

                # Check for growth
                growthfound = check_for_growth(od_float, growthmin, normalise)
            sys.stdout.flush()  # Forces prints to display immediately

            if growthfound:
                # Runs fitderiv only if growth is over growthmin
                for attemptno in range(5):
                    try:
                        fitty = fitderiv(t, np.transpose(od_float), bd=fitparams, exitearly=False, nosamples=nosamples,
                                         noruns=noruns, logs=logdata)  # Peters program
                        break
                    except KeyboardInterrupt:
                        raise KeyboardInterrupt('User aborted run')
                    except MemoryError:
                        version = platform.architecture()[0]
                        if version == '32bit':
                            raise MemoryError(
                                'Out of Memory while fitting. Try installing 64-bit python or using fewer replicates')
                        elif version == '64bit':
                            raise MemoryError('Out of memory while fitting. Try using fewer replicates')
                        else:
                            raise MemoryError(
                                'Out of memory while fitting. Unable to determine python version, try making more memory available or using fewer replicates')

                    except:
                        print('Fitting failure, retrying')
                        if attemptno == 4:
                            raise
                # Pulls stats and data from Peters routine
                gr = fitty.ds['max df']
                err = fitty.ds['max df var']
                lag = fitty.ds['lag time']
                timeofmax_gr = fitty.ds['time of max df']
                note = ["Normal"]
                fitcurve = fitty.f
                fitcurveerr = fitty.fvar
                fitdercurve = fitty.df
                fitdercurveerr = fitty.dfvar
                functime = fitty.t

                if makeplots:
                    plt.figure()
                    plt.subplot(2, 1, 1)
                    plt.plot(functime, fitty.d, 'r.')
                    plt.plot(functime, fitcurve, 'b')
                    if logdata:
                        plt.ylabel('log OD')
                    else:
                        plt.ylabel('OD')
                    plt.xlabel('Time [h]')
                    plt.subplot(2, 1, 2)
                    plt.plot(functime, fitdercurve, 'b')
                    plt.fill_between(functime, fitdercurve - np.sqrt(fitdercurveerr),
                                     fitdercurve + np.sqrt(fitdercurveerr), facecolor='blue', alpha=0.2)
                    plt.ylabel('GR [Hr$^{-1}$]')
                    plt.xlabel('Time [h]')

                    if replicatesexist:
                        picname = replicate_chosen
                    elif predefinedinput == 'BMG':
                        picname = labels.iloc[0] + str(labels.iloc[1])
                    elif len(labels) == 1:
                        picname = labels.astype('str')
                        picname = picname.str.cat()
                    else:
                        picname = labels.astype('str')
                        picname = picname.str.cat()
                    picname = picname + '.PNG'
                    picname = os.path.join(filepath, filename + ' plots', picname)
                    plt.savefig(picname)
                    if showplots:
                        plt.ion()
                        plt.show()
                    else:
                        plt.close()

            else:
                # Returns no growth if none detected
                gr = 0
                err = 0
                lag = 0
                note = ["No Growth"]
                timeofmax_gr = 0
                noofreps = 0

                fitcurve = np.zeros(datalength)
                fitcurveerr = np.zeros(datalength)
                fitdercurve = np.zeros(datalength)
                fitdercurveerr = np.zeros(datalength)
                print("No growth found! Less than " + str(growthmin) + ' change in OD detected')

            # Sticks into the output variable (allows individual debugging)
            if i == 1:
                growthrates = (
                    pd.concat([labels, pd.DataFrame([gr, err, lag, timeofmax_gr, noofreps])],
                              ignore_index=True)).transpose()
                growthcurves = (pd.DataFrame(firstline)).transpose()
                growthcurveserr = (pd.DataFrame(firstline)).transpose()
                growthcurvesder = (pd.DataFrame(firstline)).transpose()
                growthcurvesdererr = (pd.DataFrame(firstline)).transpose()
            else:
                growthratesin = (
                    pd.concat([labels, pd.DataFrame([gr, err, lag, timeofmax_gr, noofreps])],
                              ignore_index=True)).transpose()
                growthrates = pd.concat([growthrates, growthratesin], ignore_index=True)

            growthcurvesin = (pd.concat([labels, pd.DataFrame(fitcurve)], ignore_index=True)).transpose()
            growthcurveserrin = (pd.concat([labels, pd.DataFrame(fitcurveerr)], ignore_index=True)).transpose()
            growthcurvesderin = (pd.concat([labels, pd.DataFrame(fitdercurve)], ignore_index=True)).transpose()
            growthcurvesdererrin = (pd.concat([labels, pd.DataFrame(fitdercurveerr)], ignore_index=True)).transpose()

            growthcurves = pd.concat([growthcurves, growthcurvesin], ignore_index=True)
            growthcurveserr = pd.concat([growthcurveserr, growthcurveserrin], ignore_index=True)
            growthcurvesder = pd.concat([growthcurvesder, growthcurvesderin], ignore_index=True)
            growthcurvesdererr = pd.concat([growthcurvesdererr, growthcurvesdererrin], ignore_index=True)
    except KeyboardInterrupt:
        raise
    except:
        print('ERROR DURING FIT: DUMPING DATA TO OUTPUT FILE')
        raise
    finally:
        # Always runs data saving in case of error
        if 'growthrates' in locals():  # Checks that there is data to save
            if replicatesexist:
                varnames = ['Replicate Name', 'GR', 'GR Std Error', 'Lag', 'Time of max GR', 'no. of replicates']
            else:
                x = list(range(labels.shape[0]))
                for labelindex in range(labels.shape[0]):
                    x[labelindex] = 'Label'
                varnames = x + ['GR', 'GR Std Error', 'Lag', 'Time of max GR', 'no. of replicates']
            growthrates.columns = varnames

            if logdata:
                sheetnames = ['Stats', 'Fit (logn)', 'Fit Std Error (logn)', 'Derivative', 'Derivative std Error']
            else:
                sheetnames = ['Stats', 'Fit', 'Fit Std Error', 'Derivative', 'Derivative std Error']
            outputname = os.path.join(filepath, filename + ' Analysed.xlsx')
            writer = pd.ExcelWriter(outputname, engine='xlsxwriter')
            growthrates.to_excel(writer, sheet_name=sheetnames[0])
            growthcurves.to_excel(writer, sheet_name=sheetnames[1])
            growthcurveserr.to_excel(writer, sheet_name=sheetnames[2])
            growthcurvesder.to_excel(writer, sheet_name=sheetnames[3])
            growthcurvesdererr.to_excel(writer, sheet_name=sheetnames[4])
            writer.save()
    return


def sanitychecks(infile, labelcolumns, normalise, normby, time):
    # sanity check time
    timecheck = time[np.logical_not(np.isnan(time))]
    timecheck = np.gradient(timecheck)
    if any(timecheck < 0):
        print(
            '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! \nTime does not always increase along its length \n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
        print('Estimating missing values')
        timegradient = np.diff(time)
        meanstep = np.mean(timegradient)

        for i in range(len(time) - 1):
            if abs(time[i] - time[i + 1]) > meanstep * 1.5:
                time[i + 1] = time[i] + meanstep

    # sanity check data
    for i in range(1, infile.shape[0]):
        data = np.float64(infile.iloc[i, labelcolumns:].copy())
        data = data[np.logical_not(np.isnan(data))]
        if len(data) > len(timecheck):
            raise RuntimeError('Error data is longer than time')
        if len(data) > 0:
            data = normalise_traces(data, normvalue=normalise, normby=normby)
            if any(data <= 0):
                raise ArithmeticError(
                    'Error normalise value gives value <=0. Log function failed, please choose a larger value')


def remove_waterwells(indata, labelcols, deletewells=1):
    # deletewells removes from tables else forces points to zero
    cols = indata.iloc[:, 0]
    column_index = (cols == 'A') | (cols == 'H')
    cols = indata.iloc[:, 1]
    column_index = (cols == 1) | (cols == 12) | column_index
    if deletewells == 1:
        column_index = column_index == False
        data = indata.loc[column_index]
        data = data.copy()
    else:
        data = indata.copy()
        for i in range(0, len(column_index)):
            if 1 == column_index[i]:
                data.iloc[[i], 3:] = np.zeros(data.shape[1] - 3)
    return data


def cleannonreps(indata, replicol, repignore):
    # Removes wells matching a particular regex
    if repignore is None:
        return indata
    if isinstance(repignore, str):
        cols = indata.iloc[:, replicol]
        a = np.array(cols)
        for i in range(cols.size):
            try:
                if re.match(repignore, cols[i]):
                    a[i] = False
                else:
                    a[i] = True
            except:
                print('WARNING: unparsable replicate label, skipping')
                a[i] = False
                pass
        indata = (indata.loc[a]).copy()
    return indata


def multifilerepimport(filedirectory, header, skiprows, labelcols):
    # CSV import
    files = glob(os.path.join(filedirectory, '*.[cC][sS][vV]'))
    if len(files) is not 0:
        print('++++++++++ Found ' + str(len(files)) + ' files ++++++++++')
        # First need to determine max time length to use as first input
        for i in range(0, (len(files))):
            if i == 0:
                testfile = pd.read_csv(files[i], header=header, skiprows=skiprows)
                time = testfile.iloc[0, :]
                time = pd.DataFrame(time)
            else:
                testfile = pd.read_csv(files[i], header=header, skiprows=skiprows)
                if testfile.shape[1] > len(time):
                    time = testfile.iloc[0, :]
                    time = pd.DataFrame(time)

        stackfile = time.transpose()
        for i in range(0, len(files)):
            newfile = pd.read_csv(files[i], header=header, skiprows=skiprows + 1)
            stackfile = stackfile.append(newfile, ignore_index=True)

        return stackfile
    else:
        # excel import
        files = glob(os.path.join(filedirectory, '*.[xX][lL][sS][xX]'))
        print('++++++++++ Found ' + str(len(files)) + ' files ++++++++++')
        for i in range(0, (len(files))):
            if i == 0:
                testfile = pd.read_excel(files[i], header=header, skiprows=skiprows)
                time = testfile.iloc[0, :]
                time = pd.DataFrame(time)
            else:
                testfile = pd.read_excel(files[i], header=header, skiprows=skiprows)
                if testfile.shape[1] > len(time):
                    time = testfile.iloc[0, :]
                    time = pd.DataFrame(time)

        stackfile = time.transpose()

        for i in range(0, len(files)):
            newfile = pd.read_excel(files[i], header=header, skiprows=skiprows + 1)
            stackfile = stackfile.append(newfile, ignore_index=True)

        return stackfile



def normalise_traces(dataset, normvalue=0.05, normby=(4,14)):
    # Normalises line by line on points 5:15
    try:
        x = dataset.shape[1]
        for i in range(0, dataset.shape[0]):
            zeroingvalue = np.mean(dataset[i, normby[0]:normby[1]])
            zeroingvalue = normvalue - zeroingvalue
            dataset[i, :] = dataset[i, :] + zeroingvalue
        return dataset
    # For single rows
    except IndexError:
        zeroingvalue = np.mean(dataset[normby[0]:normby[1]])
        zeroingvalue = normvalue - zeroingvalue
        dataset = dataset + zeroingvalue
        return dataset
    except:
        raise


def align_replicates(dataset, normvalue=0.05, alignvalue=0.1):
    if alignvalue is not None:
        diff = check_for_growth(dataset, alignvalue, normvalue)
        invdiff = np.logical_not(diff)
        if np.any(invdiff):
            print('Error, replicate does not reach alignvalue, dropping from alignment')
            dataset = dataset[np.array(diff, bool), :]
        try:
            dataset[1]
        except:
            return dataset

        alignpoint = normvalue + alignvalue
        startindexes = np.int_(dataset[:, 0])

        # Finds where data > alignpoint for 3 consecutive points
        for i in range(0, dataset.shape[0]):
            for ii in range(0, dataset.shape[1]):
                if (dataset[i, ii] > alignpoint) & (dataset[i, ii + 1] > alignpoint) & (
                            dataset[i, ii + 2] > alignpoint):
                    x = ii
                    break
            startindexes[i] = np.int_(x)

        # Aligns data into array by dropping cells in everything but lowest lag set
        minindex = startindexes.min()
        maxindex = startindexes.max()
        maxrowlength = dataset.shape[1] - maxindex - 1
        for i in range(0, dataset.shape[0]):
            rowindex = startindexes[i]
            newdata = dataset[i, rowindex - minindex:rowindex + maxrowlength - 1]
            if i == 0:
                stack = newdata
            else:
                stack = np.vstack((stack, newdata))
        return stack
    else:
        return dataset


def check_for_growth(data, growthmin, normalise):
    try:
        test = np.array(range(0, data.shape[0]))
        for i in range(0, data.shape[0]):
            growthminnorm = growthmin + normalise
            for ii in range(0, data.shape[1] - 3):
                if (data[i, ii] > growthminnorm) & (data[i, ii + 1] > growthminnorm) & (
                            data[i, ii + 2] > growthminnorm):
                    growth = True
                    break
                else:
                    growth = False
            test[i] = growth
        return test
    # For single rows
    except IndexError:
        growthminnorm = growthmin + normalise
        for i in range(0, data.shape[0] - 3):
            if (data[i] > growthminnorm) & (data[i + 1] > growthminnorm) & (data[i + 2] > growthminnorm):
                growth = True
                break
            else:
                growth = False
        test = growth
        return test
    except:
        raise


