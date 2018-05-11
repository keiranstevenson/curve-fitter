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


class CurveFitter:
    def __init__(self, file_location=None, header=None, predefined_input=None, skip_rows=0, label_columns=3,
                 replicate_column=3,
                 replicates_exist=False, replicate_ignore=None,
                 growth_minimum=0.05, alignment_value=0.1,
                 fiting_parameters={0: [-5, 8], 1: [-6, -1], 2: [-5, 2]}, no_runs=5, no_samples=20, logdata=True,
                 make_plots=True, show_plots=False):

        if predefined_input == 'BMG':  # For BMG output from TP lab reader
            skip_rows = 6
            label_columns = 3
            replicate_column = 3
            replicate_ignore = 'Sample X*'
        elif predefined_input == 'Tecan':
            skip_rows = 63
            label_columns = 1
        ## Program parameters
        self.normalise = 10 ** (-2)

        ## datastorage
        self.od_data = None
        self.filelist = None
        self.data_dicts = None

        ## passing variables to instances of self
        waterwells = False  # no longer necessary but left
        replicate_column = replicate_column - 1  # converts to index from number

        self.header = header
        self.predefined_input = predefined_input
        self.skip_rows = skip_rows
        self.label_columns = label_columns
        self.replicate_column = replicate_column

        self.replicates_exist = replicates_exist
        self.replicate_ignore = replicate_ignore
        self.growth_minimum = growth_minimum
        self.alignment_value = alignment_value
        self.fiting_parameters = fiting_parameters
        self.no_runs = no_runs
        self.no_samples = no_samples
        self.logdata = logdata
        self.make_plots = make_plots
        self.show_plots = show_plots

        if file_location is not None:
            self.file_import(file_location)

    def file_import(self, file_location):
        file_location = os.path.normpath(file_location)
        file_location = os.path.realpath(file_location)

        if os.path.isdir(file_location):
            files = glob(os.path.join(file_location, '*.[cC][sS][vV]')) + glob(
                os.path.join(file_location, '*.[xX][lL][sS][xX]'))
            print('++++++++++Detected folder. Processing {} files++++++++++'.format(len(files)))
        elif os.path.isfile(file_location):
            files = [file_location]
            print('++++++++++ Processing file {} ++++++++++'.format(file_location))

        data_dicts = []
        # open each file in turn and append dictionaries to data_dicts
        for fileloc in files:
            try:
                input_dataframe = pd.read_csv(fileloc, header=self.header, skiprows=self.skip_rows)
            except pd.parser.CParserError:
                input_dataframe = pd.read_excel(fileloc, header=self.header, skiprows=self.skip_rows)
            except pd.parser.ParserError:
                input_dataframe = pd.read_excel(fileloc, header=self.header, skiprows=self.skip_rows)

            file_dict = dict()
            file_dict['data'] = input_dataframe.iloc[1:, :]
            file_dict['time_data'] = input_dataframe.iloc[0, :]
            location = os.path.split(fileloc)[0]
            file_dict['output_filepath'] = os.path.join(location, 'curvefitter_outputdata')
            data_dicts.append(file_dict)

        # if replicates are present and multiple files then combine all data into a single data_dict
        # this is so all replicates can be processed together regardless of starting file
        if self.replicates_exist and len(data_dicts) > 1:
            new_file_dict = dict()
            for i, data_dict in enumerate(data_dicts):
                if i == 1:
                    new_file_dict['data'] = data_dict['data']
                else:
                    new_file_dict['data'] = new_file_dict['data'].append(data_dict['data'])
            time_list = [x['time_data'] for x in data_dicts]
            time_index = np.argmax([len(x) for x in time_list])
            new_file_dict['time_data'] = time_list[time_index]
            ### needs outputdata
            data_dicts = [new_file_dict]
        self.data_dicts = data_dicts

        ## build empty values for later
        for data_dict in self.data_dicts:
            data_dict['replicates'] = None
            data_dict['growth_check'] = None
            data_dict['fit_results'] = None

    def fit_data(self):
        self.filter_data()
        self.normalise_traces()
        self.check_for_growth()
        self.set_replicates()

    def run_fitting(self):
        for data_dict in self.data_dicts:
            # name:time:data:growth
            all_datasets = []
            # each set to be fitted needs split into lists of [name,time,data,growthstatus]
            if self.replicates_exist:
                for i, replicatename in enumerate(data_dict['replicates']):
                    dataset = [replicatename]
                    time = data_dict['time']
                    time = time.loc[self.label_columns:].values
                    dataset.append(time)
                    data = data_dict['data'][data_dict.iloc[:, self.replicate_column] == replicatename]
                    if isinstance(data, pd.Series):
                        data = data.loc[self.label_columns:].values
                    else:
                        data = data.loc[:, self.label_columns:].values
                    data = self.align_replicates(data)
                    dataset.append(data)
                    growth = data_dict['growth_check'][data_dict.iloc[:, self.replicate_column] == replicatename]
                    dataset.append(growth)
                    all_datasets.append(dataset)
            else:
                for i, index in enumerate(data_dict['data'].index.values):
                    label = data_dict['data'].loc[index, :self.label_columns].values
                    dataset = [label]
                    time = data_dict['time']
                    time = time.loc[self.label_columns:].values
                    dataset.append(time)
                    data = data_dict['data'].loc[index, :]
                    data = data.loc[self.label_columns:].values
                    dataset.append(data)
                    growth = data_dict['growth_check'][i]
                    dataset.append(growth)
                    all_datasets.append(dataset)

            for dataset in all_datasets:
                label = dataset[0]
                time = dataset[1]
                data = dataset[2]
                growth = dataset[3]

                max_index_array = np.where(np.isnan(data))
                max_index = np.min(max_index_array[1])

                t = time[:, :max_index]
                od = data[:, :max_index]
                od = od[growth]
                if od.shape[0]>0:
                    # Runs fitderiv only if growth is over growth_minimum
                    for attemptno in range(5):
                        try:
                            fitty = fitderiv(t, np.transpose(od), bd=self.fiting_parameters, exitearly=False,
                                             nosamples=self.no_samples,
                                             noruns=self.no_runs, logs=self.logdata)  # Peters program
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
                else:
                    # Returns no growth if none detected
                    datalength = data.shape[1]
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
                    print('No growth found! Less than {} change in OD detected'.format(self.growth_minimum))



    def align_replicates(self, dataset):
        if self.alignment_value is not None:

            alignpoint = self.normalise + self.alignment_value
            startindexes = np.asarray(dataset[:, 1], dtype='int')

            # Finds where data > alignpoint for 3 consecutive points
            for i in range(0, dataset.shape[0]):
                for ii in range(0, dataset.shape[1]):
                    if (dataset[i, ii] > alignpoint) & (dataset[i, ii + 1] > alignpoint) & (
                            dataset[i, ii + 2] > alignpoint):
                        x = ii
                        break
                startindexes[i] = np.int(x)

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

    def set_replicates(self):
        if self.replicates_exist:
            for data_dict in self.data_dicts:
                data = data_dict['data']
                replicates = data.iloc[:, self.replicate_column]
                replicates = np.unique(replicates)
                data_dict['replicates'] = replicates

    def check_for_growth(self):
        if self.growth_minimum is not None:
            growth_min_normed = self.growth_minimum + self.normalise
            for data_dict in self.data_dicts:
                test = np.ndarray(data_dict['data'].shape[0])
                for i in data_dict['data'].index.values:
                    for ii in range(data_dict['data'].shape[1]):
                        if (data_dict['data'].loc[i, ii] > growth_min_normed) & (
                                data_dict['data'].loc[i, ii + 1] > growth_min_normed) & (
                                data_dict['data'].loc[i, ii + 2] > growth_min_normed):
                            growth = True
                            break
                        else:
                            growth = False
                    test[i] = growth
                data_dict['growth_check'] = test

    def filter_data(self):
        for data_dict in self.data_dicts:
            if self.replicate_ignore is not None:
                if isinstance(repignore, str):
                    column = data_dict['data'].loc[:, self.replicate_column]
                    filterarray = [re.match(self.replicate_ignore, x) for x in column]
                    data_dict['data'] = data_dict['data'].loc[filterarray]

    def normalise_traces(self):
        # Normalises line by line on points 5:15
        for data_dict in self.data_dicts:
            for i in data_dict['data'].index.values:
                zeroingvalue = np.nanmin(data_dict['data'].loc[i, self.replicate_column:])
                zeroingvalue = self.normalise - zeroingvalue
                data_dict[0].loc[i, self.replicate_column:] = data_dict['data'].loc[i,
                                                              self.replicate_column:] + zeroingvalue


def curvefitter(filename, header=None, predefined_input=None, skip_rows=0, label_columns=3, replicate_column=3,
                replicates_exist=False, replicate_ignore=None,
                growth_minimum=0.05, alignment_value=0.1,
                fiting_parameters={0: [-5, 8], 1: [-6, -1], 2: [-5, 2]}, no_runs=5, no_samples=20, log_data=True,
                make_plots=True, show_plots=False):
    """
    filename: filename or folder location (only if replicates_exist==1)
    predefined_input: 'BMG' sets skip_rows, label_columns, replicate_column and replicate_ignore based on standard format for BMG platereader files
    skip_rows; lines of input file to skip before retrieving data. First row assumed as time with data following immediately below
    label_columns; first n columns are labels/text and are used to populate the output
    replicate_column; column containing the strings used to match replicates_exist

    waterwells; ignores wells on the outside of 96 well plate
    replicates_exist; indicates presence of replicates_exist to be used for data sorting automatically runs normalise and align on replicates_exist to ensure most accurate GR
    replicate_ignore; regex string that defines replicates_exist to be ignored ie 'Sample
     ' for BMG files
    growth_minimum; minimum value required for growth to be counted and fitted, fitting purely flat functions consumes time for fitting and produces unreliable results
    alignment_value; aligns replicates_exist so that this value is reached at the same time for all reps
    fiting_parameters; fitparameters used by the deODouriser
    no_runs; number of fitting attempts made
    no_samples; number of samples used to calculate error
    make_plots; determines if program makes plots and saves to output folder
    show_plots; displays plots during processing
    """

    if predefined_input == 'BMG':  # For BMG output from TP lab reader
        skip_rows = 6
        label_columns = 3
        replicate_column = 3
        replicate_ignore = 'Sample X*'
    elif predefined_input == 'Tecan':
        skip_rows = 63
        label_columns = 1

    normalise = 10**(-2)
    waterwells = False  # no longer necessary but left
    replicate_column = replicate_column - 1  # converts to index from number

    # Process files before inputting
    filename = os.path.realpath(filename)
    if replicates_exist & os.path.isdir(filename) is True:
        input_dataframe = multifilerepimport(filename, header, skip_rows, label_columns)
        output_filepath = os.path.join(filename, 'curvefitter' + ' outputdata')
        filename = os.path.split(filename)[-1]
        input_dataframe = cleannonreps(input_dataframe, replicate_column, replicate_ignore)

        reps = input_dataframe.iloc[1:, replicate_column]
        unique_replicates = np.unique(reps)


        # Provide info about datashape for interation to use
        dataheight = unique_replicates.shape[0]
        datalength = input_dataframe.shape[1]
        firstline = input_dataframe.iloc[0, label_columns - 1:].copy()
        firstline = firstline.reset_index(drop=True)


        print('++++++++++ Found {} replicates ++++++++++'.format(dataheight))
        for x in unique_replicates: print(x)
        sys.stdout.flush()

    elif os.path.isdir(filename):
        files = glob(os.path.join(filename, '*.[cC][sS][vV]')) + glob(os.path.join(filename, '*.[xX][lL][sS][xX]'))
        print('++++++++++Detected folder. Processing {} files++++++++++'.format(len(files)))
        for i in range(0, len(files)):
            filename = files[i]
            print('++++++++++ Processing file {} ++++++++++'.format(filename))
            print(filename)
            # Yay recursion
            curvefitter(filename=filename, header=header, predefined_input=predefined_input, skip_rows=skip_rows,
                        label_columns=label_columns, replicate_column=replicate_column + 1, replicates_exist=replicates_exist,
                        replicate_ignore=replicate_ignore, growth_minimum=growth_minimum,
                        alignment_value=alignment_value, fiting_parameters=fiting_parameters, no_runs=no_runs,
                        no_samples=no_samples, log_data=log_data,
                        make_plots=make_plots, show_plots=show_plots)
        return

    elif replicates_exist & os.path.isfile(filename):
        try:
            input_dataframe = pd.read_csv(filename, header=header, skiprows=skip_rows)
        except pd.parser.CParserError:
            input_dataframe = pd.read_excel(filename, header=header, skiprows=skip_rows)
        except pd.parser.ParserError:
            input_dataframe = pd.read_excel(filename, header=header, skiprows=skip_rows)

        input_dataframe = cleannonreps(input_dataframe, replicate_column, replicate_ignore)
        reps = input_dataframe.iloc[1:, replicate_column]
        unique_replicates = np.unique(reps)

        dataheight = unique_replicates.shape[0]
        datalength = input_dataframe.shape[1]

        firstline = input_dataframe.iloc[0, label_columns - 1:].copy()
        firstline = firstline.reset_index(drop=True)
        print('++++++++++ Processing file {} ++++++++++'.format(filename))
        print('++++++++++ Found {} replicates ++++++++++'.format(dataheight))
        for x in unique_replicates: print(x)

        sys.stdout.flush()

        output_filepath = os.path.split(filename)[0]
        filename = os.path.split(filename)[1]
        filename = filename.split('.')[-2]
        output_filepath = os.path.join(output_filepath, 'curvefitter' + ' outputdata')

    elif os.path.isfile(filename) is True:
        try:
            input_dataframe = pd.read_csv(filename, header=header, skiprows=skip_rows)
        except pd.parser.CParserError:
            input_dataframe = pd.read_excel(filename, header=header, skiprows=skip_rows)
        except pd.parser.ParserError:
            input_dataframe = pd.read_excel(filename, header=header, skiprows=skip_rows)

        output_filepath = os.path.split(filename)[0]
        filename = os.path.split(filename)[1]
        filename = filename.split('.')[-2]
        output_filepath = os.path.join(output_filepath, 'curvefitter' + ' outputdata')

        # Gather info about raw numerical data
        dataheight = input_dataframe.shape[0] - 1  # ignore time row
        datalength = input_dataframe.shape[1]
        firstline = input_dataframe.iloc[0]

    else:
        raise ImportError('File or directory not found')

    # Checks and makes output directories
    if not os.path.isdir(output_filepath):
        os.makedirs(output_filepath)
    if make_plots:
        if not os.path.isdir(os.path.join(output_filepath, filename + ' plots', )):
            os.makedirs(os.path.join(output_filepath, filename + ' plots', ))

    # Separate time variable
    time = input_dataframe.iloc[0]
    time = time.iloc[label_columns:]
    time = np.float64(time)

    sanitychecks(input_dataframe, label_columns, normalise, time)

    try:
        for i in range(1, dataheight + 1):
            if replicates_exist:
                location = '++++++++++ Processing replicate set {} of {} ++++++++++'.format(i,dataheight)
                print(location)
                replicate_chosen = unique_replicates[i - 1]
                replicate_index = replicate_chosen == input_dataframe.iloc[:, replicate_column]
                od = input_dataframe.loc[replicate_index]

                # Defines names for results
                labels = pd.DataFrame([replicate_chosen])
                replicate_info1 = 'Fitting ' + replicate_chosen
                print(replicate_info1)
                od = (od.iloc[:, label_columns:]).copy()

                # Converts columns to float format for fitderiv
                od_float = np.array(od, dtype='float64')
                od_float = normalise_traces(od_float, normalise)
                od_float = align_replicates(od_float, normalise, alignment_value)
                noofreps = od_float.shape[0]
                replicate_info2 = 'Found {} replicates'.format(noofreps)
                print(replicate_info2)

                # Removes NaNs from analysis
                try:
                    nantest = np.isnan(od_float)
                    for ii in range(0, nantest.shape[1]):
                        if any(nantest[:, ii]):
                            x = ii - 1
                            break
                    od_float = od_float[:, :ii]
                    t = time[:ii]
                except IndexError:
                    t = time
                    od_float = odfloat[:, len(time)-1]

                if noofreps == 0:
                    growthfound = False
                else:
                    growthfound = True
            else:
                location = '++++++++++ Processing row {:03} of {:03} ++++++++++'.format(i, dataheight)
                print(location)
                od = input_dataframe.iloc[i]
                od = od[label_columns + 1:]
                noofreps = 1

                # Defines names for results
                labels = input_dataframe.iloc[i, 0:label_columns]
                labels = labels.copy()

                # Converts columns to float format for fitderiv
                od_float = np.array(od, dtype='float64')
                od_float = normalise_traces(od_float, normalise)

                # Removes NaNs from analysis
                try:
                    nantest = np.isnan(od_float)
                    for ii in range(0, len(nantest)):
                        if nantest[ii]:
                            x = ii - 1
                            break
                    od_float = od_float[:ii]
                    t = time[:ii]
                except IndexError:
                    t = time
                    od_float = odfloat[len(time)-1]

                # Check for growth
                growthfound = check_for_growth(od_float, growth_minimum, normalise)
            sys.stdout.flush()  # Forces prints to display immediately

            if growthfound:
                # Runs fitderiv only if growth is over growth_minimum
                for attemptno in range(5):
                    try:
                        fitty = fitderiv(t, np.transpose(od_float), bd=fiting_parameters, exitearly=False, nosamples=no_samples,
                                         noruns=no_runs, logs=log_data)  # Peters program
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

                if make_plots:
                    plt.figure()
                    plt.subplot(2, 1, 1)
                    plt.plot(functime, fitty.d, 'r.')
                    plt.plot(functime, fitcurve, 'b')
                    if log_data:
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

                    if replicates_exist:
                        picname = replicate_chosen
                    elif predefined_input == 'BMG':
                        picname = labels.iloc[0] + str(labels.iloc[1])
                    elif len(labels) == 1:
                        picname = labels.astype('str')
                        picname = picname.str.cat()
                    else:
                        picname = labels.astype('str')
                        picname = picname.str.cat()
                    picname = picname + '.PNG'
                    picname = os.path.join(output_filepath, filename + ' plots', picname)
                    plt.savefig(picname)
                    if show_plots:
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
                print("No growth found! Less than " + str(growth_minimum) + ' change in OD detected')

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
            if replicates_exist:
                varnames = ['Replicate Name', 'GR', 'GR Std Error', 'Lag', 'Time of max GR', 'no. of replicates']
            else:
                x = list(range(labels.shape[0]))
                for labelindex in range(labels.shape[0]):
                    x[labelindex] = 'Label'
                varnames = x + ['GR', 'GR Std Error', 'Lag', 'Time of max GR', 'no. of replicates']
            growthrates.columns = varnames

            if log_data:
                sheetnames = ['Stats', 'Fit (logn)', 'Fit Std Error (logn)', 'Derivative', 'Derivative std Error']
            else:
                sheetnames = ['Stats', 'Fit', 'Fit Std Error', 'Derivative', 'Derivative std Error']
            outputname = os.path.join(output_filepath, filename + ' Analysed.xlsx')
            writer = pd.ExcelWriter(outputname, engine='xlsxwriter')
            growthrates.to_excel(writer, sheet_name=sheetnames[0])
            growthcurves.to_excel(writer, sheet_name=sheetnames[1])
            growthcurveserr.to_excel(writer, sheet_name=sheetnames[2])
            growthcurvesder.to_excel(writer, sheet_name=sheetnames[3])
            growthcurvesdererr.to_excel(writer, sheet_name=sheetnames[4])
            writer.save()
    return


def sanitychecks(infile, labelcolumns, normalise, time):
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
            print('ERROR in line {} \n Data is longer than time, truncating data'.format(i))
            #raise RuntimeError('Error data is longer than time')
        if len(data) > 0:
            data = normalise_traces(data, normvalue=normalise)
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
        print('++++++++++ Found {} files ++++++++++'.format(len(files)))
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

def normalise_traces(dataset, normvalue=0.01):
    # Normalises line by line on points 5:15
    #normvalue = 10**(-2)
    try:
        x = dataset.shape[1]
        for i in range(0, dataset.shape[0]):
            zeroingvalue = np.nanmin(dataset[i,:])
            zeroingvalue = normvalue - zeroingvalue
            dataset[i, :] = dataset[i, :] + zeroingvalue
        return dataset
    # For single rows
    except IndexError:
        zeroingvalue = np.nanmin(dataset[:])
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


