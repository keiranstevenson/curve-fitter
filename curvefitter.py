# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 16:09:32 2016

@author: Keiran Stevenson
"""

import os
import platform
import re
from glob import glob

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

from fitderiv import fitderiv


class CurveFitter:
    def __init__(self, file_location=None, header=None, predefined_input=None, skip_rows=0, label_columns=3,
                 replicate_column=3, replicates_exist=False, replicate_ignore=None,
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
        self.data_start = int(label_columns + 1)

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

            location = os.path.join(os.path.split(fileloc)[0], 'curvefitter_outputdata')
            filename = os.path.split(fileloc)[1].split('.')[0]
            file_dict['output_filepath'] = os.path.join(location, filename + '.xlsx')
            file_dict['plot_output_filepath'] = os.path.join(location,filename+' plots')
            if os.path.exists(location) is False:
                os.mkdir(location)
            data_dicts.append(file_dict)

        # if replicates are present and multiple files then combine all data into a single data_dict
        # this is so all replicates can be processed together regardless of starting file
        if self.replicates_exist and len(data_dicts) > 1:
            new_file_dict = dict()
            for i, data_dict in enumerate(data_dicts):
                if i == 0:
                    new_file_dict['data'] = data_dict['data']
                else:
                    new_file_dict['data'] = new_file_dict['data'].append(data_dict['data'],ignore_index=True)
            time_list = [x['time_data'] for x in data_dicts]
            time_index = np.argmax([len(x) for x in time_list])
            new_file_dict['time_data'] = time_list[time_index]
            ### needs outputdata
            location = file_location
            new_file_dict['output_filepath'] = os.path.join(location, 'curvefitter_outputdata','Data.xlsx')
            new_file_dict['plot_output_filepath'] = os.path.join(location, 'curvefitter_outputdata', 'plots')

            data_dicts = [new_file_dict]
        self.data_dicts = data_dicts

        ## build empty values for later
        for data_dict in self.data_dicts:
            newfolder = data_dict['plot_output_filepath']
            if os.path.exists(newfolder) is False:
                os.mkdir(newfolder)
            data_dict['replicates'] = None
            data_dict['growth_check'] = None
            data_dict['fit_results'] = None

    def fit_data(self):
        self.filter_data()
        self.normalise_traces()
        self.check_for_growth()
        self.set_replicates()
        self.sanity_check_time()
        self.run_fitting()

    def sanity_check_time(self):
        for data_dict in self.data_dicts:
            time = np.atleast_2d(data_dict['time_data'].values)
            grad = np.gradient(time[0,self.data_start:])
            if any(grad<0):
                print('+++ERROR TIME DOES NOT INCREASE ALONG ENTIRE AXIS+++\nguessing missing time values')
                stepvalue = np.mean(grad[grad>0])
                indexes = np.where(grad<0)
                for indexval in indexes[-1]:
                    data_dict['time_data'].iloc[indexval+1+self.data_start] = data_dict['time_data'].iloc[indexval+self.data_start]+stepvalue

    def set_indexes(self, datadict):
        timelist = datadict['time_data'].values
        indexlist = []
        n = 0
        for item in timelist:
            test = isinstance(item,str)
            if type(item) == str:
                indexlist.append(item)
            else:
                newitem = float(item)
                if np.isnan(newitem):
                    indexlist.append('label_' + str(n))
                    n += 1
                else:
                    indexlist.append(newitem)
        datadict['data'].columns = indexlist

    def run_fitting(self):
        for data_dict in self.data_dicts:
            self.set_indexes(data_dict)
            # name:time:data:growth
            all_datasets = []
            # each set to be fitted needs split into lists of [name,time,data,growthstatus]
            df_stats = pd.DataFrame()
            df_fitted_curves = pd.DataFrame()
            df_fitted_curves_err = pd.DataFrame()
            df_deriv_curves = pd.DataFrame()
            df_derive_curves_err = pd.DataFrame()

            if self.replicates_exist:
                for i, replicatename in enumerate(data_dict['replicates']):


                    dataset = [pd.Series(replicatename)]
                    time = data_dict['time_data']
                    time = time.loc[self.data_start:].values
                    dataset.append(time)
                    data = data_dict['data'][data_dict['data'].iloc[:, self.replicate_column] == replicatename]
                    if isinstance(data, pd.Series):
                        data = data.iloc[self.data_start:].values
                    else:
                        data = data.iloc[:, self.data_start:].values
                    dataset.append(data)
                    testfilter = [data_dict['data'].iloc[:, self.replicate_column] == replicatename]
                    growth = np.array(data_dict['growth_check'])[testfilter]
                    dataset.append(growth)
                    all_datasets.append(dataset)
            else:
                for i, index in enumerate(data_dict['data'].index.values):
                    label = data_dict['data'].loc[index]
                    label = label.iloc[:self.label_columns]
                    dataset = [label]
                    time = data_dict['time_data']
                    time = time.iloc[self.data_start:].values
                    dataset.append(time)
                    data = data_dict['data'].loc[index, :]
                    data = data.iloc[self.data_start:].values
                    dataset.append(data)
                    growth = [data_dict['growth_check'][i]]
                    dataset.append(growth)
                    all_datasets.append(dataset)

            for dataset in all_datasets:
                label = dataset[0]
                if isinstance(label,pd.Series):
                    print('fitting {}'.format(label.values))
                else:
                    print('fitting {}'.format(label))
                time = np.atleast_2d(dataset[1])
                data = np.atleast_2d(dataset[2])
                growth = dataset[3]
                data = np.array(data, dtype=np.float64)
                time = np.array(time, dtype=np.float64)

                data = data[growth]
                if len(data)>0:
                    data = self.align_replicates(data)

                    ## NaN removal
                    data[data<0] = np.nan

                    nanfilter = np.isnan(data)
                    nanfilter = np.any(nanfilter,axis=0)
                    t = time[:,:data.shape[1]]
                    t = t[:,np.logical_not(nanfilter)]
                    t=t[0]
                    od = data[:,np.logical_not(nanfilter)]


                    ## old NaN removal
                    # max_index_array = np.where(np.isnan(data))
                    # max_index_array = max_index_array[-1]
                    # if any(max_index_array):
                    #     max_index = np.min(max_index_array)
                    # else:
                    #     max_index = data.shape[-1]
                    #
                    #
                    # t = time[:, :max_index]
                    # t = t[0]
                    # od = data[:, :max_index]
                    ##

                    if od.shape[0] > 0:
                        # Runs fitderiv only if growth is over growth_minimum
                        for attemptno in range(5):
                            try:
                                fitty = fitderiv(t.astype(np.float64), np.transpose(od), bd=self.fiting_parameters,
                                                 exitearly=False,
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
                        noofreps = od.shape[0]
                        fitcurve = fitty.f
                        fitcurveerr = fitty.fvar
                        fitdercurve = fitty.df
                        fitdercurveerr = fitty.dfvar
                        functime = fitty.t

                        if self.make_plots:
                            plt.figure()
                            plt.subplot(2, 1, 1)
                            plt.plot(functime, fitty.d[...,0], 'r.')
                            plt.plot(functime, fitcurve, 'b')
                            if self.logdata:
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

                            nametext = '_'.join([str(x) for x in label.values])
                            picname = nametext + '.PNG'
                            picname = os.path.join(data_dict['plot_output_filepath'], picname)
                            plt.savefig(picname)
                            if self.show_plots:
                                plt.ion()
                                plt.show()
                                plt.pause(0.2)
                            else:
                                plt.close()

                    else:
                        # Returns no growth if none detected
                        datalength = len(t)
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

                    labelheader = []
                    for i in range(label.size):
                        labelheader.append(str(i) + '_label')
                    label.index = labelheader
                    statnames = labelheader.copy()
                    statnames.extend(
                        ['no. of replicates', 'Growth Rate', 'Growth Rate SE', 'Lag', 'Time of max GR', 'notes on fitting'])
                    stats = label.copy()
                    stats = stats.append(pd.Series([noofreps, gr, err, lag, timeofmax_gr, note]))
                    stats.index = statnames

                    stats.name = i
                    fitcurve = label.append(pd.Series(fitcurve, index=t, name=i))
                    fitcurveerr = label.append(pd.Series(fitcurveerr, index=t, name=i))
                    fitdercurve = label.append(pd.Series(fitdercurve, index=t, name=i))
                    fitdercurveerr = label.append(pd.Series(fitdercurveerr, index=t, name=i))

                    fitcurve.name = i
                    fitcurveerr.name = i
                    fitdercurve.name = i
                    fitdercurveerr.name = i

                    df_stats = df_stats.append(stats)
                    df_fitted_curves = df_fitted_curves.append(fitcurve)
                    df_fitted_curves_err = df_fitted_curves_err.append(fitcurveerr)
                    df_deriv_curves = df_deriv_curves.append(fitdercurve)
                    df_derive_curves_err = df_derive_curves_err.append(fitdercurveerr)

            allcolumns = df_fitted_curves.columns
            allcolumns = self.sort_columns(allcolumns)

            df_fitted_curves=df_fitted_curves.reindex(allcolumns, axis=1,copy=False)
            df_fitted_curves_err=df_fitted_curves_err.reindex(allcolumns, axis=1,copy=False)
            df_deriv_curves=df_deriv_curves.reindex(allcolumns, axis=1,copy=False)
            df_derive_curves_err=df_derive_curves_err.reindex(allcolumns, axis=1,copy=False)

            outputname = data_dict['output_filepath']
            sheetnames = ['Stats', 'Fit', 'Fit std error', 'Derivative', 'Derivative std error']
            writer = pd.ExcelWriter(outputname, engine='xlsxwriter')
            df_stats.to_excel(writer, sheet_name=sheetnames[0])
            df_fitted_curves.to_excel(writer, sheet_name=sheetnames[1])
            df_fitted_curves_err.to_excel(writer, sheet_name=sheetnames[2])
            df_deriv_curves.to_excel(writer, sheet_name=sheetnames[3])
            df_derive_curves_err.to_excel(writer, sheet_name=sheetnames[4])
            writer.save()

    def sort_columns(self,a):
        d = {}
        for x in a:
            d.setdefault(type(x), []).append(x)

        # Sort each type
        d = {k: sorted(v) for k, v in d.items()}

        # The result list
        sorted_list = d[str]+d[float]
        return sorted_list

    def align_replicates(self, dataset):
        if self.alignment_value is not None and dataset.shape[0] > 1:

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
                test = [False]*data_dict['data'].shape[0]
                for n, i in enumerate(data_dict['data'].index.values):
                    for ii in range(self.data_start, data_dict['data'].shape[1]-2):
                        if (data_dict['data'].loc[i, ii] > growth_min_normed) & (
                                data_dict['data'].loc[i, ii + 1] > growth_min_normed) & (
                                data_dict['data'].loc[i, ii + 2] > growth_min_normed):
                            growth = True
                            break
                        else:
                            growth = False
                    test[n] = growth
                data_dict['growth_check'] = test

    def filter_data(self):
        for data_dict in self.data_dicts:
            if self.replicate_ignore is not None:
                if isinstance(self.replicate_ignore, str):
                    column = data_dict['data'].loc[:, self.replicate_column]
                    filterarray = []
                    for x in column:
                        if type(x) == str:
                            test = re.match(self.replicate_ignore, x)
                            if test is None:
                                filterarray.append(True)
                            else:
                                filterarray.append(False)
                        else:
                            filterarray.append(False)
                    data_dict['data'] = data_dict['data'].loc[filterarray]

    def normalise_traces(self):
        # Normalises line by line on points 5:15
        for data_dict in self.data_dicts:
            for i in data_dict['data'].index.values:
                datasub = data_dict['data'].loc[i, self.data_start:]
                zeroingvalue = np.nanmean(data_dict['data'].loc[i, self.data_start+4:14])
                zeroingvalue = 0-zeroingvalue
                data_dict['data'].loc[i, self.data_start:] = data_dict['data'].loc[i,
                                                             self.data_start:] + zeroingvalue

    def normalise_traces_to_value(self):
        # Normalises to the lowest possible value
        for data_dict in self.data_dicts:
            for i in data_dict['data'].index.values:
                datasub = data_dict['data'].loc[i, self.data_start:]
                zeroingvalue = np.nanmin(data_dict['data'].loc[i, self.data_start:])
                zeroingvalue = self.normalise - zeroingvalue
                data_dict['data'].loc[i, self.data_start:] = data_dict['data'].loc[i,
                                                             self.data_start:] + zeroingvalue


