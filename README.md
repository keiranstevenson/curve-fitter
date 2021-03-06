﻿# curve-fitter
Curve-fitter is an adapter module for the Swain _et al_ software. The Swain software is a tool for fitting of growth curves using gaussian process and extracting time derivatives, in particular the growth rate. The curve-fitter module is designed to expand the original software by processing multiple files with multiple separate curves in an automated manner.
The original fitting algorithm is provided in this package as well and can be run from the file fitderivgui. 

## Features (as of v0.3.0)
* Processing of .csv and .xlsx files in folder or individual files
* Automated selection and fitting of replicates across multiple files
* Replicates are normalised and aligned in time before fitting for accurate growth rate reporting
* Software automatically drops replicates that show no growth

## Dependencies
Scipy

# How to use curve-fitter
```Python
curvefitter(filename)
curvefitter(predefinedinput= 'BMG')
curvefitter(filename,header= None, predefinedinput= None, skiprows= 0, labelcolumns= 3, replicatecolumn= 3, replicatesexist= False, replicateignore= None, growthmin= 0.05, alignvalue= 0.1, fitparams= {0:[-5,8], 1:[-6,-1], 2:[-5,2]}, noruns= 5, nosamples= 20, makeplots = True, showplots= True):
```
### Parameters

|Parameter|Definition|Type|Default|
|---|---|---|---|
filename |Filename or folder location. If a folder is given and replicatesexist is True, then all files will be imported and any replicates with shared strings will be fitted together|string| N/A
predefinedinput|'BMG' sets skiprows, labelcols, replicols and repignore based on standard format for BMG platereader files|string or None|None
skiprows| Lines of input file to skip before retrieving data. First row assumed as time with data following immediately below| Integer| 0
labelcolumns| The number of columns at the start of the table that are not data.| integer| 3
replicatesexist| Indicates presence of replicates to be used for data sorting. This runs alignment on replicates to ensure most accurate GR. | boolean| False
replicatecolumn| Column containing the strings used to match replicates. Rows with an exact match in cells in this column will be fitted together as replicates| integer| 3
replicateignore| If true, rows with replicatecolumn = ignore will be skipped by analysis. This allows for ignoring wells that show growth but are not required| boolean| True
growthmin| Minimum value required for growth to be counted and fitted. Growth is determined as anywhere 3 consecutive points are greater than growthmin+0.01 for the curve. Fitting purely flat functions consumes time for fitting and produces unreliable results| float| 0.05
alignvalue| Aligns replicates so that this value is reached at the same time for all reps, if alignvalue=None then it is skipped. Align value must be >= growthmin. Replicates failing to meet this value will be dropped from analysis.| float or None| 0.1

|Parameters passed to Swain fitting routine|Definition (see Swain site below for better details)|Type|Default|
|---|---|---|---|
fitparams| Fit parameters used by the Swain software. Narrower parameters, fine tuned to your data are recommended for faster fitting| dictionary list of three value pairs|{0:[-5,8], 1:[-6,-1], 2:[-5,2]}
noruns| Number of fitting attempts made by the software with the best attempt selected | integer| 5
nosamples| Number of statistical samples used by the fitting program to calculate error| integer| 20
logdata| If true then data is converted to natural log before fitting occurs, such as for fitting exponential growth rates| boolean| True
makeplots| Determines if program makes plots of data+fit and derivative, and saves to output folder. | boolean| True
showplots| Displays plots during processing | boolean| False

For details on the fitting routine itself please see the references at the bottom.
While the initial parameters are designed for growth curves from a platereader, the fitting parameters can be tuned to fit a wide range of data sets. For the best fit it is recommended that you ascertain the best fitparamter by running the fitderivgui.py file and follow the instructions in the GUI to fit at least one example of your data manually. Additional info on the nature of these parameters can be found at the website in the references below.
The settings that you obtain can then be fed into the curve-fitter program to be used on the complete dataset.

### Input file format
Files should either be a csv or xlsx with the data in a row oriented format. The first line read by the program is taken as the time input.

#### example of input values and their effect on input
If replicatesexist is False

<img src=".\Images\Slide1.PNG" width="640" height="360"/>

If replicatesexist is True

<img src=".\Images\Slide2.PNG" width="640" height="360"/>

For example after filtering with replicatesexist = True, the program will find the first point at which OD>growthmin

<img src=".\Images\Slide3.PNG" width="640" height="360"/>

Program then finds where values are greater than alignvalue, and aligns by dropping values from the front of the longest row until all are aligned.

<img src=".\Images\Slide4.PNG" width="640" height="360"/>

After this fitting proceeds as described in the Swain et al algorithm.

Further examples of program calls and files are present in the examplefileimport.py script and examples folder.
# Further information
### Alignment procedure
Alignment point is determined at an index i where the values at i, i+1 and i+2 are greater than alignvalue+normalise. Where i for one replicate is greater than the i of another, data is deleted from the start of the longer replicate until both  i's are equal.
### Normalising procedure
The program finds the minimum value for a given replicate and shifts all values up or down until min(datarow) == 0.01 ie:
normaliseddata = dataset + (0.01-min(datarow))
This value was chosen as it was as close to zero as possible but without the ln function returning

# References
The original fitting routine can be found here: http://swainlab.bio.ed.ac.uk/software/fitderiv/
All references should be attributed to:
PS Swain, K Stevenson, A Leary, LF Montano-Gutierrez, IBN Clark, J Vogel, and T Pilizota. Inferring time derivatives including growth rates using Gaussian processes *Nat Commun* **7** (2016) 13766
