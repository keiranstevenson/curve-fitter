# curve-fitter
Curve-fitter is an adapter module for the Swain _et al_ software. The Swain software is a tool for fitting of growth curves using gaussian process and extracting time derivatives, in particular the growth rate. The curve-fitter module is designed to expand the original software by processing multiple files with multiple separate curves in an automated manner.
The original fitting algorithm is provided in this package as well and can be run from the file fitderivgui. 

## Features (as of v0.3.0)
* Processing of .csv and .xlsx files in folder or individual files
* Automated selection and fitting of replicates across multiple files
* Replicates are normalised and aligned in time before fitting for accurate growth rate reporting
* Software automatically drops replicates that show no growth

## Dependences
Scipy

# How to use curve-fitter
curvefitter(filename,header= None, predefinedinput= None, skiprows= 0, labelcols= 3, replicol= 3, waterwells= False, replicates= False, repignore= None, normalise= 0.05, growthmin= 0.05, alignvalue= 0.1, fitparams= {0:[-5,8], 1:[-6,-1], 2:[-5,2]}, noruns= 5, nosamples= 20, makeplots = True, showplots= True):

### Parameters

|Parameter|Definition|Type|Default|
|---|---|---|---|
filename |Filename or folder location (only if replicates==1)|string| N/A
predefinedinput|'BMG' sets skiprows, labelcols, replicol and repignore based on standard format for BMG platereader files|string or None|None
skiprows| Lines of input file to skip before retrieving data. First row assumed as time with data following immediately below| Integer| 0
labelcols| First n columns of labels/text that are used to populate the output| integer| 3
replicol| Index of column containing the values or strings used to match replicates. Rows with an exact match (including whitespace if string) in this column will be fitted together as replicates| integer| 3
waterwells| Ignores wells on the outside of a 96 well plate. Columns 1 & 12, rows A & H.| boolean| False
replicates| Indicates presence of replicates to be used for data sorting. This automatically runs normalise and align on replicates to ensure most accurate GR. | boolean| False
repignore| Regex formatted string that defines replicates to be ignored ie 'Sample *' for BMG files|regex string| None
normalise| Value that data is normalised to before fitting. Data is normalised such that the mean value of points 5-15 is equal to normalise. DO NOT USE 0(zero) or log function will fail|float| 0.05
normby| Touple of two values that define the range that is normalised to normalise value|touple of ints|(4,14)
growthmin| Minimum value required for growth to be counted and fitted. Growth is determined as anywhere 3 consecutive points are greater than growthmin+minimum value for the curve. Fitting purely flat functions consumes time for fitting and produces unreliable results| float| 0.05
alignvalue| Aligns replicates so that this value is reached at the same time for all reps, if alignvalue=None then it is skipped. Alingment point is determined as i where i, i+1 and i+2 are greater than alignvalue+normalise. Where i is ifferent for each replicate, the start of the data is removed until all i's are equal to the minimum i found.| float or None| 0.1
fitparams| Fit parameters used by the Swain software. Narrower parameters, fine tuned to your data are recommended for faster fitting| dictionary list of three value pairs|{0:[-5,8], 1:[-6,-1], 2:[-5,2]}
noruns| Number of fitting attempts made by the software with the best attempt selected | integer| 5
nosamples| Number of samples used to calculate error| integer| 20
logdata| If true then data is converted to natural log before fitting occurs, such as for fitting exponential growth rates| boolean| True
makeplots| Determines if program makes plots of data+fit and derivative, and saves to output folder. Plots made i | boolean| True
showplots| Displays plots during processing | boolean| True

For details on the fitting routine itself please see the references at the bottom.
While the initial parameters are designed for growth curves from a platereader, the fitting parameters can be tuned to fit a wide range of data sets. For the best fit it is reccomended that you ascertain the best fitparamter by running the fitderivgui.py file and follow the instructions in the GUI to fit at least one example of your data manually. Additional info on the nature of these parameters can be found at the website in the references below.
The settings that you obtain can then be fed into the curve-fitter program to be used on the complete dataset.

### Input file format
Files should either be a csv or xlsx with the data in a row oriented format. The first line read by the program is taken as the time input.

#### example of input values
<img src=".\Images\Slide1.PNG" width="600" height="500"/>

![](.\Images\Slide2.PNG =500x)

![](.\Images\Slide3.PNG)

![](.\Images\Slide4.PNG)

Final ouptu
![](.\Images\Slide5.PNG)

Further examples of program calls and files are present in the examplefileimport.py script and examples folder.
# References
The original fitting routine can be found here: http://swainlab.bio.ed.ac.uk/software/fitderiv/
All references should be attributed to:
PS Swain, K Stevenson, A Leary, LF Montano-Gutierrez, IBN Clark, J Vogel, and T Pilizota. Inferring time derivatives including growth rates using Gaussian processes *Nat Commun* **7** (2016) 13766
