# curve-fitter
___
Curve-fitter is an adapter module for the Swain _et al_ time derivative fitting an analysis software. The module is designed to allow for the processing of multiple files, curves and replicates in an automated fashion. 

### Features (as of v0.2.0-alpha)
* Processing of .csv and .xlsx files in folder or individual files
* Automated selection and fitting of replicates across multiple files
* Replicates are normalised and aligned before fitting for accurate growth rate reporting
___
## How to use curve-fitter
curvefitter(filename,header= None, predefinedinput= None, skiprows= 0, labelcols= 3, replicols= 3, waterwells= False, replicates= False, repignore= None, normalise= 0.05, growthmin= 0.05, alignvalue= 0.1, fitparams= {0:[-5,8], 1:[-6,-1], 2:[-5,2]}, noruns= 5, nosamples= 20, makeplots = True, showplots= True):

##### Parameters

|Parameter|Definition|Type|Default|
|---|---|---|---|
|filename |filename or folder location (only if replicates==1)|string| N/A
|predefinedinput|'BMG' sets skiprows, labelcols, replicols and repignore based on standard format for BMG platereader files|string or None|None
|skiprows|lines of input file to skip before retrieving data. First row assumed as time with data following immediately below| integer| 0
labelcols| first n columns of labels/text that are used to populate the output| integer| 3
replicols| column containing the strings used to match replicates| integer| 3
waterwells| ignores wells on the outside of 96 well plate| boolean| False
replicates| indicates presense of replicates to be used for data sorting automatically runs normalise and align on replicates to ensure most accurate GR| boolean| False
repignore| regex string that defines replicates to be ignored ie 'Sample *' for BMG files|regex string| None
normalise| value that data is normalised to at the start DO NOT USE 0(zero) or log function will fail|float| 0.05
growthmin| minimum value required for growth to be counted and fitted, fitting purely flat functions consumes time for fitting and produces unreliable results| float| 0.05
alignvalue| aligns replicates so that this value is reached at the same time for all reps| float| 0.1
fitparams| fitparameters used by the deODoriser| dict list of three value pairs|{0:[-5,8], 1:[-6,-1], 2:[-5,2]}
noruns| number of fitting attempts made| integer| 5
nosamples| number of samples used to calculate error| integer| 20
makeplots| determines if program makes plots and saves to output folder| boolean| True
showplots| displays plots during processing| boolean| True