# curve-fitter

Curve-fitter is an adapter module for the Swain _et al_ time derivative fitting an analysis software. The module is designed to allow for the processing of multiple files, curves and replicates in an automated fashion. The original growth

### Features (as of v0.2.0-alpha)
* Processing of .csv and .xlsx files in folder or individual files
* Automated selection and fitting of replicates across multiple files
* Replicates are normalised and aligned before fitting for accurate growth rate reporting

### Dependences
Scipy

## How to use curve-fitter
curvefitter(filename,header= None, predefinedinput= None, skiprows= 0, labelcols= 3, replicols= 3, waterwells= False, replicates= False, repignore= None, normalise= 0.05, growthmin= 0.05, alignvalue= 0.1, fitparams= {0:[-5,8], 1:[-6,-1], 2:[-5,2]}, noruns= 5, nosamples= 20, makeplots = True, showplots= True):

#### Parameters

|Parameter|Definition|Type|Default|
|---|---|---|---|
|filename |Filename or folder location (only if replicates==1)|string| N/A
|predefinedinput|'BMG' sets skiprows, labelcols, replicols and repignore based on standard format for BMG platereader files|string or None|None
|skiprows|Lines of input file to skip before retrieving data. First row assumed as time with data following immediately below| Integer| 0
labelcols| First n columns of labels/text that are used to populate the output| integer| 3
replicols| Column containing the values or strings used to match replicates. Rows with an exact match (including whitespace if string) in this column will be fitted together as replicates| integer| 3
waterwells| Ignores wells on the outside of 96 well plate. Columns 1 & 12, rows A & H | boolean| False
replicates| Indicates presense of replicates to be used for data sorting automatically runs normalise and align on replicates to ensure most accurate GR. | boolean| False
repignore| Regex string that defines replicates to be ignored ie 'Sample *' for BMG files|regex string| None
normalise| Value that data is normalised to at the start DO NOT USE 0(zero) or log function will fail|float| 0.05
growthmin| Minimum value required for growth to be counted and fitted, fitting purely flat functions consumes time for fitting and produces unreliable results| float| 0.05
alignvalue| Aligns replicates so that this value is reached at the same time for all reps| float| 0.1
fitparams| Fitparameters used by the deODoriser| dict list of three value pairs|{0:[-5,8], 1:[-6,-1], 2:[-5,2]}
noruns| Number of fitting attempts made| integer| 5
nosamples| Number of samples used to calculate error| integer| 20
makeplots| Determines if program makes plots and saves to output folder| boolean| True
showplots| Displays plots during processing| boolean| True

For details on the fitting routine itself please see the references at the bottom.
While the initial parameters are designed for growth curves from a platereader, the fitting parameters can be tuned to fit a wide range of data sets. For the best fit it is reccomended that you ascertain the best fitparamter by running the fitderivgui.py file and follow the instructions in the GUI to fit at least one example of your data manually. Additional info on the nature of these parameters can be found at the webiste in the references below.
The settings that you obtain can then be fed into the curve-fitter program to be used on the complete dataset.

#### Input file format
Files should either be a csv or xlsx with the data in a row oriented format. The first line read by the program is taken as the time input.

##### example input
For the table below:
|rowno\columnno|1|2|3|4|5|...|
|---|---|---|---|---|---|---|
1|Row|Column|Name|read1|read2
2|Row|Column|Name|0|0.1
3|A|1|Control|0.1|0.2
4|A|2|Cond1|0.1|0.2
5|A|2|Cond1|0.1|0.2
...|
```Python
from curvefitter import curvefitter

curvefitter('Example.xlsx', skiprows=1, labelcols=3, replicols=3, replicates=True)
```
Row 4 and 5 will be fitted together as replicates, row 3 will be fitted on its own. Graphs will be displated with each fitting and saved to the output folder.
```
Output structure:
\datalocation\
    inputfile.xlsx
    \inputfile outputdata\
        inputfile Analysed
        \Plots\
            Control.png
            Cond1.png
```
Further examples of program calls and files are present in the examplefileimport.py script and examples folder.
## References
The original fitting routine can be found here: http://swainlab.bio.ed.ac.uk/software/fitderiv/
All references should be attributed to:
PS Swain, K Stevenson, A Leary, LF Montano-Gutierrez, IBN Clark, J Vogel, and T Pilizota. Inferring time derivatives including growth rates using Gaussian processes *Nat Commun* **7** (2016) 13766
