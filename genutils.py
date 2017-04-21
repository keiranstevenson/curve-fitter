import numpy as np

#####
def isfloat(s):
    '''
    Tests if variable is a float

    Arguments
    --
    s: variable to be tested
    '''
    try:
        float(s)
    except ValueError:
        return False
    return True


#####
def rmwhitespace(s):
    '''
    Removes white space from a string.

    Arguments
    --
    s:
    '''
    return ''.join(s.split())


#####
def getdatestring():
    '''
    Finds the current data and time as a string
    '''
    import time
    localtime = time.localtime(time.time())
    return ( str(localtime.tm_mday).zfill(2) + str(localtime.tm_mon).zfill(2) + str(localtime.tm_year)
            + '_' + str(localtime.tm_hour).zfill(2) + str(localtime.tm_min).zfill(2) )


#####
def crosscorrelate(t, s, r, mode= 'full', title= '', ylim= False):
    '''
    Find and plot normalized cross-correlation. Returns cross-correlation and lags.

    Arguments
    --
    t : 1D array of equally spaced time points
    s : signal (1D array)
    r : response (1D array)
    title: title of plot
    '''
    from scipy.signal import correlate
    import matplotlib.pyplot as plt
    
    nodata= t.size
    dt= np.median(np.diff(t))
    ns= (s - np.mean(s))/np.std(s)
    nr= (r - np.mean(r))/np.std(r)
    corr= correlate(nr, ns, mode= mode)/nodata
    if mode == 'full':
        lags= -(np.arange(np.size(corr)) - (nodata-1))*dt
    elif mode == 'same':
        lags= -(np.arange(np.size(corr)) - (int(nodata/2)-1))*dt
    else:
        lags= np.arange(np.size(corr))*dt
    
    plt.figure()
    plt.subplot(2,1,1)
    plt.suptitle(title, fontsize= 20)
    plt.plot(t, ns, t, nr)
    plt.ylabel('normalized signals')
    plt.xlabel('t')
    plt.legend(['signal', 'response'], loc= 'upper center', bbox_to_anchor=(0.5, 1.05),
          ncol=3, fancybox=True, shadow=True)
    plt.subplot(2,1,2)
    plt.stem(lags, corr)
    plt.ylabel('normalized correlation')
    plt.xlabel('lags')
    if np.any(ylim): plt.ylim(ylim)
    plt.show()
    
    return corr, lags


#####
def estimateerrorbar(y, nopts= False):
    """
    Estimates measurement error for each data point of y by calculating the standard deviation of the nopts data points closest to that data point

    Arguments
    --
    y: data - one column for each replicate
    nopts: number of points used to estimate error bars
    """
    y= np.asarray(y)
    if y.ndim == 1:
        ebar= np.empty(len(y))
        if not nopts: nopts= np.round(0.1*len(y))
        for i in range(len(y)):
            ebar[i]= np.std(np.sort(np.abs(y[i] - y))[:nopts])
        return ebar
    else:
        print('estimateerrorbar: works for 1-d arrays only.')
        


#####
def findsmoothvariance(y, filtsig= 0.1, nopts= False):
    '''
    Estimates and then smooths the variance over replicates of data

    Arguments
    --
    y: data - one column for each replicate
    filtsig: sets the size of the Gaussian filter used to smooth the variance
    nopts: if set, uses estimateerrorbar to estimate the variance
    '''
    from scipy.ndimage import filters
    if y.ndim == 1:
        # one dimensional data 
        v= estimateerrorbar(y, nopts)**2
    else:
        # multi-dimensional data
        v= np.var(y, 1)
    # apply Gaussian filter
    vs= filters.gaussian_filter1d(v, int(len(y)*filtsig))
    return vs
    

######
def rejectionsample1D(x, y, nosamples):
    """
    Uses unadulterated rejection sampling to sample values of x from the probability distribution given by y

    Arguments
    --
    x: support of distribution
    y: histgram of distribution
    nosamples: number of samples to generate
    """
    s= np.empty(nosamples)
    ymax= max(y)
    for i in range(nosamples):
        gotone= False
        while not gotone:
            trialsample= np.random.randint(0, len(x))
            if np.random.uniform(0, ymax) <= y[trialsample]:
                gotone= True
                s[i]= x[trialsample]
    return s
    

######
def importdata(fname):
    '''
    Imports a matrix of numerical data that is stored as a .csv, .txt, or .xls file

    Arguments
    --
    fname: file to be imported
    '''
    filetype= fname.split('.')[-1]
    if filetype == 'txt' or filetype == 'csv':
        from csv import reader
        if filetype == 'txt':
            fin= reader(open(fname, 'rU'), delimiter= '\t')
        else:
            fin= reader(open(fname, 'rU'))
        d= []
        for r in fin:
            if len(r) > 0:
                d.append(r)
        d= np.array(d)
    elif filetype == 'xls' or filetype == 'xlsx':
        import xlrd
        wb= xlrd.open_workbook(fname)
        s= wb.sheets()[0]
        d= []
        for i in range(s.nrows):
            for j in range(s.ncols):
                val= s.cell(i,j).value
                if isfloat(val):
                    d.append(str(val))
                else:
                    d.append(val.encode('ascii','ignore'))
        d= np.reshape(np.array(d), (s.nrows, s.ncols))
    return d


######
def makelist(c):
    '''
    Ensures that a variable is a list

    Arguments
    --
    c: variable to be made into a list
    '''
    if type(c) is not list: 
        return [c]
    else:
        return c
    
######
def genodeint(dydt, y0, t, itype, fargs= None):
    '''
    Integrates ordinary differential equations different choices of integrator

    Arguments
    --
    dydy: systems of ODEs to be solved
    y0: initial conditions
    t: time points of interest
    itype: type of integrator - 'lsoda' or 'vode'
    fargs: passed to dydt
    '''
    r= ode(dydt, None)
    yf= [y0]
    
    if itype == 'vode':
        r.set_integrator(itype, method= 'bdf', nsteps= 100000)
    elif itype == 'lsoda':
        r.set_integrator(itype, method= 'bdf', nsteps= 10000)

    if fargs:
        r.set_initial_value(y0, t[0]).set_f_params(fargs)
    else:
        r.set_initial_value(y0, t[0])    
    for dt in np.diff(t):
        r.integrate(r.t+dt)
        if not r.successful():
            print(" Integrator error" + " for " + itype)
            yf= None
            break
        else:
            yf.append(r.y)
        
    return np.array(yf)

######
def unique_rows(a):
    '''
    Finds the unique rows in an array

    Arguments
    --
    a: array of interest
    '''
    a = np.ascontiguousarray(a)
    unique_a = np.unique(a.view([('', a.dtype)]*a.shape[1]))
    return unique_a.view(a.dtype).reshape((unique_a.shape[0], a.shape[1]))

######
def mergedicts(original, update):
    '''
    Given two dicts, merge them into a new dict

    Arguments
    --
    x: first dict
    y: second dict
    '''
    z= original.copy()
    z.update(update)
    return z

######
def dict2list(d):
    '''
    Put the values of a dictionary into a list.

    Arguments
    --
    d: dictionary
    '''
    return [d[a] for a in d.keys()]

######
def nans(shape):
    '''
    Creates an array of NaNs

    Arguments
    --
    shape: shape of array to be created
    '''
    a= np.empty(shape)
    a[:]= np.nan
    return a

######
def rmcolsofnans(a):
    '''
    Removes any columns of an array that start with a NaN

    Arguments
    --
    a: array of interest
    '''
    a= np.asarray(a)
    try:
        return a[:, ~np.isnan(a[0,:])]
    except:
        # 1D array
        return a[:, ~np.isnan(a)]


######
def plotxyerr(x, y, xerr, yerr, xlabel= 'x', ylabel= 'y', title= ''):
    '''
    Plots a noisy x versus a noisy y with errorbars shown as ellipses.

    Arguments
    --
    x: x variable (a 1D array)
    y: y variable (a 1D array)
    xerr: (symmetric) error in x (a 1D array)
    yerr: (symmetric) error in y (a 1D array)
    xlabel: label for x-axis
    ylabel: label for y-axis
    title: title of figure
    '''
    import matplotlib.pyplot as plt
    from matplotlib.patches import Ellipse
    fig= plt.figure()
    ax= fig.add_subplot(111)
    ax.plot(x, y, 'b.-')
    for i in range(len(x)):
        e= Ellipse(xy= (x[i], y[i]), width= 2*xerr[i], height= 2*yerr[i])
        ax.add_artist(e)
        e.set_facecolor('#CCCCFF')
        e.set_linewidth(0)
    plt.xlim([np.min(x-2*xerr), np.max(x+2*xerr)])
    plt.ylim([np.min(y-2*yerr), np.max(y+2*yerr)])
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.show(block= False)
