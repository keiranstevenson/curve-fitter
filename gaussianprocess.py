"""
    Encodes fitting and interpolating using Gaussian processes following Rasmussen and Williams (2006).

    The Gaussian process algorithms come from chapter 3 and the Jacobian of the negative log likelihood from chapter 5 of Rasmussen and Williams. 

    Covariance functions can either be linear, squared exponential, neural network-like, or squared exponential with a linear trend. Bounds for hyperparameters are specified in log10 space. Hyperparameters are given in log space.

    A typical workflow is:

    g= gp.nnGP({0:(-4, 2), 1: (-3, 1), 2: (-4,-2)}, x, y)
    g.findhyperparameters()
    g.predict(x)
    g.sketch('.')

    Prior functions can also be sampled. For example,

    g= gp.sqexplinGP({0: (-2,2), 1: (-2,2), 2: (-2,2), 3: (-2,2), 4: (-2,2)}, x, y)
    plot(x, g.sampleprior(3, th= [1.0, 0.1, 3.1, 1.3]))

    will plot three samples of the prior latent functions with hyperparameters 1.0, 0.1, 3.1, and 1.3. There is no need to specify the hyperparameter for measurement error: it is not used to generate prior functions.
     
"""

import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt
import genutils as gu


class gaussianprocess:
    def __init__(self, lthbounds, x, y, merrors= False):
        '''
        Creates a Gaussian process.

        Arguments
        --
        lthbounds: a dictionary of pairs of the bounds on the hyperparameters in log10 space,
        such as {0: [0,6], 1: [-3,4], 2: [-6,-4]}
        x: a 1-d array of the abscissa data
        y: a multi-dimensional array of the ordinate data
        merrors: if specified, a 1-d array of the measurement errors (as variances) 
        v'''
        self.b= gu.dict2list(lthbounds)
        self.x, self.y= x, y
        self.merrors= merrors


    def covfn(self):
        raise NotImplementedError(' No covariance function specified in class %s'
                                  % self.__class__.__name__)
    
    def d1covfn(self):
        raise NotImplementedError(' No first derivative of the covariance function specified in class %s' % self.__class__.__name__)

    def d1d2covfn(self):
        raise NotImplementedError(' No second derivative of the covariance function specified in class %s' % self.__class__.__name__)
    
    
    def kernelmatrix(self, lth, x):
        """
        Returns kernel matrix K(X,X) supplemented with measurement noise and its Cholesky decomposition.

        Arguments
        --
        lth: log of the hyperparameters
        x: abscissa values
        merrors: if specified, rescales the fitted measurement error 
        """
        k= np.empty((len(x), len(x)))
        for i in range(len(x)):
            k[i,:]= self.covfn(x[i], x, lth)[0]
        if np.any(self.merrors):
            kn= k + np.exp(lth[-1])*np.diag(self.merrors)
        else:
            kn= k + np.exp(lth[-1])*np.identity(len(x))
        L= linalg.cho_factor(kn)
        return k, L


    def nlml(self, lth):
        """
        Returns negative of log marginal likelihood.

        Arguments
        --
        lth: log of the hyperparameters 
        """
        x, y= self.x, self.y
        k, L= self.kernelmatrix(lth, x)
        al= linalg.cho_solve(L, y)
        halfdetK= np.sum(np.log(np.diagonal(L[0])))
        return 0.5*np.dot(y, al) + halfdetK + 0.5*len(y)*np.log(2*np.pi)

    
    def jacnlml(self, lth):
        """
        Returns the Jacobian of negative log marginal likelihood with respect to the hyperparameters with deriviatives being taken assuming the hyperparmaters are in log space.

        Arguments
        --
        lth: log of the hyperparameters 
        """
        x, y= self.x, self.y
        k, L= self.kernelmatrix(lth, x)
        # find derivatives of kernel matrix wrt hyperparameters
        kjac= np.empty((len(x), len(x), len(lth)))
        for i in range(len(x)):
            kjac[i, :, :-1]= self.covfn(x[i], x, lth)[1]
        if np.any(self.merrors):
            kjac[:, :, -1]= np.diag(self.merrors)*np.exp(lth[-1])
        else:
            kjac[:, :, -1]= np.identity(len(x))*np.exp(lth[-1])
        # calculate jacobian
        al= linalg.cho_solve(L, y)
        alal= np.outer(al, al)
        Kinv= linalg.cho_solve(L, np.identity(len(x)))
        return np.asarray([-0.5*np.trace(np.dot(alal - Kinv, kjac[:,:,i])) for i in range(len(lth))])

    

    def findhyperparameters(self, noruns= 1, exitearly= False, stvals= False, optmethod= 'l_bfgs_b',
                            optmessages= False, quiet= True):
        """
        Finds the best fit hyperparameters (.lth_opt) and the optimum value of negative log marginal likelihood (.nlml_opt).

        Arguments
        --
        noruns: number of attempts to find optimum hyperparameters (the best of all the runs is chosen)
        exitearly: if True, fitting stops at the first successful attempt
        stvals: an (optional) initial guess for the log hyperparameters
        optmethod: the optimization routine to be used, either 'l_bfgs_b' (default) or 'tnc'
        optmessages: if True, display messages from the optimization routine
        quiet: if True, print warning that if an optimum hyperparameter is at a bound
        """
        b= self.b
        self.hparamerr= []
        lmlml= np.empty(noruns)
        lthf= np.empty((noruns, len(b))) 
        success= np.empty(noruns)
        lth= np.empty(len(b))
        # convert b into exponential base
        b= np.array(b)*np.log(10)
        # run optimization
        for i in range(noruns):
            if np.any(stvals):
                # initial values given for hyperparameters
                lth= stvals
            else:
                # choose random initial values for hyperparameters
                lth= [np.random.uniform(b[j][0], b[j][1]) for j in range(len(b))]
            # run Gaussian process
            if optmethod == 'tnc':
                from scipy.optimize import fmin_tnc
                lthf[i,:], nf, success[i]= fmin_tnc(self.nlml, lth, fprime= self.jacnlml,
                                                  bounds= b, maxfun= 1000, messages= optmessages)
                lmlml[i]= self.nlml(lthf[i,:])
            elif optmethod == 'l_bfgs_b':
                from scipy.optimize import fmin_l_bfgs_b
                lthf[i,:], lmlml[i], dout= fmin_l_bfgs_b(self.nlml, lth, fprime= self.jacnlml,
                                                        bounds= b, disp= optmessages)
                success[i]= dout['warnflag'] + 1
            else:
                print(optmethod + ' unrecognized.')
            if success[i] != 1:
                print(' Warning: optimization failed at run ' + str(i+1))
            else:
                if exitearly: break
        # only process runs that did not converge
        if np.any(success == 1):
            lmlml= lmlml[success == 1]
            lthf= lthf[success == 1]
            # find best choice
            lthb= lthf[lmlml.argmin()]
            self.nlml_opt= lmlml.min()
            # print warning
            for i in range(len(b)):
                if (lthb[i] == b[i][1] or lthb[i] == b[i][0]):
                    if not quiet:
                        print( ' Warning: hparam[' + str(i) + '] is at a boundary')
                        print('\thparam[' + str(i) + ']= {:e}'.format(np.exp(lthb[i]))
                          + ' [{:e}'.format(np.exp(b[i][0])) + ', {:e}]'.format(np.exp(b[i][1])))
                    if lthb[i] == b[i][1]:
                        self.hparamerr.append([i, 'u'])
                    else:
                        self.hparamerr.append([i, 'l'])
            self.lth_opt= lthb
        else:
            raise gaussianprocessException('Optimization of hyperparameters failed')



    def results(self):
        '''
        Displays results from optimizing hyperparameters.
        '''
        print('log(max likelihood)= %e' % (-self.nlml_opt))
        for j, pv in enumerate(np.exp(self.lth_opt)):
            print('hparam[' + str(j) + ']= {:e}'.format(pv) + ' [{:e}'.format(10**(self.b[j][0]))
                + ', {:e}]'.format(10**(self.b[j][1])))

            
     
    def sample(self, size= 1):
        '''
        Generate samples from the Gaussian process as an array.

        Arguments
        --
        size: number of samples
        '''     
        try:
            return np.transpose(np.random.multivariate_normal(self.mnp, self.covp, size))
        except AttributeError:
            print( ' Run gp.predict() first before sampling.')



            
    def sampleprior(self, size= 1, th= False):
        '''
        Generate samples from the prior of the Gaussian process as an array.

        Arguments
        --
        size: number of samples
        th: hyperparameters to use (if not specified, the hyperparameters are chosen at random)
        '''
        x, y, b= self.x, self.y, self.b
        if th:
            # hyperparameters are given (measurement error is not necessary)
            lth= np.log(np.asarray(th))
            if len(lth) == self.noparams:
                lth= np.concatenate((lth, [1.0]))
        else:
            # sample random hyperparameters
            lth= np.log(np.power(10, [np.random.uniform(b[i][0], b[i][1]) for i in range(len(b))]))
        cov= self.kernelmatrix(lth, x)[0]
        return np.transpose(np.random.multivariate_normal(np.zeros(len(x)), cov, size))
    

         

    def predict(self, xnew, merrorsnew= False, derivs= 0, addnoise= False):
        """
        Determines the predicted mean latent function (.f) and its variance (.fvar) and potentially the predicted mean first derivative (.df) and its variance (.dfvar) and the predicted mean second derivative (.ddf) and its variance (.ddfvar) . Also .mnp is the predicted combined array of the mean latent function and its mean derivatives and .covp is the corresponding covariance matrix.

        Arguments
        --
        xnew: abscissa values for which predicted ordinate values are desired
        merrorsnew: if specified, the expected measurements errors at xnew (need not be specified if xnew= x)
        derivs: if 0, only the latent function is inferred; if 1, the latent function and the first derivative are inferred; if 2, the latent function and the first and second derivatives are inferred
        addnoise: if True, add measuremnet noise to the predicted variance
        """
        if len(self.x) == len(xnew) and (self.x == xnew).all():
            xold= True
        else:
            xold= False
        if np.any(self.merrors) and not np.any(merrorsnew) and not xold:
            print('Length of xnew is different from x.')
            raise gaussianprocessException('Measurement errors were used to find the hyperparameters and measurement errors are therefore required for any predictions.')
        elif not hasattr(self, 'lth_opt'):
            raise gaussianprocessException(' Run gp.findhyperparameters() first before making predictions.')
        else:
            # set up
            lth, x, y= self.lth_opt, self.x, self.y
            # work with an array of length 3*N: the first N values being the function,
            # the second N values being the first derivative, and the last N values being the second derivative
            Knewold= np.empty((len(xnew), len(x)))
            Knewnew= np.empty((len(xnew), len(xnew)))
            if derivs > 0:
                d1Knewold= np.empty((len(xnew), len(x)))
                d1Knewnew= np.empty((len(xnew), len(xnew)))
                d1d2Knewnew= np.empty((len(xnew), len(xnew)))
            if derivs > 1:
                d12Knewold= np.empty((len(xnew), len(x)))
                d12Knewnew= np.empty((len(xnew), len(xnew)))
                d12d2Knewnew= np.empty((len(xnew), len(xnew)))
                d12d22Knewnew= np.empty((len(xnew), len(xnew)))
            for i in range(len(xnew)):
                Knewold[i,:]= self.covfn(xnew[i], x, lth)[0]
                Knewnew[i,:]= self.covfn(xnew[i], xnew, lth)[0]
                if derivs > 0:
                    d1Knewold[i,:]= self.d1covfn(xnew[i], x, lth)[0]
                    d1Knewnew[i,:]= self.d1covfn(xnew[i], xnew, lth)[0]
                    d1d2Knewnew[i,:]= self.d1d2covfn(xnew[i], xnew, lth)[0]
                if derivs > 1:
                    d12Knewold[i,:]= self.d12covfn(xnew[i], x, lth)[0]
                    d12Knewnew[i,:]= self.d12covfn(xnew[i], xnew, lth)[0]
                    d12d2Knewnew[i,:]= self.d12d2covfn(xnew[i], xnew, lth)[0]
                    d12d22Knewnew[i,:]= self.d12d22covfn(xnew[i], xnew, lth)[0]
            if derivs == 0:
                kv= Knewold
                km= Knewnew
            elif derivs == 1:
                kv= np.bmat([[Knewold], [d1Knewold]])
                km= np.bmat([[Knewnew, np.transpose(d1Knewnew)],
                             [d1Knewnew, d1d2Knewnew]])
            elif derivs == 2:
                kv= np.bmat([[Knewold], [d1Knewold], [d12Knewold]])
                km= np.bmat([[Knewnew, np.transpose(d1Knewnew), np.transpose(d12Knewnew)],
                              [d1Knewnew, d1d2Knewnew, np.transpose(d12d2Knewnew)],
                              [d12Knewnew, d12d2Knewnew, d12d22Knewnew]])
            # find mean prediction
            k, L= self.kernelmatrix(lth, x)
            m= np.dot(kv, linalg.cho_solve(L, y))
            mnp= np.reshape(np.array(m), np.size(m))
            self.mnp= mnp
            # find variance of prediction
            covp= km - np.dot(kv, linalg.cho_solve(L, np.transpose(kv)))
            self.covp= covp
            varp= np.diag(covp)
            # for user
            self.f= mnp[:len(xnew)]
            self.fvar= varp[:len(xnew)]
            if addnoise:
                # add measurement error to the variance of the latent function
                if np.any(merrors):
                    if xold:
                        self.var += np.exp(lth[-1])*np.diag(self.merrors)
                    else:
                        self.var += merrorsnew
                else:
                    self.var += np.exp(lth[-1])*np.identity(len(xnew))
            if derivs > 0:
                self.df= mnp[len(xnew):2*len(xnew)]
                self.dfvar= varp[len(xnew):2*len(xnew)]
            if derivs > 1:
                self.ddf= mnp[2*len(xnew):]
                self.ddfvar= varp[2*len(xnew):]
            



    def sketch(self, datasymbol= 'o'):
        """
        Plots data with mean prediction plus band of twice the standard deviation.

        Arguments
        --
        datasymbol: the symbol used to mark the data points 
        """
        x, y= self.x, self.y
        f= self.f
        sd= np.sqrt(self.fvar)
        plt.figure()
        plt.plot(x, y, 'r'+datasymbol, x, f, 'b')
        plt.fill_between(x, f-2*sd, f+2*sd, facecolor='blue', alpha=0.2)
        plt.show(block= False)



####        

class lnGP(gaussianprocess):
    '''
    Gaussian process with a linear covariance function
    '''
    noparams= 2
    description= 'linear Gaussian process'
        
    def covfn(self, x, xp, lth):
        '''
        Linear covariance function.
        Returns the kernel function and Jacobian of the kernel function.

        Arguments
        --
        x: a 1-d array of abscissa
        xp: a 1-d array of alternative abscissa
        lth: the log of the hyperparameters
        '''
        th= np.exp(lth)
        k= th[0] + th[1]*x*xp
        jk= np.empty((len(xp), self.noparams))
        jk[:,0]= th[0]*np.ones(len(xp))
        jk[:,1]= th[1]*x*xp
        return k, jk

    def gradcovfn(self, x, xp, lth):
        '''
        Returns the gradient of the covariance function with respect to x.

        Arguments
        --
        x: a 1-d array of abscissa
        xp: a 1-d array of alternative abscissa
        lth: the log of the hyperparameters
        '''
        th= np.exp(lth)
        return th[1]*xp, False

    def hesscovfn(self, x, xp, lth):
        '''
        Returns the Hessian of the covariance function with respect to x and xp.

        Arguments
        --
        x: a 1-d array of abscissa
        xp: a 1-d array of alternative abscissa
        lth: the log of the hyperparameters
        '''
        th= np.exp(lth)
        return th[1]*np.ones(len(xp)), False

    
####  

class nnGP(gaussianprocess):
    '''
    Gaussian process with a neural network covariance function.
    '''
    noparams= 2
    description= 'neural network Gaussian process'

    def info(self):
        print('hparam[0] determines the initial value')
        print('hparam[1] determines the flexibility')
        print('hparam[2] determines the variance of the measurement error')
        
    def covfn(self, x, xp, lth):
        """
        Neural network covariance function from Williams.
        Returns the kernel function and Jacobian of the kernel function.

        Arguments
        --
        x: a 1-d array of abscissa
        xp: a 1-d array of alternative abscissa
        lth: the log of the hyperparameters
        """
        th= np.exp(lth)
        k= (np.arcsin(2*(th[0] + x*xp*th[1])/np.sqrt(1+2*(th[0]+x**2*th[1]))
                      /np.sqrt(1+2*(th[0]+xp**2*th[1]))))*2/np.pi
        jk= np.empty((len(xp), self.noparams))
        den= np.pi*(1+2*th[0]+2*th[1]*x**2)*(1+2*th[0]+2*th[1]*xp**2) \
          *np.sqrt(1+4*th[0]*(1+th[1]*(x-xp)**2)+2*th[1]*(x**2+xp**2))
        jk[:,0]= (4*(1+2*th[0]*(1+th[1]*(x-xp)**2) - 2*th[1]**2*x*(x-xp)**2*xp \
                     + 2*th[1]*(x**2-x*xp+xp**2)))/den*th[0]
        jk[:,1]= -(4*(2*th[0]**2*(x-xp)**2 - x*xp*(1+th[1]*(x**2+xp**2)) \
                      + th[0]*(-2*th[1]*x**3*xp+xp**2-2*x*xp*(2+th[1]*xp**2) \
                               +x**2*(1+4*th[1]*xp**2))))/den*th[1]
        return k, jk

    def d1covfn(self, x, xp, lth):
        '''
        Returns the d/dx of the covariance function.

        Arguments
        --
        x: a 1-d array of abscissa
        xp: a 1-d array of alternative abscissa
        lth: the log of the hyperparameters
        '''
        th= np.exp(lth)
        return 4*th[1]*(-2*th[0]*x+xp+2*th[0]*xp)/(1+2*th[0]+2*th[1]*x**2) \
          /np.sqrt(1+4*th[0]*(1+th[1]*(x-xp)**2)+2*th[1]*(x**2+xp**2))/np.pi, False

    def d1d2covfn(self, x, xp, lth):
        '''
        Returns the d/dx d/dxp of the covariance function.

        Arguments
        --
        x: a 1-d array of abscissa
        xp: a 1-d array of alternative abscissa
        lth: the log of the hyperparameters
        '''
        th= np.exp(lth)
        return 4*(th[1] + 4*th[0]*th[1])/(1 + 4*th[0]*(1+th[1]*(x-xp)**2) \
                                          + 2*th[1]*(x**2+xp**2))**1.5/np.pi, False

    def d12covfn(self, x, xp, lth):
        '''
        Returns d^2/dx^2 of the covariance function.

        Arguments
        --
        x: a 1-d array of abscissa
        xp: a 1-d array of alternative abscissa
        lth: the log of the hyperparameters
        '''
        th= np.exp(lth)
        return -8*th[1]*(8*th[0]**3 + th[0]**2*(6 - 8*th[1]*x*(x - 3*xp) - 16*th[1]**2*x*(x - xp)**3) + th[1]*x*xp*(3 + 6*th[1]*x**2 + 4*th[1]*xp**2) + th[0]*(1 - 2*th[1]*x*(x - 9*xp) - 8*th[1]**2*x*(x**3 - 3*x**2*xp + 3*x*xp**2 - 2*xp**3)))/(np.pi*(1 + 2*th[0] + 2*th[1]*x**2)**2*(1 + 4*th[0]*(1 + th[1]*(x - xp)**2) + 2*th[1]*(x**2 + xp**2))**1.5), False

    def d12d2covfn(self, x, xp, lth):
        '''
        Returns d^2/dx^2 d/dxp of the covariance function.

        Arguments
        --
        x: a 1-d array of abscissa
        xp: a 1-d array of alternative abscissa
        lth: the log of the hyperparameters
        '''
        th= np.exp(lth)
        return -24*(1 + 4*th[0])*th[1]**2*(x + 2*th[0]*x - 2*th[0]*xp)/(np.pi*(1 + 4*th[0]*(1 + th[1]*(x - xp)**2) + 2*th[1]*(x**2 + xp**2))**2.5), False

    def d12d22covfn(self, x, xp, lth):
        '''
        Returns d^2/dx^2 d^2/dxp^2 of the covariance function.

        Arguments
        --
        x: a 1-d array of abscissa
        xp: a 1-d array of alternative abscissa
        lth: the log of the hyperparameters
        '''
        th= np.exp(lth)
        return -48*(1 + 4*th[0])*th[1]**2*(4*th[0]**2*(-1 + 4*th[1]*(x - xp)**2) - 5*th[1]*x*xp + th[0]*(-1 + 4*th[1]*(2*x**2 - 5*x*xp + 2*xp**2)))/(np.pi*(1 + 4*th[0]*(1 + th[1]*(x - xp)**2) + 2*th[1]*(x**2 + xp**2))**3.5), False


    
####  


class sqexpGP(gaussianprocess):
    '''
    Gaussian process with a squared exponential covariance function.
    '''
    noparams= 2
    description= 'squared exponential Gaussian process'

    def info(self):
        print('hparam[0] determines the amplitude of variation')
        print('hparam[1] determines the flexibility')
        print('hparam[2] determines the variance of the measurement error')
        
    def covfn(self, x, xp, lth):
        '''
        Squared exponential covariance function.
        Returns the kernel function and Jacobian of the kernel function.

        Arguments
        --
        x: a 1-d array of abscissa
        xp: a 1-d array of alternative abscissa
        lth: the log of the hyperparameters
        '''
        th= np.exp(lth)
        xp= np.array(xp)
        e= np.exp(-th[1]/2.0*(x-xp)**2)
        k= th[0]*e 
        jk= np.empty((len(xp), self.noparams))
        jk[:,0]= e*th[0]
        jk[:,1]= -th[0]*th[1]*e/2.0*(x-xp)**2
        return k, jk

    def d1covfn(self, x, xp, lth):
        '''
        Returns d/dx of the covariance function.

        Arguments
        --
        x: a 1-d array of abscissa
        xp: a 1-d array of alternative abscissa
        lth: the log of the hyperparameters
        '''
        th= np.exp(lth)
        e= np.exp(-th[1]/2.0*(x-xp)**2)
        return -e*th[0]*th[1]*(x-xp), False

    def d1d2covfn(self, x, xp, lth):
        '''
        Returns d/dx d/dxp of the covariance function.

        Arguments
        --
        x: a 1-d array of abscissa
        xp: a 1-d array of alternative abscissa
        lth: the log of the hyperparameters
        '''
        th= np.exp(lth)
        e= np.exp(-th[1]/2.0*(x-xp)**2)
        return -e*th[0]*th[1]*(-1 + th[1]*(x-xp)**2), False

    def d12covfn(self, x, xp, lth):
        '''
        Returns d^2/dx^2 of the covariance function.

        Arguments
        --
        x: a 1-d array of abscissa
        xp: a 1-d array of alternative abscissa
        lth: the log of the hyperparameters
        '''
        th= np.exp(lth)
        e= np.exp(-th[1]/2.0*(x-xp)**2)
        return e*th[0]*th[1]*(-1 + th[1]*(x - xp)**2), False

    def d12d2covfn(self, x, xp, lth):
        '''
        Returns d^2/dx^2 d/dxp of the covariance function.

        Arguments
        --
        x: a 1-d array of abscissa
        xp: a 1-d array of alternative abscissa
        lth: the log of the hyperparameters
        '''
        th= np.exp(lth)
        e= np.exp(-th[1]/2.0*(x-xp)**2)
        return e*th[0]*th[1]**2*(-3 + th[1]*(x - xp)**2)*(x - xp), False

    def d12d22covfn(self, x, xp, lth):
        '''
        Returns d^2/dx^2 d^2/dxp^2 of the covariance function.

        Arguments
        --
        x: a 1-d array of abscissa
        xp: a 1-d array of alternative abscissa
        lth: the log of the hyperparameters
        '''
        th= np.exp(lth)
        e= np.exp(-th[1]/2.0*(x-xp)**2)
        return e*th[0]*th[1]**2*(3 - 6*th[1]*(x - xp)**2 + th[1]**2*(x - xp)**4), False

####  
    
    
class sqexplinGP(gaussianprocess):
    '''
    Gaussian process with a squared exponential covariance function with a linear trend.
    Returns the kernel function and Jacobian of the kernel function.
    
    Arguments
    --
    x: a 1-d array of abscissa
    xp: a 1-d array of alternative abscissa
    lth: the log of the hyperparameters
    '''
    noparams= 3
    description= 'squared exponential Gaussian process with a linear trend'

    def info(self):
        print('hparam[0] determines the amplitude of variation')
        print('hparam[1] determines the flexibility')
        print('hparam[2] determines the linear trend with increasing input')
        print('hparam[3] determines the variance of the measurement error')
    
    def covfn(self, x, xp, lth):
        '''
        Squared exponential covariance function with a linear trend.
        Returns the kernel function and Jacobian of the kernel function.

        Arguments
        --
        x: a 1-d array of abscissa
        xp: a 1-d array of alternative abscissa
        lth: the log of the hyperparameters
        '''
        th= np.exp(lth)
        xp= np.array(xp)
        e= np.exp(-th[1]/2.0*(x-xp)**2)
        k= th[0]*e + th[2]*x*xp
        jk= np.empty((len(xp), self.noparams))
        jk[:,0]= e*th[0]
        jk[:,1]= -th[0]*th[1]*e/2.0*(x-xp)**2
        jk[:,2]= x*xp*th[2]
        return k, jk


    def d1covfn(self, x, xp, lth):
        '''
        Returns d/dx of the covariance function.

        Arguments
        --
        x: a 1-d array of abscissa
        xp: a 1-d array of alternative abscissa
        lth: the log of the hyperparameters
        '''
        th= np.exp(lth)
        e= np.exp(-th[1]/2.0*(x-xp)**2)
        return -e*th[0]*th[1]*(x-xp) + th[2]*xp, False
    

    def d1d2covfn(self, x, xp, lth):
        '''
        Returns d/dx d/xp of the covariance function.

        Arguments
        --
        x: a 1-d array of abscissa
        xp: a 1-d array of alternative abscissa
        lth: the log of the hyperparameters
        '''
        th= np.exp(lth)
        e= np.exp(-th[1]/2.0*(x-xp)**2)
        return th[2] - e*th[0]*th[1]*(-1 + th[1]*(x-xp)**2), False

    def d12covfn(self, x, xp, lth):
        '''
        Returns d^2/dx^2 of the covariance function.

        Arguments
        --
        x: a 1-d array of abscissa
        xp: a 1-d array of alternative abscissa
        lth: the log of the hyperparameters
        '''
        th= np.exp(lth)
        e= np.exp(-th[1]/2.0*(x-xp)**2)
        return e*th[0]*th[1]*(-1 + th[1]*(x - xp)**2), False
    
    def d12d2covfn(self, x, xp, lth):
        '''
        Returns d^2/dx^2 d/dxp of the covariance function.

        Arguments
        --
        x: a 1-d array of abscissa
        xp: a 1-d array of alternative abscissa
        lth: the log of the hyperparameters
        '''
        th= np.exp(lth)
        e= np.exp(-th[1]/2.0*(x-xp)**2)
        return e*th[0]*th[1]**2*(-3 + th[1]*(x - xp)**2)*(x - xp), False

    def d12d22covfn(self, x, xp, lth):
        '''
        Returns d^2/dx^2 d^2/dxp^2 of the covariance function.

        Arguments
        --
        x: a 1-d array of abscissa
        xp: a 1-d array of alternative abscissa
        lth: the log of the hyperparameters
        '''
        th= np.exp(lth)
        e= np.exp(-th[1]/2.0*(x-xp)**2)
        return e*th[0]*th[1]**2*(3 - 6*th[1]*(x - xp)**2 + th[1]**2*(x - xp)**4), False
    
####  


class gaussianprocessException(Exception):
    __doc__= ''
    pass

