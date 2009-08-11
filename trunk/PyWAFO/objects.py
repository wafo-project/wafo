"""
Class objects for data, spectrum, covariance function and time series in WAFO

To represent data, spectra, covariance functions and time series
in WAFO, the classes WafoData, SpecData1D, CovData1D and TimeSeries are used.
Here follows a list of their attributes:

WafoData
--------
data : array_like
args : vector for 1D, list of vectors for 2D, 3D, ...
date : Date and time of creation or change.
labels : AxisLabels
children : list of WafoData objects

SpecData1D
----------
data : One sided Spectrum values
args : freguency values of freqtype
type :String: 'freq', 'dir', 'k2d', k1d', 'encdir' or 'enc'.
freqtype :'w' OR 'f' OR 'k' Frequency/wave number lag, length nf.
tr : Transformation function (default (none)).
h : Water depth (default inf).
norm : Normalization flag, Logical 1 if S is normalized, 0 if not
date : Date and time of creation or change.
v :  Ship speed, if .type = 'enc' or 'encdir'.
phi   angle of rotation of the coordinate system
           (counter-clocwise) e.g. azymuth of a ship.

CovData1D
---------
data : Covariance function values. Size [ny nx nt], all singleton dim. removed.
args : Lag of first space dimension, length nx.
h : Water depth.
tr : Transformation function.
type : 'enc', 'rot' or 'none'.
v : Ship speed, if .type='enc'
phi : Rotation of coordinate system, e.g.  direction of ship
norm : Normalization flag, Logical 1 if autocorrelation, 0 if covariance.
Rx, ... ,Rtttt :  Obvious derivatives of .R.
note : Memorandum string.
date : Date and time of creation or change.

"""
#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      pab
#
# Created:     16.09.2008
# Copyright:   (c) pab 2008
# Licence:     <your licence>
#-------------------------------------------------------------------------------
#!/usr/bin/env python


from __future__ import division
import warnings
import numpy as np
from scitools import numpytools as npt
from scitools.numpytools import (inf, pi, zeros, sqrt, where, log, exp, ceil, 
                                 floor, sin, arcsin, tanh, arctan2, newaxis, 
                                 linspace, sort)
from time import gmtime, strftime
from wafo.misc import (nextpow2, discretize, findtp, findtc, findextrema, 
        findcross, findrfc, sub_dict_select, tranproc, ecross)
from wafo.spectrum.dispersion_relation import w2k, k2w
from wafo.gaussian import Rind
#from scipy import fft
from scipy.integrate import simps, trapz
from pylab import stineman_interp
from matplotlib.mlab import psd
import scipy.signal
import scipy.interpolate as interpolate
from scipy.special import erf
from scipy.linalg import toeplitz

try:
    #import diffsumfunq
    import c_library
except:
    pass

from plotbackend import plotbackend

__all__ = ['SpecData1D','SpecData2D','WafoData', 'AxisLabels','CovData1D',
    'TimeSeries','TrGauss','LevelCrossings','CyclePairs','TurningPoints',
    'sensortypeid','sensortype']

def empty_copy(obj):
    class Empty(obj.__class__):
        def __init__(self):
            pass
    newcopy = Empty()
    newcopy.__class__ = obj.__class__
    return newcopy

def _set_seed(iseed):
    if iseed != None:
        try:
            np.random.set_state(iseed)
        except:
            np.random.seed(iseed)
def now():
    '''
    Return current date and time as a string
    '''
    return strftime("%a, %d %b %Y %H:%M:%S", gmtime())

class WafoData(object):
    '''
    Container class for data objects in WAFO

    Member variables
    ----------------
    data : array_like
    args : vector for 1D, list of vectors for 2D, 3D, ...
    labels : AxisLabels
    children : list of WafoData objects

    Member methods
    --------------
    plot :
    copy :


    Example
    -------
    >>> import numpy as np
    >>> x = np.arange(-2, 2, 0.2)

    # Plot 2 objects in one call
    >>> d2 = WafoData(np.sin(x), x, xlab='x', ylab='sin', title='sinus')
    >>> h = d2.plot()

    Plot with confidence interval
    d3 = wdata(sin(x),x)
    d3 = set(d3,'dataCI',[sin(x(:))*0.9 sin(x(:))*1.2])
    plot(d3)

    See also
    --------
    wdata/plot,
    specdata,
    covdata
    '''
    def __init__(self,data=None,args=None,**kwds):
        self.data = data
        self.args = args
        self.date = now()
        self.plotter = None
        self.children = None
        self.labels = AxisLabels(**kwds)
        self.setplotter()

    def plot(self,*args,**kwds):
        tmp = None
        if self.children!=None:
            plotbackend.hold('on')
            tmp = []
            for child in self.children:
                tmp1 = child.plot(*args,**kwds)
                if tmp1 !=None:
                    tmp.append(tmp1)
            if len(tmp)==0:
                tmp = None

        tmp2 =  self.plotter.plot(self,*args,**kwds)
        return tmp2,tmp
    
    def show(self):
        self.plotter.show()

    def copy(self):
        newcopy = empty_copy(self)
        newcopy.__dict__.update(self.__dict__)
        return newcopy


    def setplotter(self,plotmethod=None):
        '''
            Set plotter based on the data type data_1d, data_2d, data_3d or data_nd
        '''

        if isinstance(self.args,(list,tuple)): # Multidimensional data
            ndim = len(self.args)
            if ndim<2:
                warnings.warn('Unable to determine plotter-type, because len(self.args)<2.')
                print('If the data is 1D, then self.args should be a vector!')
                print('If the data is 2D, then length(self.args) should be 2.')
                print('If the data is 3D, then length(self.args) should be 3.')
                print('Unless you fix this, the plot methods will not work!')
            elif ndim==2:
                self.plotter = Plotter_2d(plotmethod)
            else:
                warnings.warn('Plotter method not implemented for ndim>2')

        else: #One dimensional data
            self.plotter = Plotter_1d(plotmethod)


class AxisLabels:
    def __init__(self,title='',xlab='',ylab='',zlab='',**kwds):
        self.title = title
        self.xlab = xlab
        self.ylab = ylab
        self.zlab = zlab
    def copy(self):
        newcopy = empty_copy(self)
        newcopy.__dict__.update(self.__dict__)
        return newcopy
        #lbkwds = self.labels.__dict__.copy()
        #labels = AxisLabels(**lbkwds)
        #return labels

    def labelfig(self):
        try:
            h1 = plotbackend.title(self.title)
            h2 = plotbackend.xlabel(self.xlab)
            h3 = plotbackend.ylabel(self.ylab)
            #h4 = plotbackend.zlabel(self.zlab)
            return h1,h2,h3
        except:
            pass

class Plotter_1d(object):
    """

    Parameters
    ----------
    plotmethod : string
        defining type of plot. Options are:
        bar : bar plot with rectangles
        barh : horizontal bar plot with rectangles
        loglog : plot with log scaling on the *x* and *y* axis
        semilogx :  plot with log scaling on the *x* axis
        semilogy :  plot with log scaling on the *y* axis
        plot : Plot lines and/or markers (default)
        stem : Stem plot
        step : stair-step plot
        scatter : scatter plot
    """

    def __init__(self,plotmethod='plot'):
        self.plotfun = None
        if plotmethod is None:
            plotmethod = 'plot'
        self.plotbackend = plotbackend
        try:
            #self.plotfun = plotbackend.__dict__[plotmethod]
            self.plotfun = getattr(plotbackend, plotmethod)
        except:
            pass
    def show(self):
        plotbackend.show()

    def plot(self,wdata,*args,**kwds):
        if isinstance(wdata.args,(list,tuple)):
            args1 = tuple((wdata.args))+(wdata.data,)+args
        else:
            args1 = tuple((wdata.args,))+(wdata.data,)+args
        h1 = self.plotfun(*args1,**kwds)
        h2 = wdata.labels.labelfig()
        return h1,h2

class Plotter_2d(Plotter_1d):
    """
    Parameters
    ----------
    plotmethod : string
        defining type of plot. Options are:
        contour (default)
        mesh
        surf


    """

    def __init__(self,plotmethod='contour'):
        if plotmethod is None:
            plotmethod = 'contour'
        super(Plotter_2d,self).__init__(plotmethod)
        #self.plotfun = plotbackend.__dict__[plotmethod]

class TrGauss(WafoData):
    """
    Non-parametric Transformation model, g.

    Information about the moments of the process can be obtained by site
    specific data, laboratory measurements or by resort to theoretical models.

    Assumption
    ----------
    The Gaussian process, Y, distributed N(0,1) is related to the
    non-Gaussian process, X, by Y = g(X), where g is a stricly increasing 
    function.

    Methods
    -------
    dist2gauss : Returns a measure of departure from the Gaussian model, i.e.,
                int (g(x)-xn)^2 dx  where int. limits is given by X.
    dat2gauss : Transform non-linear data to Gaussian scale
    gauss2dat : Transform Gaussian data to non-linear scale

    Member variables
    ----------------
    data : array-like
        Gaussian values, Y
    args : array-like
        non-Gaussian values, X
    mean, sigma : real, scalar
        mean and standard-deviation, respectively, of the non-Gaussian process. 
        Default: 
        mean = self.gauss2dat(0), 
        sigma = (self.gauss2dat(1)-self.gauss2dat(-1))/2
    
    Example
    -------
    Construct a linear transformation
    >>> import numpy as np
    >>> x = np.linspace(-5,5); y = x
    >>> g = TrGauss(x,y)
    >>> g.mean
    array([ 0.])
    >>> g.sigma
    array([ 1.])
    
    Check that the departure from a Gaussian model is zero
    >>> g.dist2gauss()<1e-16
    True
    """
    def __init__(self,*args,**kwds):
        super(TrGauss, self).__init__(*args,**kwds)
        self.labels.title = 'Transform'
        self.labels.ylab = 'g(u)'
        self.labels.xlab = 'u'
        self.mean = kwds.get('mean', None)
        self.sigma = kwds.get('sigma', None)
        if self.mean is None: 
            self.mean = self.gauss2dat(0)
        if self.sigma is None:
            self.sigma = (self.gauss2dat(1)-self.gauss2dat(-1))/2.

    def  dist2gauss(self, x=None, xnmin=-5, xnmax=5, n=513):
        """
        Return a measure of departure from the Gaussian model.

        Parameters
        ----------
        x : vector  (default sigma*linspace(xnmin,xnmax,n)+mean)
        xnmin : real, scalar
            minimum on normalized scale
        xnmax : real, scalar
            maximum on normalized scale
        n : integer, scalar
            number of evaluation points

        Returns
        -------
        t0 : real, scalar
            a measure of departure from the Gaussian model calculated as
            trapz(xn,(xn-g(x))**2.) where int. limits is given by X.
        """
        if x is None:
            xn = np.linspace(xnmin, xnmax, n)
            x = self.sigma*xn + self.mean
        else:    
            xn = (x-self.mean)/self.sigma

        y = self.dat2gauss(x)
        t0 = trapz(xn, (xn-y)**2.)
        return t0

    def gauss2dat(self, y, *yi):
        """
        Transforms Gaussian data, y, to non-linear scale.

        Parameters
        ----------
        y, y1,...,yn : array-like
            input vectors with Gaussian data values, where yi is the i'th time
            derivative of y. (n<=4)
        Returns
        -------
        x, x1,...,xn : array-like
            transformed data to a non-linear scale


        See also
        --------
        dat2gauss
        tranproc
        """
        return tranproc(self.data, self.args, y, *yi)
    
    def dat2gauss(self, x, *xi):
        """
        Transforms non-linear data, x, to Gaussian scale.

        Parameters
        ----------
        x, x1,...,xn : array-like
            input vectors with non-linear data values, where xi is the i'th time
            derivative of x. (n<=4)
        Returns
        -------
        y, y1,...,yn : array-like
            transformed data to a Gaussian scale

        See also
        --------
        gauss2dat
        tranproc.
        """
        return tranproc(self.args, self.data, x, *xi)

class LevelCrossings(WafoData):
    '''
    Container class for Level crossing data objects in WAFO

    Member variables
    ----------------
    data : array_like
    args : vector for 1D

    '''
    def __init__(self,*args,**kwds):
        super(LevelCrossings, self).__init__(*args,**kwds)
        self.labels.title = 'Level crossing spectrum'
        self.labels.xlab = 'Levels'
        self.labels.ylab = 'Count'
        self.stdev = kwds.get('stdev',None)
        self.mean = kwds.get('mean',None)
        self.setplotter(plotmethod='step')

        icmax = self.data.argmax()
        if self.data != None:
            if self.stdev is None or self.mean is None:
                logcros = where(self.data==0.0, inf, -log(self.data))
                logcmin = logcros[icmax]
                logcros = sqrt(2*abs(logcros-logcmin))
                logcros[0:icmax+1] = 2*logcros[icmax]-logcros[0:icmax+1]
                p = np.polyfit(self.args[10:-9], logcros[10:-9],1) #least square fit
                if self.stdev is None:
                    self.stdev = 1.0/p[0] #estimated standard deviation of x
                if self.mean is None:
                    self.mean = -p[1]/p[0] #self.args[icmax]
            cmax = self.data[icmax]
            x = (self.args-self.mean)/self.stdev
            y = cmax*exp(-x**2/2.0)
            self.children = [WafoData(y,self.args)]

    def sim(self,ns,alpha):
        """
        Simulates process with given irregularity factor and crossing spectrum

        Parameters
        ----------
        ns : scalar, integer
            number of sample points.
        alpha : real scalar
            irregularity factor, 0<alpha<1, small  alpha  gives
            irregular process.

        Returns
        --------
        ts : timeseries object
            with times and values of the simulated process.

        Example
        -------
        n = 10000
        S = jonswap(7)
        alpha = spec2char(S,'alpha')
        xs  = spec2sdat(S,n)
        lc  = dat2lc(xs)
        xs2 = lc2sdat(lc,n,alpha)
        Se  = dat2spec(xs2)
        plotspec(S),hold on
        plotspec(Se,'r'), hold off
        spec2char(Se,'alpha')
        lc2  = dat2lc(xs2)
        figure(gcf+1)
        subplot(211)
        lcplot(lc2)
        subplot(212)
        lcplot(lc)
        """

        # TODO % add a good example
        f = np.linspace(0,0.49999,1000)
        rho_st = 2.*sin(f*pi)**2-1.
        tmp = alpha*arcsin(sqrt((1.+rho_st)/2))
        tmp = sin(tmp)**2
        a2  = (tmp-rho_st)/(1-tmp)
        y   = np.vstack((a2+rho_st,1-a2)).min(axis=0)
        maxidx = y.argmax()
        #[maximum,maxidx]=max(y)

        rho_st = rho_st[maxidx]
        a2 = a2[maxidx]
        a1 = 2.*rho_st+a2-1.
        r0 = 1.
        r1 = -a1/(1.+a2)
        r2 = (a1**2-a2-a2**2)/(1+a2)
        sigma2 = r0+a1*r1+a2*r2
        randn = np.random.randn
        e = randn(ns)*sqrt(sigma2)
        e[:1] = 0.0
        L0 = randn(1)
        L0 = np.vstack((L0,r1*L0+sqrt(1-r2**2)*randn(1)))
        #%Simulate the process, starting in L0
        filter = scipy.signal.lfilter
        L = filter(1,[1, a1, a2],e,filter([1, a1, a2],1,L0))

        epsilon = 1.01
        min_L   = min(L)
        max_L   = max(L)
        maxi    = max(abs(np.r_[min_L, max_L]))*epsilon
        mini    = -maxi

        u = np.linspace(mini,maxi,101)
        G = (1+erf(u/sqrt(2)))/2
        G = G*(1-G)

        x = np.linspace(0,r1,100)
        factor1 = 1./sqrt(1-x**2)
        factor2 = 1./(1+x)
        integral = zeros(u.shape, dtype=float)
        for i in range(len(integral)):
            y = factor1*exp(-u[i]*u[i]*factor2)
            integral[i] = trapz(x,y)
        #end
        G = G-integral/(2*pi)
        G = G/max(G)

        Z = ((u>=0)*2-1)*sqrt(-2*log(G))

##        sumcr   = trapz(lc(:,1),lc(:,2))
##        lc(:,2) = lc(:,2)/sumcr
##        mcr     = trapz(lc(:,1),lc(:,1).*lc(:,2))
##        scr     = trapz(lc(:,1),lc(:,1).^2.*lc(:,2))
##        scr     = sqrt(scr-mcr^2)
##        g       = lc2tr(lc,mcr,scr)
##
##        f = [u u]
##        f(:,2) = tranproc(Z,fliplr(g))
##
##        process = tranproc(L,f)
##        process = [(1:length(process)) process]
##
##
##        %Check the result without reference to getrfc:
##        LCe = dat2lc(process)
##        max(lc(:,2))
##        max(LCe(:,2))
##
##        clf
##        plot(lc(:,1),lc(:,2)/max(lc(:,2)))
##        hold on
##        plot(LCe(:,1),LCe(:,2)/max(LCe(:,2)),'-.')
##        title('Relative crossing intensity')
##
##        %% Plot made by the function funplot_4, JE 970707
##        %param = [min(process(:,2)) max(process(:,2)) 100]
##        %plot(lc(:,1),lc(:,2)/max(lc(:,2)))
##        %hold on
##        %plot(levels(param),mu/max(mu),'--')
##        %hold off
##        %title('Crossing intensity')
##        %watstamp
##
##        % Temporarily
##        %funplot_4(lc,param,mu)


    def trgauss(self):
        pass

class CyclePairs(WafoData):
    '''
    Container class for Cycle Pairs data objects in WAFO

    Member variables
    ----------------
    data : array_like
    args : vector for 1D


    '''
    def __init__(self,*args,**kwds):
        super(CyclePairs, self).__init__(*args,**kwds)
        self.type_ = kwds.get('type_','max2min')
        self.stdev = kwds.get('stdev',None)
        self.mean = kwds.get('mean',None)

        self.labels.title = self.type_+ ' cycle pairs'
        self.labels.xlab = 'min'
        self.labels.ylab = 'max'

    def amplitudes(self):
        return (self.data-self.args)/2.
    
    def damage(self, beta, K=1):
        """
        Calculates the total Palmgren-Miner damage of cycle pairs.

        Parameters
        ----------
        beta : array-like, size m
            Beta-values, material parameter.                   
        K : scalar, optional
            K-value, material parameter.

        Returns
        -------
        D : ndarray, size m
            Damage.      

        Notes
        -----
        The damage is calculated according to
           D[i] = sum ( K * a**beta[i] ),  with  a = (max-min)/2

        Examples
        --------
        >>> import wafo
        >>> from matplotlib import pyplot as plt
        >>> ts = wafo.objects.mat2timeseries(wafo.data.sea())
        >>> tp = ts.turning_points()
        >>> mm = tp.cycle_pairs()
        >>> h = mm.plot('.')
        >>> bv = range(3,9)
        >>> D = mm.damage(beta=bv)
        >>> D
        array([ 138.5238799 ,  117.56050788,  108.99265423,  107.86681126,
                112.3791076 ,  122.08375071])
        >>> h = plt.plot(bv,D,'x-')

        See also
        --------
        SurvivalCycleCount
        """
        amp = np.abs(self.amplitudes())
        return np.atleast_1d([K*np.sum(amp**betai) for betai in beta])

    def level_crossings(self,type_='uM'):
        """ Return number of upcrossings from a cycle count.

        Parameters
        ----------
        type_ : int or string
            defining crossing type, options are
            0,'u'  : only upcrossings.
            1,'uM' : upcrossings and maxima (default).
            2,'umM': upcrossings, minima, and maxima.
            3,'um' :upcrossings and minima.
        Return
        ------
        lc : level crossing object
            with levels and number of upcrossings.


        Calculates the number of upcrossings from a cycle pairs, e.g.
        min2Max cycles or rainflow cycles.

        Example:
        --------
        >>> import wafo
        >>> ts = wafo.objects.mat2timeseries(wafo.data.sea())
        >>> tp = ts.turning_points()
        >>> mm = tp.cycle_pairs()
        >>> h = mm.plot('.')
        >>> lc = mm.level_crossings()
        >>> h2 = lc.plot()

        See also
        --------
        TurningPoints
        LevelCrossings
        """

        if isinstance(type_,str):
            t = dict(u=0,uM=1,umM=2,um=3)
            defnr = t.get(type_,1)
        else:
            defnr = type_

        if ((defnr<0) or (defnr>3)):
            raise ValueError('type_ must be one of (1,2,3,4).')

        index,=np.nonzero(self.args <= self.data)
        if index.size==0:
            index,=np.nonzero(self.args >= self.data)
            M = self.args[index]
            m = self.data[index]
        else:
            m = self.args[index]
            M = self.data[index]

#if isempty(index)
#  error('Error in input cc.')
#end
        ncc = len(m)
        ones = np.ones
        zeros = np.zeros
        cumsum = np.cumsum
        minima = np.vstack((m, ones(ncc), zeros(ncc), ones(ncc)))
        maxima = np.vstack((M, -ones(ncc), ones(ncc), zeros(ncc)))

        extremes = np.hstack((maxima, minima))
        index = extremes[0].argsort()
        extremes = extremes[:,index]

        ii = 0
        n = extremes.shape[1]
        extr = zeros((4,n))
        extr[:,0] = extremes[:,0]
        for i in xrange(1,n):
            if extremes[0,i]==extr[0,ii]:
                extr[1:4,ii]=extr[1:4,ii]+extremes[1:4,i]
            else:
                ii += 1
                extr[:,ii] = extremes[:,i]

        #[xx nx]=max(extr(:,1))
        nx = extr[0].argmax()+1
        levels = extr[0,0:nx]
        if defnr == 2: #% This are upcrossings + maxima
            dcount = cumsum(extr[1,0:nx]) + extr[2,0:nx]-extr[3,0:nx]
        elif defnr == 4: # % This are upcrossings + minima
            dcount = cumsum(extr[1,0:nx])
            dcount[nx-1] = dcount[nx-2]
        elif defnr == 1: #% This are only upcrossings
            dcount = cumsum(extr[1,0:nx]) - extr[3,0:nx]
        elif defnr == 3: #% This are upcrossings + minima + maxima
            dcount = cumsum(extr[1,0:nx]) + extr[2,0:nx]
        return LevelCrossings(dcount,levels,stdev=self.stdev)

def qtf(w, h=inf, g=9.81):
    """
    Return Quadratic Transfer Function

    Parameters
    ------------
    w : array-like
        angular frequencies
    h : scalar
        water depth
    g : scalar
        acceleration of gravity

    Returns
    -------
    Hs   = sum frequency effects
    Hd   = difference frequency effects
    Hdii = diagonal of Hd
    """
    w = np.atleast_1d(w)
    Nw = w.size

    kw = w2k(w,theta=0,h=h,g=g)
    kw = kw[0]
    [k1, k2] = npt.meshgrid(kw,kw)

    if h==inf: # go here for faster calculations
        Hs   = 0.25*(abs(k1)+abs(k2))
        Hd   = -0.25*abs(abs(k1)-abs(k2))
        Hdii = np.zeros(Nw)
        return Hs, Hd ,Hdii

    [w1, w2]= npt.meshgrid(w,w)



    w12  = (w1*w2)
    w1p2 = (w1+w2)
    w1m2 = (w1-w2)
    k12  = (k1*k2)
    k1p2 = (k1+k2)
    k1m2 = abs(k1-k2)
    cosh = np.cosh
    sinh = np.sinh
    if 0: # Langley
        p1 = (-2*w1p2*(k12*g**2.-w12**2.)+
            w1*(w2**4.-g**2*k2**2)+w2*(w1**4-g*2.*k1**2))/(4.*w12)
        p2= w1p2**2.*cosh((k1p2)*h)-g*(k1p2)*sinh((k1p2)*h)

        Hs = (-p1/p2*w1p2*cosh((k1p2)*h)/g-
            (k12*g**2-w12**2.)/(4*g*w12)+(w1**2+w2**2)/(4*g))

        p3 = (-2*w1m2*(k12*g**2+w12**2)-
            w1*(w2**4-g**2*k2**2)+w2*(w1**4-g**2*k1**2))/(4.*w12)
        p4= w1m2**2.*cosh(k1m2*h)-g*(k1m2)*sinh((k1m2)*h)


        Hd = (-p3/p4*(w1m2)*cosh((k1m2)*h)/g-
            (k12*g**2+w12**2)/(4*g*w12)+(w1**2.+w2**2.)/(4.*g))

    else: #  % Marthinsen & Winterstein
        tmp1 = 0.5*g*k12/w12
        tmp2 = 0.25/g*(w1**2.+w2**2.+w12)
        Hs   = (tmp1-tmp2+0.25*g*(w1*k2**2.+w2*k1**2)/
                (w12*(w1p2)))/(1.-g*(k1p2)/(w1p2)**2.*tanh((k1p2)*h))+tmp2-0.5*tmp1 #% OK

        tmp2 = 0.25/g*(w1**2+w2**2-w12) # %OK
        Hd   = (tmp1-tmp2-0.25*g*(w1*k2**2-w2*k1**2)/
            (w12*(w1m2)))/(1.-g*(k1m2)/(w1m2)**2.*tanh((k1m2)*h))+tmp2-0.5*tmp1 # % OK


    #%tmp1 = 0.5*g*kw./(w.*sqrt(g*h))
    #%tmp2 = 0.25*w.^2/g


    Cg = 0.5*g*(tanh(kw*h) +kw*h*(1.0- tanh(kw*h)**2))/w # %Wave group velocity
    Hdii = (0.5*(0.5*g*(kw/w)**2.-0.5*w**2/g+g*kw/(w*Cg))
            /(1.-g*h/Cg**2.)-0.5*kw/sinh(2*kw*h) )# % OK
    Hd.flat[0::Nw+1] = Hdii

    #%k    = find(w1==w2)
    #%Hd(k) = Hdii

    #% The NaN's occur due to division by zero. => Set the isnans to zero
    isnan = np.isnan
    Hdii = where(isnan(Hdii),0,Hdii)
    Hd = where(isnan(Hd),0,Hd)
    Hs = where(isnan(Hs),0,Hs)

    return Hs, Hd ,Hdii

class SpecData1D(WafoData):
    """ Container class for 1D spectrum data objects in WAFO

    Member variables
    ----------------
    data : array_like
    args : vector for 1D, list of vectors for 2D, 3D, ...

    type : string
        spectrum type, one of 'freq', 'k1d', 'enc' (default 'freq')
    freqtype : letter
        frequency type, one of: 'f', 'w' or 'k' (default 'w')


    Examples
    --------
    >>> import numpy as np
    >>> import wafo.spectrum.models as sm
    >>> Sj = sm.Jonswap(Hm0=3)
    >>> w = np.linspace(0,4,256)
    >>> S1 = Sj.toSpecData(w)   #Make spectrum object from numerical values
    >>> S = SpecData1D(Sj(w),w) # Alternatively do it manually

    See also
    --------
    WafoData
    CovData
    """

    def __init__(self,*args,**kwds):
        super(SpecData1D, self).__init__(*args,**kwds)
        self.name = 'WAFO Spectrum Object'
        self.type = 'freq'
        self.freqtype = 'w'
        self.angletype = ''
        self.h = np.inf
        self.tr = None
        self.phi = 0.0
        self.v = 0.0
        self.norm = False
        somekeys = ['angletype', 'phi', 'name', 'h', 'tr', 'freqtype', 'v', 'type', 'norm']

        self.__dict__.update(sub_dict_select(kwds,somekeys))

        self.setlabels()

    def tocov_matrix(self, nr=0, nt=None, dt=None):
        '''
        Computes covariance function and its derivatives, alternative version

        Parameters
        ----------
        nr : scalar integer
            number of derivatives in output, nr<=4          (default 0)
        nt : scalar integer
            number in time grid, i.e., number of time-lags.
            (default rate*(n-1)) where rate = round(1/(2*f(end)*dt)) or
                     rate = round(pi/(w(n)*dt)) depending on S.
        dt : real scalar
            time spacing for R

        Returns
        -------
        R : [R0, R1,...Rnr], shape Nt+1 x Nr+1 
            matrix with autocovariance and its derivatives, i.e., Ri (i=1:nr) 
            are column vectors with the 1'st to nr'th derivatives of R0.  

        NB! This routine requires that the spectrum grid is equidistant
           starting from zero frequency.

        Example
        -------
        >>> import wafo.spectrum.models as sm
        >>> Sj = sm.Jonswap()
        >>> S = Sj.toSpecData()
        >>> R = S.tocov_matrix(nr=3, nt=256, dt=0.1)
        >>> R[:2,:]
        array([[ 3.06075987,  0.        , -1.67750289,  0.        ],
               [ 3.05246132, -0.16662376, -1.66819445,  0.18634189]])

        See also
        --------
        cov,
        resample,
        objects
        '''

        ftype = self.freqtype # %options are 'f' and 'w' and 'k'
        freq  = self.args
        n     = len(freq)
        dt_old = self.sampling_period()
        if dt is None:
            dt = dt_old
            rate = 1
        else:
            rate = np.maximum(np.round(dt_old*1./dt),1.)


        if nt is None:
            nt = rate*(n-1)
        else: #%check if Nt is ok
            nt = np.minimum(nt,rate*(n-1))


        checkdt = 1.2*min(np.diff(freq))/2./np.pi
        if ftype in 'k':
            lagtype = 'x'
        else:
            lagtype = 't'
            if ftype in 'f':
                checkdt = checkdt*2*np.pi
        msg1 = 'The step dt = %g in computation of the density is too small.' % dt
        msg2 = 'The step dt = %g step is small, may cause numerical inaccuracies.' % dt

        if (checkdt < 2.**-16/dt):
            np.disp(msg1)
            np.disp('The computed covariance (by FFT(2^K)) may differ from the theoretical.')
            np.disp('Solution:')
            raise ValueError('use larger dt or sparser grid for spectrum.')


        #% Calculating covariances
        #%~~~~~~~~~~~~~~~~~~~~~~~~
        S2 = self.copy()
        S2.resample(dt)

        R2 = S2.tocov(nr, nt, rate=1)
        R = np.zeros((nt+1,nr+1),dtype=float)
        R[:,0] = R2.data[0:nt+1]
        fieldname = 'R' + lagtype*nr
        for ix in range(1,nr+1):
            fn = fieldname[:ix+1]
            Ri = getattr(R2,fn)
            R[:,ix] = Ri[0:nt+1]


        EPS0 = 0.0001
        cc  = R[0,0]-R[1,0]*(R[1,0]/R[0,0])
        if nt+1>=5:
            #%cc1=R(1,1)-R(3,1)*(R(3,1)/R(1,1))+R(3,2)*(R(3,2)/R(1,3))
            #%cc3=R(1,1)-R(5,1)*(R(5,1)/R(1,1))+R(5,2)*(R(5,2)/R(1,3))

            cc2 = R[0,0]-R[4, 0]*(R[4, 0]/R[0, 0])
            if (cc2<EPS0):
                warnings.warn(msg1)

        if (cc<EPS0):
            np.disp(msg2)
        return R

    def tocov(self, nr=0, nt=None, rate=None):
        '''
        Computes covariance function and its derivatives

        Parameters
        ----------
        nr : number of derivatives in output, nr<=4 (default = 0).
        nt : number in time grid, i.e., number of time-lags
              (default rate*(length(S.S)-1)).
        rate = 1,2,4,8...2**r, interpolation rate for R
               (default = 1, no interpolation)

        Returns
        -------
        R : CovData1D
            auto covariance function

        The input 'rate' gives together with the spectrum
        the t-grid-spacing: dt=pi/(S.w(end)*rate), S.w(end) is the Nyquist freq.
        This results in the t-grid: 0:dt:Nt*dt.

        What output is achieved with different S and choices of Nt,Nx and Ny:
        1) S.type='freq' or 'dir', Nt set, Nx,Ny not set: then result R(t) (one-dim)
        2) S.type='k1d' or 'k2d', Nt set, Nx,Ny not set: then result R(x) (one-dim)
        3) Any type, Nt and Nx set =>R(x,t); Nt and Ny set =>R(y,t)
        4) Any type, Nt, Nx and Ny set => R(x,y,t)
        5) Any type, Nt not set, Nx and/or Ny set => Nt set to default, goto 3) or 4)

        NB! This routine requires that the spectrum grid is equidistant
         starting from zero frequency.
        NB! If you are using a model spectrum, S, with sharp edges
         to calculate covariances then you should probably round off the sharp
         edges like this:

        Example:
        >>> import wafo.spectrum.models as sm
        >>> Sj = sm.Jonswap()
        >>> S = Sj.toSpecData()
        >>> S.data[0:40] = 0.0
        >>> S.data[100:-1] = 0.0
        >>> Nt = len(S.data)-1
        >>> R = S.tocov(nr=0,nt=Nt)

        R   = spec2cov(S,0,Nt)
        win = parzen(2*Nt+1)
        R.R = R.R.*win(Nt+1:end)
        S1  = cov2spec(R)
        R2  = spec2cov(S1)
        figure(1)
        plotspec(S),hold on, plotspec(S1,'r')
        figure(2)
        covplot(R), hold on, covplot(R2,[],[],'r')
        figure(3)
        semilogy(abs(R2.R-R.R)), hold on,
        semilogy(abs(S1.S-S.S)+1e-7,'r')

        See also
        --------
        cov2spec
        '''

        freq = self.args
        n = len(freq)

        if freq[0]>0:
            raise ValueError('Spectrum does not start at zero frequency/wave number.\n Correct it with resample, for example.')
        dw = np.abs(np.diff(freq,n=2,axis=0))
        if np.any(dw>1.0e-8):
            raise ValueError('Not equidistant frequencies/wave numbers in spectrum.\n Correct it with resample, for example.')


        if rate is None:
            rate=1 # %interpolation rate
        elif rate>16:
            rate = 16
        else: # make sure rate is a power of 2
            rate = 2**nextpow2(rate)

        if nt is None:
            nt = rate*(n-1)
        else: #check if Nt is ok
            nt = np.minimum(nt, rate*(n-1))

        S = self.copy()

        if self.freqtype in 'k':
            lagtype = 'x'
        else:
            lagtype = 't'

        dT = S.sampling_period()
        #normalize spec so that sum(specn)/(n-1)=R(0)=var(X)
        specn = S.data*freq[-1]
        if S.freqtype in 'f':
            w = freq*2*np.pi
        else:
            w = freq

        nfft = rate*2**nextpow2(2*n-2)

        Rper = np.r_[specn, np.zeros(nfft-(2*n)+2) , np.conj(specn[n-1:0:-1])] # % periodogram
        t    = np.r_[0:nt+1]*dT*(2*n-2)/nfft

        fft = np.fft.fft

        r   = fft(Rper,nfft).real/(2*n-2)
        R = CovData1D(r[0:nt+1],t,lagtype=lagtype)
        R.tr   = S.tr
        R.h    = S.h
        R.norm = S.norm

        if nr>0:
            w = np.r_[w , np.zeros(nfft-2*n+2) ,-w[n-1:0:-1] ]
            fieldname = 'R' + lagtype*nr 
            for ix in range(1,nr+1):
                Rper = -1j*w*Rper
                r    = fft(Rper,nfft).real/(2*n-2)
                setattr(R,fieldname[0:ix+1],r[0:nt+1])
        return R

    def to_t_pdf(self, u=None, pdef='Tc', paramt=None, **options):
        '''
        Density of crest/trough- period or length, version 2. 
       
        Parameters
        ----------
        u : real scalar
            reference level (default the most frequently crossed level).
        pdef : string, 'Tc', Tt', 'Lc' or 'Lt'
            'Tc',    gives half wave period, Tc (default).
            'Tt',    gives half wave period, Tt
            'Lc' and 'Lt' ditto for wave length.
        paramt : [t0, tn, nt] 
            where t0, tn and nt is the first value, last value and the number
            of points, respectively, for which the density will be computed. 
            paramt= [5, 5, 51] implies that the density is computed only for 
            T=5 and using 51 equidistant points in the interval [0,5].
        options : optional parameters
            controlling the performance of the integration. See Rind for details.

        Notes
        -----
        SPEC2TPDF2 calculates pdf of halfperiods  Tc, Tt, Lc or Lt 
        in a stationary Gaussian transform process X(t), 
        where Y(t) = g(X(t)) (Y zero-mean Gaussian with spectrum given in S). 
        The transformation, g, can be estimated using LC2TR,
        DAT2TR, HERMITETR or OCHITR.  
    
        Example
        -------
        The density of Tc is computed by:
        >>> import pylab as plb
        >>> from wafo.spectrum import models as sm
        >>> w = np.linspace(0,3,100)
        >>> Sj = sm.Jonswap()
        >>> S = Sj.toSpecData()
        >>> f = S.to_t_pdf(pdef='Tc', paramt=(0, 10, 51), speed=7) 
        >>> h = f.plot()
     
        estimated error bounds
        >>> h2 = plb.plot(f.args, f.data+f.err, 'r', f.args, f.data-f.err, 'r')  
    
        >>> plb.close('all')
    
        See also  
        --------
        Rind, spec2cov2, specnorm, dat2tr, dat2gaus, perioddef, wavedef

        '''

        opts = dict(speed=9)
        opts.update(options)
        if pdef[0] in ('l','L'):
            if self.type!='k1d':
                raise ValueError('Must be spectrum of type: k1d')
        elif pdef[0] in ('t','T'):
            if self.type!='freq':
                raise ValueError('Must be spectrum of type: freq')
        else:
            raise ValueError('pdef must be Tc,Tt or Lc, Lt')
#        if strncmpi('l',def,1)
#          spec=spec2spec(spec,'k1d')
#        elseif strncmpi('t',def,1)
#          spec=spec2spec(spec,'freq')
#        else
#          error('Unknown def')
#        end
        pdef2defnr = dict(tc=1, lc=1, tt=-1, lt=-1)
        defnr = pdef2defnr[pdef.lower()]
         
        S = self.copy()
        S.normalize()
        m, mtxt = self.moment(nr=2, even=True)
        A = sqrt(m[0]/m[1])
        
       
        if self.tr is None:
            y = np.linspace(-5,5,513)
            g = TrGauss(y,sqrt(m[0])*y)
        else:
            g = self.tr
        
        
        if u is None:
            u = g.gauss2dat(0) #% most frequently crossed level 
        
        # transform reference level into Gaussian level
        un = g.dat2gauss(u)

        #disp(['The level u for Gaussian process = ', num2str(u)])
        
        if paramt is None:
            #% z2 = u^2/2
            z  = -np.sign(defnr)*un/sqrt(2)
            expectedMaxPeriod = 2*ceil(2*pi*A*exp(z)*(0.5+erf(z)/2)) 
            paramt = [0, expectedMaxPeriod, 51]
        
        t0     = paramt[0]
        tn     = paramt[1]
        Ntime  = paramt[2]
        t      = np.linspace(0, tn/A, Ntime) #normalized times
        Nstart = max(round(t0/tn*(Ntime-1)),1)  #% index to starting point to
                                             #% evaluate
                                    
        dt = t[1] - t[0]
        nr = 2
        R = S.tocov_matrix(nr,Ntime-1,dt)
        #R  = spec2cov2(S,nr,Ntime-1,dt)
        
                        
        xc   = np.vstack((un, un))
        indI = -np.ones(4, dtype=int)
        Nd   = 2
        Nc   = 2
        XdInf = 100.e0*sqrt(-R[0,2])
        XtInf = 100.e0*sqrt(R[0,0])
        
        B_up  = np.hstack([un+XtInf, XdInf, 0])
        B_lo  = np.hstack([un,    0, -XdInf])
        #%INFIN = [1 1 0]
        #BIG   = zeros((Ntime+2,Ntime+2))
        ex    = zeros(Ntime+2, dtype=float)
        #%CC    = 2*pi*sqrt(-R(1,1)/R(1,3))*exp(un^2/(2*R(1,1)))
        #%  XcScale = log(CC)
        opts['xcscale'] = log(2*pi*sqrt(-R[0, 0]/R[0, 2]))+(un**2/(2*R[0, 0]))
        
        f = zeros(Ntime, dtype=float)
        err = zeros(Ntime, dtype=float)
        
        rind = Rind(**opts)
        #h11 = fwaitbar(0,[],sprintf('Please wait ...(start at: %s)',datestr(now)))
        for pt in xrange(Nstart,Ntime):
            Nt = pt - Nd + 1
            Ntd = Nt + Nd
            Ntdc = Ntd + Nc
            indI[1] = Nt-1
            indI[2] = Nt
            indI[3] = Ntd-1
            
            #% positive wave period  
            BIG = self._covinput(pt,R) 
          
            tmp = rind(BIG, ex[:Ntdc], B_lo, B_up, indI, xc, Nt)
            f[pt], err[pt] = tmp[:2]
            #fwaitbar(pt/Ntime,h11,sprintf('%s Ready: %d of %d',datestr(now),pt,Ntime))
        #end
        #close(h11)
        
        
        titledict = dict(tc='Density of Tc',tt='Density of Tt',lc='Density of Lc',lt='Density of Lt')
        Htxt = titledict.get(pdef.lower())
        
        if pdef[0].lower()=='l':
            xtxt = 'wave length [m]'
        else:
            xtxt = 'period [s]'
        
        Htxt = '%s_{v =%2.5g}' % (Htxt, u)
        pdf = WafoData(f/A, t*A, title=Htxt, xlab=xtxt)
        pdf.err = err/A
        pdf.u = u
        pdf.options = opts
        return pdf


    def _covinput(self,pt,R):
        """
        Return covariance matrix for Tc or Tt period problems
    
        Parameters
        ----------
        pt : scalar integer
            time
        R : array-like, shape Ntime x 3
            [R0,R1,R2] column vectors with autocovariance and its derivatives, 
            i.e., Ri (i=1:2) are vectors with the 1'st and 2'nd derivatives of R0.  
      
        The order of the variables in the covariance matrix are organized as follows: 
        For pt>1:
        ||X(t2)..X(ts),..X(tn-1)|| X'(t1) X'(tn)|| X(t1) X(tn) || 
        = [Xt                          Xd                    Xc]
    
        where 
    
        Xt = time points in the indicator function
        Xd = derivatives
        Xc=variables to condition on
    
        Computations of all covariances follows simple rules: 
            Cov(X(t),X(s))=r(t,s),
        then  Cov(X'(t),X(s))=dr(t,s)/dt.  Now for stationary X(t) we have
        a function r(tau) such that Cov(X(t),X(s))=r(s-t) (or r(t-s) will give 
        the same result).

        Consequently  
            Cov(X'(t),X(s))    = -r'(s-t)    = -sign(s-t)*r'(|s-t|)
            Cov(X'(t),X'(s))   = -r''(s-t)   = -r''(|s-t|)
            Cov(X''(t),X'(s))  =  r'''(s-t)  =  sign(s-t)*r'''(|s-t|)
            Cov(X''(t),X(s))   =  r''(s-t)   =   r''(|s-t|)
            Cov(X''(t),X''(s)) =  r''''(s-t) = r''''(|s-t|)
    
        """
        # cov(Xd)
        Sdd = -toeplitz(R[[0, pt],2])      
        # cov(Xc)
        Scc = toeplitz(R[[0, pt],0])    
        # cov(Xc,Xd)
        Scd = np.array([[0, R[pt,1]],[ -R[pt,1], 0]])
        
        if pt > 1 :
            #%cov(Xt)
            Stt = toeplitz(R[:pt-1,0]) # Cov(X(tn),X(ts))  = r(ts-tn)   = r(|ts-tn|)
            #%cov(Xc,Xt) 
            Sct = R[1:pt,0]        # Cov(X(tn),X(ts))  = r(ts-tn)   = r(|ts-tn|)
            Sct = np.vstack((Sct,Sct[::-1]))
            #%Cov(Xd,Xt)
            Sdt = -R[1:pt,1]         # Cov(X'(t1),X(ts)) = -r'(ts-t1) = r(|s-t|)
            Sdt = np.vstack((Sdt, -Sdt[::-1]))
            #N   = pt + 3
            big = np.vstack((np.hstack((Stt, Sdt.T, Sct.T)),
                             np.hstack((Sdt, Sdd, Scd.T)),
                             np.hstack((Sct, Scd, Scc))))
        else:
            #N = 4
            big =  np.vstack((np.hstack((Sdd, Scd.T)),
                                    np.hstack((Scd, Scc))))
        return big

    def to_specnorm(self):
        S = self.copy()
        S.normalize()
        return S

    def sim(self,ns=None,cases=1,dt=None,iseed=None,method='random',derivative=False):
        ''' Simulates a Gaussian process and its derivative from spectrum

        Parameters
        ----------
        ns : scalar
            number of simulated points.  (default length(S)-1=n-1).
                     If ns>n-1 it is assummed that R(k)=0 for all k>n-1
        cases : scalar
            number of replicates (default=1)
        dt : scalar
            step in grid (default dt is defined by the Nyquist freq)
        iseed : int or state
            starting state/seed number for the random number generator
            (default none is set)
        method : string
            if 'exact'  : simulation using cov2sdat
            if 'random' : random phase and amplitude simulation (default)
        derivative : bool
            if true : return derivative of simulated signal as well
            otherwise

        Returns
        -------
        xs    = a cases+1 column matrix  ( t,X1(t) X2(t) ...).
        xsder = a cases+1 column matrix  ( t,X1'(t) X2'(t) ...).

        Details
        -------
        Performs a fast and exact simulation of stationary zero mean
        Gaussian process through circulant embedding of the covariance matrix
        or by summation of sinus functions with random amplitudes and random
        phase angle.

        If the spectrum has a non-empty field .tr, then the transformation is
        applied to the simulated data, the result is a simulation of a transformed
        Gaussian process.

        Note: The method 'exact' simulation may give high frequency ripple when
        used with a small dt. In this case the method 'random' works better.

        Example:
        >>> import wafo.spectrum.models as sm
        >>> Sj = sm.Jonswap();S = Sj.toSpecData()
        >>> ns =100; dt = .2
        >>> x1 = S.sim(ns,dt=dt)

        >>> import numpy as np
        >>> import scipy.stats as st
        >>> x2 = S.sim(20000,20)
        >>> truth1 = [0,np.sqrt(S.moment(1)[0]),0., 0.]
        >>> funs = [np.mean,np.std,st.skew,st.kurtosis]
        >>> for fun,trueval in zip(funs,truth1):
        ...     res = fun(x2[:,1::],axis=0)
        ...     m = res.mean()
        ...     sa = res.std()
        ...     assert(np.abs(m-trueval)<sa)

        waveplot(x1,'r',x2,'g',1,1)

        See also
        --------
        cov2sdat, gaus2dat

        Reference
        -----------
        C.R Dietrich and G. N. Newsam (1997)
        "Fast and exact simulation of stationary
        Gaussian process through circulant embedding
        of the Covariance matrix"
        SIAM J. SCI. COMPT. Vol 18, No 4, pp. 1088-1107

        Hudspeth, R.T. and Borgman, L.E. (1979)
        "Efficient FFT simulation of Digital Time sequences"
        Journal of the Engineering Mechanics Division, ASCE, Vol. 105, No. EM2,

        '''

        fft = np.fft.fft

        S = self.copy()
        if dt is not None:
            S.resample(dt)


        ftype = S.freqtype
        freq  = S.args

        dT = S.sampling_period()
        Nt = freq.size

        if ns is None:
            ns=Nt-1

        if method in 'exact':

            #nr=0,Nt=None,dt=None
            R = S.tocov(nr=0)
            T = Nt*dT
            ix = np.flatnonzero(R.args>T)

            # Trick to avoid adding high frequency noise to the spectrum
            if ix.size>0:
                R.data[ix[0]::]=0.0

            return R.sim(ns=ns, cases=cases, iseed=iseed, derivative=derivative)

        _set_seed(iseed)

        ns = ns+np.mod(ns,2) # make sure it is even

        fi    = freq[1:-1]
        Si    = S.data[1:-1]
        if ftype in ('w','k'):
            fact = 2.*np.pi
            Si = Si*fact
            fi = fi/fact

        zeros = np.zeros

        x = zeros((ns,cases+1))

        df = 1/(ns*dT)


        # interpolate for freq.  [1:(N/2)-1]*df and create 2-sided, uncentered spectra
        f = np.arange(1,ns/2.)*df

        Fs = np.hstack((0., fi, df*ns/2.))
        Su = np.hstack((0., np.abs(Si)/2., 0.))


        Si = np.interp(f,Fs,Su)
        Su=np.hstack((0., Si, 0, Si[(ns/2)-2::-1]))
        del(Si, Fs)

        # Generate standard normal random numbers for the simulations
        randn = np.random.randn
        Zr = randn((ns/2)+1,cases)
        Zi = np.vstack((zeros((1,cases)), randn((ns/2)-1,cases), zeros((1, cases))))

        A = zeros((ns,cases), dtype=complex)
        A[0:(ns/2+1), :] = Zr - 1j*Zi
        del(Zr, Zi)
        A[(ns/2+1):ns, :] = A[ns/2-1:0:-1, :].conj()
        A[0, :] = A[0, :]*np.sqrt(2.)
        A[(ns/2), :] = A[(ns/2), :]*np.sqrt(2.)


        # Make simulated time series
        T    = (ns-1)*dT
        Ssqr = np.sqrt(Su*df/2.)

        # stochastic amplitude
        A    = A*Ssqr[:,np.newaxis]


        # Deterministic amplitude
        #A = sqrt[1]*Ssqr(:,ones(1,cases)).*exp(sqrt(-1)*atan2(imag(A),real(A)))
        del( Su, Ssqr)


        x[:, 1::] = fft(A,axis=0).real
        x[:, 0] = np.linspace(0, T, ns) #' %(0:dT:(np-1)*dT).'


        if derivative:
            xder=np.zeros(ns,cases+1)
            w = 2.*np.pi*np.hstack((0, f, 0.,-f[-1::-1]))
            A = -1j*A*w[:,newaxis]
            xder[:,1:(cases+1)] = fft(A,axis=0).real
            xder[:,0]           = x[:,0]



        if S.tr is not None:
            np.disp('   Transforming data.')
            g = S.tr
            G=np.fliplr(g) #% the invers of g
            if derivative:
                for ix in range(cases):
                    tmp=tranproc(np.hstack((x[:,ix+1], xder[:,ix+1])),G)
                    x[:,ix+1]=tmp[:,0]
                    xder[:,ix+1]=tmp[:,1]

            else:
                for ix in range(cases):
                    x[:,ix+1]=tranproc(x[:,ix+1],G)



        if derivative:
            return x,xder
        else:
            return x

# function [x2,x,svec,dvec,A]=spec2nlsdat(S,np,dt,iseed,method,truncationLimit)
    def sim_nl(self,ns=None,cases=1,dt=None,iseed=None,method='random',
        fnlimit=1.4142,reltol=1e-3,g=9.81):
        """ 
        Simulates a Randomized 2nd order non-linear wave X(t)

        Parameters
        ----------
        ns : scalar
            number of simulated points.  (default length(S)-1=n-1).
            If ns>n-1 it is assummed that R(k)=0 for all k>n-1
        cases : scalar
            number of replicates (default=1)
        dt : scalar
            step in grid (default dt is defined by the Nyquist freq)
        iseed : int or state
            starting state/seed number for the random number generator
            (default none is set)
        method : string
            'apStochastic'    : Random amplitude and phase (default)
            'aDeterministic'  : Deterministic amplitude and random phase
            'apDeterministic' : Deterministic amplitude and phase
        fnLimit : scalar
            normalized upper frequency limit of spectrum for 2'nd order
            components. The frequency is normalized with
            sqrt(gravity*tanh(kbar*waterDepth)/Amax)/(2*pi)
            (default sqrt(2), i.e., Convergence criterion [1]_).
            Other possible values are:
            sqrt(1/2)  : No bump in trough criterion
            sqrt(pi/7) : Wave steepness criterion
        reltol : scalar
            relative tolerance defining where to truncate spectrum for the
            sum and difference frequency effects


        Returns
        -------
        xs2 = a cases+1 column matrix  ( t,X1(t) X2(t) ...).
        xs1 = a cases+1 column matrix  ( t,X1'(t) X2'(t) ...).

        Details
        -------
        Performs a Fast simulation of Randomized 2nd order non-linear
        waves by summation of sinus functions with random amplitudes and
        phase angles. The extent to which the simulated result are applicable
        to real seastates are dependent on the validity of the assumptions:

        1.  Seastate is unidirectional
        2.  Surface elevation is adequately represented by 2nd order random
            wave theory
        3.  The first order component of the surface elevation is a Gaussian
            random process.

        If the spectrum does not decay rapidly enough towards zero, the
        contribution from the 2nd order wave components at the upper tail can
        be very large and unphysical. To ensure convergence of the perturbation
        series, the upper tail of the spectrum is truncated at FNLIMIT in the
        calculation of the 2nd order wave components, i.e., in the calculation
        of sum and difference frequency effects. This may also be combined with
        the elimination of second order effects from the spectrum, i.e., extract
        the linear components from the spectrum. One way to do this is to use
        SPEC2LINSPEC.

        Example
        --------
        np =100; dt = .2
        [x1, x2] = spec2nlsdat(jonswap,np,dt)
        waveplot(x1,'r',x2,'g',1,1)

        See also
        --------
          spec2linspec, spec2sdat, cov2sdat

        References
        ----------
        .. [1] Nestegaard, A  and Stokka T (1995)
                A Third Order Random Wave model.
                In proc.ISOPE conf., Vol III, pp 136-142.

        .. [2] R. S Langley (1987)
                A statistical analysis of non-linear random waves.
                Ocean Engng, Vol 14, pp 389-407

        .. [3] Marthinsen, T. and Winterstein, S.R (1992)
                'On the skewness of random surface waves'
                In proc. ISOPE Conf., San Francisco, 14-19 june.
        """


        # TODO % Check the methods: 'apdeterministic' and 'adeterministic'

        from wafo.spectrum import dispersion_relation as sm

        Hm0,Tm02 = self.characteristic(['Hm0','Tm02'])[0].tolist()
        #Hm0 = self.characteristic('Hm0')[0]
        #Tm02 = self.characteristic('Tm02')[0]


        _set_seed(iseed)
        fft = np.fft.fft

        S = self.copy()
        if dt is not None:
            S.resample(dt)


        ftype = S.freqtype
        freq  = S.args

        dT = S.sampling_period()
        Nt = freq.size

        if ns is None:
            ns = Nt-1

        ns = ns+np.mod(ns,2) # make sure it is even

        fi    = freq[1:-1]
        Si    = S.data[1:-1]
        if ftype in ('w','k'):
            fact = 2.*np.pi
            Si = Si*fact
            fi = fi/fact

        Smax = max(Si)
        waterDepth = min(abs(S.h),10.**30)

        zeros = np.zeros

        x = zeros((ns,cases+1))

        df = 1/(ns*dT)

        # interpolate for freq.  [1:(N/2)-1]*df and create 2-sided, uncentered spectra
        # ----------------------------------------------------------------------------
        f = np.arange(1,ns/2.)*df
        Fs = np.hstack((0., fi, df*ns/2.))
        w = 2.*np.pi*np.hstack((0., f, df*ns/2.))
        kw = sm.w2k(w ,0.,waterDepth,g)[0]
        Su = np.hstack((0., np.abs(Si)/2., 0.))



        Si = np.interp(f,Fs,Su)
        nmin  = (Si>Smax*reltol).argmax()
        nmax  = np.flatnonzero(Si>0).max()
        Su = np.hstack((0., Si,0, Si[(ns/2)-2::-1]))
        del(Si, Fs)

        # Generate standard normal random numbers for the simulations
        # -----------------------------------------------------------
        randn = np.random.randn
        Zr = randn((ns/2)+1,cases)
        Zi = np.vstack((zeros((1,cases)), randn((ns/2)-1,cases), zeros((1,cases))))

        A                = zeros((ns,cases),dtype=complex)
        A[0:(ns/2+1),:]  = Zr - 1j*Zi
        del(Zr, Zi)
        A[(ns/2+1):ns,:] = A[ns/2-1:0:-1,:].conj()
        A[0,:]           = A[0,:]*np.sqrt(2.)
        A[(ns/2),:]    = A[(ns/2),:]*np.sqrt(2.)


        # Make simulated time series
        # --------------------------

        T    = (ns-1)*dT
        Ssqr = np.sqrt(Su*df/2.)


        if method.startswith('apd') : # apdeterministic
            # Deterministic amplitude and phase
            A[1:(ns/2),:]    = A[1,0]
            A[(ns/2+1):ns,:] = A[1,0].conj()
            A = sqrt(2)*Ssqr[:,np.newaxis]*exp(1J*arctan2(A.imag,A.real))
        elif method.startswith('ade'): # adeterministic
            # Deterministic amplitude and random phase
            A = sqrt(2)*Ssqr[:,np.newaxis]*exp(1J*arctan2(A.imag,A.real))
        else:
            # stochastic amplitude
            A    = A*Ssqr[:,np.newaxis]
        # Deterministic amplitude
        #A     =  sqrt(2)*Ssqr(:,ones(1,cases)).*exp(sqrt(-1)*atan2(imag(A),real(A)))
        del( Su, Ssqr)


        x[:,1::] = fft(A,axis=0).real
        x[:,0] = np.linspace(0,T,ns) #' %(0:dT:(np-1)*dT).'



        x2 = x.copy()

        # If the spectrum does not decay rapidly enough towards zero, the
        # contribution from the wave components at the  upper tail can be very
        # large and unphysical.
        # To ensure convergence of the perturbation series, the upper tail of the
        # spectrum is truncated in the calculation of sum and difference
        # frequency effects.
        # Find the critical wave frequency to ensure convergence.

        sqrt = np.sqrt
        log = np.log
        pi = np.pi
        tanh = np.tanh
        numWaves = 1000. # Typical number of waves in 3 hour seastate
        kbar = sm.w2k(2.*np.pi/Tm02,0.,waterDepth)[0]
        Amax = sqrt(2*log(numWaves))*Hm0/4 #% Expected maximum amplitude for 1000 waves seastate

        fLimitUp = fnlimit*sqrt(g*tanh(kbar*waterDepth)/Amax)/(2*pi)
        fLimitLo = sqrt(g*tanh(kbar*waterDepth)*Amax/waterDepth)/(2*pi*waterDepth)

        nmax   = min(np.flatnonzero(f<=fLimitUp).max(),nmax)+1
        nmin   = max(np.flatnonzero(fLimitLo<=f).min(),nmin)+1

        #if isempty(nmax),nmax = np/2end
        #if isempty(nmin),nmin = 2end % Must always be greater than 1
        fLimitUp = df*nmax
        fLimitLo = df*nmin

        print('2nd order frequency Limits = %g,%g'% (fLimitLo, fLimitUp))



##        if nargout>3,
##        %compute the sum and frequency effects separately
##        [svec, dvec] = disufq((A.'),w,kw,min(h,10^30),g,nmin,nmax)
##        svec = svec.'
##        dvec = dvec.'
##
##        x2s  = fft(svec) % 2'nd order sum frequency component
##        x2d  = fft(dvec) % 2'nd order difference frequency component
##
##        % 1'st order + 2'nd order component.
##        x2(:,2:end) =x(:,2:end)+ real(x2s(1:np,:))+real(x2d(1:np,:))
##        else
        A = A.T
        rvec,ivec = c_library.disufq(A.real,A.imag,w,kw,waterDepth,g,nmin,nmax,cases,ns)

        svec = rvec + 1J*ivec
        svec.shape = (cases,ns)
        x2o  = fft(svec,axis=1).T # 2'nd order component


        # 1'st order + 2'nd order component.
        x2[:,1::] = x[:,1::]+ x2o[0:ns,:].real

        return x2, x



    def stats_nl(self,h=None,moments='sk',method='approximate',g=9.81):
        """
        Statistics of 2'nd order waves to the leading order.

        Parameters
        ----------
        h : scalar
            water depth (default self.h)
        moments : string (default='sk')
            composed of letters ['mvsk'] specifying which moments to compute:
                   'm' = mean,
                   'v' = variance,
                   's' = (Fisher's) skew,
                   'k' = (Fisher's) kurtosis.
        method : string
            'approximate' method due to Marthinsen & Winterstein (default)
            'eigenvalue'  method due to Kac and Siegert

        Skewness = kurtosis-3 = 0 for a Gaussian process.
        The mean, sigma, skewness and kurtosis are determined as follows:
        method == 'approximate':  due to Marthinsen and Winterstein
        mean  = 2 * int Hd(w1,w1)*S(w1) dw1
        sigma = sqrt(int S(w1) dw1)
        skew  = 6 * int int [Hs(w1,w2)+Hd(w1,w2)]*S(w1)*S(w2) dw1*dw2/m0^(3/2)
        kurt  = (4*skew/3)^2

        where Hs = sum frequency effects  and Hd = difference frequency effects

        method == 'eigenvalue'

        mean  = sum(E)
        sigma = sqrt(sum(C^2)+2*sum(E^2))
        skew  = sum((6*C^2+8*E^2).*E)/sigma^3
        kurt  = 3+48*sum((C^2+E^2).*E^2)/sigma^4

        where
        h1 = sqrt(S*dw/2)
        C  = (ctranspose(V)*[h1;h1])
        and E and V is the eigenvalues and eigenvectors, respectively, of the 2'order
        transfer matrix. S is the spectrum and dw is the frequency spacing of S.

        Example:
        --------
        #Simulate a Transformed Gaussian process:
        >>> import wafo.spectrum.models as sm
        >>> Sj = sm.Jonswap()
        >>> S = Sj.toSpecData()
        >>> me,va,sk,ku = S.stats_nl(moments='mvsk')


        Hm0=7;Tp=11
        S = jonswap([],[Hm0 Tp]); [sk, ku, me]=spec2skew(S)
        g=hermitetr([],[Hm0/4 sk ku me]);  g2=[g(:,1), g(:,2)*Hm0/4]
        ys = spec2sdat(S,15000)   % Simulated in the Gaussian world
        xs = gaus2dat(ys,g2)      % Transformed to the real world

        See also
        ---------
        hermitetr, ochitr, lc2tr, dat2tr

        References:
        -----------
        Langley, RS (1987)
        'A statistical analysis of nonlinear random waves'
        Ocean Engineering, Vol 14, No 5, pp 389-407

        Marthinsen, T. and Winterstein, S.R (1992)
        'On the skewness of random surface waves'
        In proceedings of the 2nd ISOPE Conference, San Francisco, 14-19 june.

        Winterstein, S.R, Ude, T.C. and Kleiven, G. (1994)
        'Springing and slow drift responses:
        predicted extremes and fatigue vs. simulation'
        In Proc. 7th International behaviour of Offshore structures, (BOSS)
        Vol. 3, pp.1-15
        """

        #% default options
        if h is None:
            h = self.h

        #S = ttspec(S,'w')
        w = np.ravel(self.args)
        S = np.ravel(self.data)
        if self.freqtype in ['f','w']:
            vari = 't'
            if self.freqtype=='f':
               w = 2.*pi*w
               S = S/(2.*pi)
        #m0 = self.moment(nr=0)
        m0 = simps(S,w)
        sa = np.sqrt(m0)
        Nw = w.size

        Hs, Hd,Hdii = qtf(w,h,g)

        #%return
        #%skew=6/sqrt(m0)^3*simpson(S.w,simpson(S.w,(Hs+Hd).*S1(:,ones(1,Nw))).*S1.')

        Hspd = trapz(trapz((Hs+Hd)*S[np.newaxis,:],w)*S,w)
        output = []
        if method[0] == 'a': # %approx : Marthinsen, T. and Winterstein, S.R (1992) method
            if 'm' in moments:
                output.append( 2.*trapz(Hdii*S,w))
            if 'v' in moments:
                output.append(m0)
            skew =6./sa**3*Hspd
            if 's' in moments:
                output.append(skew)
            if 'k' in moments:
                output.append( (4.*skew/3.)**2.+3.)
        else:
            raise ValueError('Unknown option!')

##        elif method[0]== 'q': #, #% quasi method
##            Fn = self.nyquist_freq()
##            dw = Fn/Nw
##            tmp1 =sqrt(S[:,np.newaxis]*S[np.newaxis,:])*dw
##            Hd = Hd*tmp1
##            Hs = Hs*tmp1
##            k = 6
##            stop = 0
##            while !stop:
##                E = eigs([Hd,Hs;Hs,Hd],[],k)
##                %stop = (length(find(abs(E)<1e-4))>0 | k>1200)
##                %stop = (any(abs(E(:))<1e-4) | k>1200)
##                stop = (any(abs(E(:))<1e-4) | k>=min(2*Nw,1200))
##                k = min(2*k,2*Nw)
##            #end
##
##
##            m02=2*sum(E.^2) % variance of 2'nd order contribution
##
##            %Hstd = 16*trapz(S.w,(Hdii.*S1).^2)
##            %Hstd = trapz(S.w,trapz(S.w,((Hs+Hd)+ 2*Hs.*Hd).*S1(:,ones(1,Nw))).*S1.')
##            ma   = 2*trapz(S.w,Hdii.*S1)
##            %m02  = Hstd-ma^2% variance of second order part
##            sa   = sqrt(m0+m02)
##            skew = 6/sa^3*Hspd
##            kurt = (4*skew/3).^2+3
##        elif method[0]== 'e': #, % Kac and Siegert eigenvalue analysis
##            Fn = self.nyquist_freq()
##            dw = Fn/Nw
##            tmp1 =sqrt(S[:,np.newaxis]*S[np.newaxis,:])*dw
##            Hd = Hd*tmp1
##            Hs = Hs*tmp1
##            k = 6
##            stop = 0
##
##
##            while (not stop):
##              [V,D] = eigs([Hd,HsHs,Hd],[],k)
##              E = diag(D)
##              %stop = (length(find(abs(E)<1e-4))>0 | k>=min(2*Nw,1200))
##              stop = (any(abs(E(:))<1e-4) | k>=min(2*Nw,1200))
##              k = min(2*k,2*Nw)
##            #end
##
##
##            h1 = sqrt(S*dw/2)
##            C  = (ctranspose(V)*[h1;h1])
##
##            E2 = E.^2
##            C2 = C.^2
##
##            ma   = sum(E)                     % mean
##            sa   = sqrt(sum(C2)+2*sum(E2))    % standard deviation
##            skew = sum((6*C2+8*E2).*E)/sa^3   % skewness
##            kurt = 3+48*sum((C2+E2).*E2)/sa^4 % kurtosis


        return output





    def moment(self,nr=2,even=True,j=0):
        ''' Calculates spectral moments from spectrum

        Parameters
        ----------
        nr   : int
            order of moments (recomended maximum 4)
        even : bool
            False for all moments,
            True for only even orders
        j : int
            0 or 1

        Returns
        -------
        m     : list of moments
        mtext : list of strings describing the elements of m, see below

        Details
        -------
        Calculates spectral moments of up to order NR by use of
        Simpson-integration.

                 /                                  /
        mj_t^i = | w^i S(w)^(j+1) dw,  or  mj_x^i = | k^i S(k)^(j+1) dk
                 /                                  /

        where k=w^2/gravity, i=0,1,...,NR

        The strings in output mtext have the same position in the list
        as the corresponding numerical value has in output m
        Notation in mtext: 'm0' is the variance,
                        'm0x' is the first-order moment in x,
                       'm0xx' is the second-order moment in x,
                       'm0t'  is the first-order moment in t,
                             etc.
        For the calculation of moments see Baxevani et al.

        Example:
        >>> import numpy as np
        >>> import wafo.spectrum.models as sm
        >>> Sj = sm.Jonswap(Hm0=3)
        >>> w = np.linspace(0,4,256)
        >>> S = SpecData1D(Sj(w),w) #Make spectrum object from numerical values
        >>> S.moment()
        ([0.56220770033914191, 0.35433180985851975], ['m0', 'm0tt'])

        References
        ----------
        Baxevani A. et al. (2001)
        Velocities for Random Surfaces
        '''
        pi= np.pi
        one_dim_spectra = ['freq','enc','k1d']
        if self.type not in one_dim_spectra:
            raise ValueError('Unknown spectrum type!')

        f = np.ravel(self.args)
        S = np.ravel(self.data)
        if self.freqtype in ['f','w']:
            vari = 't'
            if self.freqtype=='f':
               f = 2.*pi*f
               S = S/(2.*pi)
        else:
            vari = 'x'
        S1=np.abs(S)**(j+1.)
        m = [simps(S1,x=f)]
        mtxt = 'm%d' % j
        mtext = [mtxt]
        step = np.mod(even,2)+1
        df = f**step
        for i in range(step,nr+1,step):
            S1 = S1*df
            m.append(simps(S1,x=f))
            mtext.append(mtxt+vari*i)
        return m, mtext
    
    def nyquist_freq(self):
        """
        Return Nyquist frequency
        """
        return self.args[-1]
    
    def sampling_period(self):
        ''' Returns sampling interval from Nyquist frequency of spectrum

        Returns
        ---------
        dT : scalar
            sampling interval, unit:
            [m] if wave number spectrum,
            [s] otherwise

        Let wm be maximum frequency/wave number in spectrum,
        then dT=pi/wm if angular frequency, dT=1/(2*wm) if natural frequency (Hz)

        Example
        -------
        S = jonswap
        dt = spec2dt(S)

        See also
        '''

        if self.freqtype in 'f':
            wmdt = 0.5  # Nyquist to sampling interval factor
        else: # ftype == w og ftype == k
            wmdt = np.pi

        wm = self.args[-1] #Nyquist frequency
        dt = wmdt/wm #sampling interval = 1/Fs
        return dt

    def resample(self,dt=None,Nmin=0,Nmax=2**13+1,method='stineman'):
        ''' Interpolate and zero-padd spectrum to change Nyquist freq.

        Parameters
        ----------
        dt : scalar
            wanted sampling interval (default as given by S, see spec2dt)
            unit: [s] if frequency-spectrum, [m] if wave number spectrum
        Nmin : scalar
            minimum number of frequencies.
        Nmax : scalar
            minimum number of frequencies
        method : string
            interpolation method (options are 'linear', 'cubic' or 'stineman')

        To be used before simulation (e.g. spec2sdat) or evaluation of covariance
        function (spec2cov) to directly get wanted sampling interval.
        The input spectrum is interpolated and padded with zeros to reach
        the right max-frequency, w(end)=pi/dt, f(end)=1/(2*dt), or k(end)=pi/dt.
        The objective is that output frequency grid should be at least as dense
        as the input grid, have equidistant spacing and length equal to
        2^k+1 (>=Nmin). If the max frequency is changed, the number of points
        in the spectrum is maximized to 2^13+1.

        Note: Also zero-padding down to zero freq, if S does not start there.
        If empty input dt, this is the only effect.

        See also
        --------
        spec2cov, spec2sdat, covinterp, spec2dt
        '''


        ftype = self.freqtype
        w     = self.args.ravel()
        n     = w.size

        #%doInterpolate = 0

        if ftype=='f':
            Cnf2dt = 0.5 # Nyquist to sampling interval factor
        else: #% ftype == w og ftype == k
            Cnf2dt = np.pi

        wnOld = w[-1]         # Old Nyquist frequency
        dTold = Cnf2dt/wnOld # sampling interval=1/Fs


        if dt is None:
            dt = dTold

        # Find how many points that is needed
        nfft   = 2**nextpow2(max(n-1,Nmin-1))
        dttest = dTold*(n-1)/nfft

        while (dttest>dt) and (nfft<Nmax-1):
            nfft   = nfft*2
            dttest = dTold*(n-1)/nfft

        nfft = nfft+1

        wnNew         = Cnf2dt/dt #% New Nyquist frequency
        dWn           = wnNew-wnOld
        doInterpolate = dWn>0 or w[1]>0 or (nfft!=n) or dt!=dTold or any(abs(np.diff(w,axis=0))>1.0e-8)

        if doInterpolate>0:
            S1 = self.data

            dw = min(np.diff(w))

            if dWn>0:
                #% add a zero just above old max-freq, and a zero at new max-freq
                #% to get correct interpolation there
                Nz = 1  + (dWn>dw) # % Number of zeros to add
                if Nz==2:
                    w = np.hstack((w, wnOld+dw, wnNew))
                else:
                    w = np.hstack((w, wnNew))

                S1 = np.hstack((S1, np.zeros(Nz)))

            if w[0]>0:
                #% add a zero at freq 0, and, if there is space, a zero just below min-freq
                Nz = 1 + (w[0]>dw) #% Number of zeros to add
                if Nz==2:
                    w=np.hstack((0, w[0]-dw, w))
                else:
                    w=np.hstack((0, w))

                S1 = np.hstack((zeros(Nz), S1))


            #% Do a final check on spacing in order to check that the gridding is
            #% sufficiently dense:
            #np1 = S1.size
            dwMin = np.finfo(float).max
            #%wnc = min(wnNew,wnOld-1e-5)
            wnc = wnNew
            specfun = lambda xi : stineman_interp(xi,w,S1)


            x,y = discretize(specfun,0,wnc)
            dwMin = np.minimum(min(np.diff(x)),dwMin)


            newNfft = 2**nextpow2(np.ceil(wnNew/dwMin))+1
            if newNfft>nfft:
                if (nfft<=2**15+1) and (newNfft>2**15+1):
                    warnings.warn('Spectrum matrix is very large (>33k). Memory problems may occur.')

                nfft = newNfft


            self.args = np.linspace(0,wnNew,nfft)
            if method=='stineman':
                self.data = stineman_interp(self.args,w,S1)
            else:
                intfun = interpolate.interp1d(w,S1,kind=method)
                self.data = intfun(self.args)
            self.data = self.data.clip(0) # clip negative values to 0

    def normalize(self, gravity=9.81):
        '''
        Normalize a spectral density such that m0=m2=1
        
        Paramter
        --------
        gravity=9.81

        Notes
        -----
        Normalization performed such that
        INT S(freq) dfreq = 1       INT freq^2  S(freq) dfreq = 1
        where integration limits are given by  freq  and  S(freq)  is the 
        spectral density; freq can be frequency or wave number.
        The normalization is defined by
            A=sqrt(m0/m2); B=1/A/m0; freq'=freq*A; S(freq')=S(freq)*B
    
        If S is a directional spectrum then a normalized gravity (.g) is added
        to Sn, such that mxx normalizes to 1, as well as m0 and mtt.
        (See spec2mom for notation of moments)
    
        If S is complex-valued cross spectral density which has to be
        normalized, then m0, m2 (suitable spectral moments) should be given.
    
        Example:
        ------- 
        S = jonswap
        [Sn,mn4] = specnorm(S)
        mts = spec2mom(S,2)     % Should be equal to one!
        '''
        mom, mtext = self.moment(nr=4, even=True)
        m0 = mom[0] 
        m2 = mom[1] 
        m4 = mom[2]  
        
        SM0 = sqrt(m0)
        SM2 = sqrt(m2)  
        A = SM0/SM2
        B = SM2/(SM0*m0)
         
        if self.freqtype=='f':
            self.args = self.args*A/2/pi
            self.data = self.data*B*2*pi  
        elif self.freqtype=='w' :
            self.args = self.args*A
            self.data = self.data*B
            m02 = m4/gravity**2
            m20 = m02
            self.g = gravity*sqrt(m0*m20)/m2
        self.A = A
        self.norm = True
        #S.date=datestr(now)

    def bandwidth(self,factors=0):
        ''' 
        Return some spectral bandwidth and irregularity factors

        Parameters
        -----------
        factors : array-like
            Input vector 'factors' correspondence:
            0 alpha=m2/sqrt(m0*m4)                        (irregularity factor)
            1 eps2 = sqrt(m0*m2/m1^2-1)                   (narrowness factor)
            2 eps4 = sqrt(1-m2^2/(m0*m4))=sqrt(1-alpha^2) (broadness factor)
            3 Qp=(2/m0^2)int_0^inf f*S(f)^2 df            (peakedness factor)

        Returns
        --------
        bw : arraylike
            vector of bandwidth factors
            Order of output is the same as order in 'factors'

        Example:
        >>> import numpy as np
        >>> import wafo.spectrum.models as sm
        >>> Sj = sm.Jonswap(Hm0=3)
        >>> w = np.linspace(0,4,256)
        >>> S = SpecData1D(Sj(w),w) #Make spectrum object from numerical values
        >>> S.bandwidth([0,1,2,3])
        array([ 0.65354446,  0.3975428 ,  0.75688813,  2.00207912])
        '''

        if self.freqtype in 'k':
            vari = 'k'
        else:
            vari = 'w'

        m, mtxt = self.moment(nr=4, even=False)
        sqrt = np.sqrt
        fact = np.atleast_1d(factors)
        alpha = m[2]/sqrt(m[0]*m[4])
        eps2 = sqrt(m[0]*m[2]/m[1]**2.-1.)
        eps4 = sqrt(1.-m[2]**2./m[0]/m[4])
        f = self.args
        S = self.data
        Qp = 2/m[0]**2.*simps(f*S**2,x=f)
        bw = np.array([alpha,eps2,eps4,Qp])
        return bw[fact]

    def characteristic(self, fact='Hm0', T=1200, g=9.81):
        """
        Returns spectral characteristics and their covariance

        Parameters
        ----------
        fact : vector with factor integers or a string or a list of strings
            defining spectral characteristic, see description below.
        T  : scalar
            recording time (sec) (default 1200 sec = 20 min)
        g : scalar
            acceleration of gravity [m/s^2]

        Returns
        -------
        ch : vector
            of spectral characteristics
        R  : matrix
            of the corresponding covariances given T
        chtext : a list of strings
            describing the elements of ch, see example.


        Description
        ------------
        If input spectrum is of wave number type, output are factors for
        corresponding 'k1D', else output are factors for 'freq'.
        Input vector 'factors' correspondence:
        1 Hm0   = 4*sqrt(m0)                              Significant wave height
        2 Tm01  = 2*pi*m0/m1                              Mean wave period
        3 Tm02  = 2*pi*sqrt(m0/m2)                        Mean zero-crossing period
        4 Tm24  = 2*pi*sqrt(m2/m4)                        Mean period between maxima
        5 Tm_10 = 2*pi*m_1/m0                             Energy period
        6 Tp    = 2*pi/{w | max(S(w))}                    Peak period
        7 Ss    = 2*pi*Hm0/(g*Tm02^2)                     Significant wave steepness
        8 Sp    = 2*pi*Hm0/(g*Tp^2)                       Average wave steepness
        9 Ka    = abs(int S(w)*exp(i*w*Tm02) dw ) /m0     Groupiness parameter
        10 Rs    = (S(0.092)+S(0.12)+S(0.15)/(3*max(S(w))) Quality control parameter
        11 Tp1   = 2*pi*int S(w)^4 dw                      Peak Period (robust estimate for Tp)
                  ------------------
                  int w*S(w)^4 dw

        12 alpha = m2/sqrt(m0*m4)                          Irregularity factor
        13 eps2  = sqrt(m0*m2/m1^2-1)                      Narrowness factor
        14 eps4  = sqrt(1-m2^2/(m0*m4))=sqrt(1-alpha^2)    Broadness factor
        15 Qp    = (2/m0^2)int_0^inf w*S(w)^2 dw           Peakedness factor

        Order of output is same as order in 'factors'
        The covariances are computed with a Taylor expansion technique
        and is currently only available for factors 1, 2, and 3. Variances
        are also available for factors 4,5,7,12,13,14 and 15

        Quality control:
        ----------------
        Critical value for quality control parameter Rs is Rscrit = 0.02
        for surface displacement records and Rscrit=0.0001 for records of
        surface acceleration or slope. If Rs > Rscrit then probably there
        are something wrong with the lower frequency part of S.

        Ss may be used as an indicator of major malfunction, by checking that
        it is in the range of 1/20 to 1/16 which is the usual range for
        locally generated wind seas.

        Examples:
        ---------
        >>> import numpy as np
        >>> import wafo.spectrum.models as sm
        >>> Sj = sm.Jonswap(Hm0=5)
        >>> S = Sj.toSpecData() #Make spectrum ob
        >>> S.characteristic(1)
        (array([ 8.59007646]), array([[ 0.03040216]]), ['Tm01'])

        >>> [ch, R, txt] = S.characteristic([1,2,3])  # fact a vector of integers
        >>> S.characteristic('Ss')               # fact a string
        (array([ 0.04963112]), array([[  2.63624782e-06]]), ['Ss'])

        >>> S.characteristic(['Hm0','Tm02'])   # fact a list of strings
        (array([ 4.99833578,  8.03139757]), array([[ 0.05292989,  0.02511371],
               [ 0.02511371,  0.0274645 ]]), ['Hm0', 'Tm02'])

        See also
        ---------
        bandwidth,
        moment

        References
        ----------
        Krogstad, H.E., Wolf, J., Thompson, S.P., and Wyatt, L.R. (1999)
        'Methods for intercomparison of wave measurements'
        Coastal Enginering, Vol. 37, pp. 235--257

        Krogstad, H.E. (1982)
        'On the covariance of the periodogram'
        Journal of time series analysis, Vol. 3, No. 3, pp. 195--207

        Tucker, M.J. (1993)
        'Recommended standard for wave data sampling and near-real-time processing'
        Ocean Engineering, Vol.20, No.5, pp. 459--474

        Young, I.R. (1999)
        "Wind generated ocean waves"
        Elsevier Ocean Engineering Book Series, Vol. 2, pp 239
        """

        #% TODO % Need more checking on computing the variances for Tm24,alpha, eps2 and eps4
        #% TODO % Covariances between Tm24,alpha, eps2 and eps4 variables are also needed
        NaN = np.nan
        tfact = dict(Hm0=0,Tm01=1,Tm02=2,Tm24=3, Tm_10=4,Tp=5,Ss=6, Sp=7, Ka=8,
              Rs=9, Tp1=10, Alpha=11, Eps2=12, Eps4=13, Qp=14)
        tfact1 = ('Hm0','Tm01','Tm02','Tm24', 'Tm_10','Tp','Ss', 'Sp', 'Ka',
              'Rs', 'Tp1','Alpha','Eps2','Eps4','Qp')

        if isinstance(fact,str):
            fact = list((fact,))
        if isinstance(fact,(list,tuple)):
            nfact = []
            for k in fact:
                if isinstance(k,str):
                    nfact.append(tfact.get(k.capitalize(),15))
                else:
                    nfact.append(k)
        else:
            nfact = fact

        nfact = np.atleast_1d(nfact)

        if np.any((nfact>14) | (nfact<0)):
            raise ValueError('Factor outside range (0,...,14)')

        vari = self.freqtype

        f  = self.args.ravel()
        S1 = self.data.ravel()
        m, mtxt = self.moment(nr=4,even=False)

        #% moments corresponding to freq  in Hz
        for k in range(1,5):
            m[k] = m[k]/(2*np.pi)**k

        #pi = np.pi
        ind  = np.flatnonzero(f>0)
        m.append(simps(S1[ind]/f[ind],f[ind])*2.*np.pi) #  % = m_1
        m_10 = simps(S1[ind]**2/f[ind],f[ind])*(2*pi)**2/T #    % = COV(m_1,m0|T=t0)
        m_11 = simps(S1[ind]**2./f[ind]**2,f[ind])*(2*pi)**3/T  #% = COV(m_1,m_1|T=t0)

        #sqrt = np.sqrt
        #%      Hm0        Tm01        Tm02             Tm24         Tm_10
        Hm0  = 4.*sqrt(m[0])
        Tm01 = m[0]/m[1]
        Tm02 = sqrt(m[0]/m[2])
        Tm24 = sqrt(m[2]/m[4])
        Tm_10= m[5]/m[0]

        Tm12 = m[1]/m[2]

        ind = S1.argmax()
        maxS = S1[ind]
        #[maxS ind] = max(S1)
        Tp   = 2.*pi/f[ind] #                                   % peak period /length
        Ss   = 2.*pi*Hm0/g/Tm02**2 #                             % Significant wave steepness
        Sp   = 2.*pi*Hm0/g/Tp**2 #                               % Average wave steepness
        Ka   = abs(simps(S1*np.exp(1J*f*Tm02),f))/m[0] #% groupiness factor

        #% Quality control parameter
        #% critical value is approximately 0.02 for surface displacement records
        #% If Rs>0.02 then there are something wrong with the lower frequency part
        #% of S.
        Rs   = np.sum(np.interp(np.r_[0.0146, 0.0195, 0.0244]*2*pi,f,S1))/3./maxS
        Tp2  = 2*pi*simps(S1**4,f)/simps(f*S1**4,f)


        alpha1 = Tm24/Tm02 #                 % m(3)/sqrt(m(1)*m(5))
        eps2   = sqrt(Tm01/Tm12-1.)#         % sqrt(m(1)*m(3)/m(2)^2-1)
        eps4   = sqrt(1.-alpha1**2) #          % sqrt(1-m(3)^2/m(1)/m(5))
        Qp     = 2./m[0]**2*simps(f*S1**2,f)




        ch = np.r_[Hm0, Tm01, Tm02, Tm24, Tm_10, Tp, Ss, Sp, Ka, Rs, Tp2, alpha1, eps2, eps4, Qp]


        #% Select the appropriate values
        ch     = ch[nfact]
        chtxt = [tfact1[i] for i in nfact]

        #if nargout>1,
        #% covariance between the moments:
        #%COV(mi,mj |T=t0) = int f^(i+j)*S(f)^2 df/T
        mij,mijtxt = self.moment(nr=8,even=False,j=1)
        for ix,tmp in enumerate(mij):
            mij[ix] = tmp/T/((2.*pi)**(ix-1.0))


        #% and the corresponding variances for
        #%{'hm0', 'tm01', 'tm02', 'tm24', 'tm_10','tp','ss', 'sp', 'ka', 'rs', 'tp1','alpha','eps2','eps4','qp'}
        R = np.r_[4*mij[0]/m[0],
       	    mij[0]/m[1]**2.-2.*m[0]*mij[1]/m[1]**3.+m[0]**2.*mij[2]/m[1]**4.,
        	0.25*(mij[0]/(m[0]*m[2])-2.*mij[2]/m[2]**2+m[0]*mij[4]/m[2]**3),
        	0.25*(mij[4]/(m[2]*m[4])-2*mij[6]/m[4]**2+m[2]*mij[8]/m[4]**3) ,
        	m_11/m[0]**2+(m[5]/m[0]**2)**2*mij[0]-2*m[5]/m[0]**3*m_10,
        	NaN,
        	(8*pi/g)**2*(m[2]**2/(4*m[0]**3)*mij[0]+mij[4]/m[0]-m[2]/m[0]**2*mij[2]),
        	NaN*np.ones(4),
        	m[2]**2*mij[0]/(4*m[0]**3*m[4])+mij[4]/(m[0]*m[4])+mij[8]*m[2]**2/(4*m[0]*m[4]**3)-
        	m[2]*mij[2]/(m[0]**2*m[4])+m[2]**2*mij[4]/(2*m[0]**2*m[4]**2)-m[2]*mij[6]/m[0]/m[4]**2,
        	(m[2]**2*mij[0]/4+(m[0]*m[2]/m[1])**2*mij[2]+m[0]**2*mij[4]/4-m[2]**2*m[0]*mij[1]/m[1]+
                m[0]*m[2]*mij[2]/2-m[0]**2*m[2]/m[1]*mij[3])/eps2**2/m[1]**4,
        	(m[2]**2*mij[0]/(4*m[0]**2)+mij[4]+m[2]**2*mij[8]/(4*m[4]**2)-m[2]*mij[2]/m[0]+
        	m[2]**2*mij[4]/(2*m[0]*m[4])-m[2]*mij[6]/m[4])*m[2]**2/(m[0]*m[4]*eps4)**2,
        	NaN]

        #% and covariances by a taylor expansion technique:
        #% Cov(Hm0,Tm01) Cov(Hm0,Tm02) Cov(Tm01,Tm02)
        S0 = np.r_[ 2./(sqrt(m[0])*m[1])*(mij[0]-m[0]*mij[1]/m[1]),
        	1./sqrt(m[2])*(mij[0]/m[0]-mij[2]/m[2]),
        	1./(2*m[1])*sqrt(m[0]/m[2])*(mij[0]/m[0]-mij[2]/m[2]-mij[1]/m[1]+m[0]*mij[3]/(m[1]*m[2]))]

        R1  = np.ones((15,15))
        R1[:,:] = NaN
        for ix,Ri in enumerate(R):
            R1[ix,ix] = Ri



        R1[0,2:4]   = S0[:2]
        R1[1,2]     = S0[2]
        for ix in [0,1]: #%make lower triangular equal to upper triangular part
            R1[ix+1:,ix] = R1[ix,ix+1:]


        R = R[nfact]
        R1= R1[nfact,:][:,nfact]


        #% Needs further checking:
        #% Var(Tm24)= 0.25*(mij[4]/(m[2]*m[4])-2*mij[6]/m[4]**2+m[2]*mij[8]/m[4]**3) ...
        return ch, R1, chtxt

    def setlabels(self):
        ''' Set automatic title, x-,y- and z- labels on SPECDATA object

            based on type, angletype, freqtype
        '''

        N = len(self.type)
        if N==0:
            raise ValueError('Object does not appear to be initialized, it is empty!')

        labels = ['','','']
        if self.type.endswith('dir'):
            title = 'Directional Spectrum'
            if self.freqtype.startswith('w'):
                labels[0] = 'Frequency [rad/s]'
                labels[2] = 'S(w,\theta) [m^2 s / rad^2]'
            else:
                labels[0] = 'Frequency [Hz]'
                labels[2] = 'S(f,\theta) [m^2 s / rad]'

            if self.angletype.startswith('r'):
                labels[1] = 'Wave directions [rad]'
            elif self.angletype.startswith('d'):
                labels[1] = 'Wave directions [deg]'
        elif self.type.endswith('freq'):
            title = 'Spectral density'
            if self.freqtype.startswith('w'):
                labels[0] = 'Frequency [rad/s]'
                labels[1] = 'S(w) [m^2 s/ rad]'
            else:
                labels[0] = 'Frequency [Hz]'
                labels[1] = 'S(f) [m^2 s]'
        else:
            title = 'Wave Number Spectrum'
            labels[0] =  'Wave number [rad/m]'
            if self.type.endswith('k1d'):
                labels[1] = 'S(k) [m^3/ rad]'
            elif self.type.endswith('k2d'):
                labels[1] = labels[0]
                labels[2] = 'S(k1,k2) [m^4/ rad^2]'
            else:
                raise ValueError('Object does not appear to be initialized, it is empty!')
        if self.norm!=0:
            title = 'Normalized ' + title
            labels[0] = 'Normalized ' + labels[0].split('[')[0]
            if not self.type.endswith('dir'):
                labels[1] = labels[1].split('[')[0]
            labels[2] = labels[2].split('[')[0]

        self.labels.title = title
        self.labels.xlab = labels[0]
        self.labels.ylab = labels[1]
        self.labels.zlab = labels[2]
class SpecData2D(WafoData):
    """ Container class for 2D spectrum data objects in WAFO

    Member variables
    ----------------
    data : array_like
    args : vector for 1D, list of vectors for 2D, 3D, ...

    type : string
        spectrum type (default 'freq')
    freqtype : letter
        frequency type (default 'w')
    angletype : string
        angle type of directional spectrum (default 'radians')

    Examples
    --------
    >>> import numpy as np
    >>> import wafo.spectrum.models as sm
    >>> Sj = sm.Jonswap(Hm0=3)
    >>> w = np.linspace(0,4,256)
    >>> S = SpecData1D(Sj(w),w) #Make spectrum object from numerical values

    See also
    --------
    WafoData
    CovData
    """

    def __init__(self,*args,**kwds):
        super(SpecData2D, self).__init__(*args,**kwds)

        self.name='WAFO Spectrum Object'
        self.type='freq'
        self.freqtype='w'
        self.angletype=''
        self.h=np.inf
        self.tr=None
        self.phi=0.
        self.v=0.
        self.norm=0
        somekeys = ['angletype', 'phi', 'name', 'h', 'tr', 'freqtype', 'v', 'type', 'norm']

        self.__dict__.update(sub_dict_select(kwds,somekeys))

        if self.type.endswith('dir') and self.angletype=='':
            self.angletype = 'radians'

        self.setlabels()

    def toacf(self):
        pass
    def sim(self):
        pass
    def sim_nl(self):
        pass
    def rotate(self):
        pass
    def moment(self,nr=2,vari='xt',even=True):
        ''' Calculates spectral moments from spectrum

        Parameters
        ----------
        nr   : int
            order of moments (maximum 4)
        vari : string
            variables in model, optional when two-dim.spectrum,
                   string with 'x' and/or 'y' and/or 't'
        even : bool
            False for all moments,
            True for only even orders

        Returns
        -------
        m     : list of moments
        mtext : list of strings describing the elements of m, see below

        Details
        -------
        Calculates spectral moments of up to order four by use of
        Simpson-integration.

           //
        m_jkl=|| k1^j*k2^k*w^l S(w,th) dw dth
           //

        where k1=w^2/gravity*cos(th-phi),  k2=w^2/gravity*sin(th-phi)
        and phi is the angle of the rotation in S.phi. If the spectrum
        has field .g, gravity is replaced by S.g.

        The strings in output mtext have the same position in the cell array
        as the corresponding numerical value has in output m
        Notation in mtext: 'm0' is the variance,
                        'mx' is the first-order moment in x,
                       'mxx' is the second-order moment in x,
                       'mxt' is the second-order cross moment between x and t,
                     'myyyy' is the fourth-order moment in y
                             etc.
        For the calculation of moments see Baxevani et al.

        Example:
          S=demospec('dir')
          [m,mtext]=spec2mom(S,2,'xyt')

        References
        ----------
        Baxevani A. et al. (2001)
        Velocities for Random Surfaces
        '''

##% Tested on: Matlab 6.0
##% Tested on: Matlab 5.3
##% History:
##% Revised by I.R. 04.04.2001: Introducing the rotation angle phi.
##% Revised by A.B. 23.05.2001: Correcting 'mxxyy' and introducing
##% 'mxxyt','mxyyt' and 'mxytt'.
##% Revised by A.B. 21.10.2001: Correcting 'mxxyt'.
##% Revised by A.B. 21.10.2001: Adding odd-order moments.
##% By es 27.08.1999


        pi= np.pi
        two_dim_spectra = ['dir','encdir','k2d']
        if self.type not in two_dim_spectra:
            raise ValueError('Unknown 2D spectrum type!')

##        if (vari==None and nr<=1:
##            vari='x'
##        elif vari==None:
##            vari='xt'
##        else #% secure the mutual order ('xyt')
##            vari=''.join(sorted(vari.lower()))
##            Nv=len(vari)
##
##            if vari[0]=='t' and Nv>1:
##                vari = vari[1::]+ vari[0]
##
##        Nv = len(vari)
##
##        if not self.type.endswith('dir'):
##            S1 = self.tospec(self.type[:-2]+'dir')
##        else:
##            S1 = self
##        w = np.ravel(S1.args[0])
##        theta = S1.args[1]-S1.phi
##        S = S1.data
##        Sw = simps(S,x=theta)
##        m = [simps(Sw,x=w)]
##        mtext=['m0']
##
##        if nr>0:
##
##          nw=w.size
##          if strcmpi(vari(1),'x')
##            Sc=simpson(th,S1.S.*(cos(th)*ones(1,nw))).'
##            % integral S*cos(th) dth
##          end
##          if strcmpi(vari(1),'y')
##            Ss=simpson(th,S1.S.*(sin(th)*ones(1,nw))).'
##            % integral S*sin(th) dth
##            if strcmpi(vari(1),'x')
##            Sc=simpson(th,S1.S.*(cos(th)*ones(1,nw))).'
##            end
##          end
##          if ~isfield(S1,'g')
##            S1.g=gravity
##          end
##          kx=w.^2/S1.g(1) % maybe different normalization in x and y => diff. g
##          ky=w.^2/S1.g(end)
##
##          if Nv>=1
##            switch vari
##              case 'x'
##                vec = kx.*Sc
##                mtext(end+1)={'mx'}
##              case 'y'
##                vec = ky.*Ss
##                mtext(end+1)={'my'}
##              case 't'
##                vec = w.*Sw
##               mtext(end+1)={'mt'}
##            end
##          else
##            vec = [kx.*Sc ky.*Ss w*Sw]
##            mtext(end+(1:3))={'mx', 'my', 'mt'}
##          end
##          if nr>1
##          if strcmpi(vari(1),'x')
##            Sc=simpson(th,S1.S.*(cos(th)*ones(1,nw))).'
##            % integral S*cos(th) dth
##            Sc2=simpson(th,S1.S.*(cos(th).^2*ones(1,nw))).'
##            % integral S*cos(th)^2 dth
##          end
##          if strcmpi(vari(1),'y')||strcmpi(vari(2),'y')
##            Ss=simpson(th,S1.S.*(sin(th)*ones(1,nw))).'
##            % integral S*sin(th) dth
##            Ss2=simpson(th,S1.S.*(sin(th).^2*ones(1,nw))).'
##            % integral S*sin(th)^2 dth
##            if strcmpi(vari(1),'x')
##              Scs=simpson(th,S1.S.*((cos(th).*sin(th))*ones(1,nw))).'
##              % integral S*cos(th)*sin(th) dth
##            end
##          end
##          if ~isfield(S1,'g')
##            S1.g=gravity
##          end
##
##          if Nv==2
##            switch vari
##              case 'xy'
##                vec=[kx.*Sc ky.*Ss kx.^2.*Sc2 ky.^2.*Ss2 kx.*ky.*Scs]
##                mtext(end+(1:5))={'mx','my','mxx', 'myy', 'mxy'}
##              case 'xt'
##                vec=[kx.*Sc w.*Sw kx.^2.*Sc2 w.^2.*Sw kx.*w.*Sc]
##                mtext(end+(1:5))={'mx','mt','mxx', 'mtt', 'mxt'}
##              case 'yt'
##                vec=[ky.*Ss w.*Sw ky.^2.*Ss2 w.^2.*Sw ky.*w.*Ss]
##                mtext(end+(1:5))={'my','mt','myy', 'mtt', 'myt'}
##            end
##          else
##            vec=[kx.*Sc ky.*Ss w.*Sw kx.^2.*Sc2 ky.^2.*Ss2  w.^2.*Sw kx.*ky.*Scs kx.*w.*Sc ky.*w.*Ss]
##            mtext(end+(1:9))={'mx','my','mt','mxx', 'myy', 'mtt', 'mxy', 'mxt', 'myt'}
##          end
##          if nr>3
##            if strcmpi(vari(1),'x')
##              Sc3=simpson(th,S1.S.*(cos(th).^3*ones(1,nw))).'
##              % integral S*cos(th)^3 dth
##              Sc4=simpson(th,S1.S.*(cos(th).^4*ones(1,nw))).'
##              % integral S*cos(th)^4 dth
##            end
##            if strcmpi(vari(1),'y')||strcmpi(vari(2),'y')
##              Ss3=simpson(th,S1.S.*(sin(th).^3*ones(1,nw))).'
##              % integral S*sin(th)^3 dth
##              Ss4=simpson(th,S1.S.*(sin(th).^4*ones(1,nw))).'
##              % integral S*sin(th)^4 dth
##              if strcmpi(vari(1),'x')  %both x and y
##                Sc2s=simpson(th,S1.S.*((cos(th).^2.*sin(th))*ones(1,nw))).'
##                % integral S*cos(th)^2*sin(th) dth
##                Sc3s=simpson(th,S1.S.*((cos(th).^3.*sin(th))*ones(1,nw))).'
##                % integral S*cos(th)^3*sin(th) dth
##                Scs2=simpson(th,S1.S.*((cos(th).*sin(th).^2)*ones(1,nw))).'
##                % integral S*cos(th)*sin(th)^2 dth
##                Scs3=simpson(th,S1.S.*((cos(th).*sin(th).^3)*ones(1,nw))).'
##                % integral S*cos(th)*sin(th)^3 dth
##                Sc2s2=simpson(th,S1.S.*((cos(th).^2.*sin(th).^2)*ones(1,nw))).'
##                % integral S*cos(th)^2*sin(th)^2 dth
##              end
##            end
##            if Nv==2
##              switch vari
##                case 'xy'
##                  vec=[vec kx.^4.*Sc4 ky.^4.*Ss4 kx.^3.*ky.*Sc3s ...
##                        kx.^2.*ky.^2.*Sc2s2 kx.*ky.^3.*Scs3]
##                  mtext(end+(1:5))={'mxxxx','myyyy','mxxxy','mxxyy','mxyyy'}
##                case 'xt'
##                  vec=[vec kx.^4.*Sc4 w.^4.*Sw kx.^3.*w.*Sc3 ...
##                        kx.^2.*w.^2.*Sc2 kx.*w.^3.*Sc]
##                  mtext(end+(1:5))={'mxxxx','mtttt','mxxxt','mxxtt','mxttt'}
##                case 'yt'
##                  vec=[vec ky.^4.*Ss4 w.^4.*Sw ky.^3.*w.*Ss3 ...
##                        ky.^2.*w.^2.*Ss2 ky.*w.^3.*Ss]
##                  mtext(end+(1:5))={'myyyy','mtttt','myyyt','myytt','myttt'}
##              end
##            else
##              vec=[vec kx.^4.*Sc4 ky.^4.*Ss4 w.^4.*Sw kx.^3.*ky.*Sc3s ...
##                   kx.^2.*ky.^2.*Sc2s2 kx.*ky.^3.*Scs3 kx.^3.*w.*Sc3 ...
##                   kx.^2.*w.^2.*Sc2 kx.*w.^3.*Sc ky.^3.*w.*Ss3 ...
##                   ky.^2.*w.^2.*Ss2 ky.*w.^3.*Ss kx.^2.*ky.*w.*Sc2s ...
##                   kx.*ky.^2.*w.*Scs2 kx.*ky.*w.^2.*Scs]
##              mtext(end+(1:15))={'mxxxx','myyyy','mtttt','mxxxy','mxxyy',...
##              'mxyyy','mxxxt','mxxtt','mxttt','myyyt','myytt','myttt','mxxyt','mxyyt','mxytt'}
##
##            end % if Nv==2 ... else ...
##          end % if nr>3
##          end % if nr>1
##          m=[m simpson(w,vec)]
##        end % if nr>0
##      %  end %%if Nv==1... else...    to be removed
##    end % ... else two-dim spectrum



    def interp(self):
        pass
    def normalize(self):
        pass
    def bandwidth(self):
        pass
    def setlabels(self):
        ''' Set automatic title, x-,y- and z- labels on SPECDATA object

            based on type, angletype, freqtype
        '''

        N = len(self.type)
        if N==0:
            raise ValueError('Object does not appear to be initialized, it is empty!')

        labels = ['','','']
        if self.type.endswith('dir'):
            title = 'Directional Spectrum'
            if self.freqtype.startswith('w'):
                labels[0] = 'Frequency [rad/s]'
                labels[2] = 'S(w,\theta) [m**2 s / rad**2]'
            else:
                labels[0] = 'Frequency [Hz]'
                labels[2] = 'S(f,\theta) [m**2 s / rad]'

            if self.angletype.startswith('r'):
                labels[1] = 'Wave directions [rad]'
            elif self.angletype.startswith('d'):
                labels[1] = 'Wave directions [deg]'
        elif self.type.endswith('freq'):
            title = 'Spectral density'
            if self.freqtype.startswith('w'):
                labels[0] = 'Frequency [rad/s]'
                labels[1] = 'S(w) [m**2 s/ rad]'
            else:
                labels[0] = 'Frequency [Hz]'
                labels[1] = 'S(f) [m**2 s]'
        else:
            title = 'Wave Number Spectrum'
            labels[0] =  'Wave number [rad/m]'
            if self.type.endswith('k1d'):
                labels[1] = 'S(k) [m**3/ rad]'
            elif self.type.endswith('k2d'):
                labels[1] = labels[0]
                labels[2] = 'S(k1,k2) [m**4/ rad**2]'
            else:
                raise ValueError('Object does not appear to be initialized, it is empty!')
        if self.norm!=0:
            title = 'Normalized ' + title
            labels[0] = 'Normalized ' + labels[0].split('[')[0]
            if not self.type.endswith('dir'):
                labels[1] = labels[1].split('[')[0]
            labels[2] = labels[2].split('[')[0]

        self.labels.title = title
        self.labels.xlab = labels[0]
        self.labels.ylab = labels[1]
        self.labels.zlab = labels[2]

class CovData1D(WafoData):
    """ Container class for 1D covariance data objects in WAFO

    Member variables
    ----------------
    data : array_like
    args : vector for 1D, list of vectors for 2D, 3D, ...

    type : string
        spectrum type, one of 'freq', 'k1d', 'enc' (default 'freq')
    lagtype : letter
        lag type, one of: 'x', 'y' or 't' (default 't')


    Examples
    --------
    >>> import numpy as np
    >>> import wafo.spectrum.models as sm
    >>> Sj = sm.Jonswap(Hm0=3)
    >>> w = np.linspace(0,4,256)
    >>> S = SpecData1D(Sj(w),w) #Make spectrum object from numerical values

    See also
    --------
    WafoData
    CovData
    """

    def __init__(self,*args,**kwds):
        super(CovData1D, self).__init__(*args,**kwds)

        self.name='WAFO Covariance Object'
        self.type='time'
        self.lagtype='t'
        self.h=np.inf
        self.tr=None
        self.phi=0.
        self.v=0.
        self.norm=0
        somekeys = ['phi', 'name', 'h', 'tr', 'lagtype', 'v', 'type', 'norm']

        self.__dict__.update(sub_dict_select(kwds,somekeys))

        self.setlabels()
    def setlabels(self):
        ''' Set automatic title, x-,y- and z- labels

            based on type,
        '''

        N = len(self.type)
        if N==0:
            raise ValueError('Object does not appear to be initialized, it is empty!')

        labels = ['','ACF','']

        if self.lagtype.startswith('t'):
            labels[0] = 'Lag [s]'
        else:
            labels[0] = 'Lag [m]'

        if self.norm!=0:
            title = 'Auto Correlation Function '
            labels[0] = labels[0].split('[')[0]
        else:
            title = 'Auto Covariance Function '

        self.labels.title = title
        self.labels.xlab = labels[0]
        self.labels.ylab = labels[1]
        self.labels.zlab = labels[2]



##    def copy(self):
##        kwds = self.__dict__.copy()
##        wdata = CovData1D(**kwds)
##        return wdata

    def tospec(self,rate=None,method='linear',nugget=0.0,trunc=1e-5,fast=True):
        '''
        Computes spectral density from the auto covariance function

        Parameters
        ----------
        rate = scalar, int
            1,2,4,8...2^r, interpolation rate for f (default 1)

        method: string
            interpolation method 'stineman', 'linear', 'cubic'

        nugget = scalar, real
            nugget effect to ensure that round off errors do not result in
            negative spectral estimates. Good choice might be 10^-12.

        trunc : scalar, real
            truncates all spectral values where S/max(S) < trunc
                      0 <= trunc <1   This is to ensure that high frequency
                      noise is not added to the spectrum.  (default 1e-5)
        fast  : bool
             if True : zero-pad to obtain power of 2 length ACF (default)
             otherwise  no zero-padding of ACF, slower but more accurate.

        Returns
        --------
        S = SpecData1D object
            spectral density

         NB! This routine requires that the covariance is evenly spaced
             starting from zero lag. Currently only capable of 1D matrices.

        Example:
        >>> import wafo.spectrum.models as sm
        >>> import numpy as np
        >>> import scipy.signal.signaltools as st
        >>> L = 129
        >>> t = np.linspace(0,75,L)
        >>> R = np.zeros(L)
        >>> win = st.parzen(41)
        >>> R[0:21] = win[20:41]
        >>> R0 = CovData1D(R,t)
        >>> S0 = R0.tospec()

        >>> Sj = sm.Jonswap()
        >>> S = Sj.toSpecData()
        >>> R2 = S.tocov()
        >>> S1 = R2.tospec()
        >>> assert(all(abs(S1.data-S.data)<1e-4) ,'COV2SPEC')

        See also
        --------
        spec2cov
        datastructures
        '''

        dT = self.sampling_period()
        # dT = time-step between data points.

        ACF, ti = np.atleast_1d(self.data, self.args)

        if self.lagtype in 't':
            spectype = 'freq'
            ftype = 'w'
        else:
            spectype = 'k1d'
            ftype = 'k'

        if rate is None:
            rate = 1 #%interpolation rate
        else:
            rate = 2**nextpow2(rate) #%make sure rate is a power of 2


        #% add a nugget effect to ensure that round off errors
        #% do not result in negative spectral estimates
        ACF[0] = ACF[0] +nugget
        n = ACF.size
        # embedding a circulant vector and Fourier transform
        if fast:
          nfft = 2**nextpow2(2*n-2)
        else:
          nfft = 2*n-2

        nf   = nfft/2 #% number of frequencies
        fft = np.fft.fft
        ACF  = np.r_[ACF,np.zeros(nfft-2*n+2),ACF[n-1:0:-1]]

        Rper = (fft(ACF,nfft).real).clip(0) #% periodogram
        RperMax = Rper.max()
        Rper = np.where(Rper<trunc*RperMax,0,Rper)
        pi = np.pi
        S = np.abs(Rper[0:(nf+1)])*dT/pi
        w = np.linspace(0,pi/dT,nf+1)
        So = SpecData1D(S,w,type=spectype,freqtype=ftype)
        So.tr = self.tr
        So.h = self.h
        So.norm = self.norm

        if rate > 1:
            So.args = np.linspace(0, pi/dT, nf*rate)
            if method=='stineman':
                So.data = stineman_interp(So.args, w, S)
            else:
                intfun = interpolate.interp1d(w, S, kind=method)
                So.data = intfun(So.args)
            So.data = So.data.clip(0) # clip negative values to 0
        return So

    def sampling_period(self):
        ''' Returns sampling interval

        Returns
        ---------
        dT : scalar
            sampling interval, unit:
            [s] if lagtype=='t'
            [m] otherwise
        See also
        '''
        dt1 = self.args[1]-self.args[0]
        n = np.size(self.args)-1
        t = self.args[-1]-self.args[0]
        dt = t/n
        return dt

    def sim(self, ns=None, cases=1, dt=None, iseed=None, derivative=False):
        ''' 
        Simulates a Gaussian process and its derivative from ACF

        Parameters
        ----------
        ns : scalar
            number of simulated points.  (default length(S)-1=n-1).
                     If ns>n-1 it is assummed that R(k)=0 for all k>n-1
        cases : scalar
            number of replicates (default=1)
        dt : scalar
            step in grid (default dt is defined by the Nyquist freq)
        iseed : int or state
            starting state/seed number for the random number generator
            (default none is set)
        derivative : bool
            if true : return derivative of simulated signal as well
            otherwise

        Returns
        -------
        xs    = a cases+1 column matrix  ( t,X1(t) X2(t) ...).
        xsder = a cases+1 column matrix  ( t,X1'(t) X2'(t) ...).

        Details
        -------
        Performs a fast and exact simulation of stationary zero mean
        Gaussian process through circulant embedding of the covariance matrix.

        If the ACF has a non-empty field .tr, then the transformation is
        applied to the simulated data, the result is a simulation of a transformed
        Gaussian process.

        Note: The simulation may give high frequency ripple when used with a
                small dt.

        Example:
        >>> import wafo.spectrum.models as sm
        >>> Sj = sm.Jonswap()
        >>> S = Sj.toSpecData()   #Make spec
        >>> R = S.tocov()
        >>> x = R.sim(ns=1000,dt=0.2)

        See also
        --------
        spec2sdat, gaus2dat

        Reference
        -----------
        C.R Dietrich and G. N. Newsam (1997)
        "Fast and exact simulation of stationary
        Gaussian process through circulant embedding
        of the Covariance matrix"
        SIAM J. SCI. COMPT. Vol 18, No 4, pp. 1088-1107

        '''

        # TODO fix it, it does not work
        
        # Add a nugget effect to ensure that round off errors
        # do not result in negative spectral estimates
        nugget = 0 # 10**-12

        _set_seed(iseed)

        ACF = self.data.ravel()
        n = ACF.size


        I = ACF.argmax()
        if I != 0:
            raise ValueError('ACF does not have a maximum at zero lag')

        ACF.shape = (n, 1)

        dT = self.sampling_period()

        fft = np.fft.fft

        x = np.zeros((ns, cases+1))

        if derivative:
            xder = x.copy()

        #% add a nugget effect to ensure that round off errors
        #% do not result in negative spectral estimates
        ACF[0] = ACF[0] + nugget

        #% Fast and exact simulation of simulation of stationary
        #% Gaussian process throug circulant embedding of the
        #% Covariance matrix
        floatinfo = np.finfo(float)
        if (abs(ACF[-1]) > floatinfo.eps): #% assuming ACF(n+1)==0
            m2 = 2*n-1
            nfft = 2**nextpow2(max(m2, 2*ns))
            ACF = np.r_[ACF, np.zeros((nfft-m2,1)), ACF[-1:0:-1,:]]
            #disp('Warning: I am now assuming that ACF(k)=0 ')
            #disp('for k>MAXLAG.')
        else: # % ACF(n)==0
            m2 = 2*n-2
            nfft = 2**nextpow2(max(m2, 2*ns))
            ACF = np.r_[ACF, np.zeros((nfft-m2, 1)), ACF[n-1:1:-1, :]]

        #%m2=2*n-2
        S = fft(ACF,nfft,axis=0).real #% periodogram

        I = S.argmax()
        k = np.flatnonzero(S<0)
        if k.size>0:
            #disp('Warning: Not able to construct a nonnegative circulant ')
            #disp('vector from the ACF. Apply the parzen windowfunction ')
            #disp('to the ACF in order to avoid this.')
            #disp('The returned result is now only an approximation.')

            # truncating negative values to zero to ensure that
            # that this noise is not added to the simulated timeseries

            S[k] = 0.

            ix = np.flatnonzero(k>2*I)
            if ix.size>0:
##    % truncating all oscillating values above 2 times the peak
##    % frequency to zero to ensure that
##    % that high frequency noise is not added to
##    % the simulated timeseries.
                ix0 = k[ix[0]]
                S[ix0:-ix0] =0.0


        sqrt = np.sqrt
        trunc = 1e-5
        maxS = S[I]
        k = np.flatnonzero(S[I:-I]<maxS*trunc)
        if k.size>0:
            S[k+I]=0.
            #% truncating small values to zero to ensure that
            #% that high frequency noise is not added to
            #% the simulated timeseries

        cases1 = np.floor(cases/2)
        cases2 = np.ceil(cases/2)
# Generate standard normal random numbers for the simulations

        randn = np.random.randn
        epsi = randn(nfft,cases2)+1j*randn(nfft,cases2)
        Ssqr = sqrt(S/(nfft)) # %sqrt(S(wn)*dw )
        ephat = epsi*Ssqr #[:,np.newaxis]
        y = fft(ephat,nfft,axis=0)
        x[:, 1:cases+1] = np.hstack((y[2:ns+2, 0:cases2].real, y[2:ns+2, 0:cases1].imag))

        x[:, 0] = np.linspace(0,(ns-1)*dT,ns) #%(0:dT:(dT*(np-1)))'

        if derivative:
            Ssqr = Ssqr*np.r_[0:(nfft/2+1), -(nfft/2-1):0]*2*np.pi/nfft/dT
            ephat = epsi*Ssqr #[:,np.newaxis]
            y = fft(ephat,nfft,axis=0)
            xder[:, 1:(cases+1)] = np.hstack((y[2:ns+2, 0:cases2].imag -y[2:ns+2, 0:cases1].real))
            xder[:, 0] = x[:,0]

        if self.tr is not None:
            np.disp('   Transforming data.')
            g = self.tr
            if derivative:
                for ix in range(cases):
                    tmp = g.gauss2dat(x[:,ix+1], xder[:,ix+1])
                    x[:,ix+1] = tmp[0]
                    xder[:,ix+1] = tmp[1]
            else:
                for ix in range(cases):
                    x[:, ix+1] = g.gauss2dat(x[:, ix+1])

        if derivative:
            return x, xder
        else:
            return x

class TurningPoints(WafoData):
    '''
    Container class for Turning Points data objects in WAFO

    Member variables
    ----------------
    data : array_like
    args : vector for 1D


    '''
    def __init__(self,*args,**kwds):
        super(TurningPoints, self).__init__(*args,**kwds)
        self.name='WAFO TurningPoints Object'
        somekeys = ['name']
        self.__dict__.update(sub_dict_select(kwds,somekeys))

        #self.setlabels()
        if not any(self.args):
            n = len(self.data)
            self.args = xrange(0,n)
        else:
            self.args = np.ravel(self.args)
        self.data= np.ravel(self.data)

    def cycle_pairs(self,type_='min2max'):
        """ Return min2Max or Max2min cycle pairs from turning points

        Parameters
        ----------
        type_ : string
            type of cycles to return options are 'min2max' or 'max2min'

        Return
        ------
        mm : cycles object
            with min2Max or Max2min cycle pairs.

        Example
        -------
        >>> import wafo
        >>> x = wafo.data.sea()
        >>> ts = wafo.objects.mat2timeseries(x)
        >>> tp = ts.turning_points()
        >>> mM = tp.cycle_pairs()
        >>> h = mM.plot('.')


        See also
        --------
        TurningPoints
        SurvivalCycleCount
        """


        if self.data[0]>self.data[1]:
            im = 1
            iM = 0
        else:
            im = 0
            iM = 1

        # Extract min-max and max-min cycle pairs
        n = len(self.data)
        if type_.lower().startswith('min2max'):
            m = self.data[im:-1:2]
            M = self.data[im+1::2]
        else:
            type_ = 'max2min'
            M = self.data[iM:-1:2]
            m = self.data[iM+1::2]

        return CyclePairs(M,m,type=type_)

def mat2timeseries(x):
    """
    Convert 2D arrays to TimeSeries object
        assuming 1st column is time and the remaining columns contain data.
    """
    return TimeSeries(x[:,1::],x[:,0].ravel())

class TimeSeries(WafoData):
    '''
    Container class for 1D TimeSeries data objects in WAFO

    Member variables
    ----------------
    data : array_like
    args : vector for 1D, list of vectors for 2D, 3D, ...

    sensortypes : list of integers or strings
        sensor type for time series (default ['n']    : Surface elevation)
        see sensortype for more options
    position : vector of size 3
        instrument position relative to the coordinate system

    Examples
    --------
    >>> import wafo.data
    >>> x = wafo.data.sea()
    >>> ts = mat2timeseries(x)
    >>> rf = ts.tocov(lag=150)
    >>> h = rf.plot()

    '''
    def __init__(self,*args,**kwds):
        super(TimeSeries, self).__init__(*args,**kwds)
        self.name='WAFO TimeSeries Object'
        self.sensortypes=['n',]
        self.position = np.zeros(3)
        somekeys = ['sensortypes', 'position']
        self.__dict__.update(sub_dict_select(kwds,somekeys))

        #self.setlabels()
        if not any(self.args):
            n = len(self.data)
            self.args = xrange(0,n)

    def sampling_period(self):
        '''
        Returns sampling interval

        Returns
        -------
        dt : scalar
            sampling interval, unit:
            [s] if lagtype=='t'
            [m] otherwise

        See also
        '''
        dt1 = self.args[1]-self.args[0]
        N = np.size(self.args)-1
        T = self.args[-1]-self.args[0]
        dt = T/N
        return dt

    def tocov(self, lag=None, flag='biased', norm=False, dt = None):
        ''' 
        Return auto covariance function from data.

        Parameters
        ----------
        lag : scalar, int
            maximum time-lag for which the ACF is estimated. (Default lag=n-1)
        flag : string, 'biased' or 'unbiased' 
            If 'unbiased' scales the raw correlation by 1/(n-abs(k)),
            where k is the index into the result, otherwise scales the raw 
            cross-correlation by 1/n. (default) 
        norm : bool
            True if normalize output to one
        dt : scalar
            time-step between data points (default xn(2,1)-xn(1,1) or 1 Hz).

        Return
        -------
        R : CovData1D object
            with attributes:
            data : ACF vector length L+1
            args : time lags  length L+1
            stdev : estimated large lag standard deviation of the estimate
                     assuming x is a Gaussian process:
                     if R(k)=0 for all lags k>q then an approximation
                     of the variance for large samples due to Bartlett
                     var(R(k))=1/N*(R(0)^2+2*R(1)^2+2*R(2)^2+ ..+2*R(q)^2)
                     for  k>q and where  N=length(x). Special case is
                     white noise where it equals R(0)^2/N for k>0
            norm : bool
                If false indicating that R is not normalized

         Example:
         --------
         >>> import wafo.data
         >>> x = wafo.data.sea()
         >>> ts = mat2timeseries(x)
         >>> acf = ts.tocov(150)
         >>> h = acf.plot()
        '''
        n = len(self.data)
        if not lag:
            lag = n-1

        x = self.data.flatten()
        indnan = np.isnan(x)
        if any(indnan):
            x = x - np.mean(x[1-indnan]) # remove the mean pab 09.10.2000
            #indnan = find(indnan)
            Ncens = n - sum(indnan)
            x[indnan] = 0. # pab 09.10.2000 much faster for censored samples
        else:
            indnan = None
            Ncens = n
            x = x - np.mean(x)

        fft = np.fft.fft
        nfft = 2**nextpow2(n)
        Rper = abs(fft(x,nfft))**2/Ncens # Raw periodogram

        R = np.real(fft(Rper))/nfft # %ifft=fft/nfft since Rper is real!
        lags = range(0,lag+1)
        if flag.startswith('unbiased'):
            # unbiased result, i.e. divide by n-abs(lag)
            R = R[lags]*Ncens/np.arange(Ncens, Ncens-lag, -1)
        #else  % biased result, i.e. divide by n
        #  r=r(1:L+1)*Ncens/Ncens

        c0 = R[0]
        if norm:
            R = R/c0
        if dt is None:
            dt = self.sampling_period()
        t = np.linspace(0,lag*dt,lag+1)
        cumsum = np.cumsum
        acf = CovData1D(R[lags],t)
        acf.stdev=np.sqrt(np.r_[ 0, 1 ,1+2*cumsum(R[1:]**2)]/Ncens)
        acf.children = [WafoData(-2.*acf.stdev[lags],t),WafoData(2.*acf.stdev[lags],t)]
        acf.norm = norm
        return acf

    def tospec(self,*args,**kwargs):
        """ 
        Return power spectral density by Welches average periodogram method.

        Parameters
        ----------
        NFFT : int, scalar
            if len(data) < NFFT, it will be zero padded to `NFFT`
            before estimation. Must be even; a power 2 is most efficient.
        detrend : function
        Fs : real, scalar
            sampling frequency (samples per time unit).

        window : vector of length NFFT or function
            To create window vectors see numpy.blackman, numpy.hamming,
            numpy.bartlett, scipy.signal, scipy.signal.get_window etc.
        noverlap : scalar int
             gives the length of the overlap between segments.

        Returns
        -------
        S : SpecData1D
            Power Spectral Density

        Notes
        -----
        The data vector is divided into NFFT length segments.  Each segment
        is detrended by function detrend and windowed by function window.
        noverlap gives the length of the overlap between segments.  The
        absolute(fft(segment))**2 of each segment are averaged to compute Pxx,
        with a scaling to correct for power loss due to windowing.

        Reference
        ---------
        Bendat & Piersol (1986) Random Data: Analysis and Measurement
        Procedures, John Wiley & Sons
        """
        fs = 1./(2*self.sampling_period())
        S, f = psd(self.data.ravel(), Fs=fs, *args, **kwargs)
        fact = 2.0*pi
        w = fact*f
        return SpecData1D(S/fact, w)

    def turning_points(self,h=0.0,wavetype=None):
        ''' 
        Return turning points (tp) from data, optionally rainflowfiltered.

        Parameters
        ----------
        h  : scalar
            a threshold
             if  h<=0, then  tp  is a sequence of turning points (default)
             if  h>0, then all rainflow cycles with height smaller than
                      h  are removed.

        wavetype : string
            defines the type of wave. Possible options are
            'mw' 'Mw' or 'none'.
            If None all rainflow filtered min and max
            will be returned, otherwise only the rainflow filtered
            min and max, which define a wave according to the
            wave definition, will be returned.

        Returns
        -------
        tp : TurningPoints object
            with times and turning points.

        Example:
        >>> import wafo.data
        >>> x = wafo.data.sea()
        >>> x1 = x[:200,:]
        >>> ts1 = mat2timeseries(x1)
        >>> tp = ts1.turning_points(wavetype='Mw')
        >>> tph = ts1.turning_points(h=0.3,wavetype='Mw')
        >>> hs = ts1.plot()
        >>> hp = tp.plot('ro')
        >>> hph = tph.plot('k.')

        See also
        ---------
        findcross,
        findrfc
        findtp
        '''
        ind = findtp(self.data,np.max(h,0.0), wavetype)
        try:
            t = self.args[ind]
        except:
            t = ind
        return TurningPoints(self.data[ind],t)

    def trough_crest(self,v=None,wavetype=None):
        """ 
        Return trough and crest turning points

        Parameters
        -----------
        v : scalar
            reference level (default  v = mean of x).

    	wavetype : string
            defines the type of wave. Possible options are
    		'dw', 'uw', 'tw', 'cw' or None.
    		If None indices to all troughs and crests will be returned,
    		otherwise only the paired ones will be returned
    		according to the wavedefinition.

        Returns
        --------
        tc : TurningPoints object
            with trough and crest turningpoints
        """
        ind = findtc(self.data, v, wavetype)[0]
        try:
            t = self.args[ind]
        except:
            t = ind
        return TurningPoints(self.data[ind], t)
    
    def wave_periods(self, vh=None, pdef='d2d', wdef=None, index=None, rate=1):
        """ 
        Return sequence of wave periods/lengths from data.

        Parameters
        ----------
        vh : scalar
            reference level ( default v=mean(x(:,2)) ) or
            rainflow filtering height (default h=0)
        pdef : string
            defining type of waveperiod (wavelength) returned:
            Level v separated 't2c', 'c2t', 't2t' or 'c2c' -waveperiod.
            Level v 'd2d', 'u2u', 'd2u' or 'u2d' -waveperiod.
            Rain flow filtered (with height greater than h)
            'm2M', 'M2m', 'm2m' or 'M2M' -waveperiod.
            Explanation to the abbreviations:
            M=Max, m=min, d=down-crossing, u=up-crossing ,
            t=trough and c=crest.
            Thus 'd2d' means period between a down-crossing to the
            next down-crossing and 'u2c' means period between a
            u-crossing to the following crest.
        wdef : string
            defining type of wave. Possible options are
            'mw','Mw','dw', 'uw', 'tw', 'cw' or None.
            If wdef is None all troughs and crests will be used,
            otherwise only the troughs and crests which define a
            wave according to the wavedefinition are used.

        index : vector
            index sequence of one of the following :
            -level v-crossings (indices to "du" are required to
                calculate 'd2d', 'd2u', 'u2d' or 'u2u' waveperiods)
            -level v separated trough and crest turningpoints
                (indices to 'tc' are required to calculate
                't2t', 't2c', 'c2t' or 'c2c' waveperiods)
            -level v crossings and level v separated trough and
                crest turningpoints (indices to "dutc" are
                required to calculate t2u, u2c, c2d or d2t
                waveperiods)
            -rainflow filtered turningpoints with minimum rfc height h
               (indices to "mMtc" are required to calculate
               'm2m', 'm2M', 'M2m' or 'M2M' waveperiods)

        rate : scalar
            interpolation rate. If rate larger than one, then x is
            interpolated before extrating T

        Returns
        --------
        T : vector
            sequence of waveperiods (or wavelengths).
        index : vector
            of indices


        Example:
        --------
        >>> import wafo
        >>> x = wafo.data.sea()
        >>> ts = wafo.objects.mat2timeseries(x[0:400,:])
        >>> T = ts.wave_periods(vh=0.0,pdef='c2c')

        T = dat2wa(x1,0,'c2c') #% Returns crest2crest waveperiods
        subplot(121), waveplot(x1,'-',1,1),subplot(122),histgrm(T)

        See also:
        --------
        findtp,
        findtc,
        findcross, perioddef
        """

##% This is a more flexible version than the dat2hwa or tp2wa routines.
##% There is a secret option: if pdef='all' the function returns
##% all the waveperiods 'd2t', 't2u', 'u2c' and 'c2d' in sequence.
##% It is up to the user to extract the right waveperiods.
##% If the first is a down-crossing then the first is a 'd2t' waveperiod.
##% If the first is a up-crossing then the first is a 'u2c' waveperiod.
##%
##%	Example:
##%		[T ind]=dat2wa(x,0,'all') %returns all waveperiods
##%		nn = length(T)
##%		% want to extract all t2u waveperiods
##%		if x(ind(1),2)>0 % if first is down-crossing
##%			Tt2u=T(2:4:nn)
##%		else 		% first is up-crossing
##%			Tt2u=T(4:4:nn)
##%		end

        if rate>1: #% interpolate with spline
            n = np.ceil(self.data.size*rate)
            ti = linspace(self.args[0], self.args[-1], n)
            x = stineman_interp(ti, self.args, self.data)
        else:
            x = self.data
            ti = self.args


        if vh is None:
            if pdef[0] in ('m','M'):
                vh = 0
                np.disp('   The minimum rfc height, h,  is set to: %g' % vh)
            else:
                vh = x.mean()
                np.disp('   The level l is set to: %g' % vh)


        if index is None:
            if pdef in ('m2m', 'm2M', 'M2m','M2M'):
                index = findtp(x, vh, wdef)
            elif pdef in ('u2u','u2d','d2u', 'd2d'):
                index = findcross(x, vh, wdef)
            elif pdef in ('t2t','t2c','c2t', 'c2c'):
                index = findtc(x,vh,wdef)[0]
            elif pdef in ('d2t','t2u', 'u2c', 'c2d','all'):
                index, v_ind = findtc(x, vh, wdef)
                index = sort(np.r_[index, v_ind]) #% sorting crossings and tp in sequence
            else:
                raise ValueError('Unknown pdef option!')

        if (x[index[0]]>x[index[1]]): #% if first is down-crossing or max
            if pdef in  ('d2t', 'M2m', 'c2t', 'd2u' , 'M2M', 'c2c', 'd2d', 'all'):
                start = 1
            elif pdef in ('t2u', 'm2M', 't2c', 'u2d' ,'m2m', 't2t', 'u2u'):
                start = 2
            elif pdef in ('u2c'):
                start = 3
            elif pdef in ('c2d'):
                start = 4
            else:
                raise ValueError('Unknown pdef option!')
            # else first is up-crossing or min
        elif pdef in ('all', 'u2c', 'm2M', 't2c', 'u2d', 'm2m', 't2t', 'u2u'):
            start = 0
        elif pdef in ('c2d', 'M2m', 'c2t', 'd2u', 'M2M', 'c2c', 'd2d'):
            start = 1
        elif pdef in ('d2t'):
            start = 2
        elif pdef in ('t2u'):
            start = 3
        else:
            raise ValueError('Unknown pdef option!')

        # determine the steps between wanted periods
        if pdef in ('d2t', 't2u', 'u2c', 'c2d' ):
            step = 4
        elif pdef in ('all'):
            step = 1 #% secret option!
        else:
            step = 2

        #% determine the distance between min2min, t2t etc..
        if pdef in ('m2m', 't2t', 'u2u', 'M2M', 'c2c', 'd2d'):
            dist = 2
        else:
            dist = 1

        nn = len(index)
        #% New call: (pab 28.06.2001)
        if pdef[0] in ('u', 'd'):
            t0 = ecross(ti, x, index[start:(nn-dist):step], vh)
        else: # % min, Max, trough, crest or all crossings wanted
            t0 = x[index[start:(nn-dist):step]]

        if pdef[2] in ('u','d'):
            t1 = ecross(ti, x, index[(start+dist):nn:step], vh)
        else: # % min, Max, trough, crest or all crossings wanted
            t1 = x[index[(start+dist):nn:step]]

        T = t1 - t0
##        if False: #% Secret option: indices to the actual crossings used.
##            index=index.ravel()
##            ind = [index(start:(nn-dist):step) index((start+dist):nn:step)].'
##            ind = ind(:)


        return T, index

        #% Old call: kept just in case
        #%T  = x(index((start+dist):step:nn),1)-x(index(start:step:(nn-dist)),1)



    def reconstruct(self):
        pass
    def plot_wave(self):
        pass
    def plot_sp_wave(self,wave_idx,tz_idx=None,*args,**kwds):
        """
        Plot specified wave(s) from timeseries

        wave_idx : integer vector
            of indices to waves we want to plot, i.e., wave numbers.
        tz_idx : integer vector
            of indices to the beginning, middle and end of
            defining wave, i.e. for zero-downcrossing waves, indices to
            zerocrossings (default trough2trough wave)

        Examples
        --------
        Plot waves nr. 6,7,8 and waves nr. 12,13,...,17
        >>> import wafo
        >>> x = wafo.data.sea()
        >>> ts = wafo.objects.mat2timeseries(x[0:500])
        >>> ts.plot_sp_wave(np.r_[6:9,12:18])


        See also
        --------
        waveplot, findtc
        """

        if tz_idx is None:
            tc_ind,tz_idx = findtc(self.data,0,'tw') # finding trough to trough waves

##        dw, = np.nonzero(np.abs(np.diff(wave_idx))>1)
##        Nsub = dw.size+1
##        Nwp = zeros(Nsub)
##        if Nsub>1:
##            Nwp[Nsub]= wave_idx[-1]-wave_idx(dw[-1])+1
##            wave_idx = wave_idx[0:dw[-1]+1]
##            for ix in range(Nsub-2,0,-2):
##                Nwp[ix+1] = wave_idx[dw[ix+1]-1]-wave_idx[dw[ix]]+1 # # of waves pr subplot
##                wave_idx(dw[ix]+1:dw[ix+1]-1)=[]
##
##            Nwp[0]= wave_idx(dw[0]-1)-wave_idx[0]+1
##            wave_idx(1:dw[0]-1)=[]
##        else:
##            Nwp[]= wave_idx[-1]-wave_idx[0]+1
##        end
##
##        Nsub = min(6,Nsub)
##        Nfig = np.ceil(Nsub/6)
##        Nsub = min(6,np.ceil(Nsub/Nfig))
##        gcf = plotbackend.gcf
##        for iy in range(Nfig):
##            for ix in range(Nsub):
##                plotbackend.subplot(Nsub,1,mod(ix,Nsub)+1)
##                ind = np.r_[tz_ind[2*wave_idx[ix]-1]:tz_ind[2*wave_idx[ix]+2*Nwp[ix]]]
##                #% indices to wave
##                plotbackend.plot(self.args[ind],self.data[ind],*args,**kwds)
##                plotbackend.hold('on')
##                plotbackend.plot(self.args[ind[0],ind[-1]],[0, 0])
##
##
##                if Nwp[ix]==1:
##                    plotbackend.ylabel('Wave %d' % wave_idx[ix])
##                else:
##                    plotbackend.ylabel('Wave %d - %d' % (wave_idx[ix], wave_idx[ix]+Nwp[ix]-1)])
##
##            xlabel('Time [sec]')
##            #wafostamp
##            if (mod(ix-1,Nsub)+1==Nsub)&& iy!=Nfig:
##                plotbackend.figure(gcf()+1)

def sensortypeid(*sensortypes):
    ''' Return ID for sensortype name

    Parameter
    ---------
    sensortypes : list of strings defining the sensortype

    Returns
    -------
    sensorids : list of integers defining the sensortype

    Valid senor-ids and -types for time series are as follows:
        0,  'n'    : Surface elevation              (n=Eta)
        1,  'n_t'  : Vertical surface velocity
        2,  'n_tt' : Vertical surface acceleration
        3,  'n_x'  : Surface slope in x-direction
        4,  'n_y'  : Surface slope in y-direction
        5,  'n_xx' : Surface curvature in x-direction
        6,  'n_yy' : Surface curvature in y-direction
        7,  'n_xy' : Surface curvature in xy-direction
        8,  'P'    : Pressure fluctuation about static MWL pressure
        9,  'U'    : Water particle velocity in x-direction
        10, 'V'    : Water particle velocity in y-direction
        11, 'W'    : Water particle velocity in z-direction
        12, 'U_t'  : Water particle acceleration in x-direction
        13, 'V_t'  : Water particle acceleration in y-direction
        14, 'W_t'  : Water particle acceleration in z-direction
        15, 'X_p'  : Water particle displacement in x-direction from its mean position
        16, 'Y_p'  : Water particle displacement in y-direction from its mean position
        17, 'Z_p'  : Water particle displacement in z-direction from its mean position

    Example:
    >>> sensortypeid('W','v')
    [11, 10]
    >>> sensortypeid('rubbish')
    [1.#QNAN]

    See also 
    --------
    sensortype
    '''

    sensorid_table = dict(n=0, n_t=1, n_tt=2, n_x=3, n_y=4, n_xx=5,
        n_yy=6, n_xy=7, p=8, u=9, v=10, w=11, u_t=12,
        v_t=13, w_t=14, x_p=15, y_p=16, z_p=17)
    try:
        return [sensorid_table.get(name.lower(), np.NAN) for name in sensortypes]
    except:
        raise ValueError('Input must be a string!')



def sensortype(*sensorids):
    ''' 
    Return sensortype name

    Parameter
    ---------
    sensorids : vector or list of integers defining the sensortype

    Returns
    -------
    sensornames : tuple of strings defining the sensortype
        Valid senor-ids and -types for time series are as follows:
        0,  'n'    : Surface elevation              (n=Eta)
        1,  'n_t'  : Vertical surface velocity
        2,  'n_tt' : Vertical surface acceleration
        3,  'n_x'  : Surface slope in x-direction
        4,  'n_y'  : Surface slope in y-direction
        5,  'n_xx' : Surface curvature in x-direction
        6,  'n_yy' : Surface curvature in y-direction
        7,  'n_xy' : Surface curvature in xy-direction
        8,  'P'    : Pressure fluctuation about static MWL pressure
        9,  'U'    : Water particle velocity in x-direction
        10, 'V'    : Water particle velocity in y-direction
        11, 'W'    : Water particle velocity in z-direction
        12, 'U_t'  : Water particle acceleration in x-direction
        13, 'V_t'  : Water particle acceleration in y-direction
        14, 'W_t'  : Water particle acceleration in z-direction
        15, 'X_p'  : Water particle displacement in x-direction from its mean position
        16, 'Y_p'  : Water particle displacement in y-direction from its mean position
        17, 'Z_p'  : Water particle displacement in z-direction from its mean position

    Example:
    >>> sensortype(range(3))
    ('n', 'n_t', 'n_tt')

    See also 
    --------
    sensortypeid, tran
    '''
    valid_names = ('n', 'n_t', 'n_tt', 'n_x', 'n_y', 'n_xx', 'n_yy', 'n_xy',
                  'p', 'u', 'v', 'w', 'u_t', 'v_t', 'w_t', 'x_p', 'y_p', 'z_p',
                  np.NAN)
    ids = np.atleast_1d(*sensorids)
    if isinstance(ids, list):
        ids = np.hstack(ids)
    n = len(valid_names) - 1
    ids = where(((ids<0) | (n<ids)), n , ids)
    
    #try:
    return tuple(valid_names[i] for i in ids)
    #except:
    #    raise ValueError('Input must be an integer!')


def main():
    import matplotlib
    matplotlib.interactive(True)
    from wafo.spectrum import models as sm
    sensortype(range(21))
    w = np.linspace(0,3,100)
    Sj = sm.Jonswap()
    S = Sj.toSpecData()
    
    f = S.to_t_pdf(pdef='Tc', paramt=(0, 10, 51), speed=7)
    f.err
    f.plot()
    f.show()
    #pdfplot(f)
    #hold on, 
    #plot(f.x{:}, f.f+f.err,'r',f.x{:}, f.f-f.err)  estimated error bounds
    #hold off  

    
    
    #S = SpecData1D(Sj(w),w)
    R = S.tocov(nr=1)
    S1 = S.copy()
    Si = R.tospec()
    ns =5000
    dt = .2
    x1 = S.sim_nl(ns= ns,dt=dt)
    x2 = TimeSeries(x1[:,1],x1[:,0])
    R = x2.tocov(lag=100)
    R.plot()

    S.plot('ro')
    t = S.moment()
    t1 = S.bandwidth([0,1,2,3])
    S1 = S.copy()
    S1.resample(dt=0.3, method='cubic')
    S1.plot('k+')
    x = S1.sim(ns=100)
    import pylab
    pylab.clf()
    pylab.plot(x[:,0],x[:,1])
    pylab.show()

    pylab.close('all')
    np.disp('done')


if __name__ == '__main__':
    if  True: #False : #  
        import doctest
        doctest.testmod()
    else:
        main()
        import wafo
        ts = wafo.objects.mat2timeseries(wafo.data.sea())
        tp = ts.turning_points()
        mm = tp.cycle_pairs()
        lc = mm.level_crossings()
        lc.plot()
        T = ts.wave_periods(vh=0.0,pdef='c2c')
        import wafo.spectrum.models as sm
        Sj = sm.Jonswap()
        S = Sj.toSpecData()
        me,va,sk,ku = S.stats_nl(moments='mvsk')
        import wafo.data
        x = wafo.data.sea()
        ts = mat2timeseries(x)
        rf = ts.tocov(lag=150)
        rf.plot()
        #main()
        import wafo.spectrum.models as sm
        Sj = sm.Jonswap()
        S = Sj.toSpecData()

        R = S.tocov()
        x = R.sim(ns=1000,dt=0.2)
        S.characteristic(['hm0','tm02'])
        ns = 1000 
        dt = .2
        x1 = S.sim(ns,dt=dt)

        ts = TimeSeries(x1[:,1],x1[:,0])
        tp = ts.turning_points(0.0)

        x = np.arange(-2,2,0.2)

        # Plot 2 objects in one call
        d2 = WafoData(np.sin(x),x,xlab='x',ylab='sin',title='sinus')


        d0 = d2.copy()
        d0.data = d0.data*0.9
        d1 = d2.copy()
        d1.data = d1.data*1.2
        d1.children = [d0]
        d2.children = [d1]

        d2.plot()
        print 'Done'