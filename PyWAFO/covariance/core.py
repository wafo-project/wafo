from __future__ import division
import warnings
import numpy as np
import scipy.interpolate as interpolate
from pylab import stineman_interp

from wafo.wafodata import WafoData
from wafo.misc import sub_dict_select, nextpow2 #, JITImport
import wafo.spectrum as _wafospec
#_wafospec = JITImport('wafo.spectrum')


__all__ = ['CovData1D']

def _set_seed(iseed):
    if iseed != None:
        try:
            np.random.set_state(iseed)
        except:
            np.random.seed(iseed)
            
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
    >>> import wafo.spectrum as sp
    >>> Sj = sp.models.Jonswap(Hm0=3)
    >>> w = np.linspace(0,4,256)
    >>> S = sp.SpecData1D(Sj(w),w) #Make spectrum object from numerical values

    See also
    --------
    WafoData
    CovData
    """

    def __init__(self,*args,**kwds):
        super(CovData1D, self).__init__(*args,**kwds)

        self.name = 'WAFO Covariance Object'
        self.type = 'time'
        self.lagtype = 't'
        self.h = np.inf
        self.tr = None
        self.phi = 0.
        self.v = 0.
        self.norm = 0
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

        if self.norm:
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

    def tospecdata(self, rate=None, method='linear', nugget=0.0, trunc=1e-5, fast=True):
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
        >>> S0 = R0.tospecdata()

        >>> Sj = sm.Jonswap()
        >>> S = Sj.tospecdata()
        >>> R2 = S.tocovdata()
        >>> S1 = R2.tospecdata()
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
        So = _wafospec.SpecData1D(S, w, type=spectype, freqtype=ftype)
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
        ''' 
        Returns sampling interval

        Returns
        ---------
        dt : scalar
            sampling interval, unit:
            [s] if lagtype=='t'
            [m] otherwise
        '''
        dt1 = self.args[1]-self.args[0]
        n = np.size(self.args)-1
        t = self.args[-1]-self.args[0]
        dt = t/n
        if abs(dt-dt1) > 1e-10:
            warnings.warn('Data is not uniformly sampled!')
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
        >>> S = Sj.tospecdata()   #Make spec
        >>> R = S.tocovdata()
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


def test_covdata():
    import wafo.data
    x = wafo.data.sea()
    ts = mat2timeseries(x)
    rf = ts.tocovdata(lag=150)
    rf.plot()
    
def main():
    import wafo.spectrum.models as sm
    import matplotlib
    matplotlib.interactive(True)
    Sj = sm.Jonswap()
    S = Sj.tospecdata()   #Make spec
    S.plot()
    R = S.tocovdata()
    R.plot()
    #x = R.sim(ns=1000,dt=0.2)


if __name__ == '__main__':
    if  True: #False : #  
        import doctest
        doctest.testmod()
    else:
        main()