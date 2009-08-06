#-------------------------------------------------------------------------------
# Name:        kdetools
# Purpose:
#
# Author:      pab
#
# Created:     01.11.2008
# Copyright:   (c) pab2 2008
# Licence:     LGPL
#-------------------------------------------------------------------------------
#!/usr/bin/env python
import numpy as np
from scipy.special import gamma
from numpy import pi, atleast_2d
def sphere_volume(d,r=1.0):
    """
     Returns volume of  d-dimensional sphere with radius r

    Parameters
    ----------
    d : scalar or array_like
        dimension of sphere
    r : scalar or array_like
        radius of sphere (default 1)

    Reference
    ---------
    Wand,M.P. and Jones, M.C. (1995)
    'Kernel smoothing'
    Chapman and Hall, pp 105
    """
    return (r0**d)* 2.*pi**(d/2.)/(d*gamma(d/2.))
def trangood(ff,nmin=None,mini=None,maxi=None,nmax=inf):
    """
    Make sure transformation is efficient.

    Parameters
    ------------
    ff : two column array
        input transform function, [X f(X)].
    nmin : scalar, int
        minimum number of points in the good transform.
               (Default  ff.shape[0])
    mini : scalar, real
        minimum data value to transform.
              (Default  min(ff[:,0]))
    maxi : scalar, real
        maximum data value to transform.
               (Default  max(ff[:,0]))
    nmax : scalar, int
        maximum number of points in the good transform
              (default inf)
    Returns
    -------
    f    = the good transform function, [X f(X)].

    TRANGOOD interpolates ff linearly  and optionally
    extrapolate ff linearly outside the range of ff(:,1)
    with X uniformly spaced.

    See also
    ---------
    tranproc, interp1q
    """
    ff = atleast_2d(ff)
    n,d = ff.shape
    if (d!=2):
        raise ValueError('ff  must be a two column matrix.')

    if (n<2):
        raise ValueError('ff  must have at least two rows.')

    i = ff[:,0].argsort()
    f = ff[i,:]
    del i
    df    = diff(f[:,0]);
    if ( any(df<=0)):
        raise ValueError('Duplicate x-values in  ff  not allowed.')


    nf = f.shape[0]

    if maxi is None:
        maxi = f[nf,0]
    if mini is None:
        mini = f[0,0]
    if nmin is None:
        nmin = nf
    if (nmin<2):
        nmin = 2
    if (nmax<2):
        nmax = 2

    ddf = diff(df)
    xn = f[nf,0]
    x0 = f[0,0]
    L = xn-x0
    if ( (nf<nmin) or (nmax<nf) or any(abs(ddf)>10*eps*(L)) ):
##  % pab 07.01.2001: Always choose the stepsize df so that
##  % it is an exactly representable number.
##  % This is important when calculating numerical derivatives and is
##  % accomplished by the following.
        df = L/(min(nmin,nmax)-1)
        df = (df+2.)-2.
        x = arange(x0,xn,df)
        #% New call pab 11.11.2000: This is much quicker
        f = array([ x, np.interp(x,f[:,0],f[:,1])]).T


    #% f(:,1) is now uniformly spaced
    df = f[1,0]-f[0,0]

    #% Extrapolate linearly outside the range of ff
    #%----------------------------------------------
    if (mini<f[0,0]):
        f1 = df*arange(floor((mini-f[0,0])/df),-2)
        f2 = f[0,1]+f1*(f[1,1]-f[0,1])/(f[1,0]-f[0,0]);
        f  = np.vstack([[f1+f[0,0], f2],f])

    n = f.shape[0];
    if (maxi>f[n,0]):
        f1 = df*arange(1,ceil((maxi-f(n,1))/df)+1)
        f2 = f(n,2)+f1*(f(n,2)-f(n-1,2))/(f(n,1)-f(n-1,1));
        f  = np.vstack([f,[f1+f[n,0], f2]])
    return f



def tranproc(x,ff):
    """
    Transforms process X and up to four derivatives
          using the transformation f.

    Parameters
    ----------
     x = input data matrix with 1+N columns, [X X1 ... XN], where
         Xi  is the i'th time derivative of  X. 0<=N<=4.
     f = [x,f(x)], transform function, y = f(x).

    Returns
    -------
    y = output data matrix with 1+N columns, [Y Y1 ...YN], of
         transformed data, where Y = f(X) and  Yi is the i'th time
        derivative of Y = f(X).

    By the basic rules of derivation:
    Y1 = f'(X)*X1
    Y2 = f''(X)*X1^2 + f'(X)*X2
    Y3 = f'''(X)*X1^3 + f'(X)*X3 + 3*f''(X)*X1*X2
    Y4 = f''''(X)*X1^4 + f'(X)*X4 + 6*f'''(X)*X1^2*X2
      + f''(X)*(3*X2^2 + 4*X1*X3)

    The derivation of f is performed numerically with a central difference
    method with linear extrapolation towards the beginning and end of f,
    respectively.

    Example
    --------
    # Derivative of g and the transformed Gaussian model.
    >>> x = linspace(-6,6,501)';
    >>> g = hermitetr(x);
    >>> gder = tranproc([g(:,1) ones(size(g,1),1)],g);
    >>> gder(:,1) = g(:,1);
    >>> plot(g(:,1),[g(:,2),gder(:,2)])
    >>> plot(g(:,1),pdfnorm(g(:,2)).*gder(:,2),g(:,1),pdfnorm(g(:,1)))
    >>> legend('Transformed model','Gaussian model')

    See also  trangood.
    """
    x,ff = atleast_2d(x,ff)
    N    = x.shape[1]-1 # N = number of derivatives
    nmax = ceil((max(ff[:,0])-min(ff[:,0]))*10**(7./max(N,1)));
    f    = trangood(ff,ff.shape[0],min(x[:,0]),max(x[:,0]),nmax)

    n  = f.shape[0]
    y  = x.copy()
    xu = 1+(n-1)*(x[:,0]-f[0,0])/(f[n,0]-f[0,0])

    fi = floor(xu)

    i  = find(fi==n);
    fi[i] = fi[i]-1

    xu = xu-fi
    y[:,0] = f[fi,1]+(f[fi+1,1]-f[fi,1])*xu


    if N>0:
        hn = f[1,0]-f[0,0];
        if hn**N<sqrt(eps):
            disp('Numerical problems may occur for the derivatives in tranproc.')
            warnings.warn('The sampling of the transformation may be too small.')

        #% Transform X with the derivatives of  f.
        fxder = zeros((size(x,1),N));
        fder  = f
        for k in range(N): #% Derivation of f(x) using a difference method.
            n = fder.shape[0]
            #%fder = [(fder(1:n-1,1)+fder(2:n,1))/2 diff(fder(:,2))./diff(fder(:,1))];
            fder = np.vstack([(fder[0:n-1,0]+fder[1:n,0])/2, diff(fder[:,1])/hn])
            fxder[:,k] = tranproc(x[:,0],fder)

            #%(-fder(ix+2,2)+8*fder(ix+1,2) - ...
            #%	      8*fder(ix-1,2)+fder(ix-2,2))./(12*hn);
        #% Calculate the transforms of the derivatives of X.
        #% First time derivative of y: y1 = f'(x)*x1
        y[:,1]= fxder[:,0]*x[:,1]
        if N>1:
            #% Second time derivative of y:
            #%             y2 = f''(x)*x1.^2+f'(x)*x2
            y[:,2]=fxder[:,1]*x[:,1]**2. + fxder[:,0]*x[:,2]
            if N>2:
                #% Third time derivative of y:
                #%      y3 = f'''(x)*x1.^3+f'(x)*x3 +3*f''(x)*x1*x2
                y[:,3]=fxder[:,2]*x[:,1]**3 + fxder[:,1]*x[:,3] + \
                    3*fxder[:,1]*x[:,1]*x[:,2]
                if N>3:
                    #% Fourth time derivative of y:
	                #%    y4 = f''''(x)*x1.^4+f'(x)*x4
 	                #%    +6*f'''(x)*x1^2*x2+f''(x)*(3*x2^2+4x1*x3)
  	                y[:,3]=fxder[:,3]*x[:,1]**4 + fxder[:,0]*x[:,4] + \
    	               6*fxder[:,2]*x[:,1]**2*x[:,2] + \
    	                   fxder[:,2]*(3*x[:,2]**2+4.*x[:,1]*x[:,3])
   	                if N>4:
   	                    warnings.warn('Transformation of derivatives of order>4 not supported in tranproc.')
    return y



class kde(object):
    """ Representation of a kernel-density estimate using Gaussian kernels.

    Parameters
    ----------
    dataset : (# of dims, # of data)-array
        datapoints to estimate from

    Members
    -------
    d : int
        number of dimensions
    n : int
        number of datapoints

    Methods
    -------
    kde.evaluate(points) : array
        evaluate the estimated pdf on a provided set of points
    kde(points) : array
        same as kde.evaluate(points)
    kde.integrate_gaussian(mean, cov) : float
        multiply pdf with a specified Gaussian and integrate over the whole domain
    kde.integrate_box_1d(low, high) : float
        integrate pdf (1D only) between two bounds
    kde.integrate_box(low_bounds, high_bounds) : float
        integrate pdf over a rectangular space between low_bounds and high_bounds
    kde.integrate_kde(other_kde) : float
        integrate two kernel density estimates multiplied together

   Internal Methods
   ----------------
    kde.covariance_factor() : float
        computes the coefficient that multiplies the data covariance matrix to
        obtain the kernel covariance matrix. Set this method to
        kde.scotts_factor or kde.silverman_factor (or subclass to provide your
        own). The default is scotts_factor.
    """

    def __init__(self, dataset,**kwds):
        self.kernel='gauss'
        self.hs = None
        self.hsmethod=None
        self.L2 = None
        self.__dict__.update(kwds)

        self.dataset = atleast_2d(dataset)
        self.d, self.n = self.dataset.shape


        self._compute_covariance()


    def evaluate(self, points):
        """Evaluate the estimated pdf on a set of points.

        Parameters
        ----------
        points : (# of dimensions, # of points)-array
            Alternatively, a (# of dimensions,) vector can be passed in and
            treated as a single point.

        Returns
        -------
        values : (# of points,)-array
            The values at each point.

        Raises
        ------
        ValueError if the dimensionality of the input points is different than
        the dimensionality of the KDE.
        """

        points = atleast_2d(points).astype(self.dataset.dtype)

        d, m = points.shape
        if d != self.d:
            if d == 1 and m == self.d:
                # points was passed in as a row vector
                points = reshape(points, (self.d, 1))
                m = 1
            else:
                msg = "points have dimension %s, dataset has dimension %s" % (d,
                    self.d)
                raise ValueError(msg)

        result = zeros((m,), points.dtype)

        if m >= self.n:
            # there are more points than data, so loop over data
            for i in range(self.n):
                diff = self.dataset[:,i,newaxis] - points
                tdiff = dot(self.inv_cov, diff)
                energy = sum(diff*tdiff,axis=0)/2.0
                result += exp(-energy)
        else:
            # loop over points
            for i in range(m):
                diff = self.dataset - points[:,i,newaxis]
                tdiff = dot(self.inv_cov, diff)
                energy = sum(diff*tdiff,axis=0)/2.0
                result[i] = sum(exp(-energy),axis=0)

        result /= self._norm_factor

        return result

    __call__ = evaluate

##function [f, hs,lambda]= kdefun(A,options,varargin)
##%KDEFUN  Kernel Density Estimator.
##%
##% CALL:  [f, hs] = kdefun(data,options,x1,x2,...,xd)
##%
##%   f      = kernel density estimate evaluated at x1,x2,...,xd.
##%   data   = data matrix, size N x D (D = # dimensions)
##%  options = kdeoptions-structure or cellvector of named parameters with
##%            corresponding values, see kdeoptset for details.
##%   x1,x2..= vectors/matrices defining the points to evaluate the density
##%
##%  KDEFUN gives a slow, but exact kernel density estimate evaluated at x1,x2,...,xd.
##%  Notice that densities close to normality appear to be the easiest for the kernel
##%  estimator to estimate and that the degree of estimation difficulty increases with
##%  skewness, kurtosis and multimodality.
##%
##%  If D > 1 KDE calculates quantile levels by integration. An
##%  alternative is to calculate them by ranking the kernel density
##%  estimate obtained at the points DATA  i.e. use the commands
##%
##%      f    = kde(data);
##%      r    = kdefun(data,[],num2cell(data,1));
##%      f.cl = qlevels2(r,f.PL);
##%
##%  The first is probably best when estimating the pdf and the latter is the
##%  easiest and most robust for multidimensional data when only a visualization
##%  of the data is needed.
##%
##%  For faster estimates try kdebin.
##%
##% Examples:
##%     data = rndray(1,500,1);
##%     x = linspace(sqrt(eps),5,55);
##%     plotnorm((data).^(.5)) % gives a straight line => L2 = 0.5 reasonable
##%     f = kdefun(data,{'L2',.5},x);
##%     plot(x,f,x,pdfray(x,1),'r')
##%
##% See also  kde, mkernel, kdebin
##
##% Reference:
##%  B. W. Silverman (1986)
##% 'Density estimation for statistics and data analysis'
##%  Chapman and Hall , pp 100-110
##%
##%  Wand, M.P. and Jones, M.C. (1995)
##% 'Kernel smoothing'
##%  Chapman and Hall, pp 43--45
##
##
##
##
##%Tested on: matlab 5.2
##% History:
##% revised pab Feb2004
##%  -options moved into a structure
##% revised pab Dec2003
##% -removed some code
##% revised pab 27.04.2001
##% - changed call from mkernel to mkernel2 (increased speed by 10%)
##% revised pab 01.01.2001
##% - added the possibility that L2 is a cellarray of parametric
##%   or non-parametric transformations (secret option)
##% revised pab 14.12.1999
##%  - fixed a small error in example in help header
##% revised pab 28.10.1999
##%  - added L2
##% revised pab 21.10.99
##%  - added alpha to input arguments
##%  - made it fully general for d dimensions
##%  - HS may be a smoothing matrix
##% revised pab 21.09.99
##%  - adapted from kdetools by Christian Beardah
##
##  defaultoptions = kdeoptset;
##% If just 'defaults' passed in, return the default options in g
##if ((nargin==1) && (nargout <= 1) &&  isequal(A,'defaults')),
##  f = defaultoptions;
##  return
##end
##error(nargchk(1,inf, nargin))
##
##[n, d]=size(A); % Find dimensions of A,
##               % n=number of data points,
##               % d=dimension of the data.
##if (nargin<2 || isempty(options))
##  options  = defaultoptions;
##else
##  switch lower(class(options))
##   case {'char','struct'},
##    options = kdeoptset(defaultoptions,options);
##   case {'cell'}
##
##      options = kdeoptset(defaultoptions,options{:});
##   otherwise
##    error('Invalid options')
##  end
##end
##kernel   = options.kernel;
##h        = options.hs;
##alpha    = options.alpha;
##L2       = options.L2;
##hsMethod = options.hsMethod;
##
##if isempty(h)
##  h=zeros(1,d);
##end
##
##L22 = cell(1,d);
##k3=[];
##if isempty(L2)
##  L2=ones(1,d); % default no transformation
##elseif iscell(L2)   % cellarray of non-parametric and parametric transformations
##  Nl2 = length(L2);
##  if ~(Nl2==1||Nl2==d), error('Wrong size of L2'), end
##  [L22{1:d}] = deal(L2{1:min(Nl2,d)});
##  L2 = ones(1,d); % default no transformation
##  for ix=1:d,
##    if length(L22{ix})>1,
##      k3=[k3 ix];       % Non-parametric transformation
##    else
##     L2(ix) = L22{ix};  % Parameter to the Box-Cox transformation
##    end
##  end
##elseif length(L2)==1
##  L2=L2(:,ones(1,d));
##end
##
##amin=min(A);
##if any((amin(L2~=1)<=0))  ,
##  error('DATA cannot be negative or zero when L2~=1')
##end
##
##
##nv=length(varargin);
##if nv<d,
##  error('some or all of the evaluation points x1,x2,...,xd is missing')
##end
##
##xsiz = size(varargin{1}); % remember size of input
##Nx   = prod(xsiz);
##X    = zeros(Nx,d);
##for ix=1:min(nv,d),
##  if (any(varargin{ix}(:)<=0) && (L2(ix)~=1)),
##    error('xi cannot be negative or zero when L2~=1')
##  end
##  X(:,ix)=varargin{ix}(:); % make sure it is a column vector
##end
##
##
##%new call
##lX = X; %zeros(Nx,d);
##lA = A; %zeros(size(A));
##
##k1 = find(L2==0); % logaritmic transformation
##if any(k1)
##  lA(:,k1)=log(A(:,k1));
##  lX(:,k1)=log(X(:,k1));
##end
##k2=find(L2~=0 & L2~=1); % power transformation
##if any(k2)
##  lA(:,k2)=sign(L2(ones(n,1),k2)).*A(:,k2).^L2(ones(n,1),k2);
##  lX(:,k2)=sign(L2(ones(Nx,1),k2)).*X(:,k2).^L2(ones(Nx,1),k2);
##end
##% Non-parametric transformation
##for ix = k3,
##  lA(:,ix) = tranproc(A(:,ix),L22{ix});
##  lX(:,ix) = tranproc(X(:,ix),L22{ix});
##end
##
##
##hsiz=size(h);
##if (min(hsiz)==1)||(d==1)
##  if max(hsiz)==1,
##    h=h*ones(1,d);
##  else
##    h=reshape(h,[1,d]); % make sure it has the correct dimension
##  end;
##  ind=find(h<=0);
##  if any(ind)    % If no value of h has been specified by the user then
##    h(ind)=feval(hsMethod,lA(:,ind),kernel); % calculate automatic values.
##  end
##  deth = prod(h);
##else  % fully general smoothing matrix
##  deth = det(h);
##  if deth<=0
##    error('bandwidth matrix h must be positive definit')
##  end
##end
##
##if alpha>0
##  Xn   = num2cell(lA,1);
##  opt1 = kdeoptset('kernel',kernel,'hs',h,'alpha',0,'L2',1);
##  f2   = kdefun(lA,opt1,Xn{:}); % get a pilot estimate by regular KDE (alpha=0)
##  g    = exp(sum(log(f2))/n);
##
##  lambda=(f2(:)/g).^(-alpha);
##else
##  lambda=ones(n,1);
##end
##
##
##
##
##
##f=zeros(Nx,1);
##if (min(hsiz)==1)||(d==1)
##  for ix=1:n,     % Sum over all data points
##    Avec=lA(ix,:);
##    Xnn=(lX-Avec(ones(Nx,1),:))./(h(ones(Nx,1),:) *lambda(ix));
##    f = f + mkernel2(Xnn,kernel)/lambda(ix)^d;
##  end
##else % fully general
##  h1=inv(h);
##  for ix=1:n,     % Sum over all data points
##    Avec=lA(ix,:);
##    Xnn=(lX-Avec(ones(Nx,1),:))*(h1/lambda(ix));
##    f = f + mkernel2(Xnn,kernel)/lambda(ix)^d;
##  end
##end
##f=f/(n*deth);
##
##% transforming back
##if any(k1), % L2=0 i.e. logaritmic transformation
##  for ix=k1
##    f=f./X(:,ix);
##  end
##  if any(max(abs(diff(f)))>10)
##    disp('Warning: Numerical problems may have occured due to the logaritmic')
##    disp('transformation. Check the KDE for spurious spikes')
##  end
##end
##if any(k2) % L2~=0 i.e. power transformation
##  for ix=k2
##    f=f.*(X(:,ix).^(L2(ix)-1))*L2(ix)*sign(L2(ix));
##  end
##  if any(max(abs(diff(f)))>10)
##    disp('Warning: Numerical problems may have occured due to the power')
##    disp('transformation. Check the KDE for spurious spikes')
##  end
##end
##if any(k3), % non-parametric transformation
##  oneC = ones(Nx,1);
##  for ix=k3
##    gn  = L22{ix};
##    %Gn  = fliplr(L22{ix});
##    %x0  = tranproc(lX(:,ix),Gn);
##    if any(isnan(X(:,ix))),
##      error('The transformation does not have a strictly positive derivative.')
##    end
##    hg1  = tranproc([X(:,ix) oneC],gn);
##    der1 = abs(hg1(:,2)); % dg(X)/dX = 1/(dG(Y)/dY)
##    % alternative 2
##    %pp  = smooth(Gn(:,1),Gn(:,2),1,[],1);
##    %dpp = diffpp(pp);
##    %der1 = 1./abs(ppval(dpp,f.x{ix}));
##    % Alternative 3
##    %pp  = smooth(gn(:,1),gn(:,2),1,[],1);
##    %dpp = diffpp(pp);
##    %%plot(hg1(:,1),der1-abs(ppval(dpp,x0)))
##    %der1 = abs(ppval(dpp,x0));
##    if any(der1<=0),
##      error('The transformation must have a strictly positive derivative')
##    end
##    f = f.*der1;
##  end
##  if any(max(abs(diff(f)))>10)
##    disp('Warning: Numerical problems may have occured due to the power')
##    disp('transformation. Check the KDE for spurious spikes')
##  end
##end
##
##f=reshape(f,xsiz); % restore original shape
##if nargout>1
##  hs=h;
##end
##
##
##
##
##
##
##
##
##
##
##function [z,c]=mkernel(varargin)
##%MKERNEL Multivariate Kernel Function.
##%
##% CALL:  z = mkernel(x1,x2,...,xd,kernel);
##%        z = mkernel(X,kernel);
##%
##%
##%   z      = kernel function values evaluated at x1,x2,...,xd
##%   x1,x2..= input arguments, vectors or matrices with common size
##% or
##%   X      = cellarray of vector/matrices with common size
##%            (i.e. X{1}=x1, X{2}=x2....)
##%
##%   kernel = 'epanechnikov'  - Epanechnikov kernel.
##%            'epa1'          - product of 1D Epanechnikov kernel.
##%            'biweight'      - Bi-weight kernel.
##%            'biw1'          - product of 1D Bi-weight kernel.
##%            'triweight'     - Tri-weight kernel.
##%            'triangular'    - Triangular kernel.
##%            'gaussian'      - Gaussian kernel
##%            'rectangular'   - Rectangular kernel.
##%            'laplace'       - Laplace kernel.
##%            'logistic'      - Logistic kernel.
##%
##%  Note that only the first 4 letters of the kernel name is needed.
##%
##% See also  kde, kdefun, kdebin
##
##%  Reference:
##%  B. W. Silverman (1986)
##% 'Density estimation for statistics and data analysis'
##%  Chapman and Hall, pp. 43, 76
##%
##%  Wand, M. P. and Jones, M. C. (1995)
##% 'Density estimation for statistics and data analysis'
##%  Chapman and Hall, pp 31, 103,  175
##
##%Tested on: matlab 5.3
##% History:
##% Revised pab sep2005
##% -replaced reference to kdefft with kdebin
##% revised pab aug2005
##% -Fixed some bugs
##% revised pab Dec2003
##% removed some old code
##% revised pab 27.04.2001
##% - removed some old calls
##%  revised pab 01.01.2001
##%  - speeded up tri3
##%  revised pab 01.12.1999
##%   - added four weight, sphere
##%   - made comparison smarter => faster execution for d>1
##%  revised pab 26.10.1999
##%   fixed normalization fault in epan
##% by pab 21.09.99
##%  added multivariate epan, biweight and triweight
##%
##% collected all knorm,kepan ... into this file
##% adapted from kdetools CB
##
##d=length(varargin)-1;
##kstr=varargin{d+1}; % kernel string
##if iscell(varargin{1})
##  X=varargin{1};
##  d=numel(X);
##else
##  X=varargin;
##end
##
##switch lower(kstr(1:4))
##  case {'sphe','epan','biwe','triw','four'}
##    switch lower(kstr(1:4))
##      case 'sphe', r=0;  %Sphere = rect for 1D
##      case 'epan', r=1;  %Multivariate Epanechnikov kernel.
##      case 'biwe', r=2;  %Multivariate Bi-weight Kernel
##      case 'triw', r=3;  %Multi variate Tri-weight Kernel
##      case 'four', r=4;  %Multi variate Four-weight Kernel
##        % as r -> infty, b -> infty => kernel -> Gaussian distribution
##    end
##    b=1;% radius of the kernel
##    b2=b^2;
##    s=X{1}.^2;
##    k=find(s<=b2);
##    z=zeros(size(s));
##    ix=2;
##    while (any(k) && (ix<=d)),
##      s(k)=s(k)+X{ix}(k).^2;
##      k1=(s(k)<=b2);
##      k=k(k1);
##      ix=ix+1;
##    end;
##    if any(k)
##      c=2^r*prod(1:r)*vsph(d,b)/prod((d+2):2:(d+2*r)); % normalizing constant
##      %c=beta(r+1,r+1)*vsph(d,b)*(2^(2*r)); % Wand and Jones pp 31
##      % the commented c above does note yield the right scaling
##      % for d>1
##      z(k)=((1-s(k)/b2).^r)/c;
##    end
##
##  case 'rect', % 1D product Rectangular Kernel
##    z=zeros(size(X{1}));
##    k=find(abs(X{1})<=1);
##    ix=2;
##    while (any(k) && (ix<=d)),
##      k1 =(abs(X{ix}(k))<=1);
##      k=k(k1);
##      ix=ix+1;
##    end
##    if any(k)
##      z(k)=(0.5^d);
##    end
##  case {'epa1','biw1','triw1','fou1'}
##    switch lower(kstr(1:4))
##      %case 'rect', r=0;  %rectangular
##      case 'epa1', r=1;  %1D product Epanechnikov kernel.
##      case 'biw1', r=2;  %1D product Bi-weight Kernel
##      case 'tri1', r=3;  %1D product Tri-weight Kernel
##      case 'fou1', r=4;  %1D product Four-weight Kernel
##    end
##    b=1;
##    b2=b^2;
##    b21=1/b2;
##    z=zeros(size(X{1}));
##    k=find(abs(X{1})<=b);
##    ix=2;
##    while (any(k) && (ix<=d)),
##      %for ix=2:d
##      k1 =(abs(X{ix}(k))<=b);
##      k  = k(k1);
##      ix=ix+1;
##    end
##    if any(k)
##      c=2^r*prod(1:r)*vsph(1,b)/prod((1+2):2:(1+2*r)); % normalizing constant
##      z(k) = (1-X{1}(k).^2*b21).^r;
##      for ix=2:d
##        z(k)=z(k).*(1-X{ix}(k).^2*b21).^r;
##      end;
##      z(k)=z(k)/c^d;
##    end
##  case 'tria',% 1D product Triangular Kernel
##    z=zeros(size(X{1}));
##    k=find(abs(X{1})<1);
##    ix=2;
##    while (any(k) && (ix<=d)),
##      %for ix=2:d
##      k1 =(abs(X{ix}(k))<1);
##      k  = k(k1);
##      ix=ix+1;
##    end
##    if any(k)
##      z(k) = (1-abs(X{1}(k)));
##      for ix=2:d
##        z(k)=z(k).*(1-abs(X{ix}(k)));
##      end
##    end
##   case {'norm','gaus'},% multivariate gaussian  Density Function.
##     s=X{1}.^2;
##     for ix=2:d
##       s=s+X{ix}.^2;
##     end;
##     z=(2*pi)^(-d/2)*exp(-0.5*s);
##   case 'lapl' % Laplace Kernel
##     z=0.5*exp(-abs(X{1}));
##     for ix=2:d
##       z=z.*0.5*exp(-abs(X{ix}));
##     end
##   case 'logi', % Logistic Kernel
##     z1=exp(X{1});
##     z=z1./(z1+1).^2;
##     for ix=2:d
##       z1=exp(X{ix});
##       z=z.*z1./(z1+1).^2;
##     end
##
##   otherwise, error('unknown kernel')
## end
##
##


