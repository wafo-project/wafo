#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      pab
#
# Created:     30.12.2008
# Copyright:   (c) pab 2008
# Licence:     <your licence>
#-------------------------------------------------------------------------------
#!/usr/bin/env python
from __future__ import division
import numpy as np


class PPform(object):
    """The ppform of the piecewise polynomials is given in terms of coefficients
    and breaks.  The polynomial in the ith interval is
    x_{i} <= x < x_{i+1}

    S_i = sum(coefs[m,i]*(x-breaks[i])^(k-m), m=0..k)
    where k is the degree of the polynomial.

    Example
    -------
    >>> coef = np.array([[1,1]]) # unit step function
    >>> coef = np.array([[1,1],[0,1]]) # linear from 0 to 2
    >>> coef = np.array([[1,1],[1,1],[0,2]]) # linear from 0 to 2
    >>> breaks = [0,1,2]
    >>> pp = PPform(coef,breaks)
    >>> x = linspace(-1,3)
    >>> plot(x,pp(x))
    """
    def __init__(self, coeffs, breaks, fill=0.0, sort=False):
        self.coeffs = np.asarray(coeffs)
        if sort:
            self.breaks = np.sort(breaks)
        else:
            self.breaks = np.asarray(breaks)
        self.K = self.coeffs.shape[0]
        self.fill = fill
        self.a = self.breaks[0]
        self.b = self.breaks[-1]

    def __call__(self, xnew):
        saveshape = np.shape(xnew)
        xnew = np.ravel(xnew)
        res = np.empty_like(xnew)
        mask = (xnew >= self.a) & (xnew <= self.b)
        res[~mask] = self.fill
        xx = xnew.compress(mask)
        indxs = np.searchsorted(self.breaks, xx)-1
        indxs = indxs.clip(0,len(self.breaks))
        pp = self.coeffs
        diff = xx - self.breaks.take(indxs)
        V = np.vander(diff,N=self.K)
        # values = np.diag(dot(V,pp[:,indxs]))
        dot = np.dot
        values = np.array([dot(V[k,:],pp[:,indxs[k]]) for k in xrange(len(xx))])
        res[mask] = values
        res.shape = saveshape
        return res

    def derivative(self):
        """
        Return first derivative of the piecewise polynomial
        """
        import polynomial as pl
        cof = pl.polyder(self.coeffs)
        brks = self.breaks.copy()
        return PPform(cof,brks,fill=self.fill)
##        k = self.K
##        if k==1:
##            return PPform(np.zeros_like(self.coeffs[0,...]), brks, fill=self.fill)
##
##        ix  = np.arange(k-1,0,-1)
##        cof = ix[...,np.newaxis]*self.coeffs[:k-1,...]
##        return PPform(cof, brks, fill=self.fill)

    def integrate(self):
        """
        Return the indefinite integral of the piecewise polynomial
        """
        k = self.K
        ix = np.arange(k,0,-1)
        d = self.coeffs.shape[-1]
        cof = self.coeffs/ix
        if (d<2):
            return PPform(np.vstack((cof,np.zeros(d))),self.breaks,fill=self.fill)

        # evaluate each integrated polynomial at the
        # right endpoint of its interval
        xs = diff(self.breaks)
        index = np.arange(d-1)
##        if (d>1):
##            temp = xs(ones(d,1),:)
##            xs=temp(:).'
##            temp = 1+d*ones((1,d))*index+[-d:-1]*ones(1,l-1)
##            index=temp(:);

        vv = xs*coefs[0,index];
        for i in xrange(1,k):
            vv = xs*(vv + coefs[i,index])

        if (d>1):
            junk = zeros(d,l-1);
            junk[:]=vv
            last=(cumsum([ifa, junk]));
        else:
            last=cumsum([ifa,vv]);

        return PPform(coeffs,breaks,fill=self.fill)



##    def fromspline(cls, xk, cvals, order, fill=0.0):
##        N = len(xk)-1
##        sivals = np.empty((order+1,N), dtype=float)
##        for m in xrange(order,-1,-1):
##            fact = spec.gamma(m+1)
##            res = _fitpack._bspleval(xk[:-1], xk, cvals, order, m)
##            res /= fact
##            sivals[order-m,:] = res
##        return cls(sivals, xk, fill=fill)


class CSspline(object):
    """
    Calculates a smoothing spline.

    Parameters
    ----------
    x : array-like
        x-coordinates of data. (vector)
    y : array-like
        y-coordinates of data. (vector or matrix)
    p : real scalar
        smoothing parameter between 0 and 1:
        0 -> LS-straight line
        1 -> cubic spline interpolant
    xi : array-like
        the x-coordinates in which to calculate the smoothed function.
    def : bool
        if 0 regular smoothing spline (default)
        if 1 a smoothing spline with a constraint on the ends to
        ensure linear extrapolation outside the range of the data
    v : array-like
        variance of each y(i) (default  ones(length(X),1))

    Returns
    -------
    yy_or_pp : ndarray or ppform
        if xx is given, return the calculated y-coordinates of the smoothed function.
        If xx is not given, return pp-form of the spline.

    Given the approximate values

        y(i) = g(x(i))+e(i)

    of some smooth function, g, where e(i) is the error. SMOOTH tries to
    recover g from y by constructing a function, f, which  minimizes

      p * sum (Y(i) - f(X(i)))^2/d2(i)  +  (1-p) * int (f'')^2

    The call  pp = smooth(x,y,p)  gives the pp-form of the spline,
    for use with PPVAL.

    Example
    -------
    >>> x = linspace(0,1);
    >>> y = exp(x)+1e-1*randn(size(x));
    >>> pp = smooth(x,y,.9);
    >>> plot(x,y,x,smooth(x,y,.99,x,0,0.01),'g',x,ppval(pp,x),'k',x,exp(x),'r')

    See also
    --------
    lc2tr, dat2tr, ppval


    References
    ----------
    Carl de Boor (1978)
    'Practical Guide to Splines'
    Springer Verlag
    Uses EqXIV.6--9, pp 239
    """
    def __init__(self,xx,yy,p=None,xi=None,LinExtrap=0,d2=1):
        pass
        x,y = np.atleast_1d(xx,yy)
        x = x.ravel()
        dx = diff(x)
        mustSort = np.any(dx<0)
        if mustSort:
            ind = x.argsort()
            x = x[ind]
            y = y[...,ind]
            dx = np.diff(x)

        n = len(x)

        ndy = y.ndim
        szy = y.shape

        #nd = prod(szy(1:end-1));
        ny = szy[-1]


##        if n<2:
##            raise ValueError('There must be >=2 data points.')
##        elif any(dx<=0):
##            raise ValueError('Two consecutive values in x can not be equal.')
##        elif n!=ny:
##           raise ValueError('x and y must have the same length.')
##
##
##
##        dydx = np.diff(y)/dx
##
##        if (n==2) : #% straight line
##            coefs = [dydx.flatten, y[0,:]];
##        else:
##            if LinExtrap==2 and n==3:
##                p = 0;  # Force LS-fit
##
##            dx1=1./dx;
##
##            u = computeU();
##
##            zrs = zeros(1,nd);
##            if p<1:
##                ai = yi-6*(1-p)*D*diff([zrs;diff([zrs;u;zrs]).*dx1(:,idx);zrs]); % faster than yi-6*(1-p)*Q*u
##            else:
##                ai = yi;
##            #end
##            # The piecewise polynominals are written as
##            # fi=ai+bi*(x-xi)+ci*(x-xi)^2+di*(x-xi)^3
##            # where the derivatives in the knots according to Carl de Boor are:
##            #    ddfi  = 6*p*[0;u] = 2*ci;
##            #    dddfi = 2*diff([ci;0])./dx = 6*di;
##            #    dfi   = diff(ai)./dx-(ci+di.*dx).*dx = bi;
##
##            ci = [zrs;3*p*u];
##
##
##            if LinExtrap==2 && p!=0 && n>3, %Forcing linear extrapolation in the ends
##                ci([2,  end],:) = 0;
##            #% New call
##            #% fixing the coefficients so that we have continous
##            #% derivatives everywhere
##            ai(1,:) = -(ai(3,:)-ai(2,:))*dx(1)/dx(2) +ai(2,:)+ ci(3,:)*dx(1)*dx(2)/3;
##            ai(n,:) = (ai(n-1,:)-ai(n-2,:))*dx(n-1)/dx(n-2) +ai(n-1,:)+ ci(n-2,:)*dx(n-2)*dx(n-1)/3;
##          end
##
##
##          di    = (diff([ci;zrs]).*dx1(:,idx)/3);
##          bi    = (diff(ai).*dx1(:,idx)-(ci+di.*dx(:,idx)).*dx(:,idx));
##          ai(n:end,:) = [];
##          if nd>1
##            di = di.';
##            ci = ci.';
##            ai = ai.';
##          end
##          if ~any(di)
##            if ~any(ci)
##              coefs = [bi(:) ai(:)];
##            else
##              coefs = [ci(:) bi(:) ai(:)];
##            end
##          else
##            coefs = [di(:) ci(:) bi(:) ai(:)];
##          end
##        end
##
##        pp = mkpp(xi,coefs,szy(1:end-1));
##
##        if ~any(LinExtrap==[0 2])
##          % Extrapolation strategy
##          pp = extrapolate(pp);
##        end
##
##
##        if (nargin<4)||(isempty(xx)),
##          yy = pp;
##        else
##          yy = ppval(pp,xx);
##        end
##
##        %% Nested functions
##          function u = computeU()
##            if isempty(p) || p~=0
##              R = spdiags([dx(2:n-1) 2*(dx(1:n-2)+dx(2:n-1)) dx(1:n-2)],-1:1,n-2,n-2);
##            end
##            if isempty(p) || p<1
##              Q = spdiags([dx1(1:n-2) -(dx1(1:n-2)+dx1(2:n-1)) dx1(2:n-1)],0:-1:-2,n,n-2);
##              D = spdiags(d2,0,n,n);  % The variance
##              QDQ = Q.'*D*Q;
##              if isempty(p) || p<0
##                % Crude estimate
##                p = 1/(1+trace(QDQ)/(100*trace(R).^2));
##              end
##              if p==0
##                QQ = (6*(1-p))*(QDQ);
##              else
##                QQ = (6*(1-p))*(QDQ)+p*R;
##              end
##              %clear QDQ
##            else
##              QQ = R;
##            end
##            % Make sure Matlab uses symmetric matrix solver
##            u  = 2*((QQ+QQ.')\diff(dydx));  % faster than u=QQ\(Q'*yi);
##          end
##        end % function smooth
##
##    %% Subfunctions
##    function g = extrapolate(pp)
##    %EXTRAPOLATE Extrapolate a 1D PP linearly outside its basic interval
##
##
##    maxOrder = 2;
##
##    [breaks,coefs,pieces,order,dim]=unmkpp(pp);
##    if order<=maxOrder
##      g = pp;
##      return
##    end
##
##
##
##    %Add new breaks beyond each end
##    breaks2add = breaks([1,end]) + [-1,1];
##    newbreaks = [breaks2add(1),breaks, breaks2add(2)];
##
##
##    dx = newbreaks([1 end-1]) - breaks([1 end-1]);
##
##    ifl  = [1 pieces]; % index to first and last polynomial piece.
##    if dim>1 % repeat each point dim  times if necessary
##      dx = repmat(dx,dim,1);
##      ifl = repmat(dim*ifl,dim,1)+repmat((1-dim:0).',1,2);
##    end
##    dx = dx(:);
##    nx = length(dx);
##
##    % Get coefficients for the new last polynomial piece (aN)
##    % by just relocate the previous last polynomial and
##    % then set all terms of order > maxOrder to zero
##    ix  = dim+1:nx;
##    aN  = coefs(ifl(ix),:);
##    dxN = dx(ix);
##    % Relocate last polynomial using Horner's algorithm
##    aN = polyreloc(aN,dxN);
##
##    %set to zero all terms of order > maxOrder in first polynomial
##    aN(:,1:order-maxOrder) = 0;
##
##    %Get the coefficients for the new first piece (a1)
##    % by first setting all terms of order > maxOrder to zero and then
##    % relocate the polynomial.
##
##    ix = 1:dim;
##
##    %Set to zero all terms of order > maxOrder, i.e., not using them
##    a1 = coefs(ifl(ix),maxOrder+1:end);
##    % alternatively
##    % a1 = coefs(lu(ix),:);
##    % a1(:,1:order-maxOrder) = [];
##
##    dx1 = dx(ix);
##
##    % Relocate first polynomial using Horner's algorithm
##    a1 = polyreloc(a1,dx1);
##    a1 = [zeros(dim,order-maxOrder),a1];
##
##    %
##    % % fixing the coefficients so that we have continous
##    % % derivatives everywhere
##    % a1=-(ai(2)-ai(1))*dx(1)/dx(2) +ai(1)+ ci(3)*dx(1)*dx(2)/3;
##    % an=(ai(n-2)-ai(n-3))*dx(n-1)/dx(n-2) +ai(n-2)+ ci(n-2)*dx(n-2)*dx(n-1)/3;
##    % ai=[a1;ai; an];
##
##
##    newcoefs = [ a1; coefs; aN];
##
##
##    g = mkpp(newbreaks,newcoefs,dim);
##
##
##    end % function extrapolate
##
##
##
##def polyreloc(p,dx):
##    """
##    POLYRELOC Relocate polynomial.
##
##    CALL  R = POLYRELOC( P, dX)
##
##    R  = vector/matrix of relocated polynomial coefficients.
##    P  = vector/matrix of polynomial coefficients to relocate
##    dx = distance to relocate P.
##
##    POLYRELOC relocates the polynomial P by "moving" it dX
##    units to the left along the x-axis. So R is
##    relative to the point (-dX,0) as P is relative to the point (0,0).
##
##    P is a matrix of row vectors of coefficients in decreasing order.
##    """
##
##
##    r = np.atleast_2d(p);
##    n = r.shape[1]
##    #% Relocate polynomial using Horner's algorithm
##    for ii in range(n,1,-1):
##        for i in range(1,ii):
##            r[:,i] = dx*r[:,i-1]+r[:,i]
##    return r
##
##
##

def main():
    from scipy import interpolate
    import matplotlib.pyplot as plt
    t = np.arange(0,1.1,.1)
    x = np.sin(2*np.pi*t)
    y = np.cos(2*np.pi*t)
    tck1,u = interpolate.splprep([x,y],s=0)
    tck = interpolate.splmake(x, y, order=3, kind='smoothest', conds=None)
    pp = interpolate.ppform.fromspline(*tck)



if __name__=='__main__':
    main()
    import polynomial as pl
    coef = np.array([[1,1],[0,0]]) # linear from 0 to 2

    coef = np.array([[1,1],[1,1],[0,2]]) # quadratic from 0 to 1 and 1 to 2.
    dc = pl.polyder(coef,1)
    c2 = pl.polyint(dc,1)
    breaks = [0,1,2]
    pp = PPform(coef,breaks)
    pp(0.5)
    pp(1)
    pp(1.5)
    dpp = pp.derivative()
    import pylab as plb
    x = plb.linspace(-1,3)
    plb.plot(x,pp(x),x,dpp(x),'.')
    plb.show()
