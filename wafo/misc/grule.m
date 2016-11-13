function [bp,wf]=grule(n,a,b,method)
%GRULE  Computes nodes and weights for Gauss-Legendre quadrature.
%
%CALL    [bp, wf]=grule(n,a,b,method);
% 
%   bp     = base points
%   wf     = weight factors
%   n      = number of base points (integrates a (2n-1)th order
%            polynomial exactly)
%   a, b   = lower and upper integration limit (default a = -1, b = 1)
%   method = Integer defining the method used to obtain the nodes and weights
%          0 Newton Raphson to find zeros of the Legendre polynomial (Fast)
%          1 Newton Raphson to find zeros of the Legendre polynomial (Faster)
%          2 Eigenvalue of the Gauss-Legendre matrix (Slow)
% 
%   The Gauss-Legendre quadrature integrates an integral of the form
%        1                 n
%       Int (f(x)) dx  =  Sum  wf_j*f(bp_j)
%       -1                j=1
%
% Example: % integral of exp(x) from a = 0 to b = 3 is: exp(3)-exp(0)
%  [x,w] = grule(11,0,3);
%   assert(sum(exp(x).*w), 19.0855369231877, 1e-10);
% 
% See also: quadg.


% Reference 
%   method 0: copied from grule.m in NIT toolbox, see ref [1] 
%   method 1: see ref [4]
%   method 2:  Adapted from Netlib routine gaussq.f see ref [2,3]
%
% [1] Davis and Rabinowitz (1975) 'Methods of Numerical Integration', page 365,
%     Academic Press.
% [2]  Golub, G. H. and Welsch, J. H. (1969)
% 'Calculation of Gaussian Quadrature Rules'
%  Mathematics of Computation, vol 23,page 221-230,
%
% [3]. Stroud and Secrest (1966), 'gaussian quadrature formulas', 
%      prentice-hall, Englewood cliffs, n.j.
% 
% 
% [4]William H. Press, Saul Teukolsky, 
% William T. Wetterling and Brian P. Flannery (1997)
% "Numerical recipes in Fortran 77", Vol. 1, pp 145
% "Numerical recipes in Fortran 90"




% Revised GKS 5 June 92
% Revised Per A. Brodtkorb pab@marin.ntnu.no 30.03.1999
% Revised pab 12 jan 2006
%  - added alternative and slightly faster version for large n
%  - added integration limits to input

error(nargchk(1,4,nargin))
if nargin<2 || isempty(a)
  a = -1;
end
if nargin < 3 || isempty(b)
  b = 1;
end
if nargin<4 || isempty(method)
  method = 0;
end




switch method
  case 2
    muzero = 2;
    ai = zeros(n,1);
    bi = zeros(n-1,1);
    i  = (1:n-1)';
    abi = i;
    bi(i) = abi./sqrt(4*abi.^2 - 1);
    
    [bp,wf] = gausscof(ai,bi,muzero);
    
  otherwise
    m=fix((n+1)/2);

    mm = 4*m-1;
    t  = (pi/(4*n+2))*(3:4:mm);
    nn = (1-(1-1/n)/(8*n*n));
    xo = nn*cos(t);
    
    if method==1
      
      % Compute the zeros of the N+1 Legendre Polynomial
      % using the recursion relation and the Newton-Raphson method
      
      xo = xo.';
      % Legendre-Gauss Polynomials
      L=zeros(m,3);
    
      % Derivative of LGP
      Lp=zeros(m,1);
      
      releps = 1e-15;
      MAXIT = 100;
      % Compute the zeros of the N+1 Legendre Polynomial
      % using the recursion relation and the Newton-Raphson method
      
      
      
      ix = 0;
      % Iterate until new points are uniformly within epsilon of old points
      k   = 1:m;
      k0  = 1;
      kp1 = 2;
      while (any(k) && ix < MAXIT)
        ix = ix+1;
        
        
        L(k,k0)  = 1;
        L(k,kp1) = xo(k);
        
        for jj=2:n
          km1 = k0;
          k0  = kp1;
          kp1 = mod(k0,3)+1;
          L(k,kp1)=( (2*jj-1)*xo(k).*L(k,k0)-(jj-1)*L(k,km1) )/jj;
        end
      
        Lp(k) = n*( L(k,k0)-xo(k).*L(k,kp1) )./(1-xo(k).^2);
        
        
        dx    = L(k,kp1)./Lp(k);
        xo(k) = xo(k)-dx;
        k     = k((abs(dx)> releps*abs(xo(k))));
        %L(length(k)+1:end,:) = [];
      end
      if any(k)
        warning('WAFO:GRULE','Too many iterations in grule!')
      end
      bp = -xo.';
      Lp = Lp.';
      wf =2./((1-bp.^2).*(Lp.^2));
    else
      % [bp,wf]=grule(n)
      %  This function computes Gauss base points and weight factors
      %  using the algorithm given by Davis and Rabinowitz in 'Methods
      %  of Numerical Integration', page 365, Academic Press, 1975.
      
      
      e1   = n*(n+1);
      iter = 2;
      
      for j=1:iter
        pkm1 = 1;
        pk   = xo;
        for k=2:n
          t1   = xo.*pk;
          pkp1 = t1-pkm1-(t1-pkm1)/k+t1;
          pkm1 = pk;
          pk   = pkp1;
        end
        den = 1.-xo.*xo;
        d1  = n*(pkm1-xo.*pk);
        dpn = d1./den;
        d2pn = (2.*xo.*dpn-e1.*pk)./den;
        d3pn = (4*xo.*d2pn+(2-e1).*dpn)./den;
        d4pn = (6*xo.*d3pn+(6-e1).*d2pn)./den;
        u  = pk./dpn; v=d2pn./dpn;
        h  = -u.*(1+(.5*u).*(v+u.*(v.*v-u.*d3pn./(3*dpn))));
        p  = pk+h.*(dpn+(.5*h).*(d2pn+(h/3).*(d3pn+.25*h.*d4pn)));
        dp = dpn+h.*(d2pn+(.5*h).*(d3pn+h.*d4pn/3));
        h  = h-p./dp;
        xo = xo+h;
      end
      bp=-xo-h;
      fx=d1-h.*e1.*(pk+(h/2).*(dpn+(h/3).*(...
        d2pn+(h/4).*(d3pn+(.2*h).*d4pn))));
      wf=2*(1-bp.^2)./(fx.^2);
    end
    if (m+m) > n,
      bp(m) = 0;
    end
    if ~((m+m) == n),
      m = m-1;
    end
    jj  = 1:m;
    n1j = (n+1-jj);
    bp(n1j) = -bp(jj);
    wf(n1j) = wf(jj);
end



if a~=-1 || b~=1
  % Linear map from[-1,1] to [a,b]
  bp  = (b-a)/2*(bp+1)+a;
  %bp = (a*(1-bp)+b*(1+bp))/2;
  wf  = wf*(b-a)/2;
end
% end
