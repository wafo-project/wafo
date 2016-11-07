function [bp,wf]=lrule(n,alpha,method)
%LRULE Computes nodes and weights for generalized Laguerre  quadrature.
%    
% CALL    [bp, wf]=lrule(n,alpha);
%
%  bp     = base points
%  wf     = weight factors
%  n      = number of base points (integrates a (2n-1)th order
%           polynomial exactly)
%  alpha  = exponent of polynomial term
%  method = Integer defining the method used to obtain the nodes and weights
%          0 Newton Raphson to find zeros of the generalized Laguerre polynomial (Fast)
%          1 Eigenvalue of the Laguerre matrix (Slow) default
%
%  The Laguerre Quadrature integrates an integral of the form
%       inf                               n
%      Int (x^alpha exp(-x) F(x)) dx  =  Sum  wf(j)*F( bp(j) )
%        0                               j=1    
%
% Example
%  [x, w] = lrule(10,2);
%  assert(sum(x.*w), 6, 1e-10)
%  [x,w]= lrule(10,2,0);
%  assert(sum(x.*w), 6, 1e-10)
% 
%  See also quadl

% Reference 
%   method 0: see ref [4]
%   method 1:  Adapted from Netlib routine gaussq.f see ref [2,3]
%
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

% Adapted from Numerical recipes
% Revised by Per A. Brodtkorb 18.02.99 pab@marin.ntnu.no
% revised pab jan2006

if nargin < 3 || isempty(method)
  method = 1;
end
if nargin<2||isempty(alpha)
  alpha=0;
end
if alpha<=-1
  error('alpha must be greater than -1')
end

switch method
  case 0,
    
    MAXIT  = 10;
    releps = 3e-13;
%    releps = 1e-14;
    C=[9.084064e-01;...
      5.214976e-02;...
      2.579930e-03;...
      3.986126e-03];


    % Initial approximations to the roots go into z.
    anu = 4.0*n+2.0*alpha+2.0;
    rhs = (4*n-1:-4:3)*pi/anu;
    r3  = rhs.^(1/3);
    r2  = r3.^2;
    theta = r3.*(C(1)+r2.*(C(2)+r2.*(C(3)+r2.*C(4))));
    z = anu.*cos(theta).^2;
    dz = z;

    
    L  = zeros(3,length(z));
    Lp = zeros(1,length(z));
    pp = Lp;
    k0  = 1;
    kp1 = 2;
    k = 1:length(z);
    for its=1:MAXIT
      %Newton�s method carried out simultaneously on the roots.
      L(k0,k)  = 0;
      L(kp1,k) = 1;

      for jj=1:n
        %Loop up the recurrence relation to get the Laguerre
        %polynomials evaluated at z.
        km1 = k0;
        k0 = kp1;
        kp1 = mod(kp1,3)+1;

        L(kp1,k) =((2*jj-1+alpha-z(k)).*L(k0,k)-(jj-1+alpha)*L(km1,k))/jj;
      end
      %L now contains the desired Laguerre polynomials.
      %We next compute pp, the derivatives with a standard
      % relation involving the polynomials of one lower order.
      
      Lp(k) = L(k0,k);
      pp(k) = (n*L(kp1,k)-(n+alpha)*Lp(k))./z(k);
      
      dz(k) = L(kp1,k)./pp(k);
      z(k)  = z(k)-dz(k); % Newton�s formula.
      %k = find((abs(dz) > releps.*z));
      unfinished = find((abs(dz) > releps.*z));
      k = unfinished;
      if not(any(unfinished))
        break
      end
    end
    if any(unfinished)
      warning('WAFO:lrule','Too many iterations in lrule');
    end
    bp = z; %Store the root and the weight.
    wf = -exp(gammaln(alpha+n)-gammaln(n))./(pp.*n.*Lp);
  otherwise
    % this procedure uses the coefficients a(j), b(j) of the
    %      recurrence relation
    %
    %           b p (x) = (x - a ) p   (x) - b   p   (x)
    %            j j            j   j-1       j-1 j-2
    %
    %      for the various classical (normalized) orthogonal polynomials,
    %      and the zero-th moment
    %
    %           muzero = integral w(x) dx
    %
    %      of the given polynomial's weight function w(x).  since the
    %      polynomials are orthonormalized, the tridiagonal matrix is
    %      guaranteed to be symmetric.
    %
    %         the input parameter alpha is used only for laguerre and
    %      jacobi polynomials, and the parameter beta is used only for
    %      jacobi polynomials.  the laguerre and jacobi polynomials
    %      require the gamma function.



    a    = zeros(n,1);
    ii    = (1:n-1)';
    a(ii) = 2 .* ii - 1 + alpha;
    a(n)  = 2*n-1+alpha;
    b     = sqrt( ii .* (ii + alpha) );
    muzero =  gamma(alpha+1);
    
    [bp,wf] = gausscof(a,b,muzero);
    
 end
