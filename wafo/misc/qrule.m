function [bp,wf]=qrule(n,wfun,alpha,beta)
%QRULE Computes nodes and weights for Gaussian quadratures.   
%
%CALL:  [bp,wf]=qrule(n,wfun,alpha,beta)
%  
%  bp = base points (abscissas)
%  wf = weight factors
%  n  = number of base points (abscissas) (integrates a (2n-1)th order
%       polynomial exactly)
%wfun = integer defining the weight function, p(x). (default wfun = 1)     
%     1,11,21: p(x) = 1                       a =-1,   b = 1   Gauss-Legendre 
%     2,12   : p(x) = exp(-x^2)               a =-inf, b = inf Hermite
%     3,13   : p(x) = x^alpha*exp(-x)         a = 0,   b = inf Laguerre
%     4,14   : p(x) = (x-a)^alpha*(b-x)^beta  a =-1,   b = 1 Jacobi 
%     5      : p(x) = 1/sqrt((x-a)*(b-x)),    a =-1,   b = 1 Chebyshev 1'st kind
%     6      : p(x) = sqrt((x-a)*(b-x)),      a =-1,   b = 1 Chebyshev 2'nd kind
%     7      : p(x) = sqrt((x-a)/(b-x)),      a = 0,   b = 1
%     8      : p(x) = 1/sqrt(b-x),            a = 0,   b = 1
%     9      : p(x) = sqrt(b-x),              a = 0,   b = 1
%
%  The Gaussian Quadrature integrates a (2n-1)th order
%  polynomial exactly and the integral is of the form
%           b                         n
%          Int ( p(x)* F(x) ) dx  =  Sum ( wf_j* F( bp_j ) )
%           a                        j=1		          
% where p(x) is the weight function. 
% For Jacobi and Laguerre: alpha, beta >-1 (default alpha=beta=0)
% 
% Examples:
%   [bp,wf] = qrule(10);
%   assert(sum(bp.^2.*wf), 0.666666666666666, 1e-10) % integral of x^2 from a = -1 to b = 1
%   [bp,wf] = qrule(10,2);
%   assert(sum(bp.^2.*wf), 0.886226925452758, 1e-10)  % integral of exp(-x.^2)*x.^2 from a = -inf to b = inf
%   [bp,wf] = qrule(10,4,1,2);
%   assert(sum(bp.*wf), 0.266666666666668, 1e-10)   % integral of (x+1)*(1-x)^2 from  a = -1 to b = 1
% 
%  See also  gaussq

% Reference 
%   wfun 1 : copied from grule.m in NIT toolbox, see ref [2] 
%   wfun 2-4, 11: see ref [5]
%   wfun 12-14, 21: Adapted from Netlib routine gaussq.f see ref [1,3]
%   wfun 5-9: see ref [4]
%
% [1]  Golub, G. H. and Welsch, J. H. (1969)
% 'Calculation of Gaussian Quadrature Rules'
%  Mathematics of Computation, vol 23,page 221-230,
%
% [2] Davis and Rabinowitz (1975) 'Methods of Numerical Integration', page 365,
%     Academic Press.
%
% [3]. Stroud and Secrest (1966), 'gaussian quadrature formulas', 
%      prentice-hall, Englewood cliffs, n.j.
% 
% [4] Abromowitz and Stegun (1954) ''
%
% [5]William H. Press, Saul Teukolsky, 
% William T. Wetterling and Brian P. Flannery (1997)
% "Numerical recipes in Fortran 77", Vol. 1, pp 145
% "Numerical recipes in Fortran 90"


% By Bryce Gardner, Purdue University, Spring 1993.
% Modified by Per A. Brodtkorb 19.02.99 pab@marin.ntnu.no
% to compute other quadratures  than the default
% Revised pab march2006
% Restructured the contents + moved code to grule, hrule, lrule and jrule.

if nargin<4||isempty(beta),
 beta=0; 
end

if nargin<3||isempty(alpha),
  alpha=0; 
end
if alpha<=-1 || beta <=-1,
  error('alpha and beta must be greater than -1')
end

if nargin<2||isempty(wfun),
  wfun=1; 
end	


switch mod(wfun,10), %
  case 1, % Gauss-Legendre
    [bp,wf] = grule(n,-1,1, floor(wfun/10));
  case 2, % Hermite 
    [bp,wf] = hrule(n,floor(wfun/10));
  case 3, % Generalized Laguerre
    [bp,wf] = lrule(n,floor(wfun/10));
  case 4, % Gauss-Jacobi
    [bp,wf] = jrule(n,alpha,beta,floor(wfun/10));
 case 5, % p(x)=1/sqrt((x-a)*(b-x)), a=-1 and b=1 (default) 
  jj  = 1:n;
  wf = ones(1,n) * pi / n;
  bp = cos( (2*jj-1)*pi / (2*n) );

 case 6, % p(x)=sqrt((x-a)*(b-x)),   a=-1 and b=1
  jj = 1:n;
  xj = jj * pi / (n+1);
  wf = pi / (n+1) * sin( xj ).^2;
  bp = cos( xj );

 case 7, % p(x)=sqrt((x-a)/(b-x)),   a=0 and b=1
    jj = 1:n;
    xj = (jj-0.5)*pi / (2*n+1);
    bp = cos( xj ).^2;
    wf = 2*pi.*bp/(2*n+1) ;

 case 8, % p(x)=1/sqrt(b-x),         a=0 and b=1
   [bp1, wf1] = grule(2*n);
   k = find(0<=bp1);
   wf = 2*wf1(k);
   bp = 1-bp1(k).^2;
   
 case 9, % p(x)=sqrt(b-x),           a=0 and b=1
   [bp1, wf1] = grule(2*n+1);
   k = find(0<bp1);
   wf = 2*bp1(k).^2.*wf1(k);
   bp = 1-bp1(k).^2;  
otherwise, 
  error('unknown weight function')
end

