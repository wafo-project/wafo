function [bp, wf]=gausscof(a,b,mu0,def)
%GAUSSCOF Computes coefficients for Gaussian quadrature from the Jacobi matrix
%
% CALL [bp wf]=gausscof(a,b,mu0)
%
%   bp,wf = base points (abscissas) and weight factors, respectively. Length N
%   a,b   = coefficients of the recurrence relation  
%
%           b p (x) = (x - a ) p   (x) - b   p   (x)       (nonmonic)
%            j j            j   j-1       j-1 j-2
%         or
%             p (x) = (x - a ) p   (x) - b   p   (x)        (monic)
%              j            j   j-1       j-1 j-2
%
%           for the set of (normalized) orthogonal polynomials  
%
%    mu0  = zero-th moment of the given polynomial's weight function w(x).
%           (default mu0=1)  
% 
%  (A monic polynomial is one which the term of highest degree has
%  coefficient of unity)
%   Since the polynomials are orthonormalized, the tridiagonal matrix
%   is guaranteed to be symmetric.
%
% Example % Compute nodes and weights for Gauss-Legendre quadrature:
% muzero = 2; n= 5;
% ai = zeros(n,1);
% bi = zeros(n-1,1);
% ix  = (1:n-1)';
% abi = ix;
% bi(ix) = abi./sqrt(4*abi.^2 - 1);
% [bp,wf] = gausscof(ai,bi,muzero);
% sum(bp.^2.*wf) % gaussq('x.^2',-1,1)
%
% See also orthog

if nargin<4||isempty(def)
   def=0;
end
if nargin<3||isempty(mu0)
  mu0=1;
end
n = length(a); 
if def==1, % monic polynomial
  b= sqrt(b); 
end

if 1,% faster call
  A0 = diag(a);
  A0( 2 : n+1 : n*(n-1) ) = b;
  A0( n+1 : n+1 : n^2-1 ) = b;
else
  A0 = diag(a) + diag(b,1) + diag(b,-1);
end

[v, d] = eig(A0);

if 1,
  [bp, ii] = sort( diag(d).' );
  wf = (sqrt(mu0).*v(1,ii)).^2;
else % save some valuable time by not sorting
  bp = diag(d).' ;
  wf = mu0.*v(1,:).^2;
end


