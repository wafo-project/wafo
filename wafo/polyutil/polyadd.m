function r = polyadd(p, q)
%POLYADD Add polynomials.
%
%   R = POLYADD(P, Q) adds the polynomials whose coefficients are the
%   elements of the vectors P and Q.
%
%  Example
%  assert(polyadd([1,2,3],[1,2]), [1,3,5], 1e-10);
% 
%  See also POLYADD, POLYSUB, POLYMUL, POLYDIV, POLYPOW.

%   Author:      Peter J. Acklam
%   Time-stamp:  1998-06-22 20:36:21
%   E-mail:      jacklam@math.uio.no
%   WWW URL:     http://www.math.uio.no/~jacklam

   error(nargchk(2, 2, nargin));

   m = length(p);
   n = length(q);
   l = max(m, n);

   r = zeros(1, l);                     % Initialize output.
   r(l-m+1:l) = p;                      % Insert first polynomial.
   r(l-n+1:l) = r(l-n+1:l) + q;         % Add second polynomial

   r = polytrim(r);                     % Trim the result.
