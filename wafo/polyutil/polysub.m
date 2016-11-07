function r = polysub(p, q)
%POLYSUB Subtract polynomials.
%
%   R = POLYSUB(P, Q) subtract the polynomial Q from the polynomial P.
%
%   P and Q are vectors of coefficients in decreasing order.
%
% Example:
%   assert(polysub([1,2,3],[1,2]), [1,1,1], 1e-10);
%
% See also POLYADD, POLYMUL, POLYDIV, POLYPOW.

%   Author:      Peter J. Acklam
%   Time-stamp:  1998-06-22 20:32:25
%   E-mail:      jacklam@math.uio.no
%   WWW URL:     http://www.math.uio.no/~jacklam

   error(nargchk(2, 2, nargin));

   r = polyadd(p, -q);
