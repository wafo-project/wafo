function r = polypow(p, n)
%POLYPOW Polynomial raised to the Nth power.
%
%   R = POLYPOW(P, N) returns the polynomial P raised to the Nth power.
%
%   P is a vector of coefficients in decreasing order.
%
%   Example
%   assert(polypow([1,1], 2), [1,2,1], 1e-12);
%   assert(polypow([1,1], 3), [1,3,3,1], 1e-12);
%
%   See also polyadd, polysub, polymul, polydiv.

% History
% revised pab 10.01.2001
%  -fixed a bug: changed  r = conv(p, p); to  r = conv(r, p);
%   Author:      Peter J. Acklam
%   Time-stamp:  2000-08-08 19:27:25
%   E-mail:      jacklam@math.uio.no
%   WWW URL:     http://www.math.uio.no/~jacklam

   error(nargchk(2, 2, nargin));

   if any(size(n) ~= 1) || ~isreal(n) || (n < 0)
      error('N must be a non-negative integer.');
   end

   if n == 0
      r = 1;
   else
      r = p;
      for i = 2:n
         r = conv(r, p);
      end
   end
