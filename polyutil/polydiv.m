function [q, r] = polydiv(a, b)
%POLYDIV Divide polynomials.
%
%   [Q, R] = POLYDIV(A, B) returns the polynomial A devided by the
%   polynomial B which results in a quotient Q and a remainder R, i.e., 
%   
%     A/B = Q + R/B
%
%   POLYDIV is essentially a call to DECONV.
%
%   See also POLYADD, POLYSUB, POLYMUL, POLYPOW, DECONV.

%   Author:      Peter J. Acklam
%   Time-stamp:  1998-06-22 20:36:14
%   E-mail:      jacklam@math.uio.no
%   WWW URL:     http://www.math.uio.no/~jacklam

   error(nargchk(2, 2, nargin));

   [q, r] = deconv(a, b);
