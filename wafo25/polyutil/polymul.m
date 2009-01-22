function r = polymul(p, q)
%POLYMUL Multiply polynomials.
%
%   R = POLYMUL(P, Q) multiplies the polynomials whose coefficients are
%   the elements of the vectors P and Q.
%
%   POLYMUL is essentially a call to CONV.
%
%   See also POLYADD, POLYSUB, POLYDIV, POLYPOW, CONV.

%   Author:      Peter J. Acklam
%   Time-stamp:  1998-06-22 20:36:03
%   E-mail:      jacklam@math.uio.no
%   WWW URL:     http://www.math.uio.no/~jacklam

   error(nargchk(2, 2, nargin));

   r = conv(p, q);
