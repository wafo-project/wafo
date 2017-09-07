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
% Example:
%   [q, r] = polydiv([3, 6, 9, 9], [1, 2, 3]);
%   assert(q, [3, 0]);
%   assert(r, [0, 0, 0, 9]);
%
% See also POLYADD, POLYSUB, POLYMUL, POLYPOW, DECONV.

%   Author:      Peter J. Acklam
%   Time-stamp:  1998-06-22 20:36:14
%   E-mail:      jacklam@math.uio.no
%   WWW URL:     http://www.math.uio.no/~jacklam

   %error(nargchk(2, 2, nargin));
   narginchk(2,2)
   [q, r] = deconv(a, b);

   
%!test
%! [b, r] = deconv ([3, 6, 9, 9], [1, 2, 3]);
%! assert (b, [3, 0]);
%! assert (r, [0, 0, 0, 9]);

%!test
%! [b, r] = polydiv ([3, 6], [1, 2, 3]);
%! assert (b, 0);
%! assert (r, [3, 6]);

%!test
%! [b, r] = polydiv ([3, 6], [1; 2; 3]);
%! assert (b, 0);
%! assert (r, [3, 6]);

%!test
%! [b,r] = polydiv ([3; 6], [1; 2; 3]);
%! assert (b, 0);
%! assert (r, [3; 6]);

%!test
%! [b, r] = polydiv ([3; 6], [1, 2, 3]);
%! assert (b, 0);
%! assert (r, [3; 6]);

%!assert (polydiv ((1:3)',[1, 1]), [1; 1])

%!error [b, r] = polydiv ([3, 6], [1, 2; 3, 4])
%!error [b, r] = polydiv ([3, 6; 1, 2], [1, 2, 3])
