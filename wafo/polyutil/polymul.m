function r = polymul(p, q)
%POLYMUL Multiply polynomials.
%
%   R = POLYMUL(P, Q) multiplies the polynomials whose coefficients are
%   the elements of the vectors P and Q.
%
%   POLYMUL is essentially a call to CONV.
%
% Example:
%   x = ones (3,1);
%   y = ones (1,3);
%   b = 2;
%   c = 3;
%   assert (polymul (x, x), [1; 2; 3; 2; 1]);
%   assert (polymul (y, y), [1, 2, 3, 2, 1]);
%   assert (polymul (x, y), [1, 2, 3, 2, 1]);
%   assert (polymul (y, x), [1; 2; 3; 2; 1]);
%   assert (polymul (c, x), [3; 3; 3]);
%   assert (polymul (c, y), [3, 3, 3]);
%   assert (polymul (x, c), [3; 3; 3]);
%   assert (polymul (y, c), [3, 3, 3]);
%   assert (polymul (b, c), 6);
%
% See also POLYADD, POLYSUB, POLYDIV, POLYPOW, CONV.

%   Author:      Peter J. Acklam
%   Time-stamp:  1998-06-22 20:36:03
%   E-mail:      jacklam@math.uio.no
%   WWW URL:     http://www.math.uio.no/~jacklam

   %error(nargchk(2, 2, nargin));
   narginchk(2,2)
   r = conv(p, q);

%!test
%! x = ones (3,1);
%! y = ones (1,3);
%! b = 2;
%! c = 3;
%! assert (polymul (x, x), [1; 2; 3; 2; 1]);
%! assert (polymul (y, y), [1, 2, 3, 2, 1]);
%! assert (polymul (x, y), [1, 2, 3, 2, 1]);
%! assert (polymul (y, x), [1; 2; 3; 2; 1]);
%! assert (polymul (c, x), [3; 3; 3]);
%! assert (polymul (c, y), [3, 3, 3]);
%! assert (polymul (x, c), [3; 3; 3]);
%! assert (polymul (y, c), [3, 3, 3]);
%! assert (polymul (b, c), 6);

%!shared a,b
%!test
%! a = 1:10;
%! b = 1:3;
%!assert (size (polymul (a,b)), [1, numel(a)+numel(b)-1])
%!assert (size (polymul (b,a)), [1, numel(a)+numel(b)-1])

%!test
%! a = (1:10).';
%!assert (size (polymul (a,b)), [numel(a)+numel(b)-1, 1])
%!assert (size (polymul (b,a)), [numel(a)+numel(b)-1, 1])

%!test
%! a = 1:10;
%! b = (1:3).';
%!assert (size (polymul (a,b)), [1, numel(a)+numel(b)-1])
%!assert (size (polymul (b,a)), [1, numel(a)+numel(b)-1])

%!test
%! a = 1:10;
%! b = 1:3;

%!assert (polymul (a,b,"full"), polymul (a,b))
%!assert (polymul (b,a,"full"), polymul (b,a))

%!assert (polymul (a,b,"same"), [4, 10, 16, 22, 28, 34, 40, 46, 52, 47])
%!assert (polymul (b,a,"same"), [28, 34, 40])

%!assert (polymul (a,b,"valid"), [10, 16, 22, 28, 34, 40, 46, 52])
%!assert (polymul (b,a,"valid"), zeros (1,0))


%% Test input validation
%!error polymul (1)
%!error polymul (1,2,3,4)
%!error <A and B must be vectors> polymul ([1, 2; 3, 4], 3)
%!error <A and B must be vectors> polymul (3, [1, 2; 3, 4])
%!error <SHAPE argument must be> polymul (2, 3, "INVALID_SHAPE")

