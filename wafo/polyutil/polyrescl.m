function q = polyrescl( p, x, y )
%POLYRESCL Rescale polynomial.
%
%   Q = POLYRESCL( P, X, Y ) rescaled the polynomial P by a factor X in
%   x-direction and by a factor Y in y-direction.
%
%   P is a vector of coefficients in decreasing order.
%
% Example:
%   assert(polyrescl([1,1], 2), [0.5, 1], eps);
%   assert(polyrescl([1,1], 4), [0.25, 1], eps);
%   assert(polyrescl([1,1], 2, 2), [1, 2], eps);
%
% See also POLYRELOC.

%   Author:      Peter J. Acklam
%   Time-stamp:  1999-04-24 00:12:57
%   E-mail:      jacklam@math.uio.no
%   WWW URL:     http://www.math.uio.no/~jacklam

   % Check number of input arguments.
   error( nargchk( 1, 3, nargin ) );

   % Now do the job.
   if nargin > 1
      n = length( p );
      q = p.*x.^(1-n:0);
   end
   if nargin > 2
      q = y.*q;
   end
