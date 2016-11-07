function r = polyreloc( p, x, y )
%POLYRELOC Relocate polynomial.
%
%   R = POLYRELOC( P, X, Y ) relocates the polynomial by "moving" it X
%   units along the x-axis and Y units along the y-axis. So R is
%   relative to the point (X,Y) as P is relative to the point (0,0).
%
%   P is a vector of coefficients in decreasing order.
%
% Example:
%   assert(polyreloc([1,1], 2), [1, -1], eps);
%   assert(polyreloc([1,1], 4), [1, -3], eps);
%   assert(polyreloc([1,1], 2, 2), [1,1], eps);
%
% See also POLYRESCL.

%   Author:      Peter J. Acklam
%   Time-stamp:  1999-04-24 00:12:55
%   E-mail:      jacklam@math.uio.no
%   WWW URL:     http://www.math.uio.no/~jacklam

   error( nargchk( 2, 3, nargin ) );

   n = length( p );
   %
   % Move polynomial X units to the right by a polynomial version of
   % Horner's method.
   %
   f = [ 1 -x ];
  
   r = p(1);
   for i = 1:n-1
      r = conv( r, f );
      r(:,i+1) = r(:,i+1) + p(:,i+1);
   end
   
   %
   % Move polynomial Y units upwards by adding Y.
   %
   if nargin > 2
      r(n) = r(n) + y;
   end
