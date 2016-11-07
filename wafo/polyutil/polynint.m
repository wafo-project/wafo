function r = polynint( p, n )
%POLYNINT Integrate polynomial N times.
%
%   POLYNINT( P, N ) returns the polynomial P integrated N times. P is a
%   vector of coefficients in decreasing order.
%
%   POLYNINT( P, -N ) returns POLYNDER( P, N );
%
%   POLYNINT( P, 0 ) returns P.
%
%   Example
%   assert(polynint([1/2, 1, 1],-1), [1,1], 1e-10);
%   assert(polynint([1,1,1],0), [1,1,1], 1e-10);
%   assert(polynint([3,2,1],1), [1,1,1,0], 1e-10);
%   assert(polynint([12,6,2],2), [1,1,1,0,0], 1e-10);
% 
%   See also POLYNDER, POLYDER, POLYINT.

%   Author:      Peter J. Acklam
%   Time-stamp:  1998-06-22 20:35:28
%   E-mail:      jacklam@math.uio.no
%   WWW URL:     http://www.math.uio.no/~jacklam

error( nargchk( 2, 2, nargin ) );

if ~all( size(n) == 1 ) || ~isreal(n) || n ~= round(n)
   error( 'N must be real scalar.' );
end

if n > 0
   m = length( p );
   for i = 1:n
      c = m+i-1 : -1 : i;
      p = p./c;
   end
   r = zeros( 1, m+n );
   r(1:m) = p;
elseif n < 0
   r = polynder( p, -n );
else
   r = p;
end
