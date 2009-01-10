function r = polyint( p, a, b )
%POLYINT Integrate polynomial.
%
%   R = POLYINT( P ) integrates the polynomial P.
%   I = POLYINT( P, T ) integrates the polynomial from 0 to T.
%   I = POLYINT( P, A, B ) integrates the polynomial from A to B.
%
%   P is a vector of coefficients in decreasing order.
%
%   See also POLYDER, POLYNINT.

%   Author:      Peter J. Acklam
%   Time-stamp:  1998-06-22 20:35:50
%   E-mail:      jacklam@math.uio.no
%   WWW URL:     http://www.math.uio.no/~jacklam

% Check number of input arguments.
error( nargchk( 1, 3, nargin ) );

% Get number of coefficients and perform integration.
n = length(p);
r = [ p./(n:-1:1)  0 ];

% Evaluate integral.
if ( nargin == 2 )              % POLYINT( P, T )
   r = polyval( r, a );
elseif ( nargin == 3 )          % POLYINT( P, A, B )
   r = polyval( r, b ) - polyval( r, a );
end
