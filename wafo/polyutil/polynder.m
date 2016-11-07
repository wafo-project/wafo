function r = polynder( p, n )
%POLYNDER Differentiate polynomial N times.
%
%   POLYNDER( P, N ) returns the N'th derivative of the polynomial whose
%   coefficients are the elements of the vector P.
%
%   POLYNDER( P, -N ) returns POLYNINT( P, N ).
%
%   POLYNDER( P, 0 ) returns P.
%
%   Example
%   assert(polynder([3, 2, 1],-1), [1,1,1,0], 1e-10);
%   assert(polynder([1,1,1],0), [1,1,1], 1e-10);
%   assert(polynder([3,2,1],1), [6, 2], 1e-10);
%   assert(polynder([12,6,2],2), 24, 1e-10);
%
%   See also POLYNINT, POLYINT, POLYDER.

%   Author:      Peter J. Acklam
%   Time-stamp:  1998-06-22 20:35:23
%   E-mail:      jacklam@math.uio.no
%   WWW URL:     http://www.math.uio.no/~jacklam

% Check number of input arguments.
error( nargchk( 2, 2, nargin ) );

% Make sure N is real and scalar.
if any( size(n) ~= 1 ) || ~isreal(n) || n ~= round(n)
   error( 'N must be real and scalar.' );
end

if n > 0                        % Do differentiation if N > 0.
   m = length(p);
   r = p(1:m-n);
   for i = 1:n
      c = m-i : -1 : n+1-i;
      r = c.*r;
   end
elseif n < 0                    % Do integration if N < 0.
   r = polynint( p, -n );
else                            % Do nothing if N = 0.
   r = p;
end
