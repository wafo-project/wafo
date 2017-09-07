function str = hex2bin(h, n)
%HEX2BIN  Convert hexadecimal string to a binary string.
%
%   HEX2BIN(H) returns the binary representation of D as a string.
%   HEX2BIN(H, N) produces a binary representation with at least N bits.
%
%  Example
%  assert(hex2bin('d8'), '11011000')
%
%   See also BIN2DEC, DEC2BIN, DEC2HEX, DEC2BASE.

%   Author:      Peter J. Acklam
%   Time-stamp:  2000-04-20 15:55:26
%   E-mail:      jacklam@math.uio.no
%   WWW URL:     http://www.math.uio.no/~jacklam

%   error(nargchk(1, 2, nargin));
narginchk(1,2)
   if nargin < 2
      n = 1;            % Need at least one digit even for 0.
   else
      if any( size(n) ~= 1 ) || ( n < 0 )
         error( 'N must be a positive scalar.' );
      end
      n = round(n);     % Make sure n is an integer.
   end

   if ~ischar(h) ...
      || any(   ( ( h(:) < '0' ) | ( '9' < h(:) ) ) ...
             & ( ( h(:) < 'A' ) | ( 'F' < h(:) ) ) ...
             & ( ( h(:) < 'a' ) | ( 'f' < h(:) ) ) );
      error('Invalid hexadecimal string.');
   end

   s = size(h);
   v = [ 2 1 3:ndims(h) ];

   % Convert from hexadecimal to binary double.
   d = sscanf(permute( h, v ), '%1x');
   b = rem( floor( d*pow2( -3 : 0 ) ), 2 );
   b = permute(reshape(b', [ 4*s(2) s(1) s(3:end) ]), v);

   % Convert from double to char.
   str = zeros( [ s(1) max(4*s(2),n) s(3:end) ] );
   str(:,end-4*s(2)+1:end,:) = b;
   str = char( str + '0' );
