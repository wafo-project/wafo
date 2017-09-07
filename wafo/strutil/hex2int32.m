function n = hex2int32( h )
%HEX2INT32 Convert hexadecimal string to int32 type.
%   HEX2INT32(H) interprets the hexadecimal string H and returns the
%   equivalent int32 number.
%
% Example
%  assert(hex2int32('1ab'), int32(427))

%   Author:      Peter J. Acklam
%   Time-stamp:  1999-07-19 00:27:17
%   E-mail:      jacklam@math.uio.no
%   WWW URL:     http://www.math.uio.no/~jacklam

   % Check number of input arguments.
   %error( nargchk( 1, 1, nargin ) );
   narginchk(1,1)
   % Check type and size of input argument.
   if ~ischar( h ) || isempty( h )
      error( 'Argument must be a non-empty string.' );
   end
   hs = size( h );
   hd = ndims( h );
   if hs(2) > 8
      error( 'Argument can not have more than eight columns.' );
   end

   % Convert from 0-9, A-F, a-f to 0-15.
   h = h - '0';
   k = h >= 17;
   h(k) = h(k) - 7;
   k = h >= 42;
   h(k) = h(k) - 32;
   if any( h(:) < 0 ) || any( h(:) > 15 )
      error( 'String can only contain characters 0-9, A-F and a-f.' );
   end

   % Convert to the output data type.
   ns = hs;
   ns(2) = 1;
   h = permute( h, [ 2 1 3:hd ] );
   if h(2) == 8
      s = h(1,:) >= 8;
      h(1,s) = h(1,s) - 8;
   else
      s = [];
   end
   n = pow2( 4*(hs(2)-1):-4:0 ) * h(:,:);
   n(s) = n(s) - 2147483648;
   n = int32( reshape( n, ns ) );
