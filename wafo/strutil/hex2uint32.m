function n = hex2uint32( h )
%HEX2UINT32 Convert hexadecimal string to uint32 type.
%   HEX2UINT32(H) interprets the hexadecimal string H and returns the
%   equivalent uint32 number.

%   Author:      Peter J. Acklam
%   Time-stamp:  1999-07-18 23:45:34
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
   n = pow2( 4*(hs(2)-1):-4:0 ) * h(:,:);
   n = uint32( reshape( n, ns ) );
