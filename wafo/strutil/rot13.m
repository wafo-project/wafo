function varargout = rot13(varargin)
%ROT13  Rot13 (Caesar) encryption.
%   ROT13(STR1, STR2, ...) applies a rot13 (Caesar) encryption on all input
%   strings and returns the same number of strings.
%
% Example
%
%   assert(rot13('Hello World!'), 'Uryyb Jbeyq!')

%   Author:      Peter J. Acklam
%   Time-stamp:  2000-07-17 02:00:58
%   E-mail:      jacklam@math.uio.no
%   WWW URL:     http://www.math.uio.no/~jacklam

   % create mapping vectors
   i = [ 13:25 0:12 ];
   j = [ 0:'A'-1 'A'+i 'Z'+1:'a'-1 'a'+i 'z'+1:255 ];

   varargout = cell(size(varargin));
   for i = 1:nargin
      varargout{i} = char(j(1+varargin{i}));
   end
