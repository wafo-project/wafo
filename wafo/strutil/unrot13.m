function varargout = unrot13(varargin)
%UNROT13 Rot13 (Caesar) decryption.
%   UNROT13(STR1, STR2, ...) applies a rot13 (Caesar) decryption on all
%   input strings and returns the same number of strings.
%
%   Example
%   assert(unrot13('Uryyb Jbeyq!'), 'Hello World!')

%   Author:      Peter J. Acklam
%   Time-stamp:  2000-07-17 02:01:04
%   E-mail:      jacklam@math.uio.no
%   WWW URL:     http://www.math.uio.no/~jacklam

   varargout = cell(size(varargin));
   [varargout{:}] = rot13(varargin{:});
