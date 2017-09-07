function str = sizestr(array, before, between, after)
%SIZESTR Return size of array as a string.
%
%  CALL:  str = sizestr(array, os, sep, cs); 
% 
%   str   = Size string
%   array = Array
%   os    = Opening string   (default '[')
%   sep   = Separator string (default 'x')
%   cs    = Closing string   (default char(oc+2-(oc=='(')))
%
%   SIZESTR returns the size of the array ARRAY as a string, where OS
%   is the opening string written before the first size value, SEP is the
%   string written between each size value, and CS is the closing string
%   written after the last size value.
%
%   Examples:
%
%   assert(sizestr(rand(3,1,4), '[', 'x'), '[3x1x4]')
%   assert(sizestr(rand(3,1,4), '', 'x'), '[3x1x4]')
%   assert(sizestr(rand(3,1,4), ' ', '-by-'), ' 3-by-1-by-4 ')
%   assert(sizestr(rand(3,1,4), '( ', ' by '), '( 3 by 1 by 4 )')

% revised pab 24.01.2001
% - updated documentation
% - added default values for before between and after
%   Author:      Peter J. Acklam
%   Time-stamp:  2000-07-17 02:02:58
%   E-mail:      jacklam@math.uio.no
%   WWW URL:     http://www.math.uio.no/~jacklam

   %error(nargchk(1, 4, nargin));
   narginchk(1,4)
   if nargin < 2||isempty(before),  before  = '[';          end
   if nargin < 3||isempty(between), between = 'x';          end
   if nargin < 4||isempty(after), 
     after = fliplr(before);
     ind = find(~isspace(after));
     if any(ind)
       after(ind) = after(ind)+2-(after(ind)=='('); 
     end
   end

   siz = size(array);
   fmt = sprintf('%s%%d', between);
   str = sprintf(fmt, siz(2:end));
   str = sprintf('%s%d%s%s', before, siz(1), str, after);





