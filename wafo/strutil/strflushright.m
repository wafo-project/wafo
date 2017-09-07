function lines = strflushright(lines, offset)
%STRFLUSHRIGHT Flush right each line of text.
%   STRFLUSHRIGHT(LINES, OFFSET) where LINES is a cell array of lines and
%   OFFSET is a non-negative integer, does a right flush of the text in each
%   line so the last character in each line is at the specified offset.  If
%   OFFSET is omitted, 72 is used.
%
%  Example
%  lines = {'Hello world', 'Hello space'};
%  assert(strflushright(lines, 18), {'       Hello world'; ...
%                                    '       Hello space'})

%   Author:      Peter J. Acklam
%   Time-stamp:  2000-07-17 02:14:18
%   E-mail:      jacklam@math.uio.no
%   WWW URL:     http://www.math.uio.no/~jacklam

   %error(nargchk(1, 2, nargin));
   narginchk(1,2)
   if nargin < 2
      offset = 72;
   end

   spc = ' ';
   tab = sprintf('\t');
   fmt = sprintf('%%%ds', offset);
   lines = lines(:);

   for i = 1:length(lines)
      k = find(lines{i} ~= spc & lines{i} ~= tab);
      if isempty(k)
         lines{i} = '';
      else
         lines{i} = sprintf(fmt, lines{i}(k(1):k(end)));
      end
   end
