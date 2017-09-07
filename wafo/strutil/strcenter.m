function lines = strcenter(lines, offset)
%STRCENTER Center each line of text.
%   STRCENTER(LINES, OFFSET) where LINES is a cell array of lines and OFFSET
%   is a non-negative integer, does a centering of the text in each line
%   relative to the specified offset.  If OFFSET is omitted, 72 is used.
%
% Example
%
%  lines = {'Hello world', 'Hello space'};
%  assert(strcenter(lines, 25), {'       Hello world'; ...
%                                '       Hello space'})

%   Author:      Peter J. Acklam
%   Time-stamp:  2000-07-17 02:15:50
%   E-mail:      jacklam@math.uio.no
%   WWW URL:     http://www.math.uio.no/~jacklam

   %error(nargchk(1, 2, nargin));
   narginchk(1,2)
   if nargin < 2
      offset = 72;
   end

   spc = ' ';
   tab = sprintf('\t');
   lines = lines(:);
    
   for i = 1:length(lines)
      k = find(lines{i} ~= spc & lines{i} ~= tab);
      if isempty(k),
         lines{i} = '';
      else
         linlen = k(end) - k(1) + 1;
         padlen = floor((offset - linlen)/2);
         lines{i} = [spc(:,ones(padlen,1)) lines{i}(k(1):k(end))];
      end
   end
