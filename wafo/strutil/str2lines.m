function lines = str2lines(str)
%STR2LINES Convert a string into a list of lines.
%   STR2LINES(STR) converts the string STR into a list (cell array) of lines
%   by splitting STR at line boundaries.  Both DOS (CR+LF), UNIX (LF) and
%   MAC (CR) line endings are recognized.  The line endings are removed from
%   the lines.
%
%  Example
%   s = ['Hello world' char(10) 'Hello space'];
%   lines = str2lines(s);
%   assert(lines, {'Hello world'; 'Hello space'})
%
%   %To convert from a list of lines back to a string using DOS (CR+LF) line
%   %endings, use
%
%   c = cell(2, length(lines));
%   c(1,:) = lines;
%   c(2,:) = {sprintf('\n')};
%   assert(s, [c{1:end-1}])
%
%   % or
%
%   c = [lines' ; repmat({sprintf('\n')}, 1, length(lines))];
%   assert(s, [c{1:end-1}])

%   Author:      Peter J. Acklam
%   Time-stamp:  2000-07-17 01:51:06
%   E-mail:      jacklam@math.uio.no
%   WWW URL:     http://www.math.uio.no/~jacklam
%
% revised pab 18.07.2002
% -replaced code with a call to dm2unix

   % Check number of input arguments.
   %error(nargchk(1, 1, nargin));
   narginchk(1,1)
   str = dm2unix(str,'');
   
   LF = char(10); % LineFeed character

   % Find the line endings.
   k = find(str == LF);

   % Avoid empty last string in output list when string ends in a newline.
   if ~isempty(k) && k(end) == length(str)
      k = [ 0 k ];                      % add beginning
   else
      k = [ 0 k length(str)+1 ];        % add beginning and end
   end
    
   if 1,
     clines = extract(str,k(1:end-1)+1,k(2:end));
     lines = cellstr(clines);
   else
     
     nlines = length(k) - 1;              % number of lines
     lines = cell(nlines, 1);             % initialize output
     for i = 1:nlines
       lines{i} = str(k(i)+1:k(i+1)-1);
     end
   end