function s1 = strim(s,c)
% STRIM Removes the leading and trailing character of a string.
%
%  CALL: s1 = strim(s,c);
%
%   s1,s = strings or cellarray of strings
%   c    = character to be removed   (default ' ')
%
% STRIM removes the leading an trailing character as well as the leading
% and trailing blanks.
% If S is a cellarray of strings it works on every element of S.
%
% Example:
%  s  = '    Testing testing 1-2-3   ';
%  assert(strim(s), 'Testing testing 1-2-3')
%
% See also: deblankl 

% History:
% revised pab, 07.12.2000
% - added c to input
% - made it faster.
% - added the possibility that s is a cellarray
% - updated the documentation
% - added example
% Date 11-08-98
% Written by Stefan Baunack. Please send any bug
% reports or other comments to: s.baunack@ifw-dresden.de.

%error(nargchk(1,2,nargin));
narginchk(1,2)
if isempty(s), s1=s; return; end
if ~(ischar(s)||iscell(s))
  error('Input must be a string or cellarray of strings.')
end
if nargin<2||isempty(c),c = ' '; end


if iscell(s),
  s1=cell(size(s));
  for ix=1:numel(s),
    s1{ix} = strim1(s{ix},c);
  end
else % ischar
  s1  = strim1(s,c);
end

return

function s1 = strim1(s,c)
% STRIM1 Core implementation of strim

if ~(ischar(s)),
  error('Input must be a string.')
end
[r,c] = find(s ~= c(1) & s~=' ' & s ~= 0);
if isempty(c) 
  s1 = '';
else
  s1 = s(:,min(c):max(c));
end
return
% old call: much slower
%s1 = deblank(s);
%s1 = fliplr(s1);
%s1 = deblank(s1);
%s1 = fliplr(s1);
