function s1 = ddeblank(s)
%DDEBLANK Double deblank. Removes both leading and trailing blanks of string.
%
%  CALL: s1 = ddeblank(s);
%
%   s1,s = strings or cellarray of strings
%
% DDEBLANK removes leading and trailing blanks and null characters from
% the string S.  A null character is one that has an absolute value of 0.
%
% If S is a cellarray of strings it works on every element of S.
%
% Example:
%  s  = '    Testing testing 1-2-3   ';
%  cs = {s, ' deo, deo, \\n  '};
%  assert(ddeblank(s), 'Testing testing 1-2-3')
%  assert(ddeblank(cs), {'Testing testing 1-2-3', 'deo, deo, \\n'}) 
%
% See also: deblank, dewhite, ddwhite

% History:
% revised pab 30.10.2003
%  -  
% revised pab, 07.12.2000
% - renamed from strim to ddeblank
% - made it faster.
% - added the possibility that s is a cellarray
% - updated the documentation
% - added example
%
%   Author:      Peter J. Acklam
%   Time-stamp:  2000-07-17 02:05:04
%   E-mail:      jacklam@math.uio.no
%   WWW URL:     http://www.math.uio.no/~jacklam
% Date 11-08-98
% Written by Stefan Baunack. Please send any bug
% reports or other comments to: s.baunack@ifw-dresden.de.

%error(nargchk(1,1,nargin));
narginchk(1,1)
if isempty(s), 
  s1=s; 
  return; 
end
if ~(ischar(s)||iscell(s))
  error('Input must be a string or cellarray of strings.')
end

if iscell(s),
  s1=cell(size(s));
  for ix=1:numel(s),
    s1{ix} = strim1(s{ix});
  end
else % ischar
  s1  = strim1(s);
end

return

function s1 = strim1(s)
% STRIM1 Core implementation of strim

if ~(ischar(s)),
  error('Input must be a string.')
end
[r,c] = find(~isspace(s) & s ~= 0);

if isempty(c) 
  s1 = '';
elseif size(s, 1) == 1
  s1 = s(min(c):max(c));
else
  s1 = s(:,min(c):max(c));
end
return
