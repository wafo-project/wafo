function s1 = ddewhite(s)
%DDEWHITE Double dewhite. Strip both leading and trailing whitespace.
%
% CALL: s1 = ddewhite(s)
%
%    s1, s = strings or cellarray of strings
%
%   DDEWHITE removes leading and trailing white space and any null
%   characters from the string S.  A null character is one that has an
%   absolute value of 0.
%
% Example:
%
%  s  = '    Testing testing 1-2-3   ';
%  cs = {s, ' deo, deo, \\n  '};
%  assert(ddewhite(s), 'Testing testing 1-2-3')
%  assert(ddewhite(cs), {'Testing testing 1-2-3', 'deo, deo, \\n'})
%
% See also: dewhite, deblank, ddblank.

% History:
% revised pab 14.12.2000
% - added the possibillity that s is a cellarray
% - updated documentation
%   Author:      Peter J. Acklam
%   Time-stamp:  2000-07-17 02:04:23
%   E-mail:      jacklam@math.uio.no
%   WWW URL:     http://www.math.uio.no/~jacklam

%error(nargchk(1, 1, nargin));
narginchk(1,1)
if isempty(s), s1=s; return; end
if ~(ischar(s)||iscell(s))
  error('Input must be a string or cellarray of strings.')
end

if iscell(s),
  s1=cell(size(s));
  for ix=1:numel(s),
    s1{ix} = ddwhite1(s{ix});
  end
else % ischar
  s1  = ddwhite1(s);
end

return



function sout = ddwhite1(s)
% Core implementation  
if ~ischar(s)
  error('Input must be a string (char array).');
end

if isempty(s)
  sout = s;
  return;
end

[r, c] = find(~isspace(s));
if isempty(c)
  sout='';
elseif size(s, 1) == 1
  sout = s(min(c):max(c));
else
  sout = s(:,min(c):max(c));
end
