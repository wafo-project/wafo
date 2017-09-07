function s1 = dewhite(s)
%DEWHITE Dewhite. Strip trailing whitespace.
%
% CALL: s1 = dewhite(s)
%
%    s1, s = strings or cellarray of strings
%  
%   DEWHITE removes leading and trailing white space and any null
%   characters from the string S.  A null character is one that has an
%   absolute value of 0.
%
% See also: ddewhite, deblank, ddblank.

% History:
% revised pab 14.12.2000
% - added the possibillity that s is a cellarray
% - updated documentation
%   Author:      Peter J. Acklam
%   Time-stamp:  2000-07-17 02:05:44
%   E-mail:      jacklam@math.uio.no
%   WWW URL:     http://www.math.uio.no/~jacklam

%error(nargchk(1, 1, nargin));
narginchk(1,1)
if isempty(s), s1=s; return; end
if ~(ischar(s)||iscell(s))
  error('Input must be a string or cellarray of strings.');
end

if iscell(s),
  s1=cell(size(s));
  for ix=1:numel(s),
    s1{ix} = dewhite1(s{ix});
  end
else % ischar
  s1  = dewhite1(s);
end
return

function sout = dewhite1(s)
% core implementation  
if ~ischar(s)
  error( 'Input must be a string (char array).' );
end

if isempty(s)
  sout = s;
  return;
end

[r, c] = find(~isspace(s));
if isempty(c)
  sout='';
elseif size(s, 1) == 1
  sout = s(1:max(c));
else
  sout = s(:,1:max(c));
end
  
  