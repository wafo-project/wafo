function s = strrepl(s1,s2,s3,varargin)
%STRREPL  Replaces the string S2 in S1 with S3.
% 
% CALL:  S = strrepl(S1,S2,S3,flag,flag2) 
%
% S,S1,
% S2,S3 = character strings  
%  flag = 'word'  : Replaces only occurences which are
%	            words themselves and not parts of 
%                   bigger words (default).
%         'all'   : replaces all occurences
% flag2 = 'exact' : letter match (default)
%         'ignore': uppercase letters. 
%
% Note: Only the first letter of the flags are needed 
%       for a unique identification. flag and flag2 may 
%       be given in any order.
%
% Examples:
% s1 = 'How much wood would a Woodchuck chuck?';
% assert(strrepl(s1,'chuck','pack'), 'How much wood would a Woodchuck pack?')
% assert(strrepl(s1,'chuck','pack','all'), ...
%                   'How much wood would a Woodpack pack?')
% assert(strrepl(s1,'Wood', 'food','all'), ...
%                   'How much wood would a foodchuck chuck?')
% assert(strrepl(s1,'Wood', 'food','ignore','all'), ...
%                   'How much food would a foodchuck chuck?')
%
% See also: findname, strrep, findstr

% History:
% revised pab 04.10.2000
%  - updated documentation
%  - added the options 'exact'/'ignore'
%  Based on the STRREP.M program in MATLAB
%  Kirill K. Pankratov,  kirill@plume.mit.edu
%  02/05/95

flag  = 'word';  % Default for "names" vs. any string
flag2 = 'exact'; % Default exact letter match 

 % Handle input ...............................
%error(nargchk(3,5,nargin));
narginchk(3,5)
Np = length(varargin);
for ix=1:Np
  switch lower(varargin{ix}(1))
    case {'w'}, flag='word';
    case {'a'}, flag='all'; 
    case {'e'}, flag2='exact';
    case {'i'}, flag2='ignore';
  end
end

is_name = ~strcmp(flag,'all');
% Take care of trivial cases:
%   s2 > s1,  s1 is empty  or  s2 is empty
if (length(s2)>length(s1)) || isempty(s2) || isempty(s1),
  s = s1;  return;
end
s1 = s1.';
s1 = [s1(:)' ' '];

if is_name   % If only legal names (words) are to be replaced
  [names,p0,p1] = findname(s1,s2,flag2);
else         % If all coincidences are to be replaced
  %ls2 = length(s2);
  if strcmpi(flag2,'exact'),
    p0 = findstr(s1,s2);
  else
    p0 = findstr(lower(s1),lower(s2));
  end
  p1 = p0+length(s2)-1;
end

list = s3(ones(size(p1)),:);
  % Insert new string
s = insert(s1,list,p0,p1+1);
s = s(1:length(s)-1);

return



