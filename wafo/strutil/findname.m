function [names,p0,p1] = findname(s1,s2,flag)
% FINDNAME Finding name (word) within a string.
%
%  CALL: [names,P0,P1] = findname(S1,S2,flag); 
%
%  names = matrix of full names matching the shorter string 
%  P0,P1 = starting and ending indices, respectively, of 
%          all matching occurrences of the shorter string 
%          in the longer one.
%  S1,S2 = long and short string, respectively.
%  flag  = 'exact' : match letter for letter (default)
%          'ignore': uppercase letters.
%
%  The shorter string can have "wildcard" options:
%  '?' can stand for any non-blank character
%  while  '*' at the beginning or at the end
%  can stand for any number of pre- or appended
%  non-blank characters, respectively.
%  Note: Only the first letter of the flag is needed 
%        for a unique identification. 
% Example: 
%  s = 'How much wood would a Woodchuck chuck?';
%  assert(findname(s,'wo*'), strvcat('wood', 'would'))
#  assert(findname(s,'wo*','ignore'), strvcat('wood', 'would', 'Woodchuck'))  
%
% See also: getnames

% Tested on: Matlab 5.3
% History:
% revised pab 04.10.2000
%   - added flag, nargchk
%   - cleaned up documentation
%   - changed ==[] to isempty(..)
%  Kirill K. Pankratov,  kirill@plume.mit.edu
%  02/05/95

 % Handle input .............................
%error(nargchk(2,3,nargin));
narginchk(2,3)
if nargin < 3||isempty(flag),  flag='exact';  end

% Initialize output
p0 = []; p1 = []; names = [];

sz1 = size(s1);
ls2 = length(s2);
if ls2==0||max(sz1)==0, return, end

if sz1(2)<ls2 && sz1(1)==1
  ss = s1; s1 = s2; s2 = ss; ls2=max(sz1);sz1=size(s1);
end



 % Check the shorter string for wildcards ...
nn     = find(s2=='?');
lq     = length(nn);
s2(nn) = 32*ones(size(nn)); %32.01
is_wc  = s2([1 ls2])=='*';
s2     = s2(s2~='*');
ls2    = length(s2);

 % Extract all possible names ...............
if sz1(1)==1,
  [names,p0,p1] = getnames(s1);
else
  names = s1; p0 = 1:sz1(1); 
end
if ls2 == 0, return, end
szn = size(names);

 % If s2 is longer than any name, return empty
if ls2>szn(2), names=[]; p0=[]; p1=[]; return, end

if strncmpi(flag,'exact',1),
  names2 = names;
else % ignore upper case letters
  names2 = lower(names);
  s2     = lower(s2);
end

lm    = (szn(2)-ls2)*is_wc(1);
A     = s2(ones(szn(1),1),:)';
ll    = ls2-lq;
v_cnc = zeros(1,szn(1));
vl    = max(.5,sum(names2'~=' ')-ls2+1);
for jj=0:lm
  Ad = A-names2(:,jj+(1:ls2))';
  if ls2>1
    vv = sum(~Ad)==ll;
    vv = vv & sum(abs(Ad)<.1)==ll;
  else
    vv = (~Ad) | (~ll)*(abs(Ad)>.1);
  end
  nn = find(vv);
  v_cnc(nn) = (jj+1)*ones(size(nn));
end
%v_cnc = is_wc(2)*v_cnc+(~is_wc(2))*(v_cnc==vl);
if is_wc(2),
  nn = find(v_cnc);
else
  nn = find(v_cnc==vl);
end
 % Extract matching names ........
p0    = p0(nn);
p1    = p1(nn);
names = deblank(names(nn,:));
