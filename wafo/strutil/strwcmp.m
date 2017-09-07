function status = strwcmp(string1, string2, casesense)
%STRWCMP Compares two strings allowing wildcards.
%  
%  CALL: status = strwcp(str1, str2,flag);
% 
%  status    = logical array 
%  str1,str2 = string or cellarray of strings  
%  flag      = 'exact' : letter match   (default)
%              'ignore': upper case letters
%
%  STRWCMP returns TRUE if the two strings STR1 and 
%  STR2 match.  When either STR1 or STR2 is a cellarray of strings,
%  STRWCMP returns  an array the same size as STR1 and STR2 containing 1 
%  for those elements of STR1 and STR2 that match, and 0 otherwise. 
%  One of the strings may contain one or more '*' as 
%  wildcard characters. The comparison is case sensitive. If the
%  strings are matrices they are treated as row-wise vectors.
%
% Example:
%  s1 = 'How much wood would a Woodchuck chuck?';
%  assert(strwcmp(s1,'*chuck*'), true)
%  assert(strwcmp(s1,'how*'), false)  
%  assert(strwcmp(s1,'how*','i'), true)
%
% See also: strcmp, findname

% Tested on: Matlab5.X
% History
% revised pab july 2005
% -improved the wildcard filtering: Now also handles strwcmp(s1,'*chuck')
%   correctly
% revised pab 07.12.2000
% - added the possibility that string1 or string2 is a cellarray of strings
% By  Peter M. W. Nave, 100010.3276@CompuServe.com
%  1998-04-01.
%----------------------------------------------------------------------O

%error(nargchk(2,3,nargin));
narginchk(2,3)
if ~(ischar(string1)||iscell(string1)||ischar(string2)||iscell(string2))
  error('Input must be a string or cellarray of strings.');
end

if nargin<3||isempty(casesense),casesense = 'exact';end
if strncmpi(casesense, 'ignore', 1)
  string1 = lower(string1); string2 = lower(string2);
end

if iscell(string1),
  sz1 = size(string1);
  N1 = prod(sz1);
  if iscell(string2),
    sz2 = size(string2);
    N2 = prod(sz2);
    if ~(all(sz1==sz2)||N1==1||N2==1), 
      error('Cellarrays must have common size or one scalar cell.');
    end
    sz     = max(sz1,sz2);
    status = zeros(sz);
    if N1==N2,
      for ix=1:N1, status(ix) = strwcmp1(string1{ix},string2{ix});  end
    elseif N1>N2,
      for ix=1:N1, status(ix) = strwcmp1(string1{ix},string2{1});   end
    else
      for ix=1:N2, status(ix) = strwcmp1(string1{1},string2{ix});   end
    end
  else
    status = zeros(sz1);
    for ix=1:N1,   status(ix) = strwcmp1(string1{ix},string2);    end
  end
elseif iscell(string2),
  sz2    = size(string2);
  status = zeros(sz2);
  for ix=1:prod(sz2),    status(ix) = strwcmp1(string2{ix},string1);  end
else % ischar
  status = strwcmp1(string1,string2);
end
status = logical(status);
return


function status=strwcmp1(string1,string2)
% STRWCMP1 core implementation of strwcmp

status  = 1;
string1 = string1.';   string2 = string2.'; % added pab 10.12.2000
string1 = string1(:)'; string2 = string2(:)';
if ~isa(string1, 'char') || ~isa(string2, 'char')
  error('### strwcmp: both input arguments must be character arrays');
end

if any(string1 == '*')
  mask = string1;
  if any(string2 == '*')
    error('### strwcmp: only one input string may contain ''*''.');
  end
  string = string2;
else
  string = string1;
  mask = string2;
end
  
  
first = ~strncmp(mask, '*', 1);
%last = ~strncmp(mask(end), '*', 1);

% for ii = 1:inf
while true
  if isempty(mask),
    if ~isempty(string)
      status =0;
    end
    return,
  end
  if all(mask == '*'), return, end
  
  [token, mask] = strtok(mask, '*');
     
  %ind = findstr(string, token);
  ind = strfind(string,token);
  if isempty(ind),    
    status = 0;
    return,
  else
    if first
      if ind(1) == 1
        string = string(length(token) + 1:end);
      else
        status = 0;
        return
      end
      first = 1;
    elseif length(token)<length(string)+1
      string = string(ind(end) + length(token):end);
      if isempty(string), return, end
    else
      status =0;
      return
    end
  end
end
%
% End of strwcmp.
  


