function  [no,nc,ins,bracketMisMatch] = brackets(str,cho,chc,chex)
% BRACKETS  Find brackets and parts of string enclosed in brackets. 
%
%  CALL:  [no,nc,ins] = brackets(str,cho,chc,chex)
%
%  no   = indices to opening characters in str
%  nc   = indices to closing characters in str
%  ins  = mask to parts of string inside opening and closing characters,
%         including the characters themselves.
%  str  = string to be searched for opening and closing characters.
%  cho  = opening character (default '(' )
%  chc  = closing character (default char(cho+2-(cho=='(')) )
%  chex = exclusion character. Opening and closing characters preceeded 
%         by the exclusion character are excluded from the search.
%         (default no exclusion)
%
% Note: - no and nc are sorted so that  no(i) and nc(i) (i=1,2,..)
%         give the index to the opening and the corresponding closing
%         character in str, respectively.
%       - ins also gives the number (level) of enclosing characters, i.e.,
%         how many e.g. parenthesis enclose the parts found.
%       - If the number of opening and closing character do not match in
%         size, the proper number of (dummy) opening and closing characters are
%         added to the start and end of the string, respectively, before the calculation
%         starts . 
% Example:
%  s1 = '(How) much (wood (would a)) \(Woodchuck chuck\)?';
%  [no,nc,ins] = brackets(s1,[],[],'\');
% 
%  assert(s1(find(ins)), '(How)(wood (would a))')  % find level 1 and higher enclosures
%  assert(s1(find(ins==2)), '(would a)') % find level 2 enclosures only
%

% History:
% revised pab aug2006
% - added bracketMisMatch to out put
% revised pab 07.12.2000
%  - fixed a bug: Now also able to handle cases where the number of
%    opening and closing character do not match.
%  - added a warning when number of opening and
%     closing characters do not match
%
% revised pab 20.10.2000
%  - Added documentation
%  - added isempty and nargchk to handle input part.
%  - fixed a bug: ins now also includes the closing character. 
% by  Kirill K. Pankratov, kirill@plume.mit.edu, 12/09/94

% Note: ins may exclude the closing character by changing some code
% indicated in the end of this file.


% Defaults
chodflt  = '(';      % Default for open character
chexdflt = 1000.55;  % Default for exclusion character (none)
txt = 'Number of opening and closing characters do not match'; % warning text
bracketMisMatch = false;

 % Handle input ...................
%error(nargchk(1,4,nargin))
narginchk(1,4)
if nargin < 4||isempty(chex),
  chex = chexdflt;
end
if nargin < 2||isempty(cho), 
  cho  = chodflt; 
end
if nargin < 3||isempty(chc),
  chc  = cho+2-(cho=='(');
end

sz = size(str);
 % Unwind the input string
str  = str';
str  = [' ' str(:)' ' '];
lstr = length(str)-2;

 % Initialize output
ins = zeros(1,lstr);
no  = []; 
nc  = [];

 % Find opening and closing characters
numbo = find(str==cho);   % Open
numbc = find(str==chc);   % Close

if isempty(numbo)||isempty(numbc),
  return,
end

% Exclude those preceeded by exclusion character
% (such as '\' in LaTeX)
numbo = numbo(str(numbo-1)~=chex);   % Open
numbc = numbc(str(numbc-1)~=chex);   % Close

if isempty(numbo)||isempty(numbc),
  return,
end

nbo   = length(numbo);
nbc   = length(numbc);

nco   = nbc-nbo;

bracketMisMatch = nco~=0;
if bracketMisMatch, 
  % The number of opening and closing char. do not match
  % Add  nco (dummy) opening or closing char. to numbc/numbo
  warning('strutil:brackets:BracketsMisMatch',txt)
  
  if nco<0,  
    numbc = [numbc (lstr+1)*ones(1,-nco)];
    %nbc   = length(numbc);
  else
    numbo = [2*ones(1,nco) numbo];  
    nbo   = length(numbo);
  end
end


[num,ind] = sort([numbo numbc]);
a         = [ones(1,nbo) -ones(1,nbo)] ;
a         = a(ind);

a1        = cumsum(a)+(a==-1);


nn        = -(min(a1)-1); % Positive value indicate that opening and
                          % closing char. do not match
% Old call kept just in case
% if nn>0, 
%   % Closing and opening characteres do not match   
%   % add dummy opening and closing characters to a and num.
%   if ~bracketMisMatch, 
%     warning(txt),
%   end
%   num = [2*ones(1,nn) num (lstr+1)*ones(1,nn)];
%   a   = [ones(1,nn) a  -ones(1,nn)];
%   a1  = cumsum(a)+(a==-1);
%   nbo = nbo+nn;
%   nbc = nbo;
% end 

[a,ind]   = sort(a1);
num       = num(ind);
num       = reshape(num,2,nbo);
[a,ind]   = sort(num(1,:));
num       = num(:,ind)-1;

% Output
no  = num(1,:);
nc  = num(2,:);

nsiz    = sz([ 2 1 3:end]);
if nargout>2 
  ins(no) = ones(1,nbo);      % Mark opening characters
  ins(no(1)) = ins(no(1)) + max(nn,0) + max(nco,0); % add start value
  ins(nc) = -ones(1,nbo);     % Mark closing characters
  ins     = cumsum(ins);      % include only the opening characters (original call)
  ins(nc) = ins(nc)+1;        % include the closing character as well (pab)
  ins     = reshape(ins,nsiz).';
end
if prod(nsiz)~=max(nsiz)
  [I, J] = ind2sub(nsiz,no);
  no    = sub2ind(sz,J,I);
  [I, J] = ind2sub(nsiz,nc);
  nc    = sub2ind(sz,J,I);    
end
return






