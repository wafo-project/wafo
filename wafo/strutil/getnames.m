function  [names,no,nc] = getnames(string,Ltype)
%GETNAMES  Finds all possible legal names in a string.
%
%  CALL: [names,P0,P1] = getnames(str,Ltype)  
%
%  names  = list of legal names for MATLAB variables or functions.
%  P0, P1 = indices to the beginnings and ends, 
%           respectively, of these names in str.
%  str    = String to be divided into legal names.
%  Ltype  = Defines type of list. Options are
%           'cell' or 'matrix' (default).
%
% Example:
%  names = getnames('yy = foo(x1.^2,a_1*c,flag)');
%  assert(names, strvcat('yy', 'foo', 'x1', 'a_1', 'c', 'flag'))
%
% See also  extract, find, diff

% History:
% revised pab Aug 2006
% - removed redundant/empty items where nc <= no from list
% - moved some code into subfunction findnames
% - added Ltype option
% revised pab 15.05.2001
%  - improved the construction of mask
% revised pab 04.10.2000
%  - updated documentation
%  Kirill K. Pankratov,  kirill@plume.mit.edu
%  5/15/94, 12/15/94


% Handle input ...............................................

%error(nargchk(0,2,nargin));
narginchk(0,2)
if nargin == 0
  help getnames
  return
end
if nargin<2 || isempty(Ltype)
  Ltype = 'matrix';
end

 
string = string.';

[no,nc] = findnames(string);

if ~isempty(no)
  names = extract(string,no,nc,Ltype); % Extraction
  nc = nc-1;
elseif Ltype(1) =='m'
  names = '';
else
  names = cell(1,0);
end



function  [no,nc] = findnames(string)
%FINDNAMES Find starting and ending indices of all legal names of string
%
% CALL [no,nc] = findnames(string)
% 
%  no, nc = indices to the beginnings and ends, 
%           respectively, of all legal names in string.
%  string = String to be divided into legal names.
%
% Example:
%  [no,nc] = findnames('yy = foo(x1.^2,a_1*c,flag)')
%  names   = extract(string,no,nc);
%
% See also extract

 % Make string a single line ..........
string = [' ' string(:).' ' '];
no = [];
nc = [];
if isempty(deblank(string)),
  return,
end


mask = zeros(size(string));
a = find((string>='0'& string<='9')|string=='_');
if any(a),
  mask(a) = 9; % 9 is for numbers or '_'
end
a = find(isletter(string));
if any(a),
  mask(a) = 10; % 10 is for letters
end
%  Old call
%   mask = 9*(string>='0'& string<='9');       % 9 is for numbers
%   mask = mask+9*(string=='_');               % 9 is for '_'
%   mask = mask+10*(string>='A'&string<='Z'); % 10 is for letters
%   mask = mask+10*(string>='a'&string<='z');

no = find(diff(mask)>9)+1;     % Beginnings
nc = find(diff(mask)<=-9)+1;   % Ends

nn = sort([no nc+.5]);
a  = abs(diff([.5 nn]));
nn = nn(a-floor(a)>=.1);

if isempty(nn), 
 % no=[];
 % nc=[];
  return,
end

no = nn(floor(nn)==nn);
nc = floor(nn(floor(nn)~=nn));
ix = (nc <= no);
if any(ix)
  % Remove redundant/empty items.
  no(ix) = [];
  nc(ix) = [];
end
if ~isempty(no)
  no = no-1;
  nc = nc-1;
end













