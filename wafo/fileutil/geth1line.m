function H1 = geth1line(filename,format,N)
%GETH1LINE Extracts the first comment line (the H1 line) of a m-file
% 
%   CALL:  H1 = geth1line(filename,format,N);
% 
%      H1       = string containing the H1 line of the file
%      filename = string containing the name of the file
%      format   = 0 No formatting of H1 (default)
%                 1 Format H1 
%      N        = Number of characters for filename (default 14)
% 
%  GETH1LINE assumes the H1LINE is formatted as:
%  
%  <Filename> <One line comment>
%  
%  GETH1LINE is used by MKCONTENTS to simply generate a formatted
%  contents file for a directory with m-files.
%  
%  Examples:
%  assert(geth1line('geth1line',1),...
%  'geth1line     - Extracts the first comment line (the H1 line) of a m-file');
%  assert(geth1line('geth1line',1,7), ...
%   'geth1line - Extracts the first comment line (the H1 line) of a m-file');
%  assert(geth1line('geth1line',0),...
%  'GETH1LINE     Extracts the first comment line (the H1 line) of a m-file');
% 
% See also: mkcontents
  
  
% Tested on: Matlab 5.3
% History: 
% revised pab 2006
% - more robust handling of filename
% revised pab 30.10.2003
% changed name from h1line to geth1line 
% -improved the robustness of the code  
% - added ddeblank  
% revised pab 07.12.2000
%  - added examples
%  - added N to input
%  - cleaned up some code by using strtok
%  - removed newline characters from H1
% Revised pab 27.01.2000
%  - improved while loop for searching percent signs  
% Revised by pab 08.11.1999
%   - added frmt
% Revised by pab 08.10.1999
%   - updated documentation, 
%   - changed name from findline to h1line
%   - added the possibility to specify the complete file name
%   - improved formatting
% by jr 27.09.1999

%filename
if (nargin<2||isempty(format)), 
  format = 0; 
end
if (nargin<3||isempty(N)),    
  N    = 14; 
end

%filename
H1 = '';
[filepath,filename] = fileparts(filename); % strip the .m extension

fid      = fopen(fullfile(filepath,[filename, '.m']),'r');
if fid==-1,
  fid      = fopen([filename, '.m'],'r');
  if fid==-1
    disp([ 'Unable to open file ', filename ])
    return
  end
end

foundFirstComment = 0;
while ((~foundFirstComment) && (feof(fid)==0)), 
  tline = ddeblank(fgetl(fid));
  if (length(tline) > 1), % Allow that some lines may be empty or only 
			  % contain a percent sign
    k = findstr(tline,'%'); %LOOK for percent sign anywhere in the line
    foundFirstComment = ~isempty(k);
  end  
end
fclose(fid);

if (~any(isletter(tline)) || length(tline)<2||...
    isempty(k) ||~strcmp(tline(k(1)),'%')), 
  warning('No help header is found in %s',  filename)
  return
end
[fn,rr] = strtok(tline(k(1)+1:end));
[pp,fn] = fileparts(fn); 
if strcmpi(fn,filename)
  rr = ddeblank(rr);  % strip leading & trailing blanks 
else
  warning('Wrong H1LINE-style: Filename (%s) differ from first token in H1LINE.',...
    filename)
  fn = filename;
  rr = ddeblank(tline(k(1)+1:end));
end
Nbl     = max(N - length(filename) ,1); % # of blanks
if isempty(rr)
   rr =  ' ';
end
switch format
  case 0, % No formatting of H1 line
    H1 = [ fn blanks(Nbl)  rr]; % alternative
  otherwise      %  Format H1 line  
    % make sure func. name  is in lower case letters
    H1 = [ lower(fn), blanks(Nbl) , '- ', upper(rr(1)), rr(2:end) ];
end
return

function s1 = ddeblank(s)
%DDEBLANK Double deblank. Removes both leading and trailing blanks of string.
%
%  CALL: s1 = ddeblank(s);
%
%   s1,s = strings or cellarray of strings
%
% DDEBLANK removes leading and trailing blanks and null characters from
% the string S.  A null character is one that has an absolute value of 0.

% If S is a cellarray of strings it works on every element of S.
%
% Example:
%  s  = '    Testing testing 1-2-3   '
%  s1 = ddeblank(s)
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

error(nargchk(1,1,nargin));

if isempty(s), 
  s1=s; 
  return; 
end
if ~(ischar(s)||iscell(s))
  error('Input must be a string or cellarray of strings.')
end

if iscell(s),
  ssz = size(s);
  s1  = cell(ssz);
  for ix=1:prod(ssz),
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
