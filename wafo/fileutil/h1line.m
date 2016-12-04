function H1 = h1line(filename,frmt,N)
%H1LINE  Extracts the first comment line (the H1 line) of a m-file
%
% CALL:  H1 = h1line(filename,frmt,N);
%
%    H1       = string containing the H1 line of the file
%    filename = string containing the name of the file
%    frmt     = 0 No formatting of H1 (default)
%               1 Format H1 
%    N        = Number of characters for filename (default 14)
%
% Used by MKCONT to simply generate a formatted Contents file
%
% Examples:
%   h1line('h1line',1)
%   h1line('h1line',0)
%   h1line('h1line',1,7)
%
% See also: mkcont

% Tested on: Matlab 5.3
% History:
% revised pab 16.01.2001
% fixed a bug: H1 - lines with only one word is now handled correctly
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
if nargin<2||isempty(frmt), frmt = 0; end
if nargin<3||isempty(N),    N    = 14; end

%filename
H1 = [];
filename = strtok(filename,'.'); % strip the .m extension

fid      = fopen([filename, '.m'],'r');
if fid==-1,
  disp([ 'Unable to open file ', filename ])
  return
end

ix=0;
while (ix < 1) && (feof(fid)==0), 
  lin = fgetl(fid);
  %disp(lin)
  if length(deblank(lin)) > 1 %Allow that some lines may be empty or only 
                              %contain a percent sign
    k = findstr(lin,'%'); %LOOK for percent sign anywhere in the line
    if ~isempty(k) 
      ix = 10;
    end
  end  
end
fclose(fid);

if ~any(isletter(lin)) || length(lin)<2|| isempty(k) ||~strcmp(lin(k(1)),'%'), 
  disp(['No help header is found in '  filename])
  return
end

[fn,rr] = strtok(lin(k(1)+1:end));
ch_nl = [10 13];  % New line characters
inl = find(rr==ch_nl(1)|rr==ch_nl(2), 1 ); % New line ch. pointers
if ~isempty(inl)
  rr = rr(1:inl-1); % remove Newline characters from rr
end
rr      = strtrim(rr); % strip leading &
                                                % trailing blanks 
Nbl     = max(N - length(filename) ,1);         % Nr of blanks

switch frmt
  case 0, % No formatting of H1 line
    %H1 = [ filename blanks(Nbl)  rr];
    H1 = [ fn blanks(Nbl)  rr]; % alternative
  otherwise      %  Format H1 line  
    % make sure func. name  is in lower case letters
    %H1 = [ lower(filename), blanks(Nbl) , '- ', upper(rr(1)), rr(2:end) ];
    if any(rr)
      H1 = [ lower(fn), blanks(Nbl) , '- ', upper(rr(1)), rr(2:end) ]; % alternative
    else
      H1 = [ lower(fn), blanks(Nbl) ]; % alternative
    
    end
end








