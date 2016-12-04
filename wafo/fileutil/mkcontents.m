function contstr = mkcontents(h1line,wver,srt, cfile)
%MKCONTENTS Makes Contents file in current working directory.
%
%  CALL:  contstr = mkcontents(H1line,version,sort,cfile);
%
%  contstr = character array containing the H1 lines of all the files in
%            the directory.
%  H1line  = string inserted as the first line of the file specified in
%            cfile. (default: same as the old version of cfile if it
%            exists otherwise the directory name is used) 
%  version = Version number (default: version number of old cfile
%            incremented by one otherwise '1.0.0')
%  sort    = 0 Do not sort by function name (default)
%            1 Sort by H1 line
%            2 Sort by function name
%            3 Sort by the inverse function name
%  cfile   = string containing the name of the file contstr  should be
%            written to.  (Default 'Contents.m') 
%            If ischar(cfile)==0 then  no file is made.
%
% It is worth noting that any editing done to a previous 
% version of CFILE will be lost. Only the H1 line from the
% old version are copied to the new version.
% It is assumed that the 2 first lines of the old Contents file is of the
% form:
%
%     % Toolbox Description
%     % Version <Number>  <Release>   dd-mm-yyyy
%
%  This form ensures that VER is able to read the version info from the
%  Contents.m file.
%
% Example
%  mkcontents('Test Toolbox',1.1,2)
%  ver test
%
% See also: geth1line, ver

% Tested on: Matlab 5.3
% History:
% revised pab 30.10.2003 
% renamed from mkcont to mkcontents  
% revised pab 7.12.2000
%  - changed name to mkcontents
%  - added H1txt, wver to input
%  - more options for srt.
% Revised by pab 08.10.1999
%  - added sorting
%  - changed call to mkcont
% by jr 27.09.1999

if (nargin<3||isempty(srt)),
  srt   = 0;
end
if (nargin<4||isempty(cfile)),
  cfile = 'Contents.m'; 
elseif ischar(cfile)
  cfile = [strtok(cfile,'.'),'.m' ];
end


vinfo = [];
if ischar(cfile) && exist(fullfile(pwd,cfile),'file') ~= 0, 
  if isempty(vinfo),
    vinfo  = getverinfo(fullfile(pwd,cfile));
  end
  disp(['There allready exist a ', cfile , ' in this directory'])
  disp(['Copied this file to ', cfile  'old'])
  %Make backup copy before overwriting cfile
  copyfile(cfile,[cfile 'old']);
  delete(cfile)
end
if nargin<1||isempty(h1line),
  if ~isempty(vinfo)&&~isempty(vinfo.Name) , % Use old H1 line
    h1line = vinfo.Name;
  else
    [tname , dirname] = gettbname;
    % Create new H1 line based on directory name
    h1line = [ tname    ' Toolbox ' dirname]; 
  end
end

if nargin<2||isempty(wver), 
  if isempty(vinfo), 
    wver = '1.0.0';
  else
    wver = vinfo.Version;
  end
  ind  = findstr(wver,'.');
  if isempty(ind),
    wver = [wver '.1'];
  else
    vr   = str2num(wver(ind(end)+1:end))+1; % update version
    wver = [wver(1:ind(end)) num2str(vr)];
  end
elseif isnumeric(wver)
  wver = num2str(wver); % make sure it is a string
end


contstr = '';
file    = what;
Nff     = length(file.m);
if (Nff==0),
  warning('No m-files found in this directory')
  return
end
if strcmpi(computer,'pcwin')
  file.m = lower(file.m);
end


Nf     = size(char(file.m),2)-1; % Size of filenames
tmp    = cell(Nff,1);
fn     = tmp;
ind    = zeros(Nff,1);
for ix=1:Nff    
  disp(file.m{ix})
  tmp{ix} = geth1line(file.m{ix},1,Nf); % Extract a formatted H1 line
  if ~isempty(tmp{ix})&&~strcmpi(file.m{ix},'contents.m')
    ind(ix) = 1;
    fn{ix}  = fliplr(file.m{ix});
  end
end
%tmp
ind = find(ind);
switch srt
  case 1, % make a sorted character array
    contstr1 = sort(tmp(ind)); 
  case 2, % Sort by file name 
    [t I] = sort(file.m(ind));
    contstr1 = tmp(ind(I));
  case 3, % Sort by the inverse filename
    [t I] = sort(fn(ind));
    contstr1 = tmp(ind(I));
  otherwise % No sorting 
    contstr1 = tmp(ind);       % make character array
end
contstr = char(contstr1);
if ~ischar(cfile), 
  disp(h1line)
  disp(['Version ', wver,'   ',date])
  disp(' ')
  disp(' ')
  disp(contstr),
  return, 
end

disp(['Creating ' cfile  ' in ' pwd])
fid = fopen(cfile,'wt'); % open a text file
% Write Header lines
prstr='%';
fprintf(fid,'%s \n',['% ', h1line]);
fprintf(fid,'%s \n',['% Version ', wver,'   ',date]);
fprintf(fid,'%s \n',prstr);
fprintf(fid,'%s \n',prstr);

% Write contents lines
if 1,
  fprintf(fid,'%% %s \n',contstr1{:});
else
  for ix=1:size(contstr,1),
    fprintf(fid,'%s \n',['% ', deblank(contstr(ix,:)) ] );
  end
end
fclose(fid);
% Change permissions on cfile 
% only valid for Unix systems, no effect in Win32 systems
[s,msg] = dos(['chmod go+r ' cfile]);

return

  

function [tname, dirname] = gettbname
    
%wafop     = waforoot;
%ind       = find(wafop==filesep);
wafodname =  ''; %wafop(ind(end)+1:end); 
pwdstr    = pwd;

ind = findstr(pwdstr,wafodname);
 
if isempty(ind)
  % name of directory is the name of the toolbox
  [parentDir,tname] = fileparts(pwdstr);
  tname           = upper(tname);
  dirname         = '';				   
else
  tname   = upper(wafodname);
  Nwf     = length(wafodname); 
  dirname = pwdstr(ind(1)+Nwf:end) ; % name of sub directory of toolbox
end
return

function s = getverinfo(cfile)
s =  struct('Name','','Version',{},'Release',{},'Date',{});

fid = fopen([strtok(cfile,'.'),'.m'],'rt');
if fid==-1,
  disp(['Unable to open ' cfile  ])
  disp('Version unknown')
  return
end

h1line = deblank(fgetl(fid));       % H1 line
vline  = deblank(fgetl(fid));       % version line
fclose(fid);

[r,c] = find(h1line ~= '%' & ~isspace(h1line));
if ~isempty(c), 
  % remove leading percent signs and trailing blanks
  [r,c2] = find(~isspace(h1line));
  h1line = h1line(:,min(c):max(c2)); 
end 
s(1).Name = h1line;
%h1line = geth1line(cfile);

% Look for Version
k = findstr('version',lower(vline));
if isempty(k),
  disp(['Not correct format of ' cfile ])
  disp('version number unknown')
  return
end
ind          = diff(~isspace(vline));
blancstrtstp = find([0 abs(ind)] > 0.5); 
indbl        = blancstrtstp( blancstrtstp > k); 

s.Version = vline(indbl(2):indbl(3)-1);

Nbl = length(indbl);
if (Nbl>=3),
  s.Date = vline(indbl(Nbl):end);
end

if (Nbl-1>=3),
  s.Release = vline(indbl(Nbl-1):indbl(Nbl)-1);
end

return

