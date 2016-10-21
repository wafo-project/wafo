function wvers= wafoversion(dirnamein)
%WAFOVERSION Wave Analysis for Fatigue and Oceanography version.
%
%  CALL:  verstr =  wafoversion(dirname);
%
%    verstr = string containing the WAFO version number.
%   dirname = string containing the name of the directory
%             If empty WAFOVERSION  displays the current WAFO version 
%
%  WAFOVERSION extracts the version number from  the Contents.m file.
%
%   Examples:
%      wafoversion onedim      % returns the version of the ONEDIM directory
%      wafoversion             % return current WAFO  version
%
% See also   version, ver, whatsnew, info, help.

%tested on: Matlab 5.3
%History:
% revised pab 11.11.1999
%  - made it more general
% by pab 11.10.1999

if nargin<1 ||isempty(dirnamein)
  dirname=waforoot;
else
  ind=find(dirnamein== filesep);
  if length(ind)>0
    dirname=dirnamein; % Assume that the complete path is given
  else
    if isempty(ind)
      dirn=dirnamein;
    elseif ind(end)< length(dirnamein)
      dirn=dirnamein(ind(end)+1:end);
    else
      dirn=dirnamein(ind(end-1)+1:end-1);
    end
    p=path2cell(path);
    for ix=1:length(p)
      p{ix}=fliplr( p{ix});
    end
    ind2=strmatch(fliplr(dirn),p); % see if we can find directory in
                                   % search path
    if length(ind2)>1
      disp(['There are two directories with the name: ' dirn])
    end
    if length(ind2)>0
      dirname=fliplr(p{ind2(1)});
    else
      disp(['There is no directory with the name: ' dirn])
      disp('in the search path')
      disp('WAFO version returned instead')
      dirname = waforoot;
    end    
  end
end

file=fullfile(dirname ,'Contents.m');
fid=fopen(file,'rt');
if fid==-1,
  file=fullfile(dirname ,'contents.m');
  fid=fopen(file,'rt');
  if fid==-1,
    disp(['No Contents.m file in '  dirname])
    disp('Version unknown')
    wvers=[];
    return
  end
end

lin=fgetl(fid); % H1 line
lin=fgetl(fid); % version line
fclose(fid);
% Look for Version
k = findstr('version',lower(lin));
if isempty(k),
  disp(['Not correct format of Contents.m file in ', dirname])
  disp('version number unknown')
  wvers=[]; 
  return
end
ind=diff(~isspace(lin));
blancstrtstp = find([0 abs(ind)] > 0.5); 
indbl= blancstrtstp( blancstrtstp > k); 
wvers=lin(indbl(2):indbl(3)-1);

function c = path2cell(p)
%PATH2CELL Split path into cell array of strings.
%   PATH2CELL(PATH) returns a cell array of strings containing the
%   individual path elements.

seps = [0 find(p==pathsep) length(p)+1];

c = cell(length(seps)-1,1);
for i=1:length(seps)-1,
  c{i} = p(seps(i)+1:seps(i+1)-1);
  if c{i}(end) == ']',
    c{i}(end) = [];
  end
end
return
