function initwafo(opt,add,first)
%INITWAFO   Initiation of WAFO Toolbox.
%   WAFO = "Wave Analysis for Fatigue and Oceanography"
%
% CALL:  initwafo(opt,add,first)
%
% Set MATLAB-paths to WAFO toolbox.
%
% The following options may be used.
%   opt   = 'minimum' : Minimum initiation, only paths to the 
%                       necessary m-files are added.
%           'normal'  : Normal initiation, also paths to demos, data and
%                       documentation are added. (default)
%           'full'    : Full initiation, all paths are added.
%   add   = 1 : Add paths from MATLAB search path (default)
%           0 : Remove paths from MATLAB search path
%   first = 1 : Add the paths first in the search path (default)
%           0 : Add the paths last in the search path
%   
% Example:
%   initwafo;             % Initiate WAFO toolbox paths
%   initwafo([],0);       % Remove WAFO toolbox paths
%   initwafo('full',1,0); % Full init., put WAFO paths last in search path
%
% See also  wafopath, waforoot, wafoexepath.

% Tested on: Matlab 5.3
% History: 
% Revised pab March 2007
% -added help('wafolicence') and help('wafocitation')
% Revised jr 01.01.30
%   Changed 'initiation' to 'removal' (line 108)
% Revised jr 00.05.23
%   Changed demos to wdemos
% Updated by PJ 25-Feb-2000
%   Removed input option: opt='root'
%   Now keeps WAFO root-directory in search path when removing WAFO paths.
% Revised pab feb 2000
%   Enabled several subdirectories in the papers directory
% Updated by PJ 19-Jan-2000
%   Now you need to set the path to WAFO root directory 
%   before using initwafo.
% by Par Johannesson 29-Sep-1999
%   new routine 

global WAFO_WSTATS_DEFAULT_PLOTFLAG
  
% Check input and output
ni = nargin;
%no = nargout;
%error(nargchk(0,3,ni));
narginchk(0,3)

% help('wafolicence')
% disp(' ')
% help('wafocitation')
if ni == 0, opt = []; end
if ni<2, add = []; end
if ni<3, first = []; end

% Set default value
if isempty(opt)
  opt = 'normal';
end

if isempty(add)
  add = 1;
end

if isempty(first)
  first=1;
end

% Add root-path to WAFO toolbox
wafop = waforoot;

path(wafop,path); % Add WAFO path first in search path

opt = lower(opt);  % To lower case
optNr = 1*strcmp(opt,'minimum') + 2*strcmp(opt,'normal') + ...
        3*strcmp(opt,'full');

p = {wafop};
if optNr >= 1 % Add paths to WAFO routines
  p1 = wafopath('tools');
  p(end+1:end+length(p1),1) = p1;
end
if optNr >= 2 % Add paths to WAFO demos, data and documentation
  p(end+1,1) = wafopath('docs'); 
  p(end+1,1) = wafopath('wdemos'); 
  p(end+1,1) = wafopath('data');
end
if optNr >= 3 % Add paths to WAFO papers
  tmp = wafopath('papers');
  Np = length(tmp(:,1));
  for ix=1:Np,
    p(end+1,1) = tmp(ix,:);
  end
end
    
if add % Add paths
  
  rmpath(wafop) 
  if first % Add WAFO paths first in search path
    p(end+1) = {'-begin'};
   else % Add WAFO paths last in search path (first == 0)
    p(end+1) = {'-end'};
  end
  addpath(p{:});
  disp(['WAFO path: ' wafop]);
  disp(['WAFO toolbox paths set: ' opt ' initiation']);
  
  WAFO_WSTATS_DEFAULT_PLOTFLAG = 0;
  
  msgId = 'MATLAB:divideByZero';
  warning('off',msgId);
else % Remove paths (add == 0) 

  % Remove all paths ecxept the path to wafo root-directory
  rmpath(p{2:end});
  disp(['WAFO toolbox paths removed: ' opt ' removal']);
  
end
