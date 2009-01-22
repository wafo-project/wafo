function p = wafopath(part,add,first)
% WAFOPATH Adds or removes WAFO subdirectories from the search path.
%
% CALL:  wafopath(part,add,first) 
% 
%   part  = 'tools'  : All necessary WAFO routines. (default)
%           'docs'   : Documentation.
%           'wdemos'  : Demos.
%           'data'   : Data sets.
%           'papers' : Papers.
%   add   = 1 : Add paths from MATLAB search path (default)
%           0 : Remove paths from MATLAB search path 
%   first = 1 : Add the paths first in the search path (default)
%           0 : Add the paths last in the search path
%
% See also  initwafo, waforoot, wafoexepath.

%          NOTE: User no longer has to edit this file (since v1.1.14).
%          It is now automatic.

% Tested on: Matlab 5.3
% History:
% Revised pab 2007
% -added path to mex-directory
% Revised pab Feb 2007
% - Fixed a bug: Does not add class directories anymore.
% Revised ny PJ 10-Apr-2001
% - Now adds wafo/papers to path, since Contents-file added to PAPERS
% revised pab 09.12.2000
% - made sure all paths are in lower case letters
% revised jr 23.05.2000
% - changed demos -> wdemos, sim -> wsim, 
% - removed lines related to fplot
% revised pab 28.01.2000
% - changed 'papers' option: to include the actual directory of the paperscripts
% - changed 'tools'  option: only the tools which is actually installed on 
%                            the system is added to the matlab path
% Changed by Pär Johannesson 12-Jan-2000
%   Path 'exec/...' is now first in search path.
%   Searches for mex-files first, then m-files.
% Changed by Pär Johannesson 28-Sep-1999
%   Added input argument 'part'  
% by pab 11.08.99
%   new routine 

% Check input and output
ni = nargin;
no = nargout;
error(nargchk(0,3,ni));

% Set default values
if ni<1||isempty(part),  part  = 'tools'; end
if ni<2||isempty(add),   add   = 1;       end
if ni<3||isempty(first), first = 1;       end


%pref=[pathsep, waforoot];
%p = [waforoot filesep];
proot = [waforoot filesep];


if strcmp(part,'tools')
  
  p = {wafomexpath}; % Search for mex-files first
					
  w      = dir(proot);
  inddir = find(cat(1,w.isdir));
  wdir   = {w(inddir).name};
  % remove non-tools directories from dir list
  wdir(1:2) = []; % remove '.' and '..' from dir list
  ind = strmatch('exec',wdir);   wdir(ind)=[];
  ind = strmatch('docs',wdir);   wdir(ind)=[];
  ind = strmatch('wdemos',wdir); wdir(ind)=[];
  ind = strmatch('papers',wdir); wdir(ind)=[];
  ind = strmatch('data',wdir);   wdir(ind)=[];
  ind = strmatch('source',wdir); wdir(ind)=[];
  ind = strncmpi(wdir,'@',1);    wdir(ind)=[]; % remove class directories
  N   = length(wdir);
  
  % adds only those directories which is actually installed on the
  % current system
  for ix=1:N
    p(ix+1,1) = {[proot lower(wdir{ix})]};
  end
  
  
elseif strcmp(part,'docs')
  p = {[proot 'docs']};
elseif strcmp(part,'wdemos')
  p = {[proot 'wdemos']};
elseif strcmp(part,'data')
  p = {[proot 'data']};
elseif strcmp(part,'papers')
  w=dir(fullfile(proot,'papers'));
  w=w(3:end); % remove '.' and '..' from dir list
  N=length(w);
  if N<1
    disp('No paper scripts are installed')
  else
    p = {[proot 'papers']}; % PJ 10-Apr-2001
    ix=0;
    iy=1;
    while ix<N
      ix=ix+1;
      if w(ix).isdir
        iy=iy+1;
        p(iy,1) = {[proot 'papers' filesep lower(w(ix).name)]};
        %else
        % Commented by PJ 10-Apr-2001
        %	disp(['There should be no m-files in:  ' proot 'papers'])
        %	disp('only directories')
      end
    end
  end
end

if no == 0
  if add 
    if first
      addpath(p{:},'-begin');
    else
      addpath(p{:},'-end');
    end
    disp('Pathnames Successfully Set')
  else
    rmpath(p{:})
    disp('Pathnames Successfully Removed')
  end 
end
  
