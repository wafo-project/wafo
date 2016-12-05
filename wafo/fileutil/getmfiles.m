function mfiles = getmfiles(mdirs,recursive,mfiles, exclude,exclude_patterns)
%GETMFILES Return M-files found from a list of directories.
%
% CALL:  mfiles = getmfiles(mdirs,recursive) 
%  
% mfiles    = cellarray of M-files found in directories MDIRS.
% mdirs     = cellarray or character array with directories to search for
%             M-files. (default '.')
% recursive = 1 if a recursive search should be used
%             0 if only MDIRS should be searched (default)
%
% Example:
%  
%  mfiles = getmfiles(fullfile(waforoot, 'fileutil'));
%  names = {};
%  for i=1:length(mfiles),
%   [folder, names{i}] = fileparts(mfiles{i});
%  end
%  assert(names(1:3), {'Contents','bindiff','cdtomfile'});
% 
% See also: getsubdirs
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -*- Mode: Matlab -*- %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% getmfiles.m --- 
%% Author          : Per Andreas Brodtkorb
%% Revised On      : Sat Nov 01 18:31:28 2003
%% Last Modified By: Per Andreas Brodtkorb
%% Last Modified On: Wed Dec 10 11:54:16 2003
%% Update Count    : 8
%% Status          : Unknown, Use with caution!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  error(nargchk(0,3,nargin))
  if (nargin<5 || isempty(exlude_patterns))
    exlude_patterns = cell(1,0);
  end
  if (nargin<4 || isempty(exclude))
    exclude = cell(1,0);
  end
  if (nargin<3 || isempty(mfiles))
    mfiles = cell(1,0);
  end
  if (nargin<2 || isempty(recursive))
    recursive = 0;
  end
  if (nargin<1 || isempty(mdirs)),
    mdirs = {pwd};
  elseif ischar(mdirs),
    mdirs = cellstr(mdirs);
  end
  
  numDir = length(mdirs);
  for i = 1:numDir
    status = exist(mdirs{i});
    if (status == 2),     % M-file
      mfiles{end+1} = mdirs{i};
    elseif (status == 7), % Directory
      w = what(mdirs{i});
      w = w(1); %- Sometimes an array is returned...
      for j=1:length(w.m)
        mfiles{end+1} = fullfile(mdirs{i},w.m{j});
      end
      if recursive
        d = dir(mdirs{i});
        if (~isempty(d)) % pab bugfix
          d = {d([d.isdir]).name};
          d = {d{~ismember(d,{'.' '..'})}};
        end
        for j = 1:length(d)
          mfiles = getmfiles({fullfile(mdirs{i},d{j})}, recursive, mfiles);
        end
      end
    else
      fprintf('Warning: Unprocessed file %s.\n',mdirs{i});
    end
  end