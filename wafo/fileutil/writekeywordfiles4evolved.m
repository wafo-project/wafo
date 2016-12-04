function writekeywordfiles4evolved(mdirs)
%writekeywordfiles4evolved
% 
% CALL: writekeywordfiles4evolved(mdirs)
%
% mdirs     = cellarray or character array with directories to search for
%             M-files. (default '.')
% 
% Example % 
%   writekeywordfiles4evolved('signal')
%
% See also

% References
%

% Tested on: Matlab 7.1
% History:
% 

% Author: pab
% Date: 12.07.2006
% Time: 15:02:51
% Version: 1.0


if nargin<1
	mdirs = '.';	
end

filename = 'matlab';
fid = fopen(sprintf('%s.kw', filename),'w');
fprintf(fid,'Name: Matlab\n');
fprintf(fid,'Description: Matlab m-files\n');
fprintf(fid,'LineComment: \%% \n');
fprintf(fid,'Extensions: m\n');
fprintf(fid,'CaseSensitive: 1\n');
fprintf(fid,'StringChar: %s\n' ,'''');
fprintf(fid,'WholeWord: 1 \n');
fprintf(fid,'Keyword1:');

keywords = {...
    'break',
    'case',
    'catch',
    'continue',
    'else',
    'elseif',
    'end',
    'for',
    'function',
    'global',
    'if',
    'otherwise',
    'persistent',
    'return',
    'switch',
    'try',
    'while',
    };
fprintf(fid,'%s ',keywords{:});
recursive=1;
mfiles = getmfiles(mdirs,recursive);
funcNames = strrep(mfiles,'.m','');

fprintf(fid,'\nKeyword2:');
k = find(isfunction(funcNames) & ~ismember(funcNames,keywords));
if any(k)
	fprintf(fid,'%s ',funcNames{k});
end


fclose(fid)
fid = fopen(sprintf('%s.a',filename),'w');
for ix = k(:).'
[T,R] = strtok(geth1line(funcNames{ix}));
fprintf(fid,'%s %s\n',funcNames{ix},R);
end
fclose(fid)





function mfiles = getmfiles(mdirs,recursive)
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
% Examples:
%  mfiles = getmfiles;
%  mfiles = getmfiles('signal')  
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
    status = exist(mdirs{i},'file');
    if (status == 2),     % M-file
      mfiles{end+1} = mdirs{i};
    elseif (status == 7), % Directory
      w = what(mdirs{i});
      w = w(1); %- Sometimes an array is returned...
      for j=1:length(w.m)
        mfiles{end+1} = w.m{j};
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