function  direc = getsubdirs(rootdir,N,direc)
%GETSUBDIRS Recursively find subdirectories of a root directory
%
% CALL: direc = getsubdirs(rootdir,N);
%
%  rootdir  = path to root directory (Default '.', i.e., current directory).
%  direc    = cell array of directories found including the root directory.
%  N        = maximum number of recursive calls (default 10).
%
% Example
%  d = getsubdirs(waforoot, 3);
%  names = {};
%  for i=1:length(d),
%   [folder, names{i}] = fileparts(d{i});
%  end
%  assert(names(2:4), {'@data_1d', '@data_2d' '@data_3d'})
%
% See also: dir 

% History
% by pab 12.10.2003

error(nargchk(0,3,nargin))  
if (nargin<1||isempty(rootdir)),
  rootdir = {'.'};
elseif ischar(rootdir)
  rootdir = cellstr(rootdir);
end
if (nargin<2||isempty(N)),
  N = 10;
end
if nargin<3||isempty(direc)
  direc = cell(1,0);
end

for iy=1:length(rootdir)
  if ( exist(rootdir{iy}) ==7), % directory
    direc{end+1} = rootdir{iy};
    if N>0,
      d = dir(rootdir{iy});  
      if (~isempty(d)) % pab bugfix
        d = {d([d.isdir]).name};
        d = {d{~ismember(d,{'.' '..'})}};
      end
      
      for j = 1:length(d)
        direc = getsubdirs({fullfile(rootdir{iy},d{j})},N-1,direc);
      end
    end
  else
    warning(sprintf('Directory not found: %s.\n',rootdir{iy}));
  end
end
return



