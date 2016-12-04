function  direc=findsubdir(rpwd,N)
% FINDSUBDIR Recursively find subdirectories of a root directory
%
% CALL: direc = findsubdir(rpwd,N);
%
%  rpwd  = path to root directory (Default current directory).
%  direc = character array of directories found including the root directory.
%  N     = maximum number of recursive calls (default 10).
%
% See also: dir 

% History
% by pab 12.01.2000

if nargin<1||isempty(rpwd),
  rpwd = pwd;
end
if nargin<2||isempty(N)
  N=10;
end

direc = rpwd;
for iy=1:size(rpwd,1)
  w = dir(deblank(rpwd(iy,:)));
  if isempty(w), direc=[];warning([rpwd, ' not found']),end
  
  inddir = find(cat(1,w.isdir));
  if ~isempty(inddir),
    inddir(1:2)=[]; % remove '.' and '..' from dir list
  end
  if N>0 && ~isempty(inddir)
    for ix=inddir(:).',
      subdir = fullfile(rpwd,w(ix).name,'');
      direc  = strvcat(direc,findsubdir(subdir,N-1));
    end
  end
end
return


