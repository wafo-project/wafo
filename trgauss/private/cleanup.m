function cleanup(varargin)
%CLEANUP Removes files from harddisk if they exist
%
% CALL cleanup(filename1,filename2)
% 
% CLEANUP works in the same way as delete does on files, except it 
% it give no warning if the file does not exist.
%
%

error(nargchk(1,inf,nargin))
ni = nargin;
for ix =1:ni
  fn  = varargin{ix};
  if exist(fn,'file'),
    delete(fn)
  end
end