function cleanup(varargin)
%CLEANUP Removes files from harddisk if they exist
%
% CALL cleanup(filename1,filename2)
% 
% CLEANUP works in the same way as delete does on files, except it 
% it give no warning if the file does not exist.
%
% Modified by GL, 12-17-2015, replaced exist by strmatch

error(nargchk(1,inf,nargin))
ni = nargin;
inmap = ls;
for ix =1:ni
  fn  = varargin{ix};
  xfn = strmatch(fn,inmap);
  if ~isempty(xfn),
    delete(fn)
  end
end
