function pathlist = path2cell(pathstr)
%PATH2CELL Convert search path to cell array.
%
%   PATH2CELL returns a cell array where each element is a directory
%   in the search path.
%
%   PATH2CELL(MYPATH) converts MYPATH rather than the search path.
%
%   Empty directories are removed, so under UNIX, PATH2CELL('foo::bar')
%   returns { 'foo' 'bar' } rather than { 'foo' '' 'bar' }.

%   Author:      Peter J. Acklam
%   Time-stamp:  2000-02-24 15:19:09
%   E-mail:      jacklam@math.uio.no
%   WWW URL:     http://www.math.uio.no/~jacklam

   % Check number of input arguments.
   error(nargchk(0, 1, nargin));

   % Use Matlab's search path if no input path is given.
   if ~nargin
      pathstr = path;
   end

   k = find( pathstr == pathsep );         % Find path separators.
   k = [ 0 k length(pathstr)+1 ];          % Find directory boundaries.
   ndirs = length(k)-1;                    % Number of directories.
   pathlist = cell(0);                     % Initialize output argument.
   for i = 1:ndirs
      dir = pathstr( k(i)+1 : k(i+1)-1 );  % Get i'th directory.
      if ~isempty(dir)                     % If non-empty...
         pathlist{end+1,1} = dir;          % ...append to list.
      end
   end
