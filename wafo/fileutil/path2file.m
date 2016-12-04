function path2file(file, pathstr)
%PATH2FILE Write the Matlabpath to a file, one directory on each line.
%
%   PATH2FILE FILE writes Matlab's search path to the file FILE with
%   one directory on each line.
%
%   PATH2FILE(FILE, MYPATH) writes the specified path rather than the
%   search path.

%   Author:      Peter J. Acklam
%   Time-stamp:  2000-08-08 19:55:42
%   E-mail:      jacklam@math.uio.no
%   WWW URL:     http://www.math.uio.no/~jacklam

   % Check number of input arguments.
   error(nargchk(1, 2, nargin));

   % Use Matlab's search path if no input path is given.
   if nargin < 2
      pathstr = path;
   end

   dirs = path2cell(pathstr);
   n = length(dirs);

   fid = fopen(file, 'wt');
   if fid < 0
      error([file ': Unable to open file "' file '" for writing.']);
   end
   for i = 1:n
      if ~isempty(dirs{i})                      % If non-empty...
         fprintf(fid, dirs{i});                 % ...append to list.
      end
   end
   fclose(fid);
