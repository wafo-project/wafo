function vernum = matfilever(file)
%MATFILEVER Return MAT-file version information.
%
%   MATFILEVER(FILE) returns the MAT-file version information for the
%   specified file.

%   Author:      Peter J. Acklam
%   Time-stamp:  2000-08-02 16:42:44
%   E-mail:      jacklam@math.uio.no
%   WWW URL:     http://www.math.uio.no/~jacklam

   % check number of arguments
   error(nargchk(1, 1, nargin));

   % append .mat suffix if not present
   k = find(file == '.', 1);
   if isempty(k)
      file = [file '.mat'];
   end

   % see if the file exists
   if ~exist(file, 'file')
      error([file ': No such file.']);
   end

   file = which(file);

   % try to open the file for reading
   fid = fopen(file, 'r');
   if fid < 0
      error([file ': Can''t open file for reading.']);
   end

   [data, count] = fscanf(fid, '%c', 8192);
   fclose(fid);

   k = find(~data);             % find null bytes
   str = char(data(1:k(1)));
   if nargout
      vernum = str;
   else
      disp(str);
   end
