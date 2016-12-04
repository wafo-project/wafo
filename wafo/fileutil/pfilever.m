function vernum = pfilever(file)
%PFILEVER Return P-file version number.
%
%   PFILEVER(FILE) returns the P-file version number for the specified
%   file.

%   Author:      Peter J. Acklam
%   Time-stamp:  2000-05-22 18:04:44
%   E-mail:      jacklam@math.uio.no
%   WWW URL:     http://www.math.uio.no/~jacklam

   % check number of arguments
   error(nargchk(1, 1, nargin));

   file = which(file);

   % see if the file exists
   if ~exist(file, 'file')
      error([file ': No such file.']);
   end

   % try to open the file for reading
   fid = fopen(file, 'r');
   if fid < 0
      error([file ': Can''t open file for reading.']);
   end

   [data, count] = fscanf(fid, '%c', 1);
   if ~isequal(data, char(0)) && ~isequal(data, char(1))
      fclose(fid);
      error([file ': File does not appear to be a P-file.']);
   end

   [data, count] = fscanf(fid, '%s', 1);
   if ~isequal(data, 'P-file');
      fclose(fid);
      error([file ': File does not appear to be a P-file.']);
   end

   [data, count] = fscanf(fid, '%f', 1);
   if count < 1
      fclose(fid);
      error([file ': File does not appear to be a P-file.']);
   end

   vernum = data;

   [data, count] = fscanf( fid, '%c', 1 );
   if ~isequal(data, char(0))
      fclose(fid);
      error([file ': File does not appear to be a P-file.']);
   end

   fclose(fid);
