function t = filecmp(file1, file2)
%FILECMP Compare two files (binary comparison).
%
%   FILECMP(FILE1, FILE2) returns 1 if both files FILE1 and FILE2 exist and
%   are identical and 0 otherwise.

%   Author:      Peter J. Acklam
%   Time-stamp:  2000-05-30 13:53:18
%   E-mail:      jacklam@math.uio.no
%   WWW URL:     http://www.math.uio.no/~jacklam

   % Check number of input arguments.
   error(nargchk(2, 2, nargin));

   % See if first file exists.
   if ~exist(file1, 'file')
      error([file1 ': No such file.']);
   end

   % See if second file exists.
   if ~exist(file2, 'file')
      error([file2 ': No such file.']);
   end

   % Try to open first file for reading.
   fid1 = fopen(file1, 'rb');
   if fid1 < 0
      error([file1 ': Can''t open file for reading.']);
   end

   % Try to open second file for reading.
   fid2 = fopen(file2, 'rb');
   if fid2 < 0
      fclose(fid1);
      error([file2 ': Can''t open file for reading.']);
   end

   % See if the two files are actually the same file.
   if fid1 == fid2
      fclose(fid1);
      error('The two file names refer to the same file.');
   end

   % Quick exit if the file sizes differ.
   dir1 = dir(fopen(fid1));
   dir2 = dir(fopen(fid2));
   if dir1.bytes ~= dir2.bytes
      t = 0;
      return
   end

   blocksize = 128;
   t = 1;

   while t

      [data1, count1] = fread(fid1, blocksize);
      [data2, count2] = fread(fid2, blocksize);

      if count1 ~= count2 || any( data1 ~= data2 )
         t = 0;
      end

      if count1 == 0 || count2 == 0;
         break
      end

   end

   fclose(fid1);
   fclose(fid2);
