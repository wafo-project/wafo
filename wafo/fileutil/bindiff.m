function patch_data = bindiff( file1, file2 )
%BINDIFF Binary difference. Returns patchdata for BINPATCH.
%
%   PATCH_DATA = BINDIFF(FILE1, FILE2) does a binary comparison of the two
%   files FILE1 and FILE2 and returns a PATCH_DATA matrix.  PATCH_DATA is a
%   matrix with three columns for address (file offset), old_value, and
%   new_value.
%
%   The two files FILE1 and FILE2 must have the same size.
%
%   BINPATCH(FILE1, PATCH_DATA) will cause FILE1 to be identical to FILE2 if
%   the patching succeeded.
%
%   See also BINPATCH.

%   Author:      Peter J. Acklam
%   Time-stamp:  2000-05-30 13:55:26
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
      error('Files have different sizes.');
   end

   blocksize = 32*1024;
   patch_data = [];
   block_num = 0;

   while ~feof(fid1) && ~feof(fid2)

      block_num = block_num + 1;
      [data1, count1] = fread(fid1, blocksize, 'uchar');
      [data2, count2] = fread(fid2, blocksize, 'uchar');

      if count1 ~= count2
         fclose(fid1);
         fclose(fid2);
         error('Files have different sizes.');
      end

      k = find( data1 ~= data2 );
      if ~isempty(k)
         offset = k - 1 + (block_num-1)*blocksize;
         patch_data = [ patch_data ; offset data1(k) data2(k) ];
      end

   end

   fclose(fid1);
   fclose(fid2);

   patch_data = [ patch_data ; -sum(patch_data) ];
