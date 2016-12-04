function hexview(file)
%HEXVIEW View a file in hexadecimal mode.
%
%   HEXVIEW FILE displays the file FILE in hexadecimal mode with byte
%   offset and ascii strings.

%   Author:      Peter J. Acklam
%   Time-stamp:  2000-05-30 12:33:04
%   E-mail:      jacklam@math.uio.no
%   WWW URL:     http://www.math.uio.no/~jacklam

   % check number of arguments
   error(nargchk(1, 1, nargin));

   % see if the file exists
   if ~exist(file, 'file')
      error([file ': No such file.']);
   end

   % try to open the file for reading
   fid = fopen(file, 'r');
   if fid < 0
      error([file ': Can''t open file for reading.']);
   end

   blocksize = 16;      % number of bytes for each line
   dot       = '.';     % used for non-printable characters
   offset    = 0;       % file offset of data for current line

   while ~feof(fid)

      % read a block of data from the file
      [data, count] = fread(fid, blocksize);

      % print the offset and the hex values for the data
      fprintf('%08X', offset);
      fprintf(' %02X', data);

      % replace non-printable characters
      data( data < 32 | data > 126 ) = dot;

      % fill line if necessary (only used on last line if at all)
      if count < blocksize
         fprintf(repmat( '   ', 1, blocksize - count));
      end

      % print a space, the actual data and a newline
      fprintf(' ');
      fprintf('%c', data);
      fprintf('\n');

      % update offset value
      offset = offset + blocksize;

   end

   % close file
   fclose(fid);
