function splitfile(infile, outfilebase, chunksize, opt)
%SPLITFILE Split a file into fixed-size pieces.
%
%   SPLITFILE(INFILE, OUTFILEBASE, CHUNKSIZE) splits the file INFILE
%   into files of size CHUNKSIZE.  The output files are named
%   outfilebase.000, outfilebase.001, etc.
%
%   SPLITFILE('INFILE', 'OUTFILEBASE', CHUNKSIZE, '-PAUSE') where will
%   cause the program to pause after each output file has been created.
%
%   For example, split a file directly to an 1.4Mb PC floppy drive
%
%      splitfile bigfile.tar.gz a:\bigfile 1457664 -pause
%
%   The pause flag will make it possible to change floppy disks between
%   each output file is written.

%   Author:      Peter J. Acklam
%   Time-stamp:  2000-02-24 14:25:05
%   E-mail:      jacklam@math.uio.no
%   WWW URL:     http://www.math.uio.no/~jacklam

   % Check number of input arguments.
   error(nargchk(3, 4, nargin));

   pauseflag = nargin == 4;

   if isempty(infile) || ~ischar(infile)
      error('Input file name must be a non-empty string.');
   end

   if isempty(outfilebase) || ~ischar(outfilebase)
      error('Output file basename must be a non-empty string.');
   end

   if ischar(chunksize)
      chunksize = sscanf(chunksize, '%d');
   end

   if isempty(chunksize)
      error('Chunksize is empty.');
   end

   % Size of block read (and written) at a time.
   blocksize = 8192;

   % Open input file for reading.
   ifid = fopen(infile, 'r');
   if ifid < 0
      error([ 'Can''t open file "' infile '" for reading.' ]);
   end

   % Initialize output file counter and start splitting.
   outfilecount = 0;

   while ~feof(ifid)

      % Create and open output file name and initialize counter.
      outfile = sprintf('%s.%03d', outfilebase, outfilecount);
      ofid = fopen(outfile, 'wb');
      bytes_written = 0;
      fprintf('Writing "%s"', outfile);

      % Write output file.
      while (bytes_written < chunksize) && ~feof(ifid)
         readnow = min(blocksize, chunksize-bytes_written);
         data = fread(ifid, readnow);
         count = fwrite(ofid, data);
         bytes_written = bytes_written + count;
      end
      fclose(ofid);
      fprintf(', %d bytes written\n', bytes_written);

      if pauseflag && ~feof(ifid)
         disp('Press any key to continue...');
         pause;
      end

      % Increment output file counter.
      outfilecount = outfilecount + 1;

   end

   fclose(ifid);
