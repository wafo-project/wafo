function list = file2cell(file, opt)
%FILE2CELL Convert a text file to a cell array of lines.
%
%   FILE2CELL(FILENAME) returns a cell array where the i'th cell
%   contains the i'th line of the specified file.
%
%   FILE2CELL(FILENAME, OPT) might be used to remove empty cells from
%   the returned cell array.  Recognized values for OPT are
%      0 - no change
%      1 - remove leading empty cells
%      2 - remove trailing empty cells
%      3 - remove leading and trailing empty cells
%      4 - remove all empty cells
%
% Example
%
% list = file2cell('file2cell.m');
% assert(list{1},'function list = file2cell(file, opt)');
%

%   Author:      Peter J. Acklam
%   Time-stamp:  2000-05-30 12:41:48
%   E-mail:      jacklam@math.uio.no
%   WWW URL:     http://www.math.uio.no/~jacklam

   % Check number of input arguments.
   error(nargchk(1, 2, nargin));
   if nargin < 2
      opt = 0;
   end

   % See of file exists.
   if ~exist(file, 'file')
      error([file ': No such file.']);
   end

   % Try to open file for reading.
   fid = fopen(file, 'r');
   if fid < 0
      error([file ': Can''t open file for reading.']);
   end

   % Initialize output list and line counter.
   list = {};
   lineno = 0;

   % Read data from file.
   while ~feof(fid)
      lineno = lineno + 1;              % increment line counter
      list{lineno,1} = fgetl(fid);   % add line to list
   end

   % Close file
   fclose(fid);

   % Remove empty elements from output list.
   if ~isempty(opt) && ~isequal(opt, 0)
      k = find(~cellfun('isempty', list));
      if isequal(opt, 1)
         list = list(k(1):end);       % remove leading empty elements
      elseif isequal(opt, 2)
         list = list(1:k(end));
      elseif isequal(opt, 3)
         list = list(k(1):k(end));
      elseif isequal(opt, 4);
         list = list(k);              % remove all empty elements
      else
         warning('Unrecognized value for OPT, so ignored.');
      end
   end
