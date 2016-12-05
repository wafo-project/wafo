function locations = locate(varargin)
%LOCATE Locate one or more files in the search path.
%
%   LOCATE FILEGLOB prints a list of all Matlab files on the path that
%   match FILEGLOB. Private subdirectories are also searched.
%
%   LIST = LOCATE(FILEGLOB) returns a list (cell array) of the files.  The
%   list is not displayed.
%
%   Multiple filegobs may be given, e.g., LOCATE('*.m', '*.dll').
%
%   Note that the globbing is done by Matlab and not by the operating
%   system.
%
%   Examples:
%
%   loc0 = locate('*plot.m');    % Find all files ending with 'plot.m'.
%   loc1 = locate('im*');        % Find all files starting with 'im'.
%   names = {};
%   for i=1:length(loc0),
%    [folder, names{i}] = fileparts(loc0{i});
%   end
%   i0 = [strmatch('pdfplot', names)(1), strmatch('trplot', names)(1), ...
%         strmatch('nplot', names)(1)]; 
%   assert(names(i0), {'pdfplot','trplot','nplot'});
%
%   See also DIR, SYSGLOB.

%   Author:      Peter J. Acklam
%   Time-stamp:  2001-09-24 19:48:17 +0200
%   E-mail:      jacklam@math.uio.no
%   URL:         http://www.math.uio.no/~jacklam

   % check number of input arguments
   nargsin = nargin;
   error(nargchk(1, Inf, nargsin));

   dirs = path2cell;            % convert path to list
   startdir = cd;               % get starting directory

   display = 1;
   if nargout                   % if there are output arguments
      locations = {};           %   initialize output list
      display = 0;              %   and don't display results
   end

   while ~isempty(dirs)           % while there are unprocessed dirs...

      directory = dirs{1};      % get first directory
      dirs = dirs(2:end);       %   and remove it from the list

      % fprintf('%s\n', directory);
      cd(directory);            % chdir to the directory

      % get the list of file names that match the glob(s)
      found = {};
      for i = 1:nargsin
         % dirinfo = dir(fullfile(directory,varargin{i}));
         dirinfo = dir(varargin{i});
         found = { found{:} dirinfo.name };
      end

      cd(startdir);             % chdir back to starting directory

      % Append list of files found in this directory to the main list or
      % display the files if no output is wanted.
      if ~isempty(found)
         % get unique list of files
         found = unique(found);

         % prepend directory name
         for i = 1:length(found)
            found{i} = fullfile(directory, found{i});
         end
         if nargout
            locations = [ locations(:) ; found(:) ];
         else
            fprintf('%s\n', found{:});
         end
      end

      % If this directory has a private subdirectory, look there too.
      subdirectory = fullfile(directory, 'private');
      if exist(subdirectory, 'dir')
         dirs = { subdirectory , dirs{:} }';
      end

   end
