function outpathstr = trimpath(inpathstr)
%TRIMPATH Trim Matlab's search path.
%
%   NEWPATH = TRIMPATH(INPATH) trims the given path INPATH and returns
%   the trimmed version NEWPATH.  The trimming includes removing empty
%   directory strings, duplicate directories, non-existing directories
%   and trailing file separators from each directory.
%
%   If no input path is given, the current search path is used.  If no
%   output argument is specified, the new path is set as the new
%   search path.

%   Author:      Peter J. Acklam
%   Time-stamp:  2000-02-24 15:10:32
%   E-mail:      jacklam@math.uio.no
%   WWW URL:     http://www.math.uio.no/~jacklam

   % Check number of input arguments.
   error(nargchk(0, 1, nargin));

   % Use Matlab's search path if no input path is given.
   if nargin == 0
      inpathstr = path;
   end
   origpathstr = inpathstr;                % Save for later comparison.

   ps = pathsep;
   k = find( inpathstr == ps );            % Find path separators.
   k = [ 0 k length(inpathstr)+1 ];        % Find directory boundaries.
   ndirs = length(k)-1;                    % Number of directories.

   % Build list of directories for new path.  Skip directories that are
   % empty or that do not exist or that are already in the path.
   newpathlist = {};
   for i = 1:ndirs

      dir = inpathstr( k(i)+1 : k(i+1)-1 );        % Get i'th directory.

      fprintf( '''%s''', dir );

      if isempty( dir )
         fprintf( ' skipped (empty)\n' );

      elseif any(strcmp(dir, newpathlist))
         fprintf(' skipped (already in path)\n');

      else

         [t, newdir] = isdirectory(dir);
         if t
            if strcmp(dir, newdir)
               fprintf(' OK\n');
            else
               fprintf(' changed to ''%s''\n', newdir);
            end
            newpathlist{end+1} = newdir;
         else
            fprintf(' skipped (does not exist)\n');
         end

      end

   end

   % Join the list inserting a path separator between each directory.
   if isempty(newpathlist)
      newpathstr = '';
   else
      newpathstr = newpathlist{1};
      for i = 2:length(newpathlist)
         newpathstr = [ newpathstr ps newpathlist{i} ];
      end
   end

   if strcmp(newpathstr, origpathstr)
      disp('No changes were made to the path.');
   else
      disp('Changes were made to the path.');
   end

   % Return the new path or make it the new search path.
   if nargout
      outpathstr = newpathstr;     % Return path.
   else
      path(newpathstr);            % Set as new Matlab search path.
   end
