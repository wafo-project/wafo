function wd = cdtomfile(file)
%CDTOMFILE Change working directory to where a specified m-file is located.
%   CDTOMFILE FILE or CDTOMFILE(FILE) changes the current working directory
%   to the directory where FILE is located.  If FILE is a directory,
%   CDTOMFILE acts like CD and makes FILE the current working directory.
%   CDTOMFILE, by itself, prints out the current directory.
%
%   WD = CDTOMFILE returns the current directory as a string.
%
% Example
%  p0 = pwd();
%  cdtomfile('cdtomfile');
%  p1 = pwd();
%  assert(p1, fullfile(waforoot, 'fileutil'));
%  cd(p0);
%
% See also CD, PWD.

%   Author:      Peter J. Acklam
%   Time-stamp:  2000-07-19 15:11:45
%   E-mail:      jacklam@math.uio.no
%   WWW URL:     http://www.math.uio.no/~jacklam

error(nargchk(0, 1, nargin));

if ~nargin
   if nargout
      wd = cd;
   else
     switch get(0, 'FormatSpacing')
       case 'loose'
         fprintf('\n%s\n\n', cd);
       case 'compact'
         fprintf('%s\n', cd);
       otherwise
         error('Unrecognized ''FormatSpacing'' root property value.');
     end
   end
   return
end

if exist(file, 'file')                  % if argument is a file on the path
   location = which(file);              % get full path to file
   if strcmp(location, 'built-in')      % if built-in
      location = which([file '.m']);    %   append '.m' and retry
   end
   path = fileparts(location);          % extract path
   cd(path);
elseif exist(file, 'dir')
   cd(file);
else
   error('"%s" not found', file);
end

if nargout
   wd = cd;
end
