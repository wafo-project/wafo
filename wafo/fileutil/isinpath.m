function t = isinpath(varargin)
%ISINPATH True for directories in path
%
%   ISINPATH(DIR1, DIR2, ...) returns an array that contains 1's for the
%   directories that are in Matlab's search path and 0's for those that
%   are not.

%   Author:      Peter J. Acklam
%   Time-stamp:  2000-02-18 13:08:44
%   E-mail:      jacklam@math.uio.no
%   URL:         http://www.math.uio.no/~jacklam

   ispc = strncmp(computer, 'PC', 2);

   % Get path and wrap it up in path separators.
   sep = pathsep;
   pth = [ sep path sep ];

   % We do case-insensitive match on PCs, so convert to lowercase.
   if ispc
      pth = lower(pth);
   end

   t = false(1, nargin);       % Initialize output array.
   for i = 1:nargin
      dir = varargin{i};
      if ispc
         dir = lower(dir);              % Lowercase if we're on a PC.
      end
      t(i) = ~isempty(findstr(pth, [sep dir sep]));
   end
