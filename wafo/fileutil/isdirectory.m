function [t, newdir] = isdirectory(dir)
%ISDIRECTORY True if argument is a directory.
%   ISDIRECTORY(DIR) returns 1 if DIR is a directory and 0 otherwise.
%   [T, NEWDIR] = ISDIRECTORY(DIR) also returns the name of the
%   directory as returned by CD.  NEWDIR is empty if DIR is not a
%   directory.
%
%   Using EXIST(DIR, 'dir') is too unreliable.  For instance, under
%   Windows it returns false if DIR is 'c:' (where c: is an existing
%   drive) and it returns false and if DIR is a directory which is in
%   the path but does not exist.

%   Author:      Peter J. Acklam
%   Time-stamp:  1998-12-30 01:26:59
%   E-mail:      jacklam@math.uio.no
%   WWW URL:     http://www.math.uio.no/~jacklam

   cwd    = cd;
   t      = 1;
   newdir = '';

   try
      cd(dir);          % Try to chdir to the specified directory.
      newdir = cd;      % If that succeeded, get the new directory name.
   catch
      t = 0;            % chdir failed, so the directory does not exist.
   end

   if t                 % If the chdir succeeded
      cd(cwd);          %   chdir back to where we started from.
   end
