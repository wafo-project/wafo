function res = strlexcmp(a, b)
%STRLEXCMP Lexicographic comparison of two strings.
%
% CALL:  status = strlexcmp(str1,str2)
%
%   status    = integer (-1, 0, or 1)
%   str1,str2 = character arrays
%
%   STRLEXCMP returns -1, 0, or 1 depending on whether STR1
%   is stringwise less than, equal to, or greater than STR2.
%
%   This is a Matlab version of the Perl `cmp' operator.
%
%   Example
%
%   assert(strlexcmp('Hello world', 'Hello space'), 1)
%   assert(strlexcmp('Hello space','Hello world'), -1)
%   assert(strlexcmp('Hello space','Hello space '), -1)
%   assert(strlexcmp('Hello space','Hello space'), 0)

%   Author:      Peter J. Acklam
%   Time-stamp:  2000-07-17 02:06:40
%   E-mail:      jacklam@math.uio.no
%   WWW URL:     http://www.math.uio.no/~jacklam

   % Check arguments.
   %error(nargchk(2, 2, nargin));
   narginchk(2,2)
   if ~ischar(a) || ~ischar(b)
      error('Both arguments must be strings.');
   end

   % Get lengths of strings.
   na = length(a);
   nb = length(b);
   n = min(na, nb);

   % Find characters that differ.
   k = find( a(1:n) ~= b(1:n) );
   if isempty(k)
      % All characters are identical.  Compare lengths.
      if na < nb
         res = -1;
      elseif na > nb
         res = 1;
      else
         res = 0;
      end
   else
      % Compare first character that is different.
      k = k(1);
      if a(k) < b(k)
         res = -1;
      elseif a(k) > b(k)
         res = 1;
      else
         res = 0;
      end
   end

