function s = rstrrep(s1, s2, s3)
%RSTRREP Recursive string replacement.
%
%   S = RSTRREP(S1, S2, S3) replaces the first occurence of S2 in S1 with
%   S3.  Then the process is repeated until S2 no longer exists in S1.
%
%   Example:
%     s1 = 'This   is  a    good   example';
%     assert(strrep(s1, '  ', ' '), 'This  is a   good  example')
%     assert(rstrrep(s1, '  ', ' '), 'This is a good example')

%   Author:      Peter J. Acklam
%   Time-stamp:  1998-06-22 20:46:56
%   E-mail:      jacklam@math.uio.no
%   WWW URL:     http://www.math.uio.no/~jacklam

   % check number of input arguments
   %error(nargchk(3, 3, nargin));
   narginchk(3,3)
   % quick exit if no change will be made
   if ( length(s2) > length(s1) ) || isempty(s2) || isempty(s1)
      s = s1;
      return;
   end

   % get length of string to replace and replacement string
   ns2 = length(s2);
   ns3 = length(s3);

   done = 0;
   while ~done
      k = findstr(s1,s2);                  % Find every occurence of s2 in s1.
      if isempty(k)                        % If s2 does not exist in s1...
         done = 1;                         % ...we're done.
         s = s1;
      else
         ns1 = length(s1);                 % Length of old string.
         ns = ns1-ns2+ns3;                 % Length of new string.
         k = k(1);                         % First occurence of s2 in s1.
         s = char(zeros(1,ns));            % Initialise new string.
         s(1:k-1) = s1(1:k-1);             % Insert first part.
         s(k:k+ns3-1) = s3;                % Insert s3.
         s(k+ns3:ns) = s1(k+ns2:ns1);      % Insert last part.
         s1 = s;
      end
   end
