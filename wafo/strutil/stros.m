function os = stros(str)
%STROS  Return the name of the OS under which the text was created.
%
%  CALL: os = stros(str);
%
%   os  = string containing the name of the operation system under which
%         the text was created.
%   str = string vector
%
%   STROS will return either 'unix', 'dos' or 'mac' depending on what
%   kind of line-ending is used in the string STR.  If the OS was not
%   determined, and empty string is returned.
%
%   Example
%
%   assert(stros(['Hello' char(10) 'world!']), 'unix')
%   assert(stros(['Hello' char(13) char(10) 'world!']), 'dos')
%   assert(stros(['Hello' char(13) 'world!']), 'mac')
%   assert(stros(['Hello world!']), '')
%
% See also: computer

% Revised pab 13.12.2000
% - updated header
%   Author:      Peter J. Acklam
%   Time-stamp:  2000-07-17 00:36:11
%   E-mail:      jacklam@math.uio.no
%   WWW URL:     http://www.math.uio.no/~jacklam

   % Check number of input arguments.
   %error(nargchk(1, 1, nargin));
   narginchk(1,1)
   % Find DOS line endings (CR+LF)
   k = find(str(1:end-1) == 13);        % find CRs
   if ~isempty(k)
      k = k(str(k+1) == 10);            % find following LFs
   end

   % Count number of DOS line endings and remove them.
   ncrlf = length(k);                   % number of CR+LFs
   if ncrlf
      str([k k+1]) = [];
   end

   % Count number of UNIX line endings (LF) and MAC line endings (CR).
   nlf = sum(str == 10);
   ncr = sum(str == 13);

   % Now get the name of the OS.
   if ncrlf && ~nlf && ~ncr               % only CR+LFs
      os = 'dos';
   elseif ~ncrlf && nlf && ~ncr           % only LFs
      os = 'unix';
   elseif ~ncrlf && ~nlf && ncr           % only CRs
      os = 'mac';
   elseif ~ncrlf && ~nlf && ~ncr          % no CRs or LFs
      os = '';
   else         % several different kinds of line-endings are used
      warning('Several different kinds of line-endings are used.');
      disp('')
      oses = {'dos' 'unix' 'mac'};
      [dummy idx] = sort([ncrlf nlf ncr]);      % sort frequencies
      os = oses{idx(end)};                      % return most likely OS
   end
