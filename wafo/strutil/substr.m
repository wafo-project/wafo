function outstr = substr( str, offset, len )
%SUBSTR Extract a substring out of a string.
%
% CALL: substring = substr(string, offset, length ) 
% 
%  substring = character array extracted from string.
%  string    = character array.
%  offset    = starting index of substring in the string.
%  length    = length of substring (default length(string(offset:end)))
%
%   SUBSTR extracts a SUBSTRING out of STRING
%   with given LENGTH starting at the given OFFSET.  First character is
%   at offset 0.  If OFFSET is negative, starts that far from the end of
%   the string.  If LENGTH is omitted, returns everything to the end of
%   the string.  If LENGTH is negative, removes that many characters
%   from the end of the string.
%
%   SUBSTR is a Matlab version of the Perl operator with the same name.
%   However, unlike Perl's SUBSTR, no warning is produced if the
%   substring is totally outside the string.
%
% Examples:
%  s = 'How much wood would a Woodchuck chuck?';
%  assert(substr(s, 0, 1), 'H')    % Get first character
%  assert(substr(s, -1, 1), '?')   % Get last character
%  assert(substr(s,  1), 'ow much wood would a Woodchuck chuck?')     % Remove first character
%  assert(substr(s,  0, -1), 'How much wood would a Woodchuck chuck') % Remove last character
%  assert(substr(s,  1, -1), 'ow much wood would a Woodchuck chuck')  % Remove first and last character
%
% See also: extract

% History:
% Revised pab 13.12.2000
% -  updated header
%
%   Author:      Peter J. Acklam
%   Time-stamp:  2000-02-29 01:38:52
%   E-mail:      jacklam@math.uio.no
%   WWW URL:     http://www.math.uio.no/~jacklam

% Check number of input arguments.
%error( nargchk( 2, 3, nargin ) );
narginchk(2,3)
n = length(str);

% Get lower index.
lb = offset+1;               % From beginning of string.
if offset < 0
  lb = lb+n;                % From end of string.
end
lb = max( lb, 1 );

% Get upper index
if nargin == 2               % SUBSTR( STR, OFFSET )
  ub = n;
elseif nargin == 3           % SUBSTR( STR, OFFSET, LEN )
  if len < 0
    ub = n+len;
  else
    ub = lb+len-1;
  end
  ub = min( ub, n );
end

% Extract substring.
outstr = str( lb : ub );
