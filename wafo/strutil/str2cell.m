function c = str2cell( varargin )
%STR2CELL Convert a string to a cell array of lines.
%
%   C = STR2CELL( STR ) creates a cell array C where each cell contains
%   a line of the string STR.
%
%   C = STR2CELL( STR, OPTS ), where OPTS is 'L', 'T' or both, removes
%   leading and/or trailing blank lines from the string before
%   converting to a cell array.
%
%   If the string contains LFs (linefeed characters, decimal 10), the
%   input string is split at their position after all CRs (carriage
%   return characters, decimal 13) have been removed. If there are no
%   LFs, the string is split at the position of the CRs. This should
%   ensure that the string is split correctly with both UNIX (LF), DOS
%   (CR+LF) and MAC (CR) definitions of a newline.
%
%  Example
%
%  s = ['Hello world' char(10) 'Hello space'];
%  assert(str2cell(s), {'Hello world'; 'Hello space'})
%

%   Author:      Peter J. Acklam
%   Time-stamp:  1998-06-30 21:20:01
%   E-mail:      jacklam@math.uio.no
%   WWW URL:     http://www.math.uio.no/~jacklam

%error( nargchk( 1, 3, nargin ) );
narginchk(1,3)
%
% Assign default values to parameters that can be changed by command
% line options.
%
strip_lead  = 0;
strip_trail = 0;

%
% Process command line options.
%
while length( varargin ) > 1
   opt = varargin{2};
   if ~ischar( opt )
      error( 'Options must be strings.' );
   end
   switch opt
      case { 'l', 'L' }
         strip_lead = 1;
      case { 'u', 'U' }
         strip_trail = 1;
      otherwise
         error( [ 'Unknown option: ' opt ] );
   end
   varargin(2) = [];
end
str = varargin{1};

%
% Strip leading blank lines.
%
if strip_lead
   k = find( ~isspace(str) );
   if ~isempty(k)
      k = min(k);
      str = str(k:end);
   end
end

%
% Strip trailing blank lines.
%
if strip_trail
   k = find( ~isspace(str) );
   if ~isempty(k)
      k = max(k);
      str = str(1:k);
   end
end

%
% Quick exit if string is empty.
%
if isempty( str )
   c = { '' };
   return
end

%
% Find the characters that separate the lines.
%
k = find( str == 10 );                  % Find all LF chars.
l = find( str == 13 );                  % Find all CR chars.
if isempty( k )                         % If no LF chars were found
   k = l;                               %   split at CR chars.
else                                    % Or else
   if ~isempty( l )                     %   If there are CR chars
      str(l) = [];                      %      remove them.
      k = find( str == 10 );            % Find all LF chars.
   end
end

%
% Avoid empty last string in output list when string ends in a newline.
%
if ~isempty( k ) && k(end) == length(str)
   k = [ 0 k ];                         % Add beginning.
else
   k = [ 0 k length(str)+1 ];           % Add beginning and end.
end

%
% Now split the string into lines.
%
n = length(k) - 1;                      % Number of lines.
c = cell( n, 1 );                       % Initialize output.
for i = 1:n
   c{i} = str( k(i)+1 : k(i+1)-1 );     % Extract line.
end
