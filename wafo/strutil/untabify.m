function newstr = untabify( str, n )
%UNTABIFY Convert tabs to spaces.
%
%   NEWSTR = UNTABIFY( STR ) converts each tab character in STR to the
%   appropriate number of space characters.
%
%   Example
%   assert(untabify(sprintf('Hello\tworld'), 2), 'Hello  world')
%   assert(untabify(sprintf('Hello\tworld'), 0), 'Helloworld')

%   Author:      Peter J. Acklam
%   Time-stamp:  1999-03-17 15:48:11
%   E-mail:      jacklam@math.uio.no
%   WWW URL:     http://www.math.uio.no/~jacklam

%error( nargchk( 1, 2, nargin ) );
narginchk(1,2)
if ~ischar( str )
  error( 'First argument must be a string.' );
end
if nargin < 2
  n = 8;
elseif any( size(n) ~= 1 ) || ~isnumeric(n) || n ~= round(n) || n < 0;
  error( 'Second argument must be a non-negative integer.' );
end

tab = sprintf( '\t' );          % Tab character.

strlen = length(str);           % Length of input string.

newstr = '';
if strlen == 0,  return,   end

tabpos = find( str == tab );     % Position of tabs in string.
tabpos = [ tabpos strlen+1 ];    % Add "end of line mark".

ntabs = length(tabpos);          % Number of tabs.

newstr = str(1:tabpos(1)-1);
newstrlen = length(newstr);

nspc = n;
for i = 1:ntabs-1
  %nspc      = n*ceil((newstrlen+1)/n) - newstrlen;
  newstr    = [ newstr blanks(nspc) str(tabpos(i)+1:tabpos(i+1)-1) ];
  newstrlen = length(newstr);
end


