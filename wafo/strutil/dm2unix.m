function str = dm2unix(str,RC)
%DM2UNIX Convert DOS and MAC line endings of string to UNIX line endings.
%
%  CALL: unixstr = dm2unix(str,RC);
%
% unixstr = string with only UNIX line endings.
% str     = string with possible mixture of DOS and MAC line endings.
% RC      = Replacement Character for the CRs in the 
%           DOS (CR+LF) line endings. (default '');
%
%
% DM2UNIX recognizes DOS (CR+LF) and MAC (CR) line endings in the string
% and convert them to UNIX (LF) line endings. This is useful in
% manipulation of strings.
%
% Note: DOS (CR+LF) line endings are replaced with (LF). 
%       If RC = ' ' then length(str) == length(unixstr).
%       If RC = ''  then length(str) <= length(unixstr).
%
% Example
%   LF = char(10); CR = char(13); CRLF = [ CR, LF];
%   str = cat(2,'test', CRLF, 'test2',CR,'test3',LF);
%   assert(dm2unix(str), cat(2,'test', LF, 'test2',LF,'test3',LF));
%
% See also: findnl


%Tested on: Matlab 5.3
%History:
% revised pab 18.07.2002
% by pab 01.05.2001

%error(nargchk(1,2,nargin));
narginchk(1,2)
if nargin<2,
  RC = ''; 
end % replacement character

% New line characters LF and CR, respectively.
LF = char(10); CR = char(13);

% Find DOS line endings (CR+LF) and
% change to UNIX line endings by removing the CRs.
str = strrep(str,[CR LF],[RC LF]);

% Convert from MAC to UNIX line endings by converting CRs to LFs.
str = strrep(str,CR,LF);

return
