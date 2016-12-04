function  [inl,linenum] = findnl(str)
%FINDNL Find new line characters and line numbers
%
% CALL: [inl,linenum] = findnl(string);
%
%  inl     = indices to new line characters in string
%  linenum = line number mask, (i.e., vector of line numbers)
%  string  = character vector
%
% Both DOS (CR+LF), UNIX (LF) and MAC (CR) line endings are recognized.  
%
% Example:   % extract 3 first lines of a file:
%  t = freadtxt('findnl.m');
%  [inl,linenum] = findnl(t);
%  assert(dewhite(t(1:inl(1))), 'function  [inl,linenum] = findnl(str)')
%  assert(dewhite(t(linenum==2)), '%FINDNL Find new line characters and line numbers')
%
% See also: find, cumsum


% History:
% revised pab 18.07.2002
% - replace code with a call to dm2unix
% revised pab 15.05.2001
% - added example
% - New line characters now belong to the same line.
% revised pab 13.12.2000
% - Now also handles the UNIX and MAC endings correctly 


% New line characters LF and CR, respectively.
LF = char(10); 
%CR = char(13);
space = ' ';

str = dm2unix(str,space);

% Find the line endings.
inl = find(str == LF);

if nargout>1
  % Line numbering mask
  linenum      = zeros(size(str));
  linenum(inl) = ones(size(inl)); % originally a - sign here
  linenum(1)   = linenum(1)+1;
  linenum      = cumsum(linenum);
  linenum(inl) = linenum(inl)-1;  % Associate new line characters to the
                                  % previous line number
end
return






