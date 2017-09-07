function [T,R] = strltok(S,D)
%STRLTOK Find last token in string.
%
% CALL:  [T,R] = strltok(S,D)
%
%   T  = Last token
%   R  = Remainder string
%   S  = String
%   D  = Delimiter (default white space characters)
%
%   STRLTOK returns the last token in the string S delimited
%   by D.   Trailing white space is ignored.
%
%   STRLTOK otionally returns the remainder of the original
%   string.
%   If the token is not found in S then R is an empty string and T
%   is same as S.
%
%   Example
%
%   [t, r] = strltok('hello world, hello space!');
%   assert(t, 'space!')
%   assert(r, 'hello world, hello ')    
% 
%   See also: strtok, isspace.

%error(nargchk(1,2,nargin))
narginchk(1,2)
T = []; 
R = [];

len = length(S);
if len == 0
    return
end

if (nargin <2)||isempty(D)
    D = [9:13 32]; % White space characters
end

i = len;
while (any(S(i) == D))
    i = i - 1;
    if (i < 1), return, end
end
finish = i;
while (~any(S(i) == D))
    i = i - 1;
    if (i < 1), break, end
end
start = i + 1;

T = S(start:finish);

if (nargout > 1)
    R = S(1:start-1);
end




