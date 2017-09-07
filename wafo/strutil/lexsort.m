function [A,indx] = lexsort(Z)
% LEXSORT Lexical sorting of a character array 
% 
% CALL: [A,I] = lexsort(Z); 
%
%  Z = character array to be sorted in lexical order
%  A = finally sorted character array
%  I = index so that A = Z(I,:);
%
% Example
%  words = strvcat('octave', 'matlab', 'scilab', 'abc', 'Octave', 'Matlab');
%  assert(lexsort(words), strvcat('abc','Matlab','matlab','Octave','octave','scilab')
%
% See also sort, sortrows

% History
% Bryce Gardner, gardner@ecn.purdue.edu, 317-494-0231
% Revised Per A. Brodtkorb 
%  - fixed a bug: i== [] is technically incorrect, 
%    replaced isempty(i) instead (L20).
%  - added I to output

A = Z;                          % A is working matrix

if isempty(A), return, end

[nr,nc] = size(Z);

A     = abs(A)-64;              % make A-Z become 1-26
i     = find(A==(32-64));       % zero trailing spaces
A(i)  = zeros(size(i));
i     = find(A>26);             % a-z
i2    = find(A<27&A>0);         % A-Z
A(i)  = (A(i)-32)*2;
A(i2) = A(i2)*2-1;              % precedence order is AaBbCcDd...
 
val   = A*52.^((nc-1:-1:0).');   % convert each word to base 52 number
 
for level=1:nc,                 % remove leading spaces
  i = find(val<52^(nc-1));
  if isempty(i), break, end
  val(i) = val(i)*52;
end
 
[val,indx] = sort(val);         % sort base 52 encoded words
 
A = Z(indx,:);                  % reorder Z





