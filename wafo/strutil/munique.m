function [D,nr,c] = munique(R)
% MUNIQUE  Extraction of unique rows out of a character array
%
% CALL: [D,NR,C] = munique(R);
%  
%   D = matrix containing unique rows of input matrix R (size [n1 m])
%  NR = index vector so that R = D(NR,:)                (size [n 1])
%   C = vector with number of occurences of each row of D in R (size [n1 1]) 
%   R = input character array  (size [n m])
%
% Example: 
%  words = strvcat('octave', 'matlab', 'scilab', 'abc', 'octave', 'matlab');
%  [D, NR, C] = munique(words);
%  assert(D, strvcat('octave', 'matlab', 'scilab', 'abc'))
%  assert(C, [2;2;1;1])
%  assert(D(NR,:),  words)  % i.e., D(NR,:) is exactly equal to words
%
% See also: unique

% History:
%  revised pab 16.10.2000
%  - updated header information
%  - changed name form unique to munique due to naming 
%    conflict with matlab's unique
%  Kirill K. Pankratov,  kirill@plume.mit.edu
%  8/10/94

%
%  Based on the program MUNIQUE.M by Richard Aufrichtig
%  (winner of M-file contest V3.0)

if isempty(R), D=[];nr=[];c=[];return,end
  
y = R*rand(size(R,2),1);
[y, i] = sort(y);
y = find([1; diff(y)]);

nr = zeros(size(R,1),1);
nr(y) = [i(1); diff(i(y))];
nr(i) = cumsum(nr);


if nargout>2,
  y = sort(nr);
  c = find(diff([y; length(nr)+1]));
  c = [c(1); diff(c)];
end

y = zeros(size(nr));
y(nr) = ones(size(nr));
i = find(y);
y(i) = 1:length(i);
nr = y(nr);
D(nr,:) = R;








