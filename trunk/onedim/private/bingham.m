function w = bingham(n,p)
% BINGHAM  returns the N-point Bingham window in a column vector.
%
%  CALL: win = bingham(n,p)
%        
%    win = bingham window
%      n = number of points
%      p = 0..1 (default 0.2) 
%         p = 0 -> rectangular window
%         p = 1 -> hanning window
%
% See also  hanning, parzen

% Tested on: matlab 5.3
% History:
%  by Per Andreas Brodtkorb%    16.12.1997
  
if (nargin < 2)
 p=0.2;
end
if (p>1 | p<0)
  disp('Error: P must be between 0 and 1!')
  return
end
n1=floor(p/2*n);
n2=n-2*n1;
w1 = .5*(1 - cos(2*pi*(1:n1)/(p*(n+1))));
w=[w1 ones(1,n2) w1(n1:-1:1) ]';
