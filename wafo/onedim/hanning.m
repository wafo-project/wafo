function w = hanning(n,p)
% HANNING returns the N-point Hanning window in a column vector.
%
%  CALL: win = hanning(n)
%        
%    win = hanning window
%      n = number of points
%
% See also  bingham, parzen

% Tested on: matlab 5.3
% History:
%  by Per Andreas Brodtkorb%    16.09.1999
  w=bingham(n,1);