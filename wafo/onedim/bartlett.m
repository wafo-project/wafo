function [w, be] = bartlett(n,fs)
% Bartlett returns the N-point Bartlett window in a column vector.
%  
%  CALL:   [win Bw] = parzen(n,fs);
%  
%    win = Bartlett window
%     Bw = bandwidth of the window  
%      n = number of points
%     fs = sampling frequency (default 1)
%
% See also bingham, hanning

% tested on: Matlab 7.5.0
% History:  
%   by sflam 05.08.2008
if nargin<2||isempty(fs),
  fs=1;
end
  
if nargout > 1
 be=1.33/(n+1)*fs; %bandwidth in Hz if fs in Hz
end  

  
w = 2*(0:(n-1)/2)/(n-1);
if rem(n,2)
    % It's an odd length sequence
    w = [w w((n-1)/2:-1:1)]';
else
    % It's even
    w = [w w(n/2:-1:1)]';
end