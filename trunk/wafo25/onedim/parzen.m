function [w, be] = parzen(n,fs)
% PARZEN returns the N-point Parzen window in a column vector.
%  
%  CALL:   [win Bw] = parzen(n,fs);
%  
%    win = parzen window
%     Bw = bandwidth of the window  
%      n = number of points
%     fs = sampling frequency (default 1)
%
% See also: bingham, hanning

% tested on: Matlab 5.3
% History:  
%   by Per Andreas Brodtkorb 14.12.1997
if nargin<2||isempty(fs),
  fs=1;
end
  
if nargout > 1
 be=1.33/(n+1)*fs; %bandwidth in Hz if fs in Hz
end  
if 1, %old oversion
  tau=2*abs((1:n)/(n+1)-0.5) ;
  tau1=tau((tau<=0.5));
  tau2=tau((find(tau>0.5)<(n/2) ));
  tau3=fliplr(tau2);
  w=[2*(1-tau2).^3 1-6*tau1.^2+6*tau1.^3  2*(1-tau3).^3]';

end

  
