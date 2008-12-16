function ind = findextrema(x)
%FINDEXTREMA Finds indices to minima and maxima of data
%
%  CALL: ind = findextrema(x);
%
%        x  = vector with sampled values.
%
%	ind = indices to minima and maxima in the original sequence x.
%
% Example
%  t = linspace(0,7*pi,250); x = sin(t);
%  ind = findextrema(x)
%  plot(t,x,'.',t(ind),x(ind),'r.')
%
% See also findcross, crossdef


% Tested on: Matlab 5.3, 5.2 5.1

% History:
% by pab April2004  

  
ind = findcross(diff(x),0)+1;