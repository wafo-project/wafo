function y = sinc(x)
%SINC Sin(pi*x)/(pi*x) function.
%
%  CALL:  y = sinc(x);
%
%    y = sin(pi*x)/(pi*x)    if x ~= 0
%      = 1                   if x == 0
%
% Example: 
%  x =linspace(-5,5)'; 
%  plot(x,sinc(x))
%
% See also  sin

%Tested on: Matlab 5.3
%History:
% by pab 05.04.2001

y = ones(size(x));
k = find(x~=0);
if any(k);
  xk = pi*x(k);
  y(k) = sin(xk)./(xk);
end


