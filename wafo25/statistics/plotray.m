function bhat = plotray(x)
%PLOTRAY Plot data on a Rayleigh distribution paper
%
% CALL:  bhat = plotray(X) 
%
%   bhat = Parameter of the distribution estimated from the
%          plot by least squares method.
%   X = data
%
% Example:
%   R=rndray(1,1,100);
%   plotray(R);shg
%
% See also  cdfray, plotqq

% Reference: Cohen & Whittle, (1988) "Parameter Estimation in Reliability
% and Life Span Models", p. 181 ff, Marcel Dekker.

%tested on: matlab 5.1
% rewritten ms 15.06.2000

F=edf(x(:),'wdata',true);
x = F.args;
F1 = F.data;
plot(x,sqrt(-log1p(-F1)),'b.','markersize',12)
U=[ones(size(x)) x];
c=U\sqrt(-log1p(-F1));
b=1/c(2)/2^(1/2);
hold on
plot(x,U*c,'r--')
hold off
title('Rayleigh Probability Plot')
xlabel('X')
ylabel('(-log(1-F))^{1/2}')
if nargout > 0,
  bhat=b;
end
wafostamp;
