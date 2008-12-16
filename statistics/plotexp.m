function phat = plotexp(x)
%PLOTEXP Plot data on Exponential distribution paper
%
% CALL:  phat = plotexp(X)
%
%       phat = [m] Parameter (see prbexp) estimated from 
%              the plot by least squares method
%          X = data vector or matrix
%
% Example:
%   R=rndexp(2,1,100);
%   phat=plotexp(R),shg
%
% See also  cdfexp, plotweib

x = x(:);
F=edf(x,'wdata',true');
x  = F.args;
F1 = F.data;
plot(x,-log1p(-F1),'b.','markersize',12);

m = mean(x);
hold on
plot(x,1/m*x,'r--')
hold off
title(['Exponential Probability Plot, m=' num2str(m)])
xlabel('x')
ylabel('-log(1-F)')
if nargout > 0,
  phat= m;
end
wafostamp;
