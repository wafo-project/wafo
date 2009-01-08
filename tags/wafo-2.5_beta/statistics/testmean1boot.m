function [pval,cimean,cistd,Z] = testmean1boot(x,alpha,B)
%TESTMEAN1BOOT Bootstrap t-test for the mean equal to 0.
%          
%  CALL: [pval, cimean, cistd] = testmean1boot(x,alpha,B)
% 
%   pval    = pvalue, twice the one sided observed p-value
%   cimean, 
%   cisigma =  100(1-alpha)% Confidence intervals for the mean and the standard
%               deviation respectively. These are of the form 
%               [LeftLimit, PointEstimate, RightLimit]. 
%   x       = data;
%   alpha   = Confidence coefficent             (default 0.05)
%   B       = The number of bootstrap samples (default 2000)
%
% TESTMEAN1BOOT test for the mean equal to 0 and computes the confidence
% interval based on a bootstrap technique. Another name for the 
% bootstrap t is studentized bootstrap.
%
% Example
%   x = rndnorm(0,1,100,1);
%   pval = testmean1boot(x+.3) % H0: mean(x)==-0.03
%
%	See also testmean1n, testmean1r.

%       Anders Holtsberg, 27-07-95
%       Copyright (c) Anders Holtsberg

x = x(:);
if nargin<2, alpha = 0.05; end
if nargin<3, B = 2000; end

n = length(x);
m = mean(x);
s = std(x);

xB = zeros(n,B);
J = ceil(rand(n*B,1)*n);
xB(:) = x(J); 
mB = mean(xB);
sB = std(xB);
Z = (mB-m)./sB;
t = percentile(Z,[alpha/2,1-alpha/2]);
cimean = [m-t(2)*s, m, m-t(1)*s];

tt = m/s;
if tt>0
   pval = 2 * sum((mB-tt*sB)>=m)/B;
else
   pval = 2 * sum((mB-tt*sB)<=m)/B;
end

if nargout>2
   d = percentile(sB/s,[alpha/2,1-alpha/2]);
   cistd = [s/d(2), s, s/d(1)];
end

