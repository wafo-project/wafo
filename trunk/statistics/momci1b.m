function [cimean,cistd,ciskew,cikurt] = momci1b(x,alpha,B)
%MOMCI1B Moment confidence intervals using Bootstrap
%          
%  CALL: [cimean, cistd] = momci1b(x,alpha,B)
% 
%   cimean, 
%   cistd   =  100(1-alpha)% Confidence intervals for the mean, standard
%               deviation, skewness and kurtosis, respectively. These are
%               of the form  [LeftLimit, PointEstimate, RightLimit]. 
%   x       = data;
%   alpha   = Confidence coefficent             (default 0.05)
%   B       = The number of bootstrap samples (default 2000)
%
% MOMCI1B computes the confidence interval for the moments
%  based on a bootstrap technique. 
%
% Example
%   x = rndnorm(0,1,100,1);
%   [cim,cistd,cisk,ciku] = momci1b(x) 
%
%	See also testmean1boot, testmean1n, testmean1r.

% by pab nov 2007
%  Based on testmean1b by     Anders Holtsberg, 27-07-95


x = x(:);
if nargin<2, alpha = 0.05; end
if nargin<3, B = 2000; end

n = length(x);
m = mean(x);
s = std(x);

xB = zeros(n,B);
J = ceil(rand(n*B,1)*n);
xB(:) = x(J); 
clear J
mB = mean(xB);
sB = std(xB);
Z = (mB-m)./sB;
t = percentile(Z,[alpha/2,1-alpha/2]);
cimean = [m-t(2)*s, m, m-t(1)*s];

if nargout>1
   d = percentile(sB/s,[alpha/2,1-alpha/2]);
   cistd = [s/d(2), s, s/d(1)];
   if nargout>2
     sk = skew(x);
     skB = skew(xB);
     d2 = percentile(skB,[alpha/2,1-alpha/2]);
     ciskew = [d2(1),sk,d2(2)];
     ku = kurt(x);
     kuB = kurt(xB);
     d3 = percentile(kuB,[alpha/2,1-alpha/2]);
     cikurt = [d3(1),ku,d3(2)];
   end
end

