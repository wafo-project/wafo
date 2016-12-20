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
%   x = rndnorm(0,1,1000,1);
%   [cim,cistd,cisk,ciku] = momci1b(x, 0.01, 5000);
%
%   assert(cim(1)<0 && 0<cim(3))
%   assert(cistd(1)<1 && 1<cistd(3))
%   assert(cisk(1)<0 && 0<cisk(3))
%   assert(ciku(1)<3 && 3<ciku(3))
%
% See also testmean1boot, testmean1n, testmean1r.

% by pab nov 2007
%  Based on testmean1b by     Anders Holtsberg, 27-07-95
%
% This program is free software; you can redistribute it and/or modify
%     it under the terms of the GNU Lesser General Public License as published by
%     the Free Software Foundation; either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU Lesser General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.




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

