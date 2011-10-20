function [pval, cimean, cisigma] = testmean1n(x,alpha)
%TESTMEAN1N  Test for mean equals 0 using one-sample T-test
%         
%         [pval, cimean, cisigma] = testmean1n(x,alpha)
%        
%   pval   = pvalue, the probability that the mean of x is as far from 0 as
%            it is, or further away,
%   cimean, 
%   cisigma =  100(1-alpha)% Confidence intervals for the mean and the standard
%               deviation respectively. These are of the form 
%               [LeftLimit, PointEstimate, RightLimit]. 
%   x      = data;
%   alpha  = Confidence coefficent             (default 0.05)
%
% TESTMEAN1N for mean equals 0 using one-sample t-test. Assumption data 
% comes from a normal distribution with unknown mean and unknown variance. 
% The confidence intervals for the mean and standard deviation are also 
% computed.
%
% Example
%  x = rndnorm(0,1,100,1);
%  pval = testmean1n(x+.3) % H0: mean(x)==-0.03
%
%	See also testmean1boot, testmean1r

%       Anders Holtsberg, 01-11-95
%       Copyright (c) Anders Holtsberg
%
%     This program is free software; you can redistribute it and/or modify
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
n = length(x);
m = mean(x);
s = std(x);
T = m/s*sqrt(n);
pval = (1-cdft(abs(T),n-1))*2;
t = invt(1-alpha/2,n-1);
cimean = [m-t*s/sqrt(n), m, m+t*s/sqrt(n)];
cisigma = s*[sqrt((n-1)/invchi2(1-alpha/2,n-1)), 1,... 
             sqrt((n-1)/invchi2(alpha/2,n-1))];
