function [pval, cimean, cisigma] = testmean2n(x,y,alpha,vartype)
%TESTMEAN2N  Two-sample t-test for mean(x) equals mean(y)
%         
%         [pval, cidiffmean, cisigma] = testmean2n(x,y,alpha,vartype)
%         
%   pval    = pvalue, the probability that the mean of x is as far from the 
%             mean of y as it is, or further away,
% cidiffmean, 
%   cisigma = 100(1-alpha)% Confidence intervals for the difference in mean
%             and the common standard deviation respectively. These are of 
%             the form:   [LeftLimit, PointEstimate, RightLimit]. 
%   x,y     = data;
%   alpha   = Confidence coefficent             (default 0.05)
%  vartype  = 'equal'   : If variances are assumed equal (default)
%             'unequal' : if variances are assumed unequal
% 
% TESTMEAN2N test for the null hypothesis that data in the vectors x and y
% have equal means, against the alternative that the means are not equal. 
% The assumption is that x and y are independent random samples from normal
% distributions with equal but unknown variances if VARTYPE is 'equal'. 
% The confidence interval for the difference in mean and the common 
% standarddeviation are also computed.
%
% When VARTYPE == 'unequal' the test is performed with Welch correction and
% is useful when the data are normal, sample sizes are small, and the 
% variances are heterogeneous. Otherwise, use the parametric t-test for
% normal data, or the permutational t-test for skewed data. For 
% heteroscedastic data that cannot be normalized, a nonparametric
% test should be used. For more detailed recommendations type testmean2n.
% 
% Example
% sz1 = [1,10];
% sz2 = [1 15];
% x = rndnorm(24.2,20,sz1);
% y = rndnorm(23.2,10,sz2);
% [pval, cimean, cisigma] = testmean2n(x,y) 
% [pval2, cimean2, cisigma2] = testmean2n(x,y,[],'unequal') 
%
% type testmean2n % Detailed recommendations  
%
%	See also testmean2r.

% Reference
% Legendre, P. and D. Borcard. Statistical comparison of univariate tests
% of homogeneity of variances.
%
% Welch, B. L. (1947), "The generalization of "student's" problem when
% several different population variances are involved", Biometrika 34: 28-35

%       Anders Holtsberg, 18-11-93
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
y = y(:);
if nargin<3|| isempty(alpha), alpha = 0.05; end
if nargin<4 
  vartype = 'equal';
end
nx = length(x);
ny = length(y);
mx = mean(x);
my = mean(y);
m = mx - my;
isequalvariance =  strncmpi(vartype,'equal',1);
if isequalvariance
  df = nx + ny - 2;
  s = sqrt(((nx-1)*var(x)+(ny-1)*var(y))/df);
  d = s*sqrt(1/nx+1/ny);
  cisigma = s*[sqrt(df/invchi2(1-alpha/2,df)), 1,... 
    sqrt(df/invchi2(alpha/2,df))];
else
  vx = var(x)/nx;
  vy = var(y)/ny;
  d2 =(vx+vy);
  d = sqrt(d2);
  df =  max(round(d2.^2./(vx^2/(nx-1)+vy^2/(ny-1))),1);
  cisigma = [nan,nan,nan];
end
T = m/d;
pval = (1-cdft(abs(T),df))*2;
t = invt(1-alpha/2,df);
cimean = [m-t*d, m, m+t*d];

% This excerp from 
% Legendre, P. and D. Borcard. Statistical comparison of univariate tests
% of homogeneity of variances.
%
% Table II Recommendations for t-test of equality of two means.
% _________________________________________________________________________
% 1. Sample sizes equal?
%  Yes -> go to 2
%  No -> go to 6
% 2. Equal sample sizes. Distribution:
%  Normal -> go to 3
%  Skewed -> go to 5
% 3. Normal distributions. THV result:
%  Variances homogeneous -> use any one of the 3 tests (most simple: parametric t-test).
%  Variances unequal -> go to 4
% 4. Variances unequal. Sample size:
%  Small -> use the t-test with Welch correction.
%  Large -> use any one of the 3 tests (most simple: parametric t-test).
% 5. Skewed distributions. THV result:
%  Variances homogeneous -> all 3 tests are valid, but the permutational t-test is preferable
%  because it has correct type I error and the highest power.
%  Variances unequal -> normalize the data or use a nonparametric test (Wilcoxon-Mann-
%  Whitney test, median test, Kolmogorov-Smirnov two-sample test, etc.).
% 6. Unequal sample sizes. Distribution:
%  Normal -> go to 7
%  Skewed -> go to 8
% 7. Normal distributions. THV result:
%  Variances homogeneous -> use the parametric or permutational t-tests (most simple:
%  parametric t-test).
%  Variances unequal -> use any one of the 3 tests (most simple: parametric t-test). Power is
%  low when the sample sizes are strongly unequal; avoid the Welch-corrected t-test in the
%  most extreme cases of sample size inequality (lower power).
% 8. Skewed distributions. THV result
%  Variances homogeneous -> use the permutational t-test.
%  Variances unequal -> normalize the data or use a nonparametric test (Wilcoxon-Mann-
%  Whitney test, median test, Kolmogorov-Smirnov two-sample test, etc.).
% ________________________________________________________________________
% 
