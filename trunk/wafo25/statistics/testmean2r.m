function  [pval, r] = testmean2r(x,y,method)
%TESTMEAN2R Wilcoxon rank-sum test for H0: mean(x) equals mean(y).
%         
%  CALL:   [pval, ranksum] = testmean2r(x,y,method)
%         
%  pval   = pvalue, the probability that the mean of x is as far from
%           mean of y as it is, or further away,
%  ranksum = test statistic
%  x,y     = data
%  method  = 'exact'  : exact
%            'approx' : normal approximation
%
%  TESTMEAN2R test if mean of x equals  mean of y. This is known as the 
%  Wilcoxon rank-sum test or Wilcoxon-Mann-Whitney test. It assumes the 
%  data x and y are independent continiuous and symmetric and that they 
%  have the same shape and spread and differ only (possibly) in their means.
%  It is two sided. If you want a one sided alternative then divide 
%  pval by 2. The pvalue is exact and might take time to compute unless a 
%  third argument 'approx' is given, in which case a normal approximation
%  is used for the distribution of the rank-sum. 
%
% Example
% sz1 = [1,10];
% sz2 = [1 15];
% x = rndcauchy(1,24.2,sz1);
% y = rndcauchy(1,23.2,sz2);
% [pval] = testmean2r(x,y) % reject H0 if pval < alpha
%
% See also testmean2n, testmean1n, testmean1boot, testmean1r 

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


if nargin < 3
  method = 'exact';
end

x = x(:);
y = y(:);
n = length(x);
m = length(y);
[z, I] = sort([x;y]);

% --- Check for ties

J = find(diff(z)==0);
if ~isempty(J)
   if any((I(J) <= n) ~= (I(J+1) <= n))
      fprintf('\nWarning: ties in data.\n');
      fprintf('Method used is to remove all tied data.\n');
      fprintf('This is not a sophisticated method.\n');
      xx = x; yy = y;
      for i = 1:length(xx)
         j = find(yy == xx(i), 1 );
         if ~isempty(j)
            xx(i) = NaN;
            yy(j) = NaN;
         end
      end
      xx((isnan(xx))) = [];
      yy((isnan(yy))) = [];
      [pval, r] = testmean2r(xx,yy,method);
      return
   end
end

J = find(diff(z)<100000*eps*max(z));
if ~isempty(J)
   if any((I(J) <= n) ~= (I(J+1) <= n))
      fprintf('\nWarning: Nearly tie!\n');
   end
end

% --- Compute the ranks of the x's

R = 1:m+n;
r = sum(R((I<=n)));
r = r - n*(n+1)/2;
r = min(r,m*n-r);

if method(1) == 'e'

   mP = ceil((m*n+1)/2);
   P = zeros(mP,n+1);
   P(1,:) = P(1,:) + 1;
   for i=1:m
      for j = 1:n
         P(:,j+1) = (P(:,j+1)*i + [zeros(min(i,mP),1); P(1:mP-i,j)]*j)/(i+j);
      end
   end
   P = P(:,n+1);
   pval = min(2*sum(P(1:r+1)),1);

else

   % see B W Lindgren page 521
   Z = (r + 0.5 - m*n/2) / sqrt(m*n*(m+n+1)/12);
   pval = min(1, 2*cdfnorm(Z));

end
