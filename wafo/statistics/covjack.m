function [C,y] = covjack(x,fun,varargin)
%COVJACK  Jackknife estimate of the variance of a parameter estimate.
%
%	  C = covjack(X,'T')
%	  
%	  The function omputes T(X) with one observation removed at a 
%	  time  and uses the result to compute an estimate of the vari- 
%	  ance of T(X) assuming that X is a representative sample from  
%	  the underlying distribution of X. If T is multidimensional then
%	  the covariance matrix is estimated. Note that the jackknife  
%	  method does not work for some functions T that are not smooth 
%	  enough, the median being one example.

%       Anders Holtsberg, 14-12-94
%       Copyright (c) Anders Holtsberg
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




if min(size(x)) == 1
   x = x(:);
end

n = size(x,1);

xmi = x(2:n,:);
s = feval(fun,xmi,varargin{:});
y = [s(:) zeros(length(s(:)),n-1)];
for i = 2:n
   xmi = x([1:i-1,i+1:n],:);
   yy = feval(fun,xmi,varargin{:});
   y(:,i) = yy(:);
end

C = cov(y')*(n-1)^2 / n;