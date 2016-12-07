function [C,y] = covboot(x,fun,B,varargin)
%COVBOOT  Bootstrap estimate of the variance of a parameter estimate.
%
%   C = covboot(X,fun,B)
%
%   Computes the T(X) many times using resampled data and  
%   uses the result to compute an estimate of the variance  
%         of T(X) assuming that X is a representative sample from  
%         the underlying distribution of X. If T is multidimensional 
%         then the covariance matrix is estimated. An optional third 
%   input argument sets the number of resamples, default is 200.
%
%  Example 
%   X = [0.953428   0.170874  0.883878   0.862775 0.059614   0.967474;
%        0.348978   0.288828  0.790103   0.031843 0.620033   0.442350]';
%   m = mean(X); 
%   ci = covboot(X,'mean');
%
%  See also cov, covjack, and ciboot.

%       Anders Holtsberg, 11-01-95
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




if nargin < 3
   B = 200;
end
if min(size(x)) == 1
   x = x(:);
end

% Now, for functions that are known to produce columnwise 
% results it is faster to avoid the forloop! Note that
% there are functions that cannot be included, i e "cov".
if ischar(fun)
  funname = fun;
elseif isa(fun,'function_handle')
  funname = func2str(fun);
else
  funname = '';
end
colWiseFun = strcmp(funname,'mean') | ...
	     strcmp(funname,'std') | ...
	     strcmp(funname,'median') | ...
	     strcmp(funname,'percentile');
	     
[n,nx] = size(x);

xb = rndboot(x);
s = feval(fun,xb,varargin{:});
y = [s(:) zeros(length(s(:)),B-1)];
if nx == 1 && colWiseFun
   Bchunk = ceil(40000/n);
   i = 1;
   while i<B
      Bnext = min(B-i,Bchunk);
      xb = rndboot(x,Bnext);
      y(:,i+(1:Bnext)) = feval(fun,xb,varargin{:});
      i = i + Bnext;
   end
else
   for i = 2:B
      xb = rndboot(x);
      yy = feval(fun,xb,varargin{:});
      y(:,i) = yy(:);
   end
end

C = cov(y')*(B/(B-1));
