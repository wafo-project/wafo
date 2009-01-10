function [C,y] = covboot(x,fun,B,varargin)
%COVBOOT  Bootstrap estimate of the variance of a parameter estimate.
%
%	  C = covboot(X,fun,B)
%	  
%	  Computes the T(X) many times using resampled data and  
%	  uses the result to compute an estimate of the variance  
%         of T(X) assuming that X is a representative sample from  
%         the underlying distribution of X. If T is multidimensional 
%         then the covariance matrix is estimated. An optional third 
%	  input argument sets the number of resamples, default is 200.
%
%	  Example 
%    X = rand(13,2), m = mean(X), ci = covboot(X,'mean')
%
%	  See also cov, covjack, and ciboot.

%       Anders Holtsberg, 11-01-95
%       Copyright (c) Anders Holtsberg


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
