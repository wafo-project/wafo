function x = stdize(x, opt)
%STDIZE   Standardize columns to have mean 0 and standard deviation 1.
%
%	  x = stdize(x)
%	  x = stdize(x, 1)
%
%	  A 1 as second argument gives normalization with N instead 
%	  of N-1 as the one argument version does.

%       LGPL Copyright (c) Anders Holtsberg, 1998

n = size(x,1);
m = mean(x);
if nargin < 2
   s = std(x);
else
   s = std(x, opt);
end
x = (x - m(ones(n,1),:)) ./ s(ones(n,1),:);
