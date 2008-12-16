function x = center(x)
%CENTER Recenter columns to have mean 0.
%
%	  x = center(x)
%
% See also stdize

%       GPL Copyright (c) Anders Holtsberg, 1998

n = size(x,1);
m = mean(x);

x = (x - m(ones(n,1),:));
