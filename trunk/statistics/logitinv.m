function p = logitinv(z,b)
%LOGITINV Inverse logit function.
%
%	  p = logitinv(z)
%
%	  The function is equal to 1/(1+exp(-z)) if one argument is given.
%	  If lodds(X,b) is given the same thing is computed for z = X*b.
%	  If length(b) is one element too much then z = X*b(1:n-1) + b(n)
%	  is assumed.

%       Anders Holtsberg, 14-12-94
%       Copyright (c) Anders Holtsberg

if nargin>1
   X = z;
   n = length(b);
   if size(X,2) < n
      z = b(n) + X*b(1:n-1);
   else
      z = X*b;
   end
end
p = 1 ./ (1+exp(-z));
