function z = rndboot(x,B)
%RNDBOOT    Simulate a bootstrap resample from a sample.
%
%	  CALL: Z = rboot(X,B)
%
%	  Give a resample of the same size as X, which is assumed
%	  to have one independent realisation in every row.
%	  RESAMPLE(X,B) gives B columns of resamples. This form works
%	  only for X one-dimensional, ie X column vector.
%
% Example
%
% See also ciboot

%       Anders Holtsberg, 14-12-94
%       Copyright (c) Anders Holtsberg

if min(size(x)) == 1
   x = x(:);
end

if nargin > 1 && size(x,2) > 1
   if B > 1, error('X multidimensional and B > 2'), end
elseif nargin < 2
   B = 1;
end

n = size(x,1);
nn = n*B;
I = ceil(n*rand(nn,1));
if B > 1
   z = zeros(n,B);
   z(:) = x(I);
else
   z = x(I,:);
end
