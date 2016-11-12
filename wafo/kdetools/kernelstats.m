function [mu2, R, Rdd] = kernelstats(kernel)
%KERNELSTATS Return 2'nd order moment of kernel pdf
%            as well as the integral of the squared kernel
%            and integral of squared double derivative of kernel.
%  
%  CALL:  [mu2, R, Rdd] = kernelstats(kernel)
%  
% mu2    = 2'nd order moment, i.e.,int(x^2*kernel(x))
% R      = integral of squared kernel, i.e., int(kernel(x)^2)
% Rdd    = integral of squared double derivative of kernel, i.e., 
%          int( (kernel''(x))^2 ).
% kernel = string identifying the kernel, i.e., one of:
%            'epanechnikov'  - Epanechnikov kernel.
%            'biweight'      - Bi-weight kernel.
%            'triweight'     - Tri-weight kernel.  
%            'triangluar'    - Triangular kernel.
%            'gaussian'      - Gaussian kernel
%            'rectangular'   - Rectanguler kernel. 
%            'laplace'       - Laplace kernel.
%            'logistic'      - Logistic kernel.  
%
%  Note that only the first 4 letters of the kernel name is needed.
%
% Example
%  [mu2,R]=kernelstats('triweight');
%
%  assert(mu2, 0.111111111111111, eps)
%  assert(R, 0.815850815850816, eps)
%
% See also  mkernel
  
% Reference 
%  Wand,M.P. and Jones, M.C. (1995) 
% 'Kernel smoothing'
%  Chapman and Hall, pp 176.
  
%History
% by pab Dec2003
  
switch lower(kernel(1:4))
 case 'biwe', % Bi-weight kernel
  mu2 = 1/7;
  R   = 5/7;
  Rdd = 45/2;
 case {'epan' 'epa1'}, % Epanechnikov kernel
  mu2 = 1/5;
  R   = 3/5;
  Rdd = inf;
 case {'gaus','norm'}, % Gaussian kernel
  mu2 = 1;
  R   = 1/(2*sqrt(pi));
  Rdd = 3/(8*sqrt(pi));
 case  'lapl', % Laplace
  mu2 = 2; 
  R   = 1/4;
  Rdd = inf;
 case 'logi', % Logistic
  mu2 = pi^2/3;
  R=1/6;
  Rdd = 1/42;
 case {'rect','unif'}, % Rectangular
  mu2 = 1/3;
  R   = 1/2;
  Rdd = inf;
 case 'tria', % Triangular
  mu2 = 1/6;
  R   = 2/3;
  Rdd = inf;
 case 'triw', % Triweight
  mu2 = 1/9;
  R   = 350/429;
  Rdd = inf;
 otherwise
  error('Unknown kernel.')
end;