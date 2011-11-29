function h=hmns(A,kernel)
%HMNS Multivariate Normal Scale Estimate of Smoothing Parameter.
%
% CALL:  h = hmns(data,kernel)
%
%   h      = M dimensional optimal value for smoothing parameter
%            given the data and kernel.  size D x D
%   data   = data matrix, size N x D (D = # dimensions )
%   kernel = 'epanechnikov'  - Epanechnikov kernel.
%            'biweight'      - Bi-weight kernel.
%            'triweight'     - Tri-weight kernel.  
%            'gaussian'      - Gaussian kernel
%  
%  Note that only the first 4 letters of the kernel name is needed.
% 
% HMNS  only gives  a optimal value with respect to mean integrated 
% square error, when the true underlying distribution  is
% Multivariate Gaussian. This works reasonably well if the data resembles a
% Multivariate Gaussian distribution. However if the distribution is 
% asymmetric, multimodal or have long tails then HNS is maybe more 
% appropriate.
%
%  Example: 
%    data = rndnorm(0, 1,20,2)
%    h = hmns(data,'epan');
%  
% See also  hns, hste, hbcv, hboot, hos, hldpi, hlscv, hscv, hstt

% Reference:  
%  B. W. Silverman (1986) 
% 'Density estimation for statistics and data analysis'  
%  Chapman and Hall, pp 43-48, 87 
%
%  Wand,M.P. and Jones, M.C. (1995) 
% 'Kernel smoothing'
%  Chapman and Hall, pp 60--63, 86--88


%Tested on: matlab 5.3
% History:
% revised pab dec2003  
% by pab 06.12.99
% 

% TODO % implement more kernels  
  

[n d]=size(A);
if (n==1) && (d>1),
  A=A';
  n=d;
  d=1;
end
if d==1,
  h=hns(A,kernel);
  return
end

switch lower(kernel(1:4))
case {'epan' }, % Epanechnikov kernel
 a=(8*(d+4)*(2*sqrt(pi))^d/vsph(d))^(1/(4+d));
case 'biwe', % Bi-weight kernel
  a=2.7779;
  if d>2
   error('not implemented for d>2')
  end 
case 'triw', % Triweight
  a=3.12;
  if d>2
    error('not implemented for d>2')
  end 
case 'gaus', % Gaussian kernel
 a=(4/(d+2))^(1/(d+4));

 otherwise
  error('Unknown kernel.')
end;
covA = cov(A);

h=a*sqrtm(covA)*n^(-1/(d+4));

return


