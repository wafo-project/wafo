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
% Example: 
%  % data = rndnorm(0, 1, 10, 2);
%  data = [-0.0233845632050972   0.9070186193622006;...
%           0.6529594866766634   1.3689145060433903;...
%           0.4477857310723146  -0.6311953712037597;...
%          -1.9256785038579962   0.5886257667993168;...
%          -0.5290011931824666  -0.3602090880229930];
%  assert(hmns(data,'epan'), [1.8897794596118287, 0.0399142070351313;...
%                             0.0399142070351313,   1.5551985456530106], 1e-10);
%  assert(hmns(data,'biwe'), [ 2.1856344490982642,   0.0461629876760024;...
%                              0.0461629876760024,   1.7986731199125918], 1e-10);
%  assert(hmns(data,'triw'), [ 2.4547966021766752,   0.0518479864462823;...
%                              0.0518479864462822,   2.0201807603323689], 1e-10);
%  assert(hmns(data,'gauss'), [ 0.7867937827489343,   0.0166179443738084;...
%                              0.0166179443738084,   0.6474938334398618], 1e-10);
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


