function h=hns(A,kernel)
%HNS Normal Scale Estimate of Smoothing Parameter.
%
% CALL:  h = hns(data,kernel)
%
%   h      = one dimensional optimal value for smoothing parameter
%            given the data and kernel.  size 1 x D
%   data   = data matrix, size N x D (D = # dimensions )
%   kernel = 'epanechnikov'  - Epanechnikov kernel.
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
% HNS only gives an optimal value with respect to mean integrated 
% square error, when the true underlying distribution 
% is Gaussian. This works reasonably well if the data resembles a
% Gaussian distribution. However if the distribution is asymmetric,
% multimodal or have long tails then HNS may  return a to large
% smoothing parameter, i.e., the KDE may be oversmoothed and mask
% important features of the data. (=> large bias).
% One way to remedy this is to reduce H by multiplying with a constant 
% factor, e.g., 0.85. Another is to try different values for H and make a 
% visual check by eye.
%
%  Example: 
%    % data = rndnorm(0, 1,20,1)
%  data = [-0.0233845632050972   0.9070186193622006;...
%           0.6529594866766634   1.3689145060433903;...
%           0.4477857310723146  -0.6311953712037597;...
%          -1.9256785038579962   0.5886257667993168;...
%          -0.5290011931824666  -0.3602090880229930];
%  assert(hns(data,'epan'), [1.73513679136905, 1.43948322577017], 1e-10);
%  assert(hns(data,'biwe'), [2.05555487703312, 1.70530460760076], 1e-10);
%  assert(hns(data,'triw'), [2.33418149081877, 1.93645545333964], 1e-10);
%  assert(hns(data,'tria'), [1.90615281623682, 1.58135947458212], 1e-10);
%  assert(hns(data,'gaus'), [0.783780547013426, 0.650230549961770], 1e-10);
%  assert(hns(data,'rect'), [1.36382287194830, 1.13143825711994], 1e-10);
%  assert(hns(data,'lapl'), [0.579817701798895, 0.481021358025369], 1e-10);
%  assert(hns(data,'logi'), [0.438140924596874, 0.363485181466877], 1e-10);
%  
% See also  hste, hbcv, hboot, hos, hldpi, hlscv, hscv, hstt, kde

% Reference:  
%  B. W. Silverman (1986) 
% 'Density estimation for statistics and data analysis'  
%  Chapman and Hall, pp 43-48 

%  Wand,M.P. and Jones, M.C. (1995) 
% 'Kernel smoothing'
%  Chapman and Hall, pp 60--63


%Tested on: matlab 5.3
% History:
% revised pab dec2003
%  fixed a bug when D>1  
% revised pab 21.09.99
% 
% updated string comparisons
% from kdetools

if nargin<2||isempty(kernel)
 kernel='gauss';
end
[n, d] = size(A);
if (n==1) && (d>1),
  A=A.';
  n=d;
  d=1;
end
% R= int(mkernel(x)^2)
% mu2= int(x^2*mkernel(x))
[mu2,R] = kernelstats(kernel);
AMISEconstant = (8*sqrt(pi)*R/(3*mu2^2*n))^(1/5);

iqr  = abs(diff(qlevels2(A,[75 25]),1,1+(d==1))); % interquartile range
stdA = std(A);
h    = stdA*AMISEconstant;
k = find(iqr>0);
if any(k)
  % use of interquartile range guards against outliers.
  % the use of interquartile range is better if 
  % the distribution is skew or have heavy tails
  % This lessen the chance of oversmoothing.
  h(k)=min(stdA(k),iqr(k)/1.349)*AMISEconstant;
end

return


 

