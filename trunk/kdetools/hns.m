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
%    data = rndnorm(0, 1,20,1)
%    h = hns(data,'epan');
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
  %d=1;
end
% R= int(mkernel(x)^2)
% mu2= int(x^2*mkernel(x))
[mu2,R] = kernelstats(kernel);
AMISEconstant = (8*sqrt(pi)*R/(3*mu2^2*n))^(1/5);

iqr  = abs(diff(qlevels2(A,[75 25]),1,1)); % interquartile range
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


 

