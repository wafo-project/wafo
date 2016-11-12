function [h,hvec,score]=hbcv(A,kernel,hvec)
%HBCV  Biased Cross-Validation estimate of smoothing parameter.
%
% CALL: [hs,hvec,score] = hbcv(data,kernel,hvec); 
% 
%   hs     = smoothing parameter
%   hvec   = vector defining possible values of hs
%            (default linspace(0.25*h0,h0),100), h0=hos(data,kernel))
%   score  = score vector
%   data   = data vector
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
%  HBCV is a hybrid of crossvalidation and direct plug-in estimates.
%  The main attraction of HBCV is that it is more stable than HLSCV in
%  the sense that its asymptotic variance is considerably lower. However,
%  this reduction in variance comes at the expense of an increase in
%  bias, with HBCV tending to be larger than the HNS estimate.
%  Asymptotically HBCV has a relative slow convergence rate.
%
%  Example: 
%  % data = rndnorm(0, 1,5,2);
%  data = [-0.0233845632050972   0.9070186193622006;...
%           0.6529594866766634   1.3689145060433903;...
%           0.4477857310723146  -0.6311953712037597;...
%          -1.9256785038579962   0.5886257667993168;...
%          -0.5290011931824666  -0.3602090880229930];
%   [hs hvec score] = hbcv(data,'epan');
%   plot(hvec,score);
%   assert(hs,  1.39391122836191, 1e-10)
%   assert(hbcv(data,'biwe'), 1.65131708223985, 1e-10)
%   assert(hbcv(data,'triw'), 1.87515002002791, 1e-10)
%   assert(hbcv(data,'tria'), 1.53129587634976, 1e-10)
%   assert(hbcv(data,'gaus'), 0.629645172927051, 1e-10)
%   assert(hbcv(data,'rect'), 1.09561852654024, 1e-10)
%   assert(hbcv(data,'lapl'),  0.137620630736033, 1e-10)
%   assert(hbcv(data,'logi'),  0.0879944390318899, 1e-10)
%
%   close all;
%
% See also  hste, hboot, hns, hos, hldpi, hlscv, hscv, hstt, kde, kdefun  


% tested on : matlab 5.2
% history:
% revised pab aug2005
% -bug fix for kernels other than Gaussian
% revised pab dec2003  
% revised pab 20.10.1999
%   updated to matlab 5.2
% changed input arguments
% taken from kdetools     Christian C. Beardah 1995 

A=A(:);
n=length(A);

if nargin<3||isempty(hvec),
  H=hos(A,kernel);
  hvec=linspace(0.25*H,H,100);
else
  hvec=abs(hvec);
end;
steps=length(hvec);

M=A*ones(size(A'));

Y1=(M-M');
 

% R   = int(mkernel(x)^2)
% mu2 = int(x^2*mkernel(x))
[mu2,R] = kernelstats(kernel);

score = zeros(1,steps);
for i=1:steps,

  %sig = sqrt(2)*hvec(i);
  sig=hvec(i);

  Y=Y1/sig;

  term2=(Y.^4-6*Y.^2+3).*exp(-0.5*Y.^2)/sqrt(2*pi);

  Rf = sum(sum(term2-diag(diag(term2))))/(n^2*sig^5);


  score(i)=R/(n*hvec(i))+mu2^2*hvec(i)^4*Rf/4;

end;

[L,I]=min(score);

h=hvec(I);
