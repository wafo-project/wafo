function [h,hvec,score]=hboot(A,kernel,hvec)
%HBOOT  Bootstrap cross-validation estimate of smoothing parameter.
%
% CALL: [hs,hvec,score] = hboot(data,kernel,hvec); 
% 
%   hs     = smoothing parameter
%   hvec   = vector defining possible values of hs
%            (default linspace(0.25*h0,h0,100), h0=hos(data,kernel))
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
%  Example:
%  % data = rndnorm(0, 1,5,2);
%  data = [-0.0233845632050972   0.9070186193622006;...
%           0.6529594866766634   1.3689145060433903;...
%           0.4477857310723146  -0.6311953712037597;...
%          -1.9256785038579962   0.5886257667993168;...
%          -0.5290011931824666  -0.3602090880229930];
%   [hs hvec score] = hboot(data,'epan');
%   plot(hvec, score);
%   assert(hs,  1.39391122836191, 1e-10)
%   assert(hboot(data,'biwe'), 1.65131708223985, 1e-10)
%   assert(hboot(data,'triw'), 1.87515002002791, 1e-10)
%   assert(hboot(data,'tria'), 1.53129587634976, 1e-10)
%   assert(hboot(data,'gaus'), 0.629645172927051, 1e-10)
%   assert(hboot(data,'rect'), 1.09561852654024, 1e-10)
%   assert(hboot(data,'lapl'), 0.465792904029649, 1e-10)
%   assert(hboot(data,'logi'), 0.351977756127559, 1e-10)  
%
%   close all;
% 
% See also  hste, hbcv, hos, hldpi, hlscv, hscv, hstt, kde, kdefun 

% tested on : matlab 5.2
% history:
% revised pab dec2003
%  -fixed a bug in default value for hvec  
% revised pab 20.10.1999
%   updated to matlab 5.2
% changed input arguments
% taken from kdetools     Christian C. Beardah 1995 

A=A(:);
n=length(A);

if nargin<3||isempty(hvec),
  H    = hos(A,kernel);
  hvec = linspace(0.25*H,H,100);
else
  hvec=abs(hvec);
end;
  
steps=length(hvec);


M=A*ones(size(A'));

Y1=(M-M');
 
[mu2,R] = kernelstats(kernel);
STEconstant =  R /(mu2^(2)*n);

kernel2 = 'gauss';
[mu2,R] = kernelstats(kernel2);
STEconstant2 =  R /(mu2^(2)*n);

hvec = hvec*(STEconstant2/STEconstant)^(1/5);

score = zeros(1,steps);

for i=1:steps,

  Y=-Y1.^2/(4*hvec(i)^2);

  T1=exp(Y/2);

  T1=T1-diag(diag(T1));

  T2=-4*exp(Y/1.5)/sqrt(3);

  T2=T2-diag(diag(T2));

  T3=sqrt(2)*exp(Y);

  T3=T3-diag(diag(T3));

  T=T1+T2+T3;

  T=sum(sum(T))+n*sqrt(2);

  score(i)=T/(2*n^2*hvec(i)*sqrt(2*pi));

end;

[L,I]=min(score);

hvec = hvec*(STEconstant/STEconstant2)^(1/5);


h=hvec(I);
