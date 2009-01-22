function [h,hvec,score]=hlscv(A,kernel,hvec)
% HLSCV  Least Squares Cross-Validation estimate of smoothing parameter
%
% CALL: [hs,hvec,score] = hlscv(data,kernel,hvec); 
% 
%   hs     = smoothing parameter
%   hvec   = vector defining possible values of hs
%             (default linspace(0.25*h0,h0,100), h0=hos(data,kernel))
%   score  = score vector
%   data   = data vector
%   kernel = 'epanechnikov'  - Epanechnikov kernel. (default)
%            'biweight'      - Bi-weight kernel.
%            'triweight'     - Tri-weight kernel.  
%            'triangluar'    - Triangular kernel.
%            'gaussian'      - Gaussian kernel
%            'rectangular'   - Rectanguler kernel. 
%            'laplace'       - Laplace kernel.
%            'logistic'      - Logistic kernel.
%
%  Note that only the first 4 letters of the kernel name is needed.
%  Studies have shown that the theoretical and practical performance of
%  HLSCV are somewhat disappointing. In particular HLSCV is highly
%  variable and a relative slow convergence rate.
%  Example: 
%    data = rndnorm(0,1,20,1)
%    [hs hvec score] = hlscv(data,'epan');
%    plot(hvec,score) 
% See also  hste, hbcv, hboot, hos, hldpi, hscv, hstt, kde, kdefun

% Reference:  
%  B. W. Silverman (1986) 
% 'Density estimation for statistics and data analysis'  
%  Chapman and Hall pp 48--52
%
%  Wand,M.P. and Jones, M.C. (1986) 
% 'Kernel smoothing'
%  Chapman and Hall, pp 63--65

% tested on : matlab 5.2
% history:
% revised pab 20.10.1999
%   updated to matlab 5.2
% changed input arguments
% taken from kdetools     Christian C. Beardah 1995 




A=A(:);
n=length(A);

if nargin<2||isempty(kernel),
  kernel='epan';
end;

inc = 512;

if nargin<3||isempty(hvec),
  H=hos(A,kernel);
  hvec=linspace(0.25*H,H,100);
else
  hvec=abs(hvec);
end;
  
steps=length(hvec);

score = zeros(1,steps);

M=A*ones(size(A'));

Y1=(M-M');

kopt = kdeoptset('kernel',kernel,'inc',inc);
for i=1:steps,
  kopt.hs = hvec(i);
  f = kdebin(A,kopt);

  %delta=f.x{1}(2)-f.x{1}(1);
  %L1=delta*sum(f.f.^2);
  L1 = trapz(f.x{1},f.f.^2);

  Y=Y1/hvec(i);

  L2=sum(sum(mkernel(Y,kernel)))/(hvec(i)*(n-1))-n*mkernel(0,kernel)/(hvec(i)*(n-1));

  score(i)=L1-2*L2/n;

end;

[L,I]=min(score);

h=hvec(I);
