function R=mkernelrnd(kernel,N,D)
% MKERNELRND Random numbers from the Multivariate Kernel Function.
%  
% CALL:  R = mkernelrnd(kernel,N,D);
%
%   R      = random samples from kernel, size N X D
%   kernel = 'epa1'          - product of 1D Epanechnikov kernel. 
%            'gaussian'      - Gaussian kernel
%  N,D     = size of random sample
%  Note that only the first 4 letters of the kernel name is needed.
%
% See also  kde, kdefun, kdebin
  
%  Reference:  
%  B. W. Silverman (1986) 
% 'Density estimation for statistics and data analysis'  
%  Chapman and Hall, pp. 143 
%  
%  Wand, M. P. and Jones, M. C. (1995) 
% 'Density estimation for statistics and data analysis'  
%  Chapman and Hall, 

%Tested on: matlab 5.3
% History:
% revised pab  sep2005
% replaced reference to kdefft with kdebin
% revised pab Dec 2003
%  -removed   
% Revised by gl 13.07.2000
%    changed call to wnormrnd
% By pab 01.12.1999

% TODO % Add more kernels  
%            'epanechnikov'  - Epanechnikov kernel. 
%            'biweight'      - Bi-weight kernel.
%            'biw1'          - product of 1D Bi-weight kernel.
%            'triweight'     - Tri-weight kernel.  
%            'triangular'    - Triangular kernel.
%            'rectangular'   - Rectangular kernel. 
%            'laplace'       - Laplace kernel.
%            'logistic'      - Logistic kernel. 
switch lower(kernel(1:4))
case 'epa1', % 1D product kernel
   V1=-1+2*rand(N,D);V2=-1+2*rand(N,D);R=-1+2*rand(N,D);
   k=find(abs(R)>=abs(V2) & abs(R) >= abs(V1));
   if any(k)
    R(k)=V2(k);
   end  
case {'gaus', 'norm'}
   R=rndnorm(0,1,N,D);
otherwise , error('Not known or not implemented')
end
