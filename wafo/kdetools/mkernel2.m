function [z,c]=mkernel2(X,kstr)
%MKERNEL2 Multivariate Kernel Function, alternative version.
%  
% CALL:  z = mkernel2(X,kernel);
%         
%   z   = kernel function values evaluated at X  
%   X   = matrix  size N x D (D = # dimensions)
%
%   kernel = 'epanechnikov'  - Epanechnikov kernel. 
%            'epa1'          - product of 1D Epanechnikov kernel. 
%            'biweight'      - Bi-weight kernel.
%            'biw1'          - product of 1D Bi-weight kernel.
%            'triweight'     - Tri-weight kernel.  
%            'triangular'    - Triangular kernel.
%            'gaussian'      - Gaussian kernel
%            'rectangular'   - Rectangular kernel. 
%            'laplace'       - Laplace kernel.
%            'logistic'      - Logistic kernel.  
%  
%  Note that only the first 4 letters of the kernel name is needed.
%
% See also  kde, kdefun, kdebin
  
%  Reference:  
%  B. W. Silverman (1986) 
% 'Density estimation for statistics and data analysis'  
%  Chapman and Hall, pp. 43, 76 
%  
%  Wand, M. P. and Jones, M. C. (1995) 
% 'Density estimation for statistics and data analysis'  
%  Chapman and Hall, pp 31, 103,  175  

%Tested on: matlab 5.3
% History:
% Revised pab Dec2003
% removed some code  
% revised pab 27.04.2001
% - renamed from mkernel to mkernel2
% - removed some old calls
% - improved speed for all kernels
%  revised pab 01.01.2001
%  - speeded up tri3
%  revised pab 01.12.1999
%   - added four weight, sphere
%   - made comparison smarter => faster execution for d>1
%  revised pab 26.10.1999
%   fixed normalization fault in epan
% by pab 21.09.99  
%  added multivariate epan, biweight and triweight
%
% collected all knorm,kepan ... into this file 
% adapted from kdetools CB 

%error(nargchk(2,2,nargin))
narginchk(2,2)
[n,d]=size(X);
               % n=number of evaluation points,
               % d=dimension of the data.  
z = zeros(n,1);
switch lower(kstr(1:4)) 
  case {'sphe','epan','biwe','triw','four'}
    switch lower(kstr(1:4))
      case 'sphe', r=0;  %Sphere = rect for 1D
      case 'epan', r=1;  %Multivariate Epanechnikov kernel.
      case 'biwe', r=2;  %Multivariate Bi-weight Kernel
      case 'triw', r=3;  %Multi variate Tri-weight Kernel 
      case 'four', r=4;  %Multi variate Four-weight Kernel
        % as r -> infty, b -> infty => kernel -> Gaussian distribution
    end 
    b  = 1;    % radius of the kernel
    b2 = b^2;  % radius squared
    s  = sum(X.^2,2);
    k  = find(s<=b2);
    
    if any(k)
      c=2^r*prod(1:r)*vsph(d,b)/prod((d+2):2:(d+2*r)); % normalizing constant
      %c=beta(r+1,r+1)*vsph(d,b)*(2^(2*r)); % Wand and Jones pp 31
      % the commented c above does note yield the right scaling
      % for d>1   
      z(k)=((1-s(k)/b2).^r)/c;
    end
    
  case 'rect', % 1D product Rectangular Kernel
    k=find(all(abs(X)<=1,2));
    if any(k) 
      z(k)=(0.5^d);
    end
  case {'epa1','biw1','triw1','fou1'}
    switch lower(kstr(1:4))
      %case 'rect', r=0;  %rectangular
      case 'epa1', r=1;  %1D product Epanechnikov kernel.
      case 'biw1', r=2;  %1D product Bi-weight Kernel
      case 'tri1', r=3;  %1D product Tri-weight Kernel 
      case 'fou1', r=4;  %1D product Four-weight Kernel
    end 
    b   = 1;
    b2  = b^2;
    b21 = 1/b2;
    k   = find(all(abs(X)<=b,2));
    if any(k)
      c=2^r*prod(1:r)*vsph(1,b)/prod((1+2):2:(1+2*r)); % normalizing constant
      z(k) = prod(1-X(k,:).^2*b21,2)/c^d;
    end
  case 'tria',% 1D product Triangular Kernel 
    k  = find(all(abs(X)<=1,2));
    if any(k)
      z(k) = prod(1-abs(X(k,:)),2);
    end    
   case {'norm','gaus'},% multivariate gaussian  Density Function.
     z = (2*pi)^(-d/2)*exp(-0.5*sum(X.^2,2));   
   case 'lapl' % Laplace Kernel 
     z=0.5^d*exp(-sum(abs(X),2));
   case 'logi', % Logistic Kernel 
     s = exp(X);
     z = prod(s./(s+1).^2,2);    
   otherwise, error('unknown kernel')
 end

return








