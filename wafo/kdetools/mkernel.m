function [z,c]=mkernel(varargin)
%MKERNEL Multivariate Kernel Function.
%  
% CALL:  z = mkernel(x1,x2,...,xd,kernel);
%        z = mkernel(X,kernel);
%         
%
%   z      = kernel function values evaluated at x1,x2,...,xd  
%   x1,x2..= input arguments, vectors or matrices with common size
% or 
%   X      = cellarray of vector/matrices with common size
%            (i.e. X{1}=x1, X{2}=x2....)
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
% Revised pab sep2005
% -replaced reference to kdefft with kdebin
% revised pab aug2005
% -Fixed some bugs
% revised pab Dec2003
% removed some old code  
% revised pab 27.04.2001
% - removed some old calls
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

d=length(varargin)-1;
kstr=varargin{d+1}; % kernel string
if iscell(varargin{1})
  X=varargin{1};
  d=numel(X);
else
  X=varargin;
end

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
    b=1;% radius of the kernel
    b2=b^2;
    s=X{1}.^2;
    k=find(s<=b2);
    z=zeros(size(s));
    ix=2;
    while (any(k) && (ix<=d)), 
      s(k)=s(k)+X{ix}(k).^2;
      k1=(s(k)<=b2);
      k=k(k1);
      ix=ix+1; 
    end;
    if any(k)
      c=2^r*prod(1:r)*vsph(d,b)/prod((d+2):2:(d+2*r)); % normalizing constant
      %c=beta(r+1,r+1)*vsph(d,b)*(2^(2*r)); % Wand and Jones pp 31
      % the commented c above does note yield the right scaling
      % for d>1   
      z(k)=((1-s(k)/b2).^r)/c;
    end
    
  case 'rect', % 1D product Rectangular Kernel
    z=zeros(size(X{1}));
    k=find(abs(X{1})<=1);
    ix=2;
    while (any(k) && (ix<=d)),
      k1 =(abs(X{ix}(k))<=1);
      k=k(k1);
      ix=ix+1; 
    end
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
    b=1;
    b2=b^2;
    b21=1/b2;
    z=zeros(size(X{1}));
    k=find(abs(X{1})<=b);
    ix=2; 
    while (any(k) && (ix<=d)),
      %for ix=2:d
      k1 =(abs(X{ix}(k))<=b);
      k  = k(k1);
      ix=ix+1; 
    end
    if any(k)
      c=2^r*prod(1:r)*vsph(1,b)/prod((1+2):2:(1+2*r)); % normalizing constant
      z(k) = (1-X{1}(k).^2*b21).^r;
      for ix=2:d
        z(k)=z(k).*(1-X{ix}(k).^2*b21).^r;
      end;   
      z(k)=z(k)/c^d;
    end
  case 'tria',% 1D product Triangular Kernel 
    z=zeros(size(X{1}));
    k=find(abs(X{1})<1);
    ix=2;
    while (any(k) && (ix<=d)),
      %for ix=2:d
      k1 =(abs(X{ix}(k))<1);
      k  = k(k1); 
      ix=ix+1;
    end
    if any(k)  
      z(k) = (1-abs(X{1}(k)));
      for ix=2:d
        z(k)=z(k).*(1-abs(X{ix}(k)));
      end
    end
   case {'norm','gaus'},% multivariate gaussian  Density Function.
     s=X{1}.^2;
     for ix=2:d  
       s=s+X{ix}.^2;
     end;
     z=(2*pi)^(-d/2)*exp(-0.5*s);   
   case 'lapl' % Laplace Kernel 
     z=0.5*exp(-abs(X{1}));
     for ix=2:d
       z=z.*0.5*exp(-abs(X{ix}));
     end
   case 'logi', % Logistic Kernel 
     z1=exp(X{1});
     z=z1./(z1+1).^2;
     for ix=2:d
       z1=exp(X{ix});
       z=z.*z1./(z1+1).^2;
     end
    
   otherwise, error('unknown kernel')
 end


return