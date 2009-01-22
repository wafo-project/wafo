function [area,epsi,a,b] = simpson(x,f,dim)
% SIMPSON Numerical integration with the Simpson method
%
% CALL: [area,epsi,a,b] = simpson(x,f,dim)
%       [area,epsi,a,b] = simpson(f,dim)
%
%    area = result, approximative integral
%    epsi = estimate of the absolute error
%    a,b  = integration limits
%      x  = argument of function (vector or matrix, see below)
%      f  = function (vector or matrix, see below)
%     dim = dimension to integrate 
%
% Evaluates an approximation of the integral of the
% vector function F(x) wrt X from A to B using the Simpson method. 
%
% X and F must be vectors of the same length,
% OR X must be a vector and F an array whose first non-singleton
% dimension is length(X), then SIMPSON operates along this dimension.
% OR, if X is an array with the same dimension as F the integration
% is performed column-wise. If DIM is given SIMPSON integrates across
% dimension DIM of F. The length of X must be the same as size(F,DIM)).
%
% Notice that if N = length(F) is an even number,
% SIMPSON is used on the N-1 first values and the Trapezoidal rule 
% on the two last values. If X is not equidistant
% spaced use TRAPZ instead.
%
% Example:%
%          x = linspace(0,4,201);
%          f = exp(-x.^2);
%          [area,epsi] = simpson(x,f)
%
% See also  trapz

% Tested on: Matlab 5.3
% History:
% revised pab feb2007
% - fixed a bug when x is vector and F is matrix.
% revised by pab 14.05.2000
%  - added dim
%  - added trapezoidal rule when lengt(f) is even
%  - changed the documentation according to the new changes.
% revised by es 28.09.1999 (documentation, help)
% by pab 1997, updated 16.06.1998


% if 1 % NEW CALL
 
  perm = []; nshifts = 0;
  if nargin == 3,                 % simpson(x,f,dim)
    perm = [dim:max(ndims(f),dim) 1:dim-1];  
    f = permute(f,perm);
  elseif nargin==2 && length(f)==1 % simpson(f,dim)
    dim = f; f = x; x=[];
    perm = [dim:max(ndims(f),dim) 1:dim-1];
    f = permute(f,perm);
  else                            % simpson(x,f) or simpson(f)
    if nargin < 2, f = x; x=[]; end
    [f,nshifts] = shiftdim(f);
  end
  
  fsiz = size(f);
  m = fsiz(1);
  if m> 5 ,
    n=3;
  else
    n=1;
  end
  if isempty(x) % assume unit spacing
    a=1;b=m;bmn=m-n;
  else
    if isempty(perm),
      x = shiftdim(x);
      xsiz=size(x);
    else
      xsiz = size(x);
      if prod(xsiz)==prod(fsiz)
        x = permute(x,perm);
      end
    end
    if prod(xsiz) == m
      %make sure x is a column vector
      x=x(:);
      %xsiz = [m,1];
    elseif xsiz(1) ~= m
      error('length(x) must equal length of first non-singleton dim of f.');
    end
    a =x(1,:);b=x(m,:);
    bmn=x(m-n,:); 
  end
    
  if m<3, 
    error('The vector must have more than 3 elements!')
  end
  
  if (mod(m, 2) == 0), % checking if m is even
    
    if n==1 % Use trapezoidal  from bmn to b
      area = (f(m,:)+f(m-1,:))/2.*(b-bmn); 
    else    % Use 4-points Newton-Cotes rule
      area = (f(m,:)+3*(f(m-1,:)+f(m-2,:))+ f(m-3,:))/8.*(b-bmn); 
    end
    m = m-n;     
    hn=2*(bmn-a)/(m-1);
  else
    area = 0;
    hn=2*(b-a)/(m-1);
  end
  
  % Midpoint approximation
  Un = hn.*sum(f(2:2:(m-1),:)); 
  
  % Trapeziodal approximation
  Tn = hn/2.*(f(1,:) + 2*sum(f(3:2:(m-2),:)) + f(m,:)); 
  
  % Simpson's approximation
  area = area+(Tn+2*Un)/3; 
  
  % Asymptotically the simpson approximation is a better estimate 
  % than both Tn and Un. If this is the case, 
  % then an estimate of the absolute error is
  if nargout>1,  epsi = abs(Tn-Un)/4; end
  
  fsiz(1) = 1;
  area = reshape(area,[ones(1,nshifts),fsiz]);
  if ~isempty(perm), area = ipermute(area,perm); end


% else % Old call
%   if nargin<2,
%     f = x;
%     [f,nshift] = shiftdim(f);
%     [nx m]=size(f);
%     fsiz(1)=1;fsiz(2)=1;
%     if nx<m
%       tmp=m;m=nx;nx=tmp; 
%       f=f.';
%       % fsiz(1)=m;  fsiz(2)=1;
%     else
%       %fsiz(1)=1;fsiz(2)=m;
%     end
%     switch m
%       case 1,  a=1;b=nx; bm1=nx-1;
%       case 2,  a=f(1,1);b=f(nx,1);bm1=f(nx-1,1); f=f(:,2);
%       otherwise, error('Wrong dimension of input! dim must be 2xN, 1xN, Nx2 or Nx1 ') 
%     end
%     m=1;
%   end  
%   
%   if nargin==2,
%     fsiz=size(f);
%     [n m]=size(f);
%     [nx mx]=size(x);
%     
%     %make sure x is a column vector
%     if nx*mx==n || nx*mx==m  , x=x(:);nx=max(nx,mx);end 
%     
%     a=x(1,:);b=x(nx,:);bm1=x(nx-1,:); 
%     if nx~=n
%       tmp=m;m=n;n=tmp; 
%       f=f.';
%       fsiz(2)=1;
%       if (nx~=n) || (ndims(f)>2)
%         error('Fx must have dimensions which equals the  size of x ')
%       end
%     else
%       fsiz(1)=1;
%     end
%   end  
%   
% 
%   if nx<3, 
%     error('The vector must have more than 3 elements!')
%   end
%   
%   if (mod(nx, 2) == 0), % checking if n is even
%     f(nx,:)=[];
%     nx=nx-1; 
%     
%     disp('Warning: Changed the integration limits ')
%     disp(['         from ' num2str(b) ' to ' num2str(bm1) ])
%     disp('         so that the length of f(x) is odd.')
%     
%     b=bm1;
%   end
%   
%   hn=2*(b-a)/(nx-1);
%   
%   
%   % Midpoint approximation
%   Un=hn.*sum(f(2:2:(nx-1),1:m)); 
%   
%   % Trapeziodal approximation
%   Tn=hn/2.*(f(1,1:m) + 2*sum(f(3:2:(nx-2),1:m)) + f(nx,1:m)); 
%   
%   % Simpson's approximation
%   area=(Tn+2*Un)/3; 
%   
%   % Asymptotically the simpson approximation is a better estimate 
%   % than both Tn and Un. If this is the case, 
%   % then an estimate of the absolute error is
%   
%   epsi = abs(Tn-Un)./4; 
%   %fsiz
%   %area
%   area=squeeze(reshape(area,fsiz));
%   epsi=squeeze(reshape(epsi,fsiz));
%   
  
%end