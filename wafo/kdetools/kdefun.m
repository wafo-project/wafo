function [f, hs,lambda]= kdefun(A,options,varargin)
%KDEFUN  Kernel Density Estimator.
%
% CALL:  [f, hs] = kdefun(data,options,x1,x2,...,xd)
%
%   f      = kernel density estimate evaluated at x1,x2,...,xd.
%   data   = data matrix, size N x D (D = # dimensions)
%  options = kdeoptions-structure or cellvector of named parameters with
%            corresponding values, see kdeoptset for details.    
%   x1,x2..= vectors/matrices defining the points to evaluate the density 
%  
%  KDEFUN gives a slow, but exact kernel density estimate evaluated at x1,x2,...,xd.
%  Notice that densities close to normality appear to be the easiest for the kernel
%  estimator to estimate and that the degree of estimation difficulty increases with 
%  skewness, kurtosis and multimodality.
%
%  If D > 1 KDE calculates quantile levels by integration. An
%  alternative is to calculate them by ranking the kernel density
%  estimate obtained at the points DATA  i.e. use the commands  
%
%      f    = kde(data);
%      r    = kdefun(data,[],num2cell(data,1));
%      f.cl = qlevels2(r,f.PL); 
%
%  The first is probably best when estimating the pdf and the latter is the
%  easiest and most robust for multidimensional data when only a visualization 
%  of the data is needed.
%
%  For faster estimates try kdebin.
%
% Examples: 
%   data = rndray(1,500,1);
%   x = linspace(sqrt(eps),5,55);
%   plotnorm((data).^(.5)); % gives a straight line => L2 = 0.5 reasonable
%   f = kdefun(data,{'L2',.5},x);
%   plot(x,f,x,pdfray(x,1),'r');
%
%   close all;
%
% See also  kde, mkernel, kdebin

% Reference:  
%  B. W. Silverman (1986) 
% 'Density estimation for statistics and data analysis'  
%  Chapman and Hall , pp 100-110
%
%  Wand, M.P. and Jones, M.C. (1995) 
% 'Kernel smoothing'
%  Chapman and Hall, pp 43--45




%Tested on: matlab 5.2
% History:
% revised pab Feb2004
%  -options moved into a structure  
% revised pab Dec2003
% -removed some code  
% revised pab 27.04.2001
% - changed call from mkernel to mkernel2 (increased speed by 10%)
% revised pab 01.01.2001
% - added the possibility that L2 is a cellarray of parametric
%   or non-parametric transformations (secret option)
% revised pab 14.12.1999
%  - fixed a small error in example in help header
% revised pab 28.10.1999
%  - added L2
% revised pab 21.10.99 
%  - added alpha to input arguments
%  - made it fully general for d dimensions
%  - HS may be a smoothing matrix
% revised pab 21.09.99  
%  - adapted from kdetools by Christian Beardah

  defaultoptions = kdeoptset;  
% If just 'defaults' passed in, return the default options in g
if ((nargin==1) && (nargout <= 1) &&  isequal(A,'defaults')),
  f = defaultoptions;
  return
end   
error(nargchk(1,inf, nargin))

[n, d]=size(A); % Find dimensions of A, 
               % n=number of data points,
               % d=dimension of the data.  
if (nargin<2 || isempty(options))
  options  = defaultoptions;
else
  switch lower(class(options))
   case {'char','struct'},
    options = kdeoptset(defaultoptions,options);
   case {'cell'}
   
      options = kdeoptset(defaultoptions,options{:});
   otherwise
    error('Invalid options')
  end
end
kernel   = options.kernel;
h        = options.hs;
alpha    = options.alpha;
L2       = options.L2;
hsMethod = options.hsMethod;

if isempty(h) 
  h=zeros(1,d);
end

L22 = cell(1,d);
k3=[];
if isempty(L2)
  L2=ones(1,d); % default no transformation
elseif iscell(L2)   % cellarray of non-parametric and parametric transformations
  Nl2 = length(L2);
  if ~(Nl2==1||Nl2==d), error('Wrong size of L2'), end
  [L22{1:d}] = deal(L2{1:min(Nl2,d)});
  L2 = ones(1,d); % default no transformation
  for ix=1:d,
    if length(L22{ix})>1,
      k3=[k3 ix];       % Non-parametric transformation
    else 
     L2(ix) = L22{ix};  % Parameter to the Box-Cox transformation
    end
  end
elseif length(L2)==1
  L2=L2(:,ones(1,d));
end

amin=min(A);
if any((amin(L2~=1)<=0))  ,
  error('DATA cannot be negative or zero when L2~=1')
end


nv=length(varargin);
if nv<d,
  error('some or all of the evaluation points x1,x2,...,xd is missing')
end

xsiz = size(varargin{1}); % remember size of input
Nx   = prod(xsiz);
X    = zeros(Nx,d);
for ix=1:min(nv,d),
  if (any(varargin{ix}(:)<=0) && (L2(ix)~=1)),
    error('xi cannot be negative or zero when L2~=1')
  end
  X(:,ix)=varargin{ix}(:); % make sure it is a column vector
end


%new call
lX = X; %zeros(Nx,d);
lA = A; %zeros(size(A));  

k1 = find(L2==0); % logaritmic transformation
if any(k1)
  lA(:,k1)=log(A(:,k1));
  lX(:,k1)=log(X(:,k1));
end
k2=find(L2~=0 & L2~=1); % power transformation
if any(k2)
  lA(:,k2)=sign(L2(ones(n,1),k2)).*A(:,k2).^L2(ones(n,1),k2);
  lX(:,k2)=sign(L2(ones(Nx,1),k2)).*X(:,k2).^L2(ones(Nx,1),k2);
end
% Non-parametric transformation
for ix = k3,
  lA(:,ix) = tranproc(A(:,ix),L22{ix});
  lX(:,ix) = tranproc(X(:,ix),L22{ix});
end


hsiz=size(h);
if (min(hsiz)==1)||(d==1)
  if max(hsiz)==1,
    h=h*ones(1,d);
  else
    h=reshape(h,[1,d]); % make sure it has the correct dimension
  end;
  ind=find(h<=0);
  if any(ind)    % If no value of h has been specified by the user then 
    h(ind)=feval(hsMethod,lA(:,ind),kernel); % calculate automatic values.  
  end
  deth = prod(h);
else  % fully general smoothing matrix
  deth = det(h);
  if deth<=0
    error('bandwidth matrix h must be positive definit')
  end
end

if alpha>0
  Xn   = num2cell(lA,1);
  opt1 = kdeoptset('kernel',kernel,'hs',h,'alpha',0,'L2',1);
  f2   = kdefun(lA,opt1,Xn{:}); % get a pilot estimate by regular KDE (alpha=0)
  g    = exp(sum(log(f2))/n);
  
  lambda=(f2(:)/g).^(-alpha);
else
  lambda=ones(n,1);
end





f=zeros(Nx,1);
if (min(hsiz)==1)||(d==1)
  for ix=1:n,     % Sum over all data points
    Avec=lA(ix,:);
    Xnn=(lX-Avec(ones(Nx,1),:))./(h(ones(Nx,1),:) *lambda(ix));
    f = f + mkernel2(Xnn,kernel)/lambda(ix)^d;
  end
else % fully general
  h1=inv(h);
  for ix=1:n,     % Sum over all data points
    Avec=lA(ix,:);
    Xnn=(lX-Avec(ones(Nx,1),:))*(h1/lambda(ix));
    f = f + mkernel2(Xnn,kernel)/lambda(ix)^d;
  end    
end
f=f/(n*deth);

% transforming back
if any(k1), % L2=0 i.e. logaritmic transformation
  for ix=k1
    f=f./X(:,ix);
  end
  if any(max(abs(diff(f)))>10)
    disp('Warning: Numerical problems may have occured due to the logaritmic')
    disp('transformation. Check the KDE for spurious spikes')
  end
end
if any(k2) % L2~=0 i.e. power transformation
  for ix=k2
    f=f.*(X(:,ix).^(L2(ix)-1))*L2(ix)*sign(L2(ix));
  end
  if any(max(abs(diff(f)))>10)
    disp('Warning: Numerical problems may have occured due to the power')
    disp('transformation. Check the KDE for spurious spikes')
  end
end
if any(k3), % non-parametric transformation
  oneC = ones(Nx,1);
  for ix=k3
    gn  = L22{ix}; 
    %Gn  = fliplr(L22{ix});
    %x0  = tranproc(lX(:,ix),Gn);
    if any(isnan(X(:,ix))),
      error('The transformation does not have a strictly positive derivative.')
    end
    hg1  = tranproc([X(:,ix) oneC],gn);
    der1 = abs(hg1(:,2)); % dg(X)/dX = 1/(dG(Y)/dY)
    % alternative 2
    %pp  = cssmooth(Gn(:,1),Gn(:,2),1,[],1);
    %dpp = diffpp(pp);
    %der1 = 1./abs(ppval(dpp,f.x{ix}));
    % Alternative 3
    %pp  = cssmooth(gn(:,1),gn(:,2),1,[],1);
    %dpp = diffpp(pp);
    %%plot(hg1(:,1),der1-abs(ppval(dpp,x0)))
    %der1 = abs(ppval(dpp,x0));
    if any(der1<=0), 
      error('The transformation must have a strictly positive derivative')
    end
    f = f.*der1;
  end
  if any(max(abs(diff(f)))>10)
    disp('Warning: Numerical problems may have occured due to the power')
    disp('transformation. Check the KDE for spurious spikes')
  end
end	
  
f=reshape(f,xsiz); % restore original shape
if nargout>1
  hs=h;
end

return






