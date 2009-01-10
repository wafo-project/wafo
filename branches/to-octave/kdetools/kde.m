function f = kde(A,options,varargin)
%KDE Kernel Density Estimator.
%
% CALL:  f = kde(data,options,x1,x2,...,xd)
%
%   f      = kernel density estimate evaluated at meshgrid(x1,x2,..).
%   data   = data matrix, size N x D (D = # dimensions)
%  options = kdeoptions-structure or cellvector of named parameters with
%            corresponding values, see kdeoptset for details.  
%   x1,x2..= vectors defining the points to evaluate the density 
%            (default depending on data)
%  
%  KDE gives a slow, but exact kernel density estimate evaluated at meshgrid(x1,x2,..).
%  Notice that densities close to normality appear to be the easiest for the kernel
%  estimator to estimate and that the degree of estimation difficulty increases with 
%  skewness, kurtosis and multimodality.
%
%  If D > 1 KDE calculates quantile levels by integration. An
%  alternative is to calculate them by ranking the kernel density
%  estimate obtained at the DATA points i.e. use the commands  
%
%      f    = kde(data);
%      tmp  = num2cell(data,1);
%      r    = kdefun(data,[],tmp{:});
%      f.cl = qlevels2(r,f.pl); 
%
%  The first is probably best when estimating the pdf and the latter is the
%  easiest and most robust for multidimensional data when only a visualization
%  of the data is needed.
%
%  For faster estimates try kdebin.
%
% Examples: 
%     data = rndray(1,500,1);
%     x = linspace(sqrt(eps),5,55);
%     plotnorm((data).^(.5)) % gives a straight line => L2 = 0.5 reasonable
%     f = kde(data,{'L2',.5},x);
%     pdfplot(f)
%
% See also  kdefun, mkernel, kdebin, kdeoptset

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
% revised pab July 2008 
% fixed a bug for alpha>0 and Nd=2 or Nd=3
% revised pab Feb2005
% -changed input options into a struct.  
% revised pab 01.01.2001
% - added the possibility that L2 is a cellarray of parametric
%   or non-parametric transformations (secret option)
% revised pab 12.12.1999
% - small modification of example in help header
% revised pab 28.10.1999
%  - added L2
% revised pab 21.10.99 
%  - added alpha to input arguments
%  - made it fully general for d dimensions
%  - HS may be a smoothing matrix
% revised pab 21.09.99  
% adapted from kdetools by Cristian Beardah

  
defaultoptions = kdeoptset;  
% If just 'defaults' passed in, return the default options in g
if ((nargin==1) && (nargout <= 1) &&  isequal(A,'defaults')),
  f = defaultoptions;
  return
end   
error(nargchk(1,inf, nargin))
[n d]=size(A); % Find dimensions of A, 
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
%hsMethod = options.hsMethod;

if isempty(h) 
  h=zeros(1,d);
end

L22 = cell(1,d);
k3=[];
if isempty(L2)
  L2=ones(1,d); % default no transformation
elseif iscell(L2)   % cellarray of non-parametric and parametric transformations
  Nl2 = length(L2);
  if ~(Nl2==1 || Nl2==d), error('Wrong size of L2'), end
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

f=createpdf(d);

nv=length(varargin);
for ix=1:min(nv,d),
  if (any(varargin{ix}<=0) && (L2(ix)~=1)),
    error('xi cannot be negative or zero when L2~=1')
  end
  f.x{ix}=varargin{ix}(:); % make sure it is a column vector.
end
amin=min(A);
if any((amin(L2~=1)<=0)),
  error('DATA cannot be negative or zero when L2~=1')
end

if nv<d
  amax=max(A);
  xyzrange=amax-amin;
  inc = options.inc;
  for ix=nv+1:d
    if ~isempty(k3) && any(k3==ix)
      lo = max(tranproc(1.25*tranproc(amin(ix),L22{ix})....
        -tranproc(amax(ix),L22{ix})/4,fliplr(L22{ix})),sqrt(eps));
    else
      switch L2(ix)
        case 1, lo=amin(ix)-xyzrange(ix)/4;
        case 0, lo=max(sqrt(eps),exp(1.25*log(amin(ix))-log(amax(ix))/4 ));
        otherwise,
          lo=max(sqrt(eps),min(amin(ix)-eps,...
            sign(L2(ix))*(sign(L2(ix))*(1.25*amin(ix)^L2(ix)-amax(ix)^L2(ix)/4 )).^(1/L2(ix))));
      end
    end
    f.x{ix}=transpose(linspace(min(lo,amin(ix)),amax(ix)+xyzrange(ix)/4,inc));
  end
end


X=cell(d,1);
switch d 
 case 1,     [X{1}] = deal(f.x{1});
 case {2 3} ,[X{:}] = meshgrid(f.x{:});
 otherwise ,  
 disp('Dimension of data large, this will take a while.')
 [X{:}]    = ndgrid(f.x{:});
end

[f.f, hs, lambda]=kdefun(A,options,X{:});	
options.hs = hs;
switch lower(kernel(1:4))
  case 'epan', tstr = 'Epanechnikov';%  - Epanechnikov kernel. (default)
  case 'biwe', tstr = 'Biweight'; %     - Bi-weight kernel.
  case 'triw', tstr = 'Triweight';%     - Tri-weight kernel.  
  case 'tria', tstr = 'Triangular';%    - Triangular kernel.
  case 'gaus', tstr = 'Gaussian';%      - Gaussian kernel
  case 'rect', tstr = 'Rectangular';%   - Rectanguler kernel. 
  case 'lapl', tstr = 'Laplace';%       - Laplace kernel.
  case 'logi', tstr = 'Logistic';%      - Logistic kernel.
  otherwise , tstr=[];
end

if alpha>0,
  f.title=['Adaptive Kernel density estimate ( ',tstr,' )'];
  [As ,ind]=unique(A,'rows');
  Ai=num2cell(As,1);method='linear';
  
  switch d
   case 1, lambda=interp1(Ai{1},lambda(ind),X{:},method);
   case 2, lambda=griddata(Ai{:},lambda(ind),X{:},method);
   case 3, lambda=griddata3(Ai{:},lambda(ind),X{:},method);
   otherwise ,
       
     Xflat = zeros(numel(X{1}),d);
     for ix=1:length(X)
         Xflat(:,ix) = X{ix}(:);
     end
     lambda = griddatan(As,lambda(ind),Xflat,method);
     lambda = reshape(lambda,size(X{1}));
  end
  f.lambda=lambda;
else
  f.title=['Kernel density estimate ( ',tstr,' )'];
end
f.note   = f.title;
f.options = options;
%f.kernel = tstr;
%f.alpha  = alpha;
%f.l2     = L2;
if  d>1
 [ql PL] = qlevels(f.f);
 f.cl    = ql;
 f.pl    = PL;
end


