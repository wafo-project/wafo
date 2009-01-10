function f = kdebin(A,options,xlo,xup)
%KDEBIN Binned Kernel Density Estimator.
%
% CALL:  f = kdebin(data,options,xlo,xup)
%
%   f      = pdf structure containing:
%            kernel density estimate evaluated in meshgrid(x1,x2,..).
%   data   = data matrix, size N x D (D = # dimensions )
%  options = kdeoptions-structure or cellvector of named parameters with
%            corresponding values, see kdeoptset for details.
% xlo,xup  = vectors specifying the range of f. Must include the range of the
%            data. (default min(data)-range(data)/4, max(data)-range(data)/4)
%            If a single value of xlo or xup is given then the range is
%            the is the same in all directions.
%
%  KDEBIN gives very fast and accurate kernel density estimate.
%  Notice that densities close to normality appear to be the easiest for
%  the kernel estimator to estimate and that the degree of estimation
%  difficulty increases with skewness, kurtosis and multimodality.
%
%  If D>1 KDEBIN calculates quantile levels by integration. An alternative is 
%  to calculate them by ranking the kernel density estimate obtained at the 
%  points given in DATA i.e. use the commands 
%
%   f    = kdebin( data , kernel)
%   tmp  = num2cell(data,1);
%   r    = evalpdf(f,tmp{:})
%   f.cl = qlevels2(r,f.pl); 
%
%  The first is probably best when estimating the pdf and the latter is the
%  easiest and most robust for 2D data when only a contour visualization of 
%  the data is needed.  
%
% Examples: 
%  data = rndray(1,500,1);
%                               %Box-Cox transform data before estimation
%  f = kdebin(data,{'L2',.5,'inc',64}); 
%  pdfplot(f)
%                               %Non-parametric transformation
%  g   = cdf2tr(edf(data,'wdata',false),mean(data),std(data));
%  opt = kdeoptset('L2',{g},'inc',64);  
%  f1  = kdebin(data,opt);
%  hold on, pdfplot(f1,'r'), hold off
%
% See also  kde, kdeoptset, mkernel

% Reference:  
%  Wand,M.P. and Jones, M.C. (1995) 
% 'Kernel smoothing'
%  Chapman and Hall, pp 182-192
%
%  B. W. Silverman (1986) 
% 'Density estimation for statistics and data analysis'  
%  Chapman and Hall pp 61--66



%Tested on: matlab 5.2, 5.3
% History:
% revised pab Feb 2005
%  -moved options into a options structure  
% revised pab Dec2003
% -removed the binning to a separate function, gridcount.  
% revised pab 05.08.2001
% - fixed a bug in the binning of c.
% - made the binning of c even faster using sparse and binc
% revised pab 27.04.2001
% -added call to mkernel2
% revised pab 19.12.2000
% - added the possibility that L2 is a cellarray of parametric
%   or non-parametric transformations (secret option)
% - fixed a bug in the calculation of xlo and xup
% revised pab 05.01.2000
%  - fixed a bug in back transformation
% revised pab 17.12.99
%  - fixed a bug in back transformation
% revised pab 09.12.99
%  -added alpha,L2
% revised pab 5.11.99
%  - fixed a bug: changed fftshift to ifftshift
% revised pab 21.10.99
%  - added the possibility that Hs is a smoothing matrix 
% revised pab 15.10.99
%  updated documentation
% revised pab 21.09.99  
%  - made it fully general for d dimensions
%  - improoved gridding
% adapted from kdfft1 and kdfft2 from kdetools by Christian C. Beardah 1994

defaultoptions = kdeoptset;  
% If just 'defaults' passed in, return the default options in g
if ((nargin==1) && (nargout <= 1) &&  isequal(A,'defaults')),
  f = defaultoptions;
  return
end  
error(nargchk(1,4, nargin))

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
inc      = options.inc;
L2       = options.L2;
hsMethod = options.hsMethod;
if isempty(h),      h      = zeros(1,d);    end



L22 = cell(1,d);
k3  = [];
if isempty(L2)
  L2 = ones(1,d);     % default no transformation
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




amin = min(A);
if any((amin(L2~=1)<=0))  ,
  error('DATA cannot be negative or zero when L2~=1')
end

%new call
lA = A;  

k1 = find(L2==0); % logaritmic transformation
if any(k1)
  lA(:,k1)=log(A(:,k1));
end
k2 = find(L2~=0 & L2~=1); % power transformation
if any(k2)
  lA(:,k2)=sign(L2(ones(n,1),k2)).*A(:,k2).^L2(ones(n,1),k2);
end
% Non-parametric transformation
for ix = k3,
  lA(:,ix) = tranproc(A(:,ix),L22{ix});
end

amax = max(lA);
amin = min(lA);
xyzrange = amax-amin;



if nargin<3||isempty(xlo)
  xlo=amin-xyzrange/4;
  if any(k1)
      xlo(k1)=max(min(.5*log(eps),amin(k1)-eps),xlo(k1));
  end      
  if any(k2)     
    %xlo(k2)=max(min(sign(L2(k2)).*(eps).^(L2(k2)/2),amin(k2)-eps),xlo(k2) );  
    ki=find(xlo(k2)<=0);
    xlo(k2(ki)) = max(sign(L2(k2(ki))).*(eps).^(L2(k2(ki))/2),amin(k2(ki))-eps); 
    %xlo(k2(ki)) = amin(k2(ki))/2;
  end
  for ix = k3,
    xlo(ix) = max(tranproc(sqrt(eps),L22{ix}),amin(ix)-eps);
  end
  xlo = min(xlo,amin); % make sure xlo<=amin pab 19.12.2000
else
  if length(xlo)<d
    xlo=xlo(1)*ones(1,d);
  end
  if any(k1), xlo(k1)=log(xlo(k1));  end
  if any(k2), xlo(k2)=sign(L2(k2)).*xlo(k2).^L2(k2);  end
  for ix = k3,
    xlo(ix) = tranproc(xlo(ix),L22{ix});
  end
  xlo = min(xlo,amin-eps); % make sure the range of the data is included
end

if nargin<4||isempty(xup)
  xup=amax+xyzrange/4;
else
  if length(xup)<d
    xup=xup(1)*ones(1,d);
  end
  if any(k1), xup(k1) = log(xup(k1));  end
  if any(k2), xup(k2) = sign(L2(k2)).*xup(k2).^L2(k2);  end
  for ix = k3,
    xup(ix) = tranproc(xup(ix),L22{ix});
  end
  xup = max(xup,amax+eps); % make sure the range of the data is included
end


f = createpdf(d);
X = zeros(inc,d);
for ix=1:d
  X(:,ix)=transpose(linspace(xlo(ix),xup(ix),inc));
  f.x{ix}=X(:,ix);
end

hsiz=size(h);
if (min(hsiz)==1)||(d==1)
  if max(hsiz)==1,
    h=h*ones(1,d);
  else
    h = reshape(h,[1 d]); % make sure it has the correct shape
  end

  ind=find(h<=0);
  if any(ind)      % If no value of h has been specified by the user then 
    h(ind)=feval(hsMethod,lA(:,ind),kernel); % calculate optimal values.  
  end
  deth = prod(h);
  hvec = h;
  HG   = 0;
else
  deth=det(h);
  if deth<=0,
    error('bandwidth matrix h must be positive definit')
  end
  hvec=diag(h).';
  h1=inv(h);
  HG=1;
end
options.hs = h;

% The kernel must be symmetric and compactly supported on [-tau tau]
% if the kernel has infinite support then the kernel must have 
% the effective support in [-tau tau], i.e., be negligible outside the range
tau=1;
switch lower(kernel(1:4))
  case 'epan', tstr= 'Epanechnikov';%  - Epanechnikov kernel. (default)
  case 'biwe', tstr= 'Biweight'; %     - Bi-weight kernel.
  case 'triw', tstr= 'Triweight';%     - Tri-weight kernel.  
  case 'tria', tstr= 'Triangular';%    - Triangular kernel.
  case 'gaus', tstr= 'Gaussian';%      - Gaussian kernel
    tau=4;
  case 'rect', tstr= 'Rectangular';%   - Rectangular kernel. 
  case 'lapl', tstr= 'Laplace';%       - Laplace kernel.
    tau=7;
  case 'logi', tstr= 'Logistic';%      - Logistic kernel.
    tau=7;
end

L1 = max(floor(tau*hvec.*(inc-1)./(xyzrange)));
L  = min(L1,inc-1);
if d<2
  fsiz=[inc,1];
  nfft=[2*inc,1];
else
  fsiz=inc*ones(1,d);
  nfft=2*inc.*ones(1,d);
end

dx = (xup-xlo)./(inc-1);
X1 = cell(d,1);
if HG,
    % new call
    X1=num2cell((-L:L)'*dx,1);
    % X1=num2cell((0:L)'*dx,1);
    if d<=3,
      [X1{:}]=meshgrid(X1{:});
    else  
      disp('Dimension of data large, this will take a while.')
      [X1{:}]=ndgrid(X1{:});
    end
    
    for ix=1:d
      X1{ix}=X1{ix}(:);
    end
    X1=num2cell([X1{:}]*h1,1);
    for ix=1:d
      X1{ix}=reshape(X1{ix},(2*L+1)*ones(1,d));
    end  
else
  switch d 
    case 1, X1{1}=transpose(dx/h*(-L:L));  
    case {2,3} 
      X1      = num2cell((-L:L)'*(dx./h),1);
      [X1{:}] = meshgrid(X1{:});
    otherwise ,  
      disp('Dimension of data large, this will take a while.')
      % new call
      X1      = num2cell((-L:L)'*(dx./h),1);
      [X1{:}] = ndgrid(X1{:});
  end
end
% Obtain the kernel weights.

if 1 %HG, %
  kw = zeros(nfft);
  indk(1:d) = {(inc-L+1):(inc+L+1)};
  kw(indk{:}) = mkernel(X1{:},kernel)/(n*deth);
  % Apply 'ifftshift' to the kernel weights, kw.
  kw = ifftshift(kw);
else
  kw1 = zeros(floor(nfft/2)+1);  
  indk(1:d)  = {1:(L+1)};
  kw1(indk{:}) = mkernel(X1{:},kernel)/(n*deth);
  kw = fftce(kw1); % circulant embedding
  clear kw1 
end

% Find the binned kernel weights, c.
c = gridcount(lA,X);

% Perform the convolution.
z = real(ifftn(fftn(c,nfft).*fftn(kw)));
clear kw

% New call pab 19.04.2001
%indk = repmat({1:inc}, d, 1); % initialize subscripts
indk(1:d) = {1:inc};  % alternatively
f.f = z(indk{:}).*(z(indk{:})>0);
clear z indk

%f.kernel=tstr;
%f.hs=h;

if (alpha>0), % adaptive kde
  f.f=f.f(:); 
  indc=transpose(find(c>0));

  if 0 % clipping to make sure lambda do not get too large
    minf=sqrt(eps); % make sure we sum over f.f(xi)>0
    ind=find(f.f(:)>minf); 
  else
    ind=indc;
  end
  % this gives smaller lambda than kdefun! Why???
  g=exp(sum(log(f.f(ind)))/length(ind(:))); % geometric mean of f
  lambda=ones(fsiz);
 
  lambda(ind)=(f.f(ind)/g).^(-alpha); %  local bandwidth factor 
  
 
  Nx=inc^d;
  f.f=zeros(Nx,1);
  switch d 
    case 1, X1{1}=transpose(dx*(1:inc));
    case {2,3}
      X1=num2cell((1:inc)'*dx,1);
      [X1{:}]=meshgrid(X1{:});
    otherwise ,
      %disp('Dimension of data large, this will take a while.')
      % new call
      X1=num2cell((1:inc)'*dx,1);
      [X1{:}]=ndgrid(X1{:});
  end
  for ix=1:d
    X1{ix}=X1{ix}(:);
  end
  lX=[X1{:}];
  
  if ~HG ,%(min(hsiz)==1)|(d==1)
    lX=lX./h(ones(Nx,1),:); 
    for ix=indc,%Sum over all data points
      Avec=lX(ix,:);
      Xnn=(lX-Avec(ones(Nx,1),:))/(lambda(ix));  
      f.f=f.f+mkernel2(Xnn,kernel)*c(ix)/(lambda(ix)^d);% pab 27.04.01
    end
  else % adaptive kde % fully general
    for ix=indc,%  Sum over all data points
     Avec=lX(ix,:);
     Xnn=(lX-Avec(ones(Nx,1),:))*(h1/lambda(ix));
     f.f=f.f+mkernel2(Xnn,kernel)*c(ix)/(lambda(ix)^d);% pab 27.04.01
    end    
  end
  f.f=reshape(f.f,fsiz)/(n*deth);
  %f.alpha=alpha;
  f.lambda=lambda;
  f.title=['Adaptive Binned Kernel density estimate ( ',tstr,' )'];
else
  f.title=['Binned Kernel density estimate ( ',tstr,' )'];
end  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        transforming back         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if any(k1)||any(k2)||any(k3),
  X1=cell(d,1);
  switch d
    case 1,X1{1}=f.x{1};
    case {2 3}, [X1{:}]=meshgrid(f.x{:});
    otherwise, [X1{:}]=ndgrid(f.x{:});
  end
  if any(k1), % L2=0 i.e. logaritmic transformation
    for ix=k1
      f.f=f.f./exp(X1{ix});
      f.x{ix}=exp(f.x{ix});
    end
    if any(max(abs(diff(f.f)))>10)
      disp('Warning: Numerical problems may have occured due to the logaritmic')
      disp('transformation. Check the KDE for spurious spikes')
    end
  end
  if any(k2) % L2~=0 i.e. power transformation
    for ix=k2
      f.f=f.f.*((sign(L2(ix)).*X1{ix}.^(1/L2(ix)))....
	  .^(L2(ix)-1))*L2(ix)*sign(L2(ix));
      f.x{ix}=(sign(L2(ix)).*f.x{ix}).^(1/L2(ix));
    end
    if any(max(abs(diff(f.f)))>10)
      disp('Warning: Numerical problems may have occured due to the power')
      disp('transformation. Check the KDE for spurious spikes')
    end
  end
  if any(k3), % non-parametric transformation
    oneC = ones(inc,1);
    oneD = ones(1,d);
    for ix=k3
      gn  = L22{ix}; 
      Gn  = fliplr(L22{ix});
      x0  = tranproc(f.x{ix},Gn);
      if any(isnan(x0)),
	error('The transformation does not have a strictly positive derivative.')
      end
      hg1  = tranproc([x0 oneC],gn);
      der1 = abs(hg1(:,2)); % dg(X)/dX = 1/(dG(Y)/dY)
      % alternative 2
      %pp  = smooth(Gn(:,1),Gn(:,2),1,[],1);
      %dpp = diffpp(pp);
      %der1 = 1./abs(ppval(dpp,f.x{ix}));
      % Alternative 3
      %pp  = smooth(gn(:,1),gn(:,2),1,[],1);
      %dpp = diffpp(pp);
      %%plot(hg1(:,1),der1-abs(ppval(dpp,x0)))
      %der1 = abs(ppval(dpp,x0));
      if any(der1<=0), 
        error(['The transformation must have a strictly positive derivative'])
      end
	
      switch d,
        case 1,  f.f = f.f.*der1;
        case {2,3},
          oneD(ix)  = inc;
          oneD(1:2) = oneD(2:-1:1);
          f.f = f.f.*repmat(reshape(der1,oneD),inc+1-oneD);
          oneD(1:2) = oneD(2:-1:1);
          oneD(ix)  = 1;
        otherwise
          oneD(ix) = inc;
          f.f = f.f.*repmat(reshape(der1,oneD),inc+1-oneD);
          oneD(ix) = 1;
      end
      
      f.x{ix} = x0;
    end
    if any(max(abs(diff(f.f)))>10)
      disp('Warning: Numerical problems may have occured due to the power')
      disp('transformation. Check the KDE for spurious spikes')
    end
  end
  if 1
  tmp=f;
 % pdfplot(f)
  for ix=[k1 k2 k3]
    f.x{ix}=transpose(linspace(f.x{ix}(1),f.x{ix}(end),inc));
  end
  switch d
    case 1,X1{1}=f.x{1};
    case {2 3}, [X1{:}] = meshgrid(f.x{:});
    otherwise,  [X1{:}] = ndgrid(f.x{:});
  end
  % interpolating to obtain equidistant grid
  if 1,
    f.f = evalpdf(tmp,X1{:},'linear');
  else
    tmp.f = log(tmp.f+eps);
    f.f=max(exp(evalpdf(tmp,X1{:},'spline'))-eps,0);
  end
  clear tmp X1
  end
  %who
end
f.note   = f.title;
%f.kernel = tstr;
%f.alpha  = alpha;
%f.l2     = L2;
f.options = options;
f.n      = n;
if d>1
  [ql PL] = qlevels(f.f);
  f.cl = ql;
  f.pl = PL;
end







