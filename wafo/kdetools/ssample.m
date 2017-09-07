function s=ssample(A,m,varargin)
%SSAMPLE  Random sampling from a smoothed empirical distribution
%
% CALL: s = ssample(data,m,options)
%  
%   s      = sampled selection from data,  size m x D
%   data   = data matrix, size N x D (D = # dimensions)
%   m      = sampling size 
%  options = kdeoptions-structure or cellvector of named parameters with
%            corresponding values, see kdeoptset for details.  
% 
%  SSAMPLE(DATA,M) selects a random sample of M data points from the
%  multivariate data-set in the matrix DATA.
%
% Example:
%     data = rndnorm(0,1,500,1);
%     s = ssample(data,100,'kernel','gauss');
%     f = kdebin(s);  
%     pdfplot(f); hold on;
%     x = linspace(-5,5);
%     plot(x,pdfnorm(x),'r');
%
%     close all;
%
% See also sample
  
%  Reference:  
%  B. W. Silverman (1986) 
% 'Density estimation for statistics and data analysis'  
%  Chapman and Hall, pp. 143 
%  
%  Wand, M. P. and Jones, M. C. (1995) 
% 'Density estimation for statistics and data analysis'  
%  Chapman and Hall, 

% History:
% revised pab dec2003  
% Revised by gl 13.07.2000
%    changed ind generation to avoid stats
% by pab 10.12.1999

% TODO % Needs further checking 
  
   defaultoptions = kdeoptset('kernel','gauss');  
% If just 'defaults' passed in, return the default options in g
if ((nargin==1) && (nargout <= 1) &&  isequal(A,'defaults')),
  s = defaultoptions;
  return
end   

%error(nargchk(2,inf, nargin))
narginchk(2,inf)
[n d]=size(A);

if (nargin<3)
  options  = defaultoptions;
else
  switch lower(class(varargin{1}))
   case {'char','struct'},
    options = kdeoptset(defaultoptions,varargin{:});
   case {'cell'}
   
      options = kdeoptset(defaultoptions,varargin{1}{:});
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

k3 = [];
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

%new call
lA = A;  

k1=find(L2==0); % logaritmic transformation
if any(k1)
  lA(:,k1)=log(A(:,k1));
end
k2=find(L2~=0 & L2~=1); % power transformation
if any(k2)
  lA(:,k2)=sign(L2(ones(n,1),k2)).*A(:,k2).^L2(ones(n,1),k2);
end
for ix = k3,
  lA(:,ix) = tranproc(A(:,ix),L22{ix});
  %lX(:,ix) = tranproc(X(:,ix),L22{ix});
end

hsiz=size(h);
if (min(hsiz)==1)||(d==1)
  if max(hsiz)==1,
    h=h*ones(1,d);
  else
    h=reshape(h,[1 d]); % make sure it has the correct shape
  end

  ind=find(h<=0);
  if any(ind)      % If no value of h has been specified by the user then 
    h(ind)=feval(hsMethod,lA(:,ind),kernel); % calculate automatic values.  
  end
  hmat=diag(h);
else
  deth=det(h);
  if deth<=0,
    error('bandwidth matrix h must be positive definit')
  end
  hmat=h;
end

ma=mean(lA);
sa=cov(lA); 

E=mkernelrnd(kernel,m,d);
% ind = unidrnd(n,m,1);  
% new generation of random sample
ind = floor(n*rand(m,1))+1;
if alpha>0,
  tmp=num2cell(lA,1);
  opt1 = kdeoptset('kernel',kernel,'hs',h,'alpha',alpha,'L2',L2);
  [f, hs,lambda]= kdefun(lA,opt1,tmp{:});
 E=E.*lambda(ind,ones(1,d));
end

 switch lower(kernel(1:4))
  case 'epa1', sk=eye(d)/5;
  case {'norm','gaus'}, sk=eye(d);
  end

s=ma(ones(m,1),:)+(lA(ind,:)-ma(ones(m,1),:)+ E*hmat)/sqrtm(eye(d)+hmat^2*sk/sa);

% transforming back
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if any(k1), % L2=0 i.e. logaritmic transformation
 s(:,k1)=exp(s(:,k1));
end
if any(k2) % L2~=0 i.e. power transformation
  for ix=k2
    s(:,ix)=sign(L2(ix)).*(s(:,ix)).^(1/L2(ix));
  end
end
