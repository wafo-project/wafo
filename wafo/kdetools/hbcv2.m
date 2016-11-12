function h=hbcv2(A,kernel)
%HBCV2 Biased Cross-Validation smoothing parameter for 2D data.
%          
% CALL: hs = hbcv2(data,kernel); 
% 
%   hs     = smoothing parameter
%   data   = two column data matrix
%   kernel = 'gaussian'      - Gaussian kernel
%            'epanechnikov'  - Epanechnikov kernel.
%            'biweight'      - Bi-weight kernel.
%            'triweight'     - Tri-weight kernel.  
%            'triangluar'    - Triangular kernel. 
%            'rectangular'   - Rectanguler kernel. 
%            'laplace'       - Laplace kernel.
%            'logistic'      - Logistic kernel.
%
%  Note that only the first 4 letters of the kernel name is needed.
%  
%  HBCV is a hybrid of crossvalidation and direct plug-in estimates.
%  The main attraction of HBCV is that it is more stable than HLSCV in
%  the sense that its asymptotic variance is considerably lower. However,
%  this reduction in variance comes at the expense of an increase in
%  bias, with HBCV tending to be larger than the HNS estimate.
%  Asymptotically HBCV has a relative slow convergence rate.
%
%  Example: 
%  % data = rndnorm(0, 1,5,2);
%  data = [-0.0233845632050972   0.9070186193622006;...
%           0.6529594866766634   1.3689145060433903;...
%           0.4477857310723146  -0.6311953712037597;...
%          -1.9256785038579962   0.5886257667993168;...
%          -0.5290011931824666  -0.3602090880229930];
%  assert(hbcv2(data,'gaus'), [0.482681614203333 1.133242297848955], 1e-10);
%  assert(hbcv2(data,'epan'), [1.05495723163671 2.55257402255233], 1e-10);
%  assert(hbcv2(data,'biwe'), [1.24293653366899 2.98949541257312], 1e-10);
%  assert(hbcv2(data,'triw'), [1.43099550062034, 3.39546227138374], 1e-10);
%  assert(hbcv2(data,'tria'), [ 1.15255694491467, 2.76211028087912], 1e-10);
%  assert(hbcv2(data,'rect'), [0.837614274875370, 1.997573853072259], 1e-10);
%  assert(hbcv2(data,'lapl'), [0.362634943481836, 0.873248947994519], 1e-10);
%  assert(hbcv2(data,'logi'), [ 0.266350664349379, 0.652154211531075], 1e-10);
%
% See also  hste, hboot, hns, hos, hldpi, hlscv, hscv, hstt, kde, kdefun   


% tested on : matlab 5.2
% history:
% revised pab aug2005
% -added ad hoc support for other kernels than Gaussian
% Revised pab dec 2003
%  -fixed some bugs  
%  -added todo comments 
% revised pab 20.10.1999
%   updated to matlab 5.2
% changed input arguments
% taken from kdetools     Christian C. Beardah 1993-94

% TODO % Add support for other kernels than Gaussian  


% Transform A so rows have st. dev.=1:

[n d]=size(A);
if d~=2,
  error('Wrong shape of data')
end
oldA = A;
A = (oldA-ones(n,1)*mean(oldA))*diag(1./std(oldA)); % centre data about
                                                    % its mean
maxit=20;

tol=1e-2;

gradf=zeros(2,1);

delta=1e-6;

%i=0;

vx(1)=hste(A(:,1),kernel);
vy(1)=hste(A(:,2),kernel);

x0=vx(1);
y0=vy(1);

%fn(1)=bcv2(A,vx(1),vy(1));



% R= int(mkernel(x)^2)
% mu2= int(x^2*mkernel(x))
kernel2 = 'gauss';
% values for gaussian kernel
[mu2ns,Rns] = kernelstats(kernel2);
Rns = Rns^d;
STEconstant2 = Rns /(mu2ns^(2)*n);

[mu2,R] = kernelstats(kernel);
R = R^d;

STEconstant = R /(mu2^(2)*n);


Nm = 31;
M = zeros(Nm,2);

i=1; % i is the iteration counter.

while i <= maxit,

  gradf(1)=(bcv2(A,x0+delta,y0)-bcv2(A,x0-delta,y0))/(2*delta);
  gradf(2)=(bcv2(A,x0,y0+delta)-bcv2(A,x0,y0-delta))/(2*delta);

  gradf=gradf/(i*norm(gradf));

  r0=[x0 y0]';

  for j=1:Nm,
    M(j,:)=(r0+(j-16)*gradf/15)';
  end;

  L=(M(:,1)>0) & (M(:,2)>0);

  L=remzero([M,L]); % check this wrong???

  M=L(:,1:2);
  Nm = size(M,1);
  t = zeros(Nm,1);
  for j=1:Nm,
    t(j,1)=bcv2(A,M(j,1),M(j,2));
  end;

  [N,I]=min(t);

  r=M(I,:)';

  %vx(i+1)=r(1); % Store iterate.
  %vy(i+1)=r(2); 

  %fn(i+1)=N;

  if (abs(norm(r-r0)) < tol),

    h=r';

% Transform values of h back:

    h=(diag(std(oldA))*h'*(STEconstant/STEconstant2)^(1/6))';

    return;
  end;

  i=i+1; % Update the iteration counter.

  x0=r(1);
  y0=r(2);

  clear t
 
end; % of while statement.

disp('Maximum number of iterations reached - terminating execution.');
return

function z=bcv2(A,h1,h2)
% BCV2        Biased cross-validation criterion function.
%          
%             BCV2(A,H1,H2)
% 
%             Christian C. Beardah 1993-94 

n=length(A);

A1=A(:,1);

A2=A(:,2);

M1=A1*ones(size(A1'));

T1=(M1-M1')/h1;

M2=A2*ones(size(A2'));

T2=(M2-M2')/h2;

z=1/(4*pi*n*h1*h2);

T=T1.^2+T2.^2;

T=(T.^2-8*T+8).*mkernel(T1,'gauss').*mkernel(T2,'gauss');

T=T-diag(diag(T));
z=z+sum(sum(T))/(4*n*(n-1)*h1*h2);

return

function X=remzero(X)
%REMZERO Removes rows containing zeroes from a matrix or vector.
%
%  CALL Y = remzero(X)
%  
% Example: 
%   x = [1 2 0 4 0 6].'
%   x = remzero(x); 
%

%History
%   
% by            Christian C. Beardah 1993

[r c]=size(X);

if r>1 && c>1,
  [i, j]=find(X==0);
  X(i,:)=[];
else
  X=X(X~=0);
end;

 return

 