function h=hstt(A,kernel,inc,maxit,releps,abseps)
%HSTT Scott-Tapia-Thompson estimate of smoothing parameter.
%
% CALL: hs = hstt(data,kernel)
%
%       hs = one dimensional value for smoothing parameter
%            given the data and kernel.  size 1 x D
%   data   = data matrix, size N x D (D = # dimensions )
%   kernel = 'epanechnikov'  - Epanechnikov kernel. (default)
%            'biweight'      - Bi-weight kernel.
%            'triweight'     - Tri-weight kernel.  
%            'triangular'    - Triangular kernel.
%            'gaussian'      - Gaussian kernel
%            'rectangular'   - Rectangular kernel. 
%            'laplace'       - Laplace kernel.
%            'logistic'      - Logistic kernel.  
%
% HSTT returns Scott-Tapia-Thompson (STT) estimate of smoothing
% parameter. This is a Solve-The-Equation rule (STE).
% Simulation studies shows that the STT estimate of HS
% is a good choice under a variety of models. A comparison with
% likelihood cross-validation (LCV) indicates that LCV performs slightly
% better for short tailed densities.
% However, STT method in contrast to LCV is insensitive to outliers.
% 
%  Example: 
%   x  = rndnorm(0,1,50,1);
%   hs = hstt(x,'gauss');
%
% See also  hste, hbcv, hboot, hos, hldpi, hlscv, hscv, kde, kdebin 

%tested on: matlab 5.2
% history
% Revised pab dec2003  
% added inc and maxit, releps and abseps as inputs  
% revised pab 16.10.1999
% updated string comparison to matlab 5.x
%
% taken from kdetools by  Christian C. Beardah 1995 

% Reference:  
%  B. W. Silverman (1986) 
% 'Density estimation for statistics and data analysis'  
%  Chapman and Hall, pp 57--61 

if nargin<2||isempty(kernel)
  kernel='gauss';
end
if nargin<3||isempty(inc)
  inc = 128/2;
end
if nargin<4 || isempty(maxit)
  maxit=100;
end
if nargin<5 || isempty(releps)
  releps = 0.01;
end
if nargin<6 || isempty(abseps)
  abseps = 0.0;
end
[n, d] = size(A);
if (n==1) && (d>1),
  A=A.';
  n=d;
  d=1;
end 

[mu2,R] = kernelstats(kernel);

STEconstant = R /(mu2^(2)*n);


h = hns(A,kernel);

% This iteration can cycle 
% Don't allow more than maxit iterations.

kopt = kdeoptset('kernel',kernel,'inc',inc);


for dim = 1:d,
  count = 1;
  h_old = 0;
  h1 = h(dim);
  
  while ((abs((h_old-h1))>max(releps*h1,abseps)) && (count<maxit)),

    h_old = h1;
    
    kopt.hs = h1;
    f = kdebin(A(:,dim),kopt);

    delta=f.x{1}(2)-f.x{1}(1);

    % Estimate psi4=R(f'') using simple finite differences and quadrature.

    ix=2:(inc-1);
    z = ((f.f(ix+1)-2*f.f(ix)+f.f(ix-1))/delta^2).^2;
    
    psi4 = delta*sum(z(:));

    h1 = (STEconstant/psi4)^(1/5);

    count = count+1;
  end;
  
  h(dim) = h1;
  
  if count>= maxit
    disp('The obtained value did not converge.')
  end
end