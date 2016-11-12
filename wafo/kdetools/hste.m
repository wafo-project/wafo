function h=hste(A,kernel,h,inc,maxit,releps,abseps)
%HSTE 2-Stage Solve the Equation estimate of smoothing parameter.
%
% CALL:  hs = hste(data,kernel,h0)
% 
%       hs = one dimensional value for smoothing parameter
%            given the data and kernel.  size 1 x D
%   data   = data matrix, size N x D (D = # dimensions )
%   kernel = 'gaussian'  - Gaussian kernel (default)
%             ( currently the only supported kernel)
%       h0 = initial starting guess for hs (default h0=hns(A,kernel))
%
%  Example: 
%  % data  = rndnorm(0,1,5,2);
%  data = [-0.0233845632050972   0.9070186193622006;...
%           0.6529594866766634   1.3689145060433903;...
%           0.4477857310723146  -0.6311953712037597;...
%          -1.9256785038579962   0.5886257667993168;...
%          -0.5290011931824666  -0.3602090880229930];
%   assert(hste(data,'gauss'), [0.564411435057012   0.431314424407532], 1e-10);
%
% See also  hbcv, hboot, hos, hldpi, hlscv, hscv, hstt, kde, kdefun

% Reference:  
%  B. W. Silverman (1986) 
% 'Density estimation for statistics and data analysis'  
%  Chapman and Hall, pp 57--61
%
%  Wand,M.P. and Jones, M.C. (1986) 
% 'Kernel smoothing'
%  Chapman and Hall, pp 74--75
  
% tested on: matlab 5.2
% revised pab aug 2005
% - All kernels supported
% Revised pab dec 2003
% added todo comments  
% added inc, maxit,abseps,releps as inputs  
% revised pab 16.10.1999
% added h0 as input
%  the gridding is made much faster
% taken from kdetools   Christian C. Beardah 1995 

% TODO % NB: this routine can be made faster:
% TODO % replace the iteration in the end with a Newton Raphson scheme


if nargin<2||isempty(kernel)
 kernel='gauss';
end
if nargin<3||isempty(h)
  h=hns(A,kernel);
end
if nargin<4 || isempty(inc)
  inc=128;
end
if nargin<5 || isempty(maxit)
  maxit = 100;
end
if nargin<6 || isempty(releps)
  releps = 0.01;
end
if nargin<7 || isempty(abseps)
  abseps = 0.0;
end
[n, d] = size(A);
if (n==1) && (d>1),
  A=A.';
  n=d;
  d=1;
end


% R   = int(mkernel(x)^2)
% mu2 = int(x^2*mkernel(x))
[mu2,R] = kernelstats(kernel);

STEconstant = R /(mu2^(2)*n);


nfft = inc*2;


xmin   = min(A);    % Find the minimum value of A.
xmax   = max(A);    % Find the maximum value of A.
xrange = xmax-xmin; % Find the range of A.

sigmaA = std(A);
iqr = abs(diff(qlevels2(A,[75 25]),1,1+(d==1))); % interquartile range
k = find(iqr>0);
if any(k)
  sigmaA(k) = min(sigmaA(k), iqr(k)/1.349);
end

% xa holds the x 'axis' vector, defining a grid of x values where 
% the k.d. function will be evaluated.

ax1 = xmin-xrange/8;
bx1 = xmax+xrange/8;

kernel2 = 'gaus'; % kernel
[mu2,R] = kernelstats(kernel2);
STEconstant2 = R /(mu2^(2)*n);
for dim = 1:d
  s = sigmaA(dim);
  ax = ax1(dim);
  bx = bx1(dim);
  
  xa = linspace(ax,bx,inc).'; 
  xn = linspace(0,bx-ax,inc);
  
  c = gridcount(A(:,dim),xa);
  
  %deltax = (bx-ax)/(inc-1);
  %xn     = (0:(inc-1))*deltax;
  %binx   = floor((A(:,dim)-ax)/deltax)+1;
  % Obtain grid counts
  %c = full(sparse(binx,1,(xa(binx+1)-A(:,dim)),inc,1)+...
	%   sparse(binx+1,1,(A(:,dim)-xa(binx)),inc,1))/deltax;


  % Step 1
  psi6NS = -15/(16*sqrt(pi)*s^7);
  psi8NS = 105/(32*sqrt(pi)*s^9);

  % Step 2
  [k40,k60] = deriv(0,kernel2);
  g1 = (-2*k40/(mu2*psi6NS*n))^(1/7);
  g2 = (-2*k60/(mu2*psi8NS*n))^(1/9);

  % Estimate psi6 given g2.
  [kw4,kw6] = deriv(xn/g2,kernel2); % kernel weights.
  kw   = [kw6,0,kw6(inc:-1:2)].';             % Apply 'fftshift' to kw.
  z    = real(ifft(fft(c,nfft).*fft(kw)));     % convolution.
  psi6 = sum(c.*z(1:inc))/(n*(n-1)*g2^7);

  % Estimate psi4 given g1.
  kw4  = deriv(xn/g1,kernel2); % kernel weights.
  kw   = [kw4,0,kw4(inc:-1:2)]';            % Apply 'fftshift' to kw.
  z    = real(ifft(fft(c,nfft).*fft(kw)));    % convolution.
  psi4 = sum(c.*z(1:inc))/(n*(n-1)*g1^5);



  %
  h1    = h(dim);
  h_old = 0;
  count = 0;
  
  while ((abs(h_old-h1)>max(releps*h1,abseps)) && (count < maxit)),
    count = count+1;
    % save old value
    h_old = h1;
  
    % Step 3
    gamma=((2*k40*mu2*psi4*h1^5)/(-psi6*R))^(1/7);

    % Now estimate psi4 given gamma.
    kw4 = deriv(xn/gamma,kernel2); %kernel weights. 
    kw  = [kw4,0,kw4(inc:-1:2)]'; % Apply 'fftshift' to kw.
    z   = real(ifft(fft(c,nfft).*fft(kw))); % convolution.

    psi4Gamma  = sum(c.*z(1:inc))/(n*(n-1)*gamma^5);
  
    % Step 4
    h1 = (STEconstant2/psi4Gamma)^(1/5);
  end;  
  % Kernel other than Gaussian scale bandwidth
  h1  = h1*(STEconstant/STEconstant2)^(1/5);
  

  if count>= maxit
    disp('The obtained value did not converge.')
  end
  h(dim) = h1;
end % for dim loop