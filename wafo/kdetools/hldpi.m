function h=hldpi(A,kernel,L,inc)
%HLDPI L-stage Direct Plug-In estimate of smoothing parameter.
%
% CALL: hs = hldpi(data,kernel,L)
%
%       hs = one dimensional value for smoothing parameter
%            given the data and kernel.  size 1 x D
%   data   = data matrix, size N x D (D = # dimensions )
%   kernel = 'epanechnikov'  - Epanechnikov kernel.
%            'biweight'      - Bi-weight kernel.
%            'triweight'     - Tri-weight kernel.
%            'triangluar'    - Triangular kernel.
%            'gaussian'      - Gaussian kernel
%            'rectangular'   - Rectanguler kernel.
%            'laplace'       - Laplace kernel.
%            'logistic'      - Logistic kernel.
%        L = 0,1,2,3,...   (default 2)
%
%  Note that only the first 4 letters of the kernel name is needed.
%
%  Example:
%  % data = rndnorm(0, 1,5,2);
%  x = [-0.0233845632050972   0.9070186193622006;...
%        0.6529594866766634   1.3689145060433903;...
%        0.4477857310723146  -0.6311953712037597;...
%       -1.9256785038579962   0.5886257667993168;...
%       -0.5290011931824666  -0.3602090880229930];
%  assert(hldpi(x,'gaus',1), [0.740583842868078, 0.672480449056617], 1e-10);
%  assert(hldpi(x,'epan',1), [1.63950773944363, 1.48874014937057], 1e-10);
%  assert(hldpi(x,'biwe',1), [1.94226653858676, 1.76365753403174], 1e-10);
%  assert(hldpi(x,'triw',1), [2.20553713027085, 2.00271800966064], 1e-10);
%  assert(hldpi(x,'tria',1), [1.80109851299778, 1.63547118733416], 1e-10);
%  assert(hldpi(x,'rect',1), [1.28865814206224, 1.17015435105684], 1e-10);
%  assert(hldpi(x,'lapl',1), [0.547862055772363, 0.497481176283924], 1e-10);
%  assert(hldpi(x,'logi',1), [0.413993548184052, 0.375923090775485], 1e-10);
%
% See also  hste, hbcv, hboot, hos, hlscv, hscv, hstt, kde, kdefun

%  Wand,M.P. and Jones, M.C. (1995)
% 'Kernel smoothing'
%  Chapman and Hall, pp 67--74


%tested on: matlab 5.2
% revised pab Agu 2005
% -fixed a bug for kernels other than Gaussian
% -made it more general: L>3 is now possible
% revised pab nov2004
% - kernel2 = 'gaus'
% -added nargchk
% revised pab Dec2003
% revised pab 16.10.1999
%  updated to matlab 5.x changed name from hdpi -> hldpi
%  changed the gridding -> hldpi is now faster
% taken from kdetools             Christian C. Beardah 1994


error(nargchk(1,4,nargin))
if nargin<2||isempty(kernel)
  kernel='gauss';
end
if nargin<3||isempty(L),
  L=2;
else
  L=abs(L);
end;
if nargin<4||isempty(inc)
  inc=128;
end

nfft = inc*2;
[n, d] = size(A);
if (n==1) && (d>1),
  A=A.';
  n=d;
  d=1;
end

xmin = min(A);    % Find the minimum value of A.
xmax = max(A);    % Find the maximum value of A.
xrange = xmax-xmin; % Find the range of A.

sigmaA = std(A);
iqr = abs(diff(qlevels2(A,[75 25]),1,1+(d==1))); % interquartile range
k = find(iqr>0);
if any(k)
  sigmaA(k) = min(sigmaA(k), iqr(k)/1.349);
end

% R= int(mkernel(x)^2)
% mu2= int(x^2*mkernel(x))
[mu2, R] = kernelstats(kernel);

STEconstant = R /(mu2^(2)*n);

% xa holds the x 'axis' vector, defining a grid of x values where
% the k.d. function will be evaluated and plotted.

ax1 = xmin-xrange/8;
bx1 = xmax+xrange/8;

kernel2 = 'gaus';
[mu2] = kernelstats(kernel2); % Bug fix: Do not delete this: pab aug2005

h = zeros(1,d);
for dim=1:d
  s = sigmaA(dim);
  ax = ax1(dim);
  bx = bx1(dim);

  xa = linspace(ax,bx,inc).';
  xn = linspace(0,bx-ax,inc);

  c = gridcount(A(:,dim),xa);

  %c=zeros(inc,1);
  %deltax =(bx-ax)/(inc-1);
  %xn     = (0:(inc-1))*deltax;
  %binx   = floor((A(:,dim)-ax)/deltax)+1;


  % Obtain the grid counts.
  %c = full(sparse(binx,1,(xa(binx+1)-A(:,dim)),inc,1)+...
  %         sparse(binx+1,1,(A(:,dim)-xa(binx)),inc,1))/deltax;

  r   = 2*L+4;
  rd2 = L+2;

  % Eq. 3.7 in Wand and Jones (1995)
  PSI_r =  (-1)^(rd2)*prod(rd2+1:r)/(sqrt(pi)*(2*s)^(r+1));
  PSI   = PSI_r;
  if L>0
    % High order derivatives of the Gaussian kernel
    [Kd{1:L}] = deriv(0,kernel2);
    
    % L-stage iterations to estimate PSI_4
    for ix = L:-1:1
      gi = (-2*Kd{ix}/(mu2*PSI*n))^(1/(2*ix+5));

      % Obtain the kernel weights.
      [KW0{1:ix}] = deriv(xn/gi,kernel2);
      kw=[KW0{ix},0,KW0{ix}(inc:-1:2)].'; % Apply 'fftshift' to kw.

      % Perform the convolution.
      z=real(ifft(fft(c,nfft).*fft(kw)));

      PSI = sum(c.*z(1:inc))/(n^2*gi^(2*ix+3));
    end
  end
  h(dim)=(STEconstant/PSI)^(1/5);
end