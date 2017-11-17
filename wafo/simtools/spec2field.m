function W = spec2field(Spec,options,varargin)
%SPEC2FIELD Spectral simulation of space-time Gaussian field 
% 
%CALL: W = spec2field(Spec,options)
%
%   W     = Gaussian wave structure W with fields 
%       .Z      = matrix of size [N N N] with singleton dim removed
%       .x (.y) = space coordinates along x-axis (y-axis if Nv>1)
%       .t      = time coordinates
%
%   Spec  = a directional frequency spectral density structure in 
%           angular frequency ('w') and directional ('theta') form
%           Alt. in wave number ('k', 'k2') form
%           If Spec.tr exists, W is a transformed Gaussian field
%
%   options = struct with fields 
%       .Nt    = giving  Nt  time points.  (default length(S)-1=n-1).
%                If Nt>n-1 it is assummed that S.S(k)=0 for all k>n-1
%       .Nu    = giving  Nu  space points along x-axis (defult = Nt)
%       .Nv    = giving  Nv  space points along y-axis (defult = Nt)
%       .dt    = step in grid (default defined by the Nyquist freq) 
%       .du    = step in grid (default defined by the Nyquist freq)
%       .dv    = step in grid (default defined by the Nyquist freq)
%       .iseed  = method or starting seed number for random number generator 
%                (default 'shuffle')
%       .plotflag = 0 (no plotting)
%
% Based on WAFO-routines seasim and spec2ldat3D
%
% Example: Generate 20 seconds of Gaussian wavefield
%           over a rectangle with length = 1024 m and
%           width = 512 m. Make movie with seasim.
%
%   S = demospec('dir')
%   opt = simoptset('Nt',20,'dt',1,'Nu',1024,'du',1,'Nv',512,'dv',1)
%   W = spec2field(S,opt)
%   M = seamovie(W,1);
%
% See also: spec2ldat,spec2sdat,cov2sdat, gaus2dat

% Tested on Matlab 8.1, 8.6, 9.2
% History
% spec2field created by GL 2017 from spec2ldat by removing X, Y variation
% Modified March 2015 by GL to give correct slope and curvature 
% Modified by GL 2014
% spec2ldat created by GL and FL 2010 (for Stochastic models article)

tic
% Initialization
if ~isfield(Spec,'g')
  Spec.g=gravity;
end
if ~isfield(Spec,'h')
    h=Inf;
    
else
    depth=Spec.h;
    Spec.h=Inf;
end
if strcmpi(Spec.type(end-2:end),'dir')
  spec=spec2spec(Spec,'k2d');
  Sdir=Spec;
elseif strcmpi(Spec.type(end-2:end),'k2d')
  spec=Spec;
  Sdir=spec2spec(Spec,'dir');
else
    error('Spectrum must be two-dimensional')
end
if isfield(Spec,'f')
  Spec = ttspec(Spec,'w');
end

if nargin<2,
    opt=simoptset;
else
    opt=options;
end
if nargin>=3,
    opt=simoptset(opt,varargin{:});
end

if isfield(opt,'iseed') 
    iseed=opt.iseed;
else 
    iseed=[];
end

if verLessThan('matlab','7.12'),
    if isempty(iseed) || strcmp(iseed,'shuffle'),
        rand('seed',double(int32(sum(100*clock))));
    elseif isnumeric(iseed)
        rand('seed',int32(iseed));
    else
        rand('seed',double(int32(sum(100*clock))));
    end
else
    if isempty(iseed)
        iseed='shuffle';
    end
    rng(iseed); 
end

if isempty(opt.Nu), 
    Nx=2^nextpow2(length(Sdir.w));
else Nx=opt.Nu;
end
if isempty(opt.Nv), 
    Ny=2^nextpow2(length(Sdir.w));
else Ny=opt.Nv;
end
if isempty(opt.Nt), 
    Nt=2^nextpow2(length(Sdir.w));
else Nt=opt.Nt;
end

if isempty(opt.du)
    dx=pi/(Sdir.w(end)^2/Sdir.g)/2;
else
    dx=opt.du;
end
if isempty(opt.dv),
    dy=dx;
    disp(['dy not specified, set equal to dx = ' num2str(dy)])
else dy=opt.dv;
end
if isempty(opt.dt)
    dt=pi/(Sdir.w(end)^2/Sdir.g)/2;
else
    dt=opt.dt;
end

% End initialization

t=(0:Nt-1)'*dt;
x=(0:Nx-1)'*dx;  
y=(0:Ny-1)'*dy; 
W=[]; 

% Checking space/time - spectrum compatibility
if spec.k(end) > pi/dx,
    dx=pi/spec.k(end);
    Umax=Nx*dx;
    disp('Too large x-space step. Aliasing may occur')
    disp(['dx changed to ' num2str(dx)]) 
    disp(['Max x changed to ' num2str(Umax)])
end
if spec.k2(end) > pi/dy,
    dy=pi/spec.k2(end);
    Vmax=Ny*dy;
    disp('Too large y-space step. Aliasing may occur')
    disp(['dy changed to ' num2str(dy)]) 
    disp(['Max y changed to ' num2str(Vmax)])
end
% Checking space/time - spectrum compativility

nfftx=2^nextpow2(max(Nx,length(Spec.w)));
nffty=2^nextpow2(max(Ny,length(Spec.w)));

dk1=2*pi/nfftx/dx; % necessary wave number lag to satisfy space lag
dk2=2*pi/nffty/dy; % necessary wave number lag to satisfy space lag
S=spec; % spec2spec(spec,'k2d');

if S.k(1)>=0
    S.S=[zeros(size(S.S,1),size(S.S,2)-1) S.S];
    S.k=[-S.k(end:-1:2) S.k];
end
% add  zeros just above old max-freq, and a zero at new max-freq
% to get non-NaN interpolation 
S.S=[zeros(2,size(S.S,2)+4); zeros(size(S.S,1),2) S.S ...
   zeros(size(S.S,1),2);zeros(2,size(S.S,2)+4)];
dk1old=S.k(2)-S.k(1);
dk2old=S.k2(2)-S.k2(1);
S.k=[min(dk1*(-nfftx/2), S.k(1)-2*dk1old) S.k(1)-dk1old S.k...
  S.k(end)+dk1old max(dk1*(nfftx/2-1),S.k(end)+2*dk1old)];
S.k2=[min(dk2*(-nffty/2),S.k2(1)-2*dk2old); S.k2(1)-dk2old;...
  S.k2; S.k2(end)+dk2old; max(dk2*(nffty/2-1),S.k2(end)+2*dk2old)];

% Interpolate in spectrum to get right frequency grid to satisfy input
S.S=interp2(S.k,S.k2,S.S,dk1*(-nfftx/2:nfftx/2-1),dk2*(-nffty/2:nffty/2-1)');
Sdiscr=S.S*dk1*dk2;

S.k=dk1*(-nfftx/2:nfftx/2-1);
S.k2=dk2*(-nffty/2:nffty/2-1)';

if sum(Sdiscr(:))<0.95*spec2mom(spec,0)
[sum(Sdiscr(:)) spec2mom(spec,0)]
    disp('WARNING: Too small dx and/or dy, or too small Nx and/or Ny')
    disp('         Information in the spectrum is lost')
    disp('         May be a good idea to interrupt')
end

% Simulation
rng(opt.iseed);
Z0 = sqrt(Sdiscr).*randn(nffty,nfftx)+1i*sqrt(Sdiscr).*randn(nffty,nfftx);
clear Sdiscr
[K,K2]=meshgrid((-nfftx/2:nfftx/2-1)*dk1,(-nffty/2:nffty/2-1)*dk2);
Kabs=sqrt(K.^2+K2.^2);
Theta=atan2(K2,K);
Omega = k2w((-nfftx/2:nfftx/2-1)*dk1,(-nffty/2:nffty/2-1)*dk2,S.h,S.g);

Z0=fftshift(Z0);
Omega_ = fftshift(Omega); 
clear Omega Theta K K2 Kabs H1u H1v H2u H2v Hu Hv
W.Z = zeros(nffty,nfftx,Nt);
nt = length(t);
for j=1:nt
    W.Z(:,:,j)=real(fft2(Z0.*exp(-1i*Omega_*t(j))));
    waitbar(j/nt);%,hand)
end

W.Z=squeeze(W.Z(1:Ny,1:Nx,:));
clear Z0 

if Nx>1
  W.y=x;
end
if Ny>1
  W.x=y;
end
  W.t=t;
toc
