function L=spec2lseries(Spec,PP,options,varargin)
%SPEC2LSERIES Spectral simulation of time series in 3D Lagrangian sea 
%
%CALL: L=spec2lseries(Spec,Points,options)
%
%   L      = struct with  n  time series
%   
%   Spec   = a directional frequency spectral density structure in 
%           angular frequency ('w') and directional ('theta') form
%           Alt. in wave number ('k', 'k2') form
%   Points = [x1, ..., xn; y1, ..., yn] array with coordinates of 
%               measurement points
%   options = struct with fields 
%       .Nt    = giving  Nt  time points.  (default length(S)-1=n-1).
%                If Nt>n-1 it is assummed that S.S(k)=0 for all k>n-1
%       .Nu    = giving  Nu  space points along x-axis (default = Nt)
%       .Nv    = giving  Nv  space points along y-axis (default = Nt)
%       .dt    = step in grid (default dt is defined by the Nyquist freq) 
%       .du    = step in grid (default du is defined by the Nyquist freq)
%       .dv    = step in grid (default dv is defined by the Nyquist freq)
%       .lalpha = alpha value for modified Lagrange (default = 0)
%       .iseed  = starting seed number for the random number generator 
%                (default 'shuffle')
%       .plotflag = 0 (no plotting)
%
% Version corresponding to Applied Ocean Research 2009 with respect 
% to .lalpha 

% Tested on Matlab 8.1, 8.6
% Created GL April 2015 based on spec2ldat3D and ldat2lwav3D

tic
% Initialization
if ~isfield(Spec,'g')
  Spec.g=gravity;
end
if ~isfield(Spec,'h')
    Spec.h=Inf;
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

if nargin<3,
    opt=simoptset;
else
    opt=options;
end
if nargin>=4,
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
    Nu=2^nextpow2(length(Sdir.w));
else Nu=opt.Nu;
end
if isempty(opt.Nv), 
    Nv=2^nextpow2(length(Sdir.w));
else Nv=opt.Nv;
end
if isempty(opt.Nt), 
    Nt=2^nextpow2(length(Sdir.w));
else Nt=opt.Nt;
end

if isempty(opt.du)
    du=pi/(Sdir.w(end)^2/Sdir.g)/2;
else
    du=opt.du;
end
if isempty(opt.dv),
    dv=du;
    disp(['dv not specified, set equal to du = ' num2str(dv)])
else dv=opt.dv;
end
if isempty(opt.dt)
    dt=pi/(Sdir.w(end)^2/Sdir.g)/2;
else
    dt=opt.dt;
end
if isempty(opt.lalpha)
    alpha=0;
else
    alpha=opt.lalpha;
end
if nargin<2 || isempty(PP),
    PP = [opt.Nu*opt.du/2;opt.Nv*opt.dv/2];
end
% End initialization

% Checking space/time - spectrum compatibility
if spec.k(end) > pi/du,
    du=pi/spec.k(end);
    Umax=Nu*du;
    disp('Too large u-space step. Aliasing may occur')
    disp(['du changed to ' num2str(du)]) 
    disp(['Max u changed to ' num2str(Umax)])
end
if spec.k2(end) > pi/dv,
    dv=pi/spec.k2(end);
    Vmax=Nv*dv;
    disp('Too large v-space step. Aliasing may occur')
    disp(['dv changed to ' num2str(dv)]) 
    disp(['Max v changed to ' num2str(Vmax)])
end
% Checking space/time - spectrum compativility

nfftu=2^nextpow2(max(Nu,length(Spec.w)));
nfftv=2^nextpow2(max(Nv,length(Spec.w)));

t=(0:Nt-1)'*dt;
u=(0:Nu-1)'*du;  
v=(0:Nv-1)'*dv; 
W=struct('Z',[],'t',t,'u',u,'v',v);
W.Z=zeros(nfftv,nfftu);
X=W; Y=W;

dk1=2*pi/nfftu/du; % necessary wave number lag to satisfy input space lag
dk2=2*pi/nfftv/dv; % necessary wave number lag to satisfy input space lag
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
S.k=[min(dk1*(-nfftu/2), S.k(1)-2*dk1old) S.k(1)-dk1old S.k...
  S.k(end)+dk1old max(dk1*(nfftu/2-1),S.k(end)+2*dk1old)];
S.k2=[min(dk2*(-nfftv/2),S.k2(1)-2*dk2old); S.k2(1)-dk2old;...
  S.k2; S.k2(end)+dk2old; max(dk2*(nfftv/2-1),S.k2(end)+2*dk2old)];

% Interpolate in spectrum to get right frequency grid to satisfy input
S.S=interp2(S.k,S.k2,S.S,dk1*(-nfftu/2:nfftu/2-1),dk2*(-nfftv/2:nfftv/2-1)');
Sdiscr=S.S*dk1*dk2;

S.k=dk1*(-nfftu/2:nfftu/2-1);
S.k2=dk2*(-nfftv/2:nfftv/2-1)';

if sum(Sdiscr(:))<0.95*spec2mom(spec,0)
    disp('WARNING: Too small du and/or du, or too small Nu and/or Nv')
    disp('         Information in the spectrum is lost')
    disp('         May be a good idea to interrupt')
end
toc
% Simulation
rng(opt.iseed);
Z0 = sqrt(Sdiscr).*randn(nfftv,nfftu)+1i*sqrt(Sdiscr).*randn(nfftv,nfftu);
clear Sdiscr
[K,K2]=meshgrid((-nfftu/2:nfftu/2-1)*dk1,(-nfftv/2:nfftv/2-1)*dk2);
Kabs=sqrt(K.^2+K2.^2);
Theta=atan2(K2,K);
Omega = k2w((-nfftu/2:nfftu/2-1)*dk1,(-nfftv/2:nfftv/2-1)*dk2,S.h,S.g);
H1=1i*(tanh(depth*Kabs)).^(-1);
H1(Kabs==0)=0;

H2=alpha*Omega.^(-2);
H2(Omega==0)=0;
H1u=H1.*cos(Theta);
H1v=H1.*sin(Theta);

% Alt 1
%H2u=-H2.*(cos(Theta)).^2.*sign(sin(Theta));
%H2v=-H2.*(sin(Theta)).^2.*sign(sin(Theta)).*sign(cos(Theta));

% Alt 2
%H2u=H2.*(cos(Theta/2)).^2;
%H2v=H2.*(sin(Theta/2)).^2.*sign(sin(Theta/2));

% Alt 3
H2u=H2.*cos(Theta).^2.*abs(cos(Theta));
H2v=H2.*cos(Theta).^2.*sin(Theta).*sign(cos(Theta));
Hu=H1u+H2u;
Hv=H1v+H2v;

Z0u=Z0.*Hu;
Z0v=Z0.*Hv;

Z0=fftshift(Z0);
Z0u=fftshift(Z0u);
Z0v=fftshift(Z0v);
Omega_ = fftshift(Omega); 
clear Omega Theta K K2 Kabs H1u H1v H2u H2v Hu Hv

[mp,np] = size(PP); 
if mp==1 && np==2, 
    PP = PP';
elseif mp>=2,
    PP=PP(1:2,:); % Points shall be defined as an  2 x np array
end
    
L = struct('t',W.t,'Z',[],'P',PP);
L.Z = cell(1,np);
    
[xI,yI]=meshgrid(PP(1,:),PP(2,:));
%hand = fwaitbar(0,'Loop over t values...');
nt = length(t);

for j=1:nt
    waitbar(j/nt);% ,hand)
    W.Z(:,:)=real(fft2(Z0.*exp(-1i*Omega_*t(j))));
    X.Z(:,:)=-real(fft2(Z0u.*exp(-1i*Omega_*t(j))));
    Y.Z(:,:)=-real(fft2(Z0v.*exp(-1i*Omega_*t(j))));
    xx=squeeze(X.Z(1:Nv,1:Nu,:));
    yy=squeeze(Y.Z(1:Nv,1:Nu,:));
    xx=repmat(X.u',Nv,1)+xx;
    yy=repmat(Y.v,1,Nu)+yy;
    ww=squeeze(W.Z(1:Nv,1:Nu,:));
    nxy=numel(ww);
 
    ww=reshape(ww,nxy,1);
    xx=reshape(xx,nxy,1);
    yy=reshape(yy,nxy,1);
%plot(xx,yy,'.')
%axis 'equal'
    Vq=griddata(xx,yy,ww,xI,yI);
    for p=1:np,
        L.Z{p}(j)=Vq(p,p);
    end    
end

toc

clear Z0 Z0u Z0v 


