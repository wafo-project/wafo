function [Y,Mv]=seasim(spec,Nx,Ny,Nt,dx,dy,dt,fftdim,plotflag,use_waitbar)
%SEASIM Spectral simulation of a Gaussian sea, 2D (x,t) or 3D (x,y,t)
%
% CALL:  [Y, Mv]=seasim(spec,Nx,Ny,Nt,dx,dy,dt,fftdim,plotflag)
%
%       Y    = output struct with simulated path
%       Mv   = movie of output (run with movie(Mv))
%       spec = spectrum struct, must be two-dimensional spectrum
%       Nx,Ny,Nt = number of points in grid. Put Ny=0 for 2D simulation
%                  (default 2^7)
%       dx,dy,dt = steps in grid (default dx,dt by the Nyquist freq, dy=dx)
%       fftdim   = 1 or 2 gives one- or two dimensional FFT. (default 2) 
%                  If 1 then FFT over t and looping over x and y.
%                  If 2 then FFT2 over x and y, looping over t. 
%       plotflag = 0,1 or 2. If 0 no plot, if 2 movie (takes a couple of min.) 
%                  if 1 single (first) frame of movie (default 0)
%                   
% The output struct contains:
%          .Z matrix of size [Ny Nx Nt], but with singleton dimensions removed 
%          .x, .t (and .y if Ny>1)  
% For the removal of singleton dimensions, see SQUEEZE.
%  
% The output movie can also be created separately by calling SEAMOVIE.
% Then can also different plot styles be chosen for the movie. The style
% used here is a surf-plot.  
%  
% Limitations: memory demanding! fftdim=2 is in general better.
% When fftdim=1 and NX*NY is large (the limit depends on spectrum, dt etc.) 
% then a slower version with more looping is used. This version can be 
% forced by putting fftdim to something greater than 2, then fftdim will equal
% the number of extra loops. Try this after an interrupt by 'Out of Memory'. 
%
% If the spectrum has a non-empty field .tr, then the transformation is 
% applied to the simulated data, the result is a simulation of a transformed
% Gaussian sea.
%
% Example: % Simulate using demospec, plot the first frame 
%  Nx = 2^7; Ny = 2^7; Nt = 100; dx = 10; dy = 10;dt = 1;
%  S = demospec('dir');
%  plotflag = 1;
%  use_waitbar = 0;
%  Y=seasim(S,Nx,Ny,Nt,dx,dy,dt,2,plotflag, use_waitbar);  
%
%  close all;
%
% See also  seamovie, spec2sdat
  
% Tested on Matlab 5.3 
% revised pab jan 2007
% -Now correctly takes S.phi into account in the simulations when fftdim==1.
% revised pab April 2006
% -Improved the stepwise procedure when fftdim~=2
% -fixed a bug thanks to Simon Streltsov: 
%   Variable Nxyp was mistyped as Nyxp some places. 
% revised pab June 2005
% added fwaitbar and removed disp statements
% revised jr 310301, added some status info, displayed when pg is running
% revised jr 180101, line 77, g -> spec.g
% revised es 200600, removed orient landscape  
% revised es 130600, corrected two errors, 
%                    added more plotting alternatives depending on dimension  
% revised es 060600, added transformation of data  
% By es 23.05.00  


% TODO % Make the simulation take phi into account correctly for fftdim==2
% Initialization
if ~isfield(spec,'g')
  spec.g=gravity;
end
if strcmpi(spec.type(end-2:end),'k2d')
  spec=spec2spec(spec,'dir');
elseif ~strcmpi(spec.type(end-2:end),'dir')
  error('Spectrum must be two-dimensional')
end
if isfield(spec,'f')
  spec = ttspec(spec,'w');
  %spec.w=spec.f*2*pi;
  %spec.S=spec.S/2/pi;
  %spec=rmfield(spec,'f');
end
if nargin<2||isempty(Nx)
  Nx=2^nextpow2(length(spec.w));
end
if nargin<3||isempty(Ny)
  Ny=2^nextpow2(length(spec.w));
end
if Ny==0
  Ny=1;
end
if nargin<4||isempty(Nt)
  Nt=2^nextpow2(length(spec.w));
end
if nargin<5||isempty(dx)
  dx=pi/(spec.w(end)^2/spec.g)/2;
end
if nargin<6||isempty(dy)
  dy=dx;
end
if nargin<7||isempty(dt)
  dt=pi/spec.w(end); %directly by the Nyquist frequency
end
if nargin<8||isempty(fftdim)
  fftdim=1;
end
if nargin<9||isempty(plotflag)
  plotflag=0;
end
if nargin<10||isempty(use_waitbar)
  use_waitbar=1;
end

t=(0:Nt-1)'*dt;
x=(0:Nx-1)'*dx;  
y=(0:Ny-1)'*dy;  

if fftdim~=2
  spec = specinterp(spec,dt,Nt);
  nf = length(spec.w); %number of frequencies
  np = length(spec.theta); % number of directions
  Nxy  = Nx*Ny;
  Nxyp = 2^10;   % limit for Nxy to do in one step
  [Xi,Yi] = meshgrid(x,y);
  %Normalize to get variance back after discretization
  Sdiscr  = (spec.S)*spec.w(end)/nf*diff(spec.theta([1,end]))/np; % size np x nf
  Sdiscr  = (randn(np,nf)+i*randn(np,nf)).*sqrt(Sdiscr);
  [k1,k2] = w2k(spec.w,spec.theta+spec.phi,spec.h,spec.g);
  Sxy = zeros(nf,Nxy);
  if use_waitbar
	h = fwaitbar(0,[],' Integration ... ');
  end
  for ix=1:Nxy,
    Sxy(:,ix) = simpson(0:np-1,Sdiscr.*exp(i*(Xi(ix)*k1+Yi(ix)*k2))).';
	if use_waitbar
		fwaitbar(ix/Nxy,h)
	end
  end 
  if use_waitbar
	close(h)
  end
  clear Sdiscr k1 k2
  nfft=2^nextpow2(2*nf-2);
  
  % size(Z) is Nxy x nfft
  Y.Z=zeros(Nxy,nfft);
  if (Nxy<Nxyp) || fftdim>2
    Sxy = [Sxy; zeros(nfft-(2*nf)+2,Nxy); conj(Sxy(nf-1:-1:2,:))];
    Y.Z = real(ifft(Sxy,nfft,1)).'*nfft/(2*nf-2)*nf;
  else %do it stepwise
    if fftdim>2 % new Nxyp to get fftdim steps
      Nxyp=ceil(Nxy/fftdim);
    end
    nr=ceil(Nxy/Nxyp);
    
    h = fwaitbar(0,[],'  ... ');
    
    ixy  = 1:Nxyp;
    if2   = nfft-nf +2 + 1 : nfft;
    %Sxy1 = [Sxy(:,ixy); zeros(nfft-(2*nf)+2,Nxyp); conj(Sxy(nf-1:-1:2,ixy))];
    Sxy1 = zeros(nfft,Nxyp);
    for j=1:nr-1    
      Sxy1(1:nf,:) = Sxy(:,ixy);
      Sxy1(if2,:)  = conj(Sxy(nf-1:-1:2,ixy));
      Y.Z(ixy,:)   = real(ifft(Sxy1,nfft,1)).'*nfft/(2*nf-2)*nf;
      ixy = ixy + Nxyp      ;
      fwaitbar(j/nr,h)
    end
    ixy1 = 1:rem(Nxy,Nxyp);
    ixy  = ixy(ixy1);
    Sxy1(1:nf,ixy1) = Sxy(:,ixy);
    Sxy1(if2,ixy1)  = conj(Sxy(nf-1:-1:2,ixy));
    Y.Z(ixy,:)      = real(ifft(Sxy1(:,ixy1),nfft,1)).'* nfft/(2*nf-2)*nf;
    clear Sxy1
    close(h)
  end
  clear Sxy
  Y.Z = squeeze(reshape(Y.Z(1:Nxy,1:Nt),Ny,Nx,Nt)); 
else  % if fftdim ...
  nfftx=2^nextpow2(max(Nx,length(spec.w)));
  nffty=2^nextpow2(max(Ny,length(spec.w)));
  if nffty<2
    nffty=2;
  end
  dk1=2*pi/nfftx/dx; % necessary wave number lag to satisfy input space lag
  dk2=2*pi/nffty/dy; % necessary wave number lag to satisfy input space lag
  S=spec2spec(spec,'k2d');
  if S.k(1)>=0
    S.S=[zeros(size(S.S,1),size(S.S,2)-1) S.S];
    S.k=[-S.k(end:-1:2) S.k];
  end
  % add  zeros just above old max-freq, and a zero at new max-freq
  % to get non-NaN interpolation 
  S.S=[zeros(2,size(S.S,1)+4); zeros(size(S.S,1),2) S.S ...
       zeros(size(S.S,1),2);zeros(2,size(S.S,1)+4)];
  dk1old=S.k(2)-S.k(1);
  dk2old=S.k2(2)-S.k2(1);
  S.k=[min(dk1*(-nfftx/2), S.k(1)-2*dk1old) S.k(1)-dk1old S.k S.k(end)+dk1old max(dk1*(nfftx/2-1),S.k(end)+2*dk1old)];
  S.k2=[min(dk2*(-nffty/2),S.k2(1)-2*dk2old); S.k2(1)-dk2old; S.k2; S.k2(end)+dk2old; max(dk2*(nffty/2-1),S.k2(end)+2*dk2old)];
  
  % Interpolate in spectrum to get right frequency grid to satisfy input
  %disp(' Interpolating in spectrum')
  S.S=interp2(S.k,S.k2,S.S,dk1*(-nfftx/2:nfftx/2-1),dk2*(-nffty/2:nffty/2-1)');
  Sdiscr=S.S*dk1*dk2;
  
  if sum(Sdiscr(:))<0.9*spec2mom(spec,0)
    disp('WARNING: Too small dx and/or dy, or too small NX and/or NY')
    disp('         Information in the spectrum is lost')
    disp('         May be a good idea to interrupt')
  end
  
  Sdiscr = fftshift(Sdiscr);
  % Simulation
  Z0 = sqrt(Sdiscr).*randn(nffty,nfftx)+i*sqrt(Sdiscr).*randn(nffty,nfftx);
  W  = k2w((-nfftx/2:nfftx/2-1)*dk1,(-nffty/2:nffty/2-1)*dk2,S.h,S.g);
  W(:,nfftx/2+1:nfftx)=-W(:,nfftx/2+1:nfftx); 
  W   = fftshift(W); 
  Y.Z = zeros(nffty,nfftx,Nt);
  %disp(' Loop over t values')
  h = fwaitbar(0,[],'Loop over t values...');
  nt = length(t);
  for j=1:nt
    Y.Z(:,:,j)=real(fft2(Z0.*exp(-i*W*t(j))));
    % Old call pab feb2007
   %  Y.Z(:,:,j)=real(fft2(Z0.*exp(i*W*t(j))));
     fwaitbar(j/nt,h)
  end
  close(h)
  Y.Z=squeeze(Y.Z(1:Ny,1:Nx,:));
  clear Z0 W Sdiscr
  
end % if fftdim ...

if Nx>1
  Y.x=x;
end
if Ny>1
  Y.y=y;
end
if Nt>1
  Y.t=t;
end

% Transformation of data, if given
if isfield(spec,'tr') && ~isempty(spec.tr)
  disp('   Transforming data.')
  G=fliplr(spec.tr); % the inverse of the transformation
  Y.Z(:)=tranproc(Y.Z(:),G); 
end
  
% Plotting if requested
Mv=[];
if plotflag==2
  Mv=seamovie(Y,1);
elseif plotflag==1
  figure(gcf)
  if min(Nx,Ny)>2
    colormap('winter')
    surfl(Y.x,Y.y,Y.Z(:,:,1),[-30, 45]);
    shading interp;
    view(-37.5,20)
    axis([Y.x(1) Y.x(end) Y.y(1) Y.y(end) 7*min(Y.Z(:)) 7*max(Y.Z(:))])
    set(gca,'xtick',[]);
    set(gca,'ytick',[]);
    axis('square');
    axis('off');
  elseif Nt>1
    if Nx>1
      contour(Y.t,Y.x,Y.Z,[0 0])
    else
      contour(Y.t,Y.y,Y.Z,[0 0])
    end      
    xlabel('[s]')
    ylabel('[m]')
  else
    if Nx>1
      plot(Y.x,Y.Z)
    elseif Ny>1
      plot(Y.y,Y.Z)
    end
    xlabel('[m]')
  end
end
