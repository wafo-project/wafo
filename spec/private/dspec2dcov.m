function R=dspec2dcov(S,nr,Nt,rate,Nx,Ny,dx,dy)
%DSPEC2DCOV Return covariance function given a directional spectrum 
%            (used in spec2cov.)
%
%  CALL: R=dspec2dcov(S,nr,Nt,rate,Nx,Ny,dx,dy)
%
%   DCOV = covariance structure (see datastructures)
%      S = spectrum struct, type 'dir' or 'freq' (see datastructures)
%     nr = 0,1,2, number of derivatives 
%     Nt = number of time lags 
%   rate = 1,2,4,8...2^r, interpolation rate for R
%  Nx,Ny = number of space lags, (default 10 resp. 0)
%  dx,dy = space lag step  (default such that a wave length is covered)
%

% NB! requires simpson

% Tested on: Matlab 5.3
% History:
% revised ir 05.04.01, theta modified by a rotation angle phi.
% revised by es 24.05.00, fixing shape of output from w2k, since w2k changed  
% revised by es 10.11.1999, nr>2  
% revised by es 13.10.1999
% by pab 16.04.1999

% TODO % Needs more work  
%initialize
if strcmpi(S.type(end-2:end),'k2d')
  S=spec2spec(S,'dir');
elseif strcmpi(S.type(end-2:end),'k1d')
  S=spec2spec(S,'freq');
end  
if nargin<5|isempty(Nx)
  Nx=10;
end
if nargin<6 |isempty(Ny)
  Ny=0;
end

g=gravity;
if isfield(S,'g')
  g=S.g;
end
m=spec2mom(S,2,'xyt',1);
if length(m)==2 % result when S of type 'freq'
  m=spec2mom(S,4,'t',1);  %m=[m0 m002 m004]
  m=[m(1) m(3)/g^2 m(3)/g^2]; %m=[m0 m200 m020];
end
if nargin<7|isempty(dx)
  dx=4*pi*sqrt(m(1)/m(2))/Nx; % divide the average wave into Nx/2 
end

if nargin<8|isempty(dy),
  if Ny>0
    dy=dx;
    if isinf(dy) % in case of Nx==0
      dy=4*pi*sqrt(m(1)/m(3))/Ny;
      dx=dy;
    end
    if isinf(dy) % in case of both Nx and Ny ==0
      dx=1;dy=1; % just to get some value
    end
  else
    dy=1;
    if isinf(dx)
      dx=dy;
    end
  end
end

vari='t'; 
if Nx>0
  vari=strcat(vari,'x');
end
if Ny>0
  vari=strcat(vari,'y');
end
R=createcov(nr,vari);
R.tr=S.tr;
R.h=S.h;
R.norm=S.norm;
if isfield(S,'note')
  R.note=S.note;
end

warning off
comment=1;
if comment,
  disp('   Calculating covariances.')
end
comment=0;
if isfield(S,'f')
  S.w=S.f*2*pi;
  S.S=S.S/2/pi;
end
if ~isfield(S,'g')
  S.g=gravity;
end
if ~isfield(S,'theta')
  S.theta=0;
  S.S=S.S(:)';
end

nf= length(S.w); %number of frequencies
np= length(S.theta); % number of angles

if Nx==Ny & dx==dy,
  symmetry=1;%only able to exploit the symmetry S(x,y,w)=conj(S(y,x,w)) 
         %if Nx=Ny,dx==dy
else
  symmetry=0; 
end
if symmetry,  % exploit the symmetry S(x,y,w)=conj(S(y,x,w)) 
  xi=((-Nx):0)'*dx; 
else     
  xi=((-Nx):Nx)'*dx;  
end
yi=((-Ny):Ny)'*dy;

[Xi Yi]=meshgrid(xi,yi);
NXYi=length(Xi(:));

if ~isfield(S,'phi') | isempty(S.phi), S.phi=0; end

theta=S.theta(:)-S.phi;     %make sure it is a column
S.theta=S.theta-S.phi;
Yi=Yi(:);Xi=Xi(:);

if comment,
  disp('   ... initiate') 
end

DS2=(S.S)*S.w(end); % Normalizing
Sxy=zeros(nf,NXYi);
Sxyx=Sxy;Sxyxx=Sxy;Sxyxxx=Sxy;Sxyxxxx=Sxy;
if np>1
  [kCtheta,kStheta]=w2k(S.w,S.theta,S.h,S.g); %size np x nf
  for ix=1:NXYi,
    Stheta=DS2.*exp(i*(Xi(ix).*kCtheta +Yi(ix).*kStheta));
    Sxy(:,ix)=simpson(theta,Stheta).';
    if nr>0
      Sxyx(:,ix)=simpson(theta,i*kCtheta.*Stheta).';
      if nr>1
	Sxyxx(:,ix)=simpson(theta,-kCtheta.^2.*Stheta).';
	if nr>2
	  Sxyxxx(:,ix)=simpson(theta,(-i*kCtheta.^3).*Stheta).';
	  if nr>3
	    Sxyxxxx(:,ix)=simpson(theta,(kCtheta.^4).*Stheta).';
	  end
	end
      end
    end  
  end
  clear kStheta
else %theta=0
  [kCtheta]=w2k(S.w,S.theta,S.h,S.g);
  kCtheta=kCtheta(:)'; %size 1 x nf
  Stheta=(ones(size(Xi))*DS2).*exp(i*(Xi*kCtheta));
  Sxy=Stheta.';
  if nr>0
    Sxyx=(i*(ones(NXYi,1)*kCtheta).*Stheta).';
    if nr>1
      Sxyxx=(-(ones(NXYi,1)*kCtheta.^2).*Stheta).';
      if nr>2
	Sxyxxx=(-i*(ones(NXYi,1)*kCtheta.^3).*Stheta).';
	if nr>3
	  Sxyxxxx=((ones(NXYi,1)*kCtheta.^4).*Stheta).';
	end
      end
    end
  end  
end

clear DS2 kCtheta  Xi Yi theta Stheta 
if symmetry,   % exploit the symmetry in calculating Sxy
  Sxy=[Sxy conj(Sxy(1:nf,end-2*Ny-1:-1:1))];
  if nr>0
    Sxyx=[Sxyx -conj(Sxyx(1:nf,end-2*Ny-1:-1:1))];
    if nr>1
      Sxyxx=[Sxyxx conj(Sxyxx(1:nf,end-2*Ny-1:-1:1))];
      if nr>2
	Sxyxxx=[Sxyxxx -conj(Sxyxxx(1:nf,end-2*Ny-1:-1:1))];
	if nr>3
	  Sxyxxxx=[Sxyxxxx conj(Sxyxxxx(1:nf,end-2*Ny-1:-1:1))];
	end
      end
    end
  end
  xi=((-Nx):Nx)'*dx; 
  NXYi=length(xi)*length(yi);
end
% size(Sxy)=[nf NXYi]

if comment,
  disp('   ... FFT')
end
nfft=rate*2^nextpow2(2*nf-2);
rat=nfft/(2*nf-2);
Sxy=[Sxy; zeros(nfft-(2*nf)+2,NXYi); conj(Sxy(nf-1:-1:2,:))];
cov=real(ifft(Sxy,nfft,1)).'*rat; % size(cov) is NXYi  x nfft
% size(cov) will now become 2*Ny+1 x 2*Nx+1 x Nt+1
R.R=rshape(cov,Nt,Nx,Ny);         % for rshape: see end of this file

if Nx>0, R.x=xi; end, if Ny>0, R.y=yi; end, R.t=(0:Nt)'*pi/(S.w(nf))/rat;

if nr>0
  w=S.w(:);
  w=[w ; zeros(nfft-2*nf+2,1) ;-w(nf-1:-1:2) ]*ones(1,NXYi);
  Sxy=i*w.*Sxy;
  cov=real(ifft(Sxy,nfft,1)).'*rat;
  R.Rt=rshape(cov,Nt,Nx,Ny);                               %r_t
  if Nx>0
    Sxyx=[Sxyx; zeros(nfft-(2*nf)+2,NXYi); conj(Sxyx(nf-1:-1:2,:))];
    cov=real(ifft(Sxyx,nfft,1)).'*rat;
    R.Rx=rshape(cov,Nt,Nx,Ny);                             %r_x
  end
  if nr>1
    Sxy=i*w.*Sxy;
    cov=real(ifft(Sxy,nfft,1)).'*rat;
    R.Rtt=rshape(cov,Nt,Nx,Ny);                            %r_tt
    if Nx>0
      Sxyxx=[Sxyxx; zeros(nfft-(2*nf)+2,NXYi); conj(Sxyxx(nf-1:-1:2,:))];
      cov=real(ifft(Sxyxx,nfft,1)).'*rat;
      R.Rxx=rshape(cov,Nt,Nx,Ny);                          %r_xx
      cov=real(ifft(i*w.*Sxyx,nfft,1)).'*rat;
      R.Rxt=rshape(cov,Nt,Nx,Ny);                          %r_xt
    end
    if nr>2
      Sxy=i*w.*Sxy;
      cov=real(ifft(Sxy,nfft,1)).'*rat;
      R.Rttt=rshape(cov,Nt,Nx,Ny);                         %r_ttt
      if Nx>0
	Sxyxxx=[Sxyxxx; zeros(nfft-(2*nf)+2,NXYi); conj(Sxyxxx(nf-1:-1:2,:))];
	cov=real(ifft(Sxyxxx,nfft,1)).'*rat;
	R.Rxxx=rshape(cov,Nt,Nx,Ny);                       %r_xxx
	cov=real(ifft(i*w.*Sxyxx,nfft,1)).'*rat;
	R.Rxxt=rshape(cov,Nt,Nx,Ny);                       %r_xxt
	cov=real(ifft(-w.^2.*Sxyx,nfft,1)).'*rat;
	R.Rxtt=rshape(cov,Nt,Nx,Ny);                       %r_xtt	
      end
       if nr>3
	Sxy=i*w.*Sxy;
	cov=real(ifft(Sxy,nfft,1)).'*rat;
	R.Rtttt=rshape(cov,Nt,Nx,Ny);                      %r_tttt
	if Nx>0
	  cov=real(ifft([Sxyxxxx; zeros(nfft-(2*nf)+2,NXYi);...
		   conj(Sxyxxxx(nf-1:-1:2,:))],nfft,1)).'*rat;  
	  R.Rxxxx=rshape(cov,Nt,Nx,Ny);                    %r_xxxx
	  cov=real(ifft(-i*w.^3.*Sxyx,nfft,1)).'*rat;
	  R.Rxttt=rshape(cov,Nt,Nx,Ny);                    %r_xttt
	  cov=real(ifft(-w.^2.*Sxyxx,nfft,1)).'*rat;
	  R.Rxxtt=rshape(cov,Nt,Nx,Ny);                    %r_xxtt
	  cov=real(ifft(i*w.*Sxyxxx,nfft,1)).'*rat;
	  R.Rxxxt=rshape(cov,Nt,Nx,Ny);                    %r_xxxt
	end
      end
    end
  end
end
return %end dspec2dcov


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cov=rshape(cov,Nt,Nx,Ny)
% RSHAPE Change size of matrix, special for covariance function routine
%  
%  CALL:   Cout=rshape(Cin,Nt,Nx,Ny)
%
%     Cout = matrix, size [2*Nx+1 2*Ny+1 Nt+1], singleton dimensions removed
%     Cin  = matrix, size=[m n] where m>=(2*Nx+1)*(2*Ny+1), n>=Nt+1
%           rows below (2*Nx+1)*(2*Ny+1) and columns beyond Nt+1 are neglected
%     Nt,Nx,Ny = integers defining  size of output matrix  

if Nx>0
  if Ny>0
    cov=reshape(cov(1:(2*Nx+1)*(2*Ny+1),1:Nt+1),2*Ny+1,2*Nx+1,Nt+1);
  else
    cov=reshape(cov(1:2*Nx+1,1:Nt+1),2*Nx+1,Nt+1);
  end
elseif Ny>0
  cov=reshape(cov(1:2*Ny+1,1:Nt+1),2*Ny+1,Nt+1);
else
  cov=reshape(cov(1:Nt+1),Nt+1);
end
return  %end rshape
