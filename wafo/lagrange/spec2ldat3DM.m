function [W,X,Y,W2,X2,Y2] = spec2ldat3DM(Spec,order,options,varargin)
%SPEC2LDAT3DM Particle trajectory simulation according to Marc Prevosto
%             2D or 3D, first or second order Lagrange waves
%
%CALL: [W,X,Y,W2,X2,Y2] = spec2ldat3DM(spec,order,options) 
%
%   I   - Spec     : one-sided spectral structure with fields 
%         'S'      : spectral density values
%         'freq'   : 1D frequency spectrum over 'w'/'f'
%         'D'      : spreading structure from  D = spreading
%   I   - order    : 1, first order (default), or 2, second order, waves 
%   I   - options  : struct with fields 
%       .Nt     = minimum number of time points
%       .Nu/.Nv = giving  Nu/Nv  space points along x/y-axis (default = 100)
%       .dt     = approximate time-step  
%       .du/.dv = step in grid (default du=dv=1)
%       .iseed  = starting seed number for the random number generator 
%                (default 'shuffle')
%       h       = water depth (field in Spec, default = Inf)
%
%       z0      = secret parameter = particle depth between 0 and -h
%                 can be set inside program, default = 0
%       trf,tls = secret parameters for to truncate spectrum to avoid 
%                 instability and increase speed, 
%                 can be set inside program, default = 0.01,0.001
%
%   O   - W,X,Y    = first order components
%   O   - W2,X2,Y2 = second order components

% Tested on Matlab 2012, 32 bit
% Tested on Matlab 8.1, 8.6
% History 
%  Author : Marc Prevosto - 25 August 1998
%  Adapted to WAFO convention by Georg Lindgren 2015
%  and expanded to full 3D fields 

tic
if nargin<2,
    order = 1;
else 
    order = 1 + (order>1);
end
    
if nargin<3,
    opt=simoptset('Nu',100,'du',1,'Nv',100,'dv',1,'iseed','shuffle');
else
    opt=options;
end
if nargin>=4,
    opt=simoptset(options,varargin{:});
end

if isfield(opt,'iseed') 
    iseed=opt.iseed;
else 
    iseed=[];
end
if isempty(iseed)
    iseed='shuffle';
end

rng(iseed); 
Nu=opt.Nu;
U=(0:Nu-1)*opt.du;
if ~isempty(opt.Nv) && opt.Nv>0,
    Nv=opt.Nv;
    V=(0:Nv-1)*opt.dv;
else V=0; Nv=1;
end

S=Spec; h=S.h;
alpha=opt.lalpha;

if isfield(S,'type'),
    spectype=S.type;
else
    disp('Unknown spec-type')
end

Nt=opt.Nt; dt=opt.dt;
%Nt=4*(length(S.S)+1);
%Nt=8*(length(S.S)+1);

nf=1/4*2^nextpow2(4*ceil(Nt/4));
Df=1/Nt/dt;

if isfield(S,'f'),
    xf = S.f';
    sp = S.S';
elseif isfield(S,'w'),
    xf = S.w'/2/pi;
    sp = S.S'*2*pi;
else 
    error('Unknown freq-type')
end

if xf(1)==0,
    xf=xf(2:end);
    sp=sp(2:end);
end

if isfield(S,'D'),
    D = S.D;
else
    D.S = 1;
    D.theta = 0;
end

theta = D.theta;
ndir = length(theta);
if ndir == 1,
    Nv = 1;
    V=0;
    x_theta = 1; %sqrt(D.S'*2*pi/ndir);
else
    x_theta = sqrt(D.S'*2*pi/ndir);
end

inc = diff(theta) ;
if abs(min(inc)-max(inc))>1.e6*eps,
    error('Increment d''angle non constant') ;
end

x_f = (1:nf)*Df;
sp = interp1(xf,sp,x_f,[],0); 
[sum(sp)*Df spec2mom(S)];
Ntt=4*(length(sp)+1);
[x_f(end) 1/2/dt];

W=struct('Z',zeros(Nv,Nu,Ntt),'u',U','v',V','t',[]);
X=W; Y=W; 
if order == 2,
    W2 = W; X2 = W; Y2 = W;
    X2.drift = W2.Z; Y2.drift = W2.Z;
end

t=[]; 
z=0; % Particle depth
nfreq = length(sp) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Truncation of spectrum
% to avoid non convergence of second order and decrease computing time
trf = 0.01 ; % frequency truncation
tls = 0.001; % low value truncation

if trf
    % spectrum high-low frequency truncation
    cs = cumsum(sp) ; cs = cs/cs(end) ;
    il = find(cs<trf, 1, 'last' ) ; ih = find(cs>1-trf, 1 ) ;
    sp(1:il+1) = 0*sp(1:il+1) ;
    sp(ih-1:end) = 0*sp(ih-1:end) ;
end
 
if tls
    % spectrum low value truncation
    if ndir>1
        stheta = size(theta) ; 
        if min(stheta)==1, 
            x_theta = x_theta(ones(1,nfreq),:)' ;
        end
        spd = x_theta.^2.*sp(ones(1,ndir),:) ; spd_t = sum(sum(spd)) ;
        [spds,ind] = sort(spd(:)) ; 
        spds = cumsum(spds) ; il = find(spds<spd_t*tls, 1, 'last' ) ;
        x_theta(ind(1:il)) = 0 ; 
    else
        cs = cumsum(sp) ;
        [spds,ind] = sort(sp) ; spds = cumsum(spds) ; 
        il = find(spds<cs*tls, 1, 'last' ) ;
        sp(ind(1:il)) = 0 ; 
    end
end
% End of truncation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate random elements and set some more parameters
df = x_f(1) ;
x_e = sqrt(sp * (df/4)); %(sp*df)/2 pour passage du one-sided au spectre
                         % et puis /2 pour la variance des composantes
              
x_e = x_e(ones(1,ndir),:) ;

cx_theta=x_theta; 

if size(cx_theta,2)==nfreq,
    x_e = cx_theta.*x_e ;
else
    x_e = cx_theta(ones(1,nfreq),:)'.*x_e ;
end

dB = randn(ndir,nfreq)+randn(ndir,nfreq)*1i ;

x_w = (2*pi)*x_f ;
x_k = 2*pi./disper2(2*pi./x_w,h)' ;
npts = 4*(nfreq+1) ;
%npts = 8*(nfreq+1) ;
cosd = zeros(ndir) ;
t = (0:npts-1)/(npts*df) ;

% End Generate random elements and parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate fields

%WB=waitbar(0,'Looping over  u');

[x_eta2r,x_eta2i,coshkz,sinhkz,coshk,sinhk,expk,cost,sint] =...
    partkern(x_w, x_k, x_e*0, theta, cosd, h, [], z) ;

x_temp = -1i*(expk.*(coshkz./sinhk)); % MP notation
H1 = x_temp; %-i pour tenir compte du signe dans la ifft
H2 = -alpha*x_w.^(-2);
H2 = H2'; 
% Set directional effects according to Stochastic Models:

Tu = cos(theta).*abs(cos(theta));
Tv = cos(theta).^2.*sign(cos(theta));
H2u = H2*Tu';
H2v = H2*Tv';
% Define full directional transfer functions 
Hu = repmat(H1,1,ndir)+H2u;
Hv = repmat(H1,1,ndir)+H2v;    

for ku=1:Nu,
    waitbar(ku/Nu);% ,WB)
    for kv=1:Nv,
        facuv=exp(-1i*(cos(theta)*x_k*U(ku)+sin(theta)*x_k*V(kv)));
        dBuv=facuv.*dB;
        
        % Generate random field elements
        if order==2,
            [x_eta2r,x_eta2i,coshkz,sinhkz,coshk,sinhk,expk,cost,sint] =...
                partkern(x_w, x_k, x_e, theta, cosd, h, dBuv, z) ;
        end  
        
        dBxe = (dBuv.*x_e).' ;

        % First order fields
        xyz1 = zeros(nfreq,3) ;
        xyz1(:,1) = (dBxe.*Hu)*cost ;
        xyz1(:,2) = (dBxe.*Hv)*sint ;
        xyz1(:,3) = (dBxe*ones(ndir,1)).*(expk.*(sinhkz./sinhk)) ;

        xyz1 = [zeros(1,3);xyz1;zeros(2*nfreq+3,3);conj(xyz1(nfreq:-1:1,:))] ;
%        xyz1 = [zeros(1,3);xyz1;zeros(6*nfreq+7,3);conj(xyz1(nfreq:-1:1,:))] ;
        xyz1 = real(ifft(xyz1))*npts ; 
  
        W.Z(kv,ku,:) = flipud(xyz1(:,3));       
        X.Z(kv,ku,:) = flipud(xyz1(:,1));
        Y.Z(kv,ku,:) = flipud(xyz1(:,2));
        
        if order==2,
        % Second order fields
            deriv = x_eta2r(1,1:2) ; x_eta2r(1,1:2) = zeros(1,2) ;
            x_eta2 = x_eta2r+1i*x_eta2i ;
            x_eta2(:,[1,2]) = -x_eta2(:,[1,2]) ; 
            % -1 to take into account sign in ifft

            nl = 2*nfreq+1 ; 
            x_eta2 = [x_eta2(1:nl,:);zeros(3,3);conj(x_eta2(nl:-1:2,:))] ;
%            x_eta2 = [x_eta2(1:nl,:);zeros(4*nfreq+7,3);conj(x_eta2(nl:-1:2,:))] ;
            eta2 = real(ifft(x_eta2))*npts ;

            W2.Z(kv,ku,:) = flipud(eta2(:,3));       
            X2.Z(kv,ku,:) = flipud(eta2(:,1));
            Y2.Z(kv,ku,:) = flipud(eta2(:,2));
            X2.drift(kv,ku,:)=-deriv(1)*t';
            Y2.drift(kv,ku,:)=-deriv(2)*t';
        end
    end
end

t = (0:npts-1)/(npts*df) ;
W.t=t';
X.t=t';
Y.t=t';
W2.t=t';
X2.t=t';
Y2.t=t';

if ndir==1,
    W = rmfield(W,'v');
    W.Z = squeeze(W.Z(1,:,:));
    X = rmfield(X,'v');
    X.Z = squeeze(X.Z(1,:,:));
%    Y = rmfield(Y,'v');
%    Y.Z = squeeze(Y.Z(1,:,:));
    Y = [];
    if order==2,
        W2 = rmfield(W2,'v');
        W2.Z = squeeze(W2.Z(1,:,:));
        X2 = rmfield(X2,'v');
        X2.Z = squeeze(X2.Z(1,:,:));
        X2.drift=squeeze(X2.drift(1,:,:));
%        Y2 = rmfield(Y2,'v');
%        Y2.Z = squeeze(Y2.Z(1,:,:));
%        Y2.drift=squeeze(Y2.drift(1,:,:));
        Y2 = [];
    end
end
%delete(WB)
toc
