function [s2,s1,L0]=spec2slcomp(S,options,varargin)
%SPEC2SLCOMP Compares 2nd order Stokes and 1st order Lagrange time waves 
%           
%CALL: [s2,s1,L0] = spec2slcomp(S,options,varargin)
% 
%   s2    = 2nd order Stokes wave structure s2.Z,s2.t
%   s1    = Gaussian wave structure s1.Z,s1.t
%   L0    = 1st order (smoothed) Lagrange wave structure L.Z,L.t 
%   spec  = S   a frequency spectral density structure in 
%               angular frequency ('w') or frequency ('f') form 
%   options = struct with fields 
%       .Nt = giving  Nt  time points.  (default length(S)-1=n-1).
%             If Nt>n-1 it is assummed that  S.S(k)=0  for  k>n-1
%       .dt = step in grid (default dt is defined by the Nyquist freq) 
%       .u  = [-u1 u1 Nu] gives  u = linspace(-u1,u1,Nu) grid
%             Nu should be an odd integer
%             Generated waves are time waves observed at  u = 0  
%  
% The routine is a combination of 
%   spec2nldat (Stokes) and spec2ldat/ldat2lwav (Lagrange)
%
% See also  spec2nldat spec2linspec, spec2ldat

% Tested on 8.1, 8.6
% History:
% by GL May 2015

% Variables controlling the truncation of the spectrum for sum and
% difference frequency effects  
reltol2ndorder     = 1e-3; %

ftype = freqtype(S); %options are 'f' and 'w' and 'k'
n     = length(S.(ftype));

numWaves = 1000; % Typical number of waves in 3 hour seastate
%C = pi/7; % Wave steepness criterion
%C = 1/2 ; % No bump in trough
C = 2;    % Convergence criterion as given in Nestegaard and Stokka (1995)
if (nargin<6 || isempty(truncationLimit)), 
  truncationLimit = sqrt(C);
end

if nargin<2,
    opt=simoptset;
else
    opt=options;
end
if nargin>=3, opt = simoptset(opt,varargin{:}); end

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

Nu=opt.Nu;
Nt=opt.Nt;
du=opt.du;
dt=opt.dt;
np=Nt;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U=opt.u;
if isempty(U),
    [w0,x0]=spec2ldat(S);
    sx=std(x0.Z(:));
    U=[-round(10*sx) round(10*sx) 401];
end

if ~isempty(U),
    Nu=U(3);
    u1=U(1);
    u2=U(2);
    u=linspace(u1,u2,Nu);
    u=u;
    du=(u2-u1)/Nu;
end

t=(0:Nt-1)*dt;
w.Z=zeros(Nu,Nt);
w.u=u';
w.t=t;
x.Z=zeros(Nu,Nt);
x.u=u';
x.t=t;

mom=spec2mom(S);
w.meanperiod=2*pi*sqrt(mom(1)/mom(2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

S=specinterp(S,dt);

np = np + mod(np,2); % make sure np is even    

fs    = S.(ftype);
Si    = S.S(2:end-1);
h     = S.h;
if isempty(h), h = inf;end

switch ftype
  case 'f'   
  case {'w','k'}
    Si = Si*2*pi;
    fs = fs/2/pi;
  otherwise
    error('Not implemented for wavenumber spectra')
end
dT = 1/(2*fs(end)); % dT
df = 1/(np*dT);

% interpolate for freq.  [1:(N/2)-1]*df and create 2-sided, uncentered spectra
% ----------------------------------------------------------------------------
f = (1:(np/2)-1)'*df;

fs(1)   = []; 
fs(end) = [];
Fs      = [0; fs(:); (np/2)*df];
Su      = [0; abs(Si(:))/2; 0];

Smax = max(Su);
%Si = interp1q(Fs,Su,f);
Si = interp1(Fs,Su,f, 'linear');

tmp  = find(Si>Smax*reltol2ndorder);
Hm0  = spec2char(S,'Hm0');
Tm02 = spec2char(S,'Tm02');
waterDepth = abs(S.h);
kbar = w2k(2*pi/Tm02,0,waterDepth);
Amax = sqrt(2*log(numWaves))*Hm0/4; % Expected maximum amplitude for 1000 waves seastate
fLimitUp = truncationLimit*sqrt(gravity*tanh(kbar*waterDepth)/Amax)/(2*pi);
fLimitLo = sqrt(gravity*tanh(kbar*waterDepth)*Amax/waterDepth^3)/(2*pi);

nmax   = min(max(find(f<=fLimitUp)),max(find(Si>0)))+1;
nmin   = max(min(find(fLimitLo<=f)),min(tmp))+1;

if isempty(nmax),nmax = np/2;end
if isempty(nmin),nmin = 2;end % Must always be greater than 1
fLimitUp = df*nmax;
fLimitLo = df*nmin;

%disp(sprintf('2nd order frequency Limits = %g,%g',fLimitLo, fLimitUp))
Su = [0; Si; 0; Si((np/2)-1:-1:1)];
clear Si Fs

T       = (np-1)*dT;
s1       = zeros(np,2);
s1(:,1)  = linspace(0,T,np)'; %(0:dT:(np-1)*dT).';
s2      = s1;

ww  = 2*pi*[0; f; np/2*df];

g  = gravity;
kw = w2k(ww ,[],h,g);
Ku = kw*u;

% Generate standard normal random numbers for the simulations
% -----------------------------------------------------------
Zr = randn((np/2)+1,1);
Zi = [zeros(1,1); randn((np/2)-1,1); zeros(1,1)];

%A = zeros(np,1);
%A(1:(np/2+1),:)  = Zr - sqrt(-1)*Zi; clear Zr Zi

A1 = Zr - sqrt(-1)*Zi; clear Zr Zi
B1 = exp(1i*Ku).*repmat(A1,1,Nu);

%A((np/2+2):np,:) = conj(A(np/2:-1:2,:));
B2 = conj(B1(np/2:-1:2,:));
A=[B1;B2];
A(1,:)           = A(1,:)*sqrt(2);
A((np/2)+1,:)    = A((np/2)+1,:)*sqrt(2);

icoth = 1i*(tanh(kw*S.h)).^(-1);
H = [0; icoth(2:end)]; % + (lalpha*ones(size(kappa')) + 1i*lbeta*kappa')./(omega').^2;
C1 = B1.*repmat(H,1,Nu);
C2 = conj(C1(np/2:-1:2,:));
C=[C1;C2];
C(1,:)           = C(1,:)*sqrt(2);
C((np/2)+1,:)    = C((np/2)+1,:)*sqrt(2);

% Make simulated time series
% --------------------------

Ssqr = sqrt(Su*df/2);
A = A.*repmat(Ssqr,1,Nu);
C = C.*repmat(Ssqr,1,Nu);

clear Su Ssqr

%x(:,2:end) = real(fft(A));
A0=A(:,round((Nu-1)/2)+1);
s1(:,2) = real(fft(A0));
W = real(fft(A));
X = real(fft(C));
w.Z=W';
x.Z=X';

[L,L0] = ldat2lwav(w,x,'time',0,5,0);

%if nargout>3, 
   %compute the sum and frequency effects separately
  [svec, dvec] = disufq((A0.'),ww,kw,min(h,10^30),g,nmin,nmax);
  svec = svec.';
  dvec = dvec.';
  
  x2s  = fft(svec); % 2'nd order sum frequency component 
  x2d  = fft(dvec); % 2'nd order difference frequency component
  
  % 1'st order + 2'nd order component.
  s2(:,2) =s1(:,2)+ real(x2s(1:np,:))+real(x2d(1:np,:)); 
  
%else
%  svec = disufq((A0.'),ww,kw,min(h,10^30),g,nmin,nmax).';
%  x2o  = fft(svec); % 2'nd order component 
% 1'st order + 2'nd order component.
%x2(:,2:end)=x(:,2:end)+ real(x2o(1:np,:)); 
%end

if opt.plotflag ~= 0,
    figure(2)
    clf
    plot(s1(:,1),s1(:,2),'b')
    hold on
    plot(s2(:,1),s2(:,2),'r')
    plot(L0.t,L0.Z,'k')
end
return




