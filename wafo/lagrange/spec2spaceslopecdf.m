function [Fu,Fd]=spec2spaceslopecdf(S,y,lev,opt,varargin)
%SPEC2SPACESLOPECDF Computes cdf for slope at crossings of space waves 
%
%CALL: [Fu,Fd] = spec2spaceslopepdf(S,y,type,levels,options,varargin)
%
%   Fu, Fd  = cdf for Lagrange space wave slopes 
%             at up- and down-crossings of specified levels
%             For time waves, use spec2timeslopecdf
%
%   y       = cdf calculated at y
%   levels  = vector of relative levels (default = [-1 0 1 2]
%   options = struct with fields (plus some more)
%       .lp     = if .p and .lbeta exist, the output  x  is modified 
%                 Lagrange with extra -beta/(-i*omega)^p in the transfer
%                 function
%       .lalpha = alpha value for modified Lagrange
%       .lbeta  = beta value for modified Lagrange
%
% Example:
%   S=jonswap; mom=spec2mom(S);
%   opt=simoptset('du',0.125,'Nt',512,'dt',0.25);
%   levels=[0 1 2];
%   y=linspace(-2,2,1001);
%   [Fu,Fd]=spec2spaceslopecdf(S,y,levels,opt)
%   clf
%   subplot(211)
%   plot(Fu.x,Fu.f{1},Fu.x,Fu.f{2},Fu.x,Fu.f{3}); hold on
%   plot(Fd.x,Fd.f{1},'-.',Fd.x,Fd.f{2},'-.',Fd.x,Fd.f{3},'-.')
%   [Fu,Fd]=spec2spaceslopecdf(S,y,levels,opt,'lalpha',1)
%   title('Slope CDF at space up- and downcrossings, symmetric waves')
%   axis([-2 2 0 1])
%   subplot(212)
%   plot(Fu.x,Fu.f{1},Fu.x,Fu.f{2},Fu.x,Fu.f{3}); hold on
%   plot(Fd.x,Fd.f{1},'-.',Fd.x,Fd.f{2},'-.',Fd.x,Fd.f{3},'-.')
%   title('Slope CDF at space up- and downcrossings, asymmetric waves')
%   axis([-2 2 0 1])
%
% See also: spec2spaceslopepdf, spec2timeslopecdf, spec2ldat, ldat2lwav, wav2slope

% Used in Lindgren & Åberg, JOMAE, 131 (2009) 031602-:

% Tested on Matlab 8.1
% History
% GL adapted spec2spaceslopepdf to give cdf Dec 2014 for use in WafoL

m=spec2mom(S);
if nargin<4,
    opt=simoptset;
end

if nargin>=5,  opt  = simoptset(opt,varargin{:});
end

alpha=opt.lalpha;
beta=opt.lbeta;

if nargin<2 || isempty(y),
    y=(-2:0.01:2)*sqrt(m(2));
end

mom=spec2mom(S);
if nargin<3 || isempty(lev) 
    levels=[-1 0 1 2]*sqrt(mom(1));
else
    levels=lev*sqrt(mom(1));
end
Nindex=length(levels);
fu = struct('x',y);
fd = struct('x',y);
Fu = struct('x',y,'levels',levels,'note','CDF for space slopes at upcrossings');
Fd = struct('x',y,'levels',levels,'note','CDF for space slopes at downcrossings');
for ind=1:Nindex,
    fu.f{ind}=[];
    fd.f{ind}=[];
end

[rww,~,~]=spec2lcov(S,0,0,1,alpha,beta);
[rwwuu,rwxuu,rxxuu]=spec2lcov(S,0,0,5,alpha,beta);
[~,rwxu0,~]=spec2lcov(S,0,0,3,alpha,beta);
rwx0u.R=-rwxu0.R;
    
for ind=1:Nindex,
    v=levels(ind);
    a=sqrt(rwwuu.R)/(1+v*rwx0u.R/rww.R);
    b=sqrt(rxxuu.R-rwx0u.R^2/rww.R-rwxuu.R^2/rwwuu.R)/(1+v*rwx0u.R/rww.R);
    c=rwxuu.R/sqrt(rwwuu.R)/(1+v*rwx0u.R/rww.R);
    [~,fup]=slopedensity(y,a,b,c);
    [~,fdo]=slopedensity(-y,a,b,-c);
    fu.f{ind}=fup;
    fd.f{ind}=fdo;
    Fu.x=fu.x;
    Fd.x=fd.x+y(2)-y(1);
    Fu.f{ind}=cumsum(fu.f{ind})*(y(2)-y(1));
    Fd.f{ind}=cumsum(fd.f{ind})*(y(2)-y(1));
end
Fu.levels=levels;
Fd.levels=levels;

function [f0,f]=slopedensity(x,a,b,c)
%
% SLOPEDENSITY computes the density of aR/(1+bU+cR), R Rayleigh, U Normal
%
% CALL: [f0,f] = slopedensity(x,a,b,c)
%
%       f0  = pdf of aR/(1+bU)
%       f   = pdf of aR/(1+bU+cR)
%       x   = vector of x-values for which pdf should be computed

if nargin<4,
    c=0;
end

d=1-c*x/a;

z=zeros(1,length(d));
for k=1:length(d),
  if (d(k)>0) || (d(k)<0),
    z(k)=x(k)/d(k);
  else
    z(k)=0;
  end
end  

sigma2=(a^2+b^2*z.^2);
sigma=sqrt(sigma2);

f1=a*sign(z)./sigma.^3.*z.*exp(-z.^2/2./sigma2);
f1=f1.*(b.^2+a.^2./sigma2).*cdfnorm(a/b*sign(z)./sigma);
f2=1/sqrt(2*pi)*(a.^2*b./sigma2.^2).*z.*exp(-0.5/b^2);
f0=f1+f2;

f=zeros(1,length(d));
for k=1:length(d),
  if (d(k)>0) || (d(k)<0),
    f(k)=f0(k)/d(k)^2;
  else
    f(k)=-99;
  end
end

for k=1:length(d),
  if f(k)==-99,
    f(k)=(f(k-1)+f(k+1))/2;
  end
end

