function [fu,fd]=spec2spaceslopepdf(S,y,levels,opt,varargin)
%SPEC2SPACESLOPEPDF Computes pdf for slope at crossings of space waves 
%
%CALL: [fu,fd] = spec2spaceslopepdf(S,y,levels,options,varargin)
%
%   fu, fd  = pdf for Lagrange space wave slopes 
%                       at up- and down-crossings of specified levels
%                       For time waves, use spec2timeslopecdf
%
%   y       = pdf calculated at  y
%   levels  = vector of relative levels (default = [-2 -1 0 1 2])
%   options = struct with fields (plus some more)
%       .lp     = if .p and .lbeta exist, the output  x  is modified 
%                 Lagrange with extra -beta/(-i*omega)^p in the transfer
%                 function
%       .lalpha = alpha value for modified Lagrange
%       .lbeta  = beta value for modified Lagrange
%       .ltheta = if exist produces [w,x] to be transformed, theta = 
%                 angle for the transformation
%                 .lp, .lbeta and .ltheta are empty ([]) by default
%       .plotflag - 'off', no plotting (default)
%                 - 'on' 
%
% Example: See example for  spec2spaceslopepdf
%
% Ref: Lindgren & Åberg, JOMAE, 131 (2009) 031602-1
% See also: spec2ldat, spec2timeslopecdf, ldat2lwav, wav2slope

% Tested on Matlab 8.1, 8.6
% History:
% Original GL december 2007
% Modified GL Nov 2014 to include code for  slopedensity

m=spec2mom(S);
if nargin<4,
    opt=simoptset;
end

if nargin>=5,  opt  = simoptset(opt,varargin{:});
end

plotta=strcmp(opt.plotflag,'on');

alpha=opt.lalpha;
beta=opt.lbeta;

if nargin<2 || isempty(y),
    y=(-2:0.01:2)*sqrt(m(2));
end

if nargin<3 || isempty(levels) 
    mom=spec2mom(S);
    levels=[-2 -1 0 1 2]*sqrt(mom(1));
end

Nindex=length(levels);
fu = struct('x',y,'levels',levels);
fd = struct('x',y,'levels',levels);
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
    if plotta, 
        figure(10)
        clf
%        plot(fu.x', cumsum(fu.f{ind})'*(y(2)-y(1)))
%        hold on
%        plot(fd.x', cumsum(fd.f{ind})'*(y(2)-y(1)),'r')
        plot(fu.x, fu.f{ind})
        hold on
        plot(fd.x, fd.f{ind},'r')
    end
end

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



