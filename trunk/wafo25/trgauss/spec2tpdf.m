function [f] = spec2tpdf(spec,utc,def,paramt,h,nit,speed,plotflag)
%SPEC2TPDF Density of crest/trough- period or length.      
%  
%  CALL:  F   = spec2tpdf(S,u,def,paramt,h,nit,speed,plotflag);
%
%        F    = density structure.
%        S    = spectral density structure
%        u    = reference level (default the most frequently crossed level).
%       def   = 'Tc',    gives half wave period, Tc (default).
%               'Tt',    gives half wave period, Tt
%               'Lc' and 'Lt' ditto for wave length.
%     paramt  = [t0 tn Nt] where t0, tn and Nt is the first value, last value 
%               and the number of points, respectively, for which
%               the density will be computed. paramt= [5 5 51] implies
%               that the density is computed only for T=5 and
%               using 51 equidistant points in the interval [0,5].
%          h  = amplitude condition: density for waves with crests above h.
%               (note  h >= 0), (default 0.) 
%               if abs(h) less 1% of standard deviation of process X then h=0.
%        nit  = 0,...,9. Dimension of numerical integration  (default 2).
%               -1,-2,-3,... different important sampling type integrations.
%       speed = defines accuraccy of calculations by choosing different 
%               parameters, possible values: 1,2...,9 (9 fastest, default 5).
%     plotflag= if 0 then do not plot, else plot (default 0).
%       []    = default values are used.
%
% SPEC2TPDF calculates pdf of halfperiods Tc, Tt, Lc or Lt  such that Ac>h or At>h,
% in a stationary Gaussian transform process X(t), where Y(t) = g(X(t))
% (Y zero-mean Gaussian with spectrum given in S). The tr. g, can be estimated
% using lc2tr, dat2tr, hermitetr or ochitr.
%
% Example:% For the directional sea compute density of encountered Tc in
%         % the direction  pi/4 from the principal wave direction (0) at
%         % points 0.,0.1,...,10.
%         
%    Hm0=7; Tp=11;  S = jonswap(6*pi/Tp,[Hm0 Tp]); 
%    D=spreading(linspace(-pi,pi,51),'cos2s',[],[],S.w);               
%    Sdir=mkdspec(S,D,1); 
%    Senc=spec2spec(Sdir,'enc',pi/8,10); paramt=[0 10 101];
%    f = spec2tpdf(Senc,[],'Tc',paramt,[],-1,[],1);
%    hold on; f1 = spec2tpdf(Sdir,[],'Tc',paramt,[],-1,[],1);
%
% See also  spec2cov, specnorm, datastructures, wavedef, wafomenu

% Tested on : matlab 5.3
% History: by I. Rychlik 01.10.1998 with name wave_th1.m
% revised by Per A. Brodtkorb 19.09.1999
% revised by I.R. 30.09.1999, bugs removing.
% continued by s.v.i 10.11.1999 by adding crests level in the period densities
% an then calculation of crests distribution.
% changed name and introduced possibility of computation of upper and lower
% bounds by I.R. 17.12.1999.
% revised by es 28.01.2000  Adjusting for directional spectrum and wave length.
% revised by es 28.01.2000 help text
% revised by IR removing error in transformation 29 VI 2000
% revised by IR adopting to Matlab 6.0 -  6 II 2001 
% revised pab 30nov2003
% revised pab May 2007
% Revised pab Aug2007
% -replaced some code with call to writecov

startTime = clock;
if nargin<3||isempty(def)
  def='tc';
end
if strcmpi('l',(def(1)))
  spec=spec2spec(spec,'k1d');
elseif strcmpi('t',(def(1)))
  spec=spec2spec(spec,'freq');
else
  error('Unknown def')
end
switch lower(def)
  case  {'tc','lc'},    defnr = 1; % 'tc' or 'lc'
  case  {'tt','lt'},    defnr =-1; % 'tt' or 'lt'
  otherwise 
    error('Unknown def')
end			    

[S, xl4, L0, L2]=specnorm(spec);
A=sqrt(L0/L2);
SCIS=0;
if nargin<6||isempty(nit)
  nit=2;
elseif nit<0
 SCIS=min(abs(nit),11);
 nit=1;
end

if isfield(spec,'tr')
  g=spec.tr;
else
  g=[];
end
if isempty(g)
  g=[(-5:0.02:5)', (-5:0.02:5)'];
  g(:,1)=sqrt(L0)*g(:,1);
end


if nargin<2||isempty(utc)
  utc_d=gaus2dat([0, 0],g); % most frequent crossed level 
  utc=utc_d(1,2);
end

% transform reference level into Gaussian level
uu=dat2gaus([0., utc],g);
u=uu(2);
disp(['The level u for Gaussian process = ' num2str(u)])



if nargin<7||isempty(speed)
  speed=5;
end			    

    
if nargin<8||isempty(plotflag)
  plotflag=0;
end			    


if nargin<4||isempty(paramt)
  Ntime = 51;
  t0    = 0.;
  tn    = 2*ceil(2*pi*sqrt(L0/L2)*exp(u^2/2)*(0.5+erf(-sign(defnr)*u/sqrt(2))/2));
else
  t0     = paramt(1);
  tn     = paramt(2);
  Ntime  = paramt(3);
end

t  = levels([0, tn/A, Ntime]); % normalized times

N0 = 1+round(t0/tn*(Ntime-1)); % the starting point to evaluate
%if Ntime>101
%  disp('nr. of wavelengths limited to 101.')
%end
 


px=gaus2dat([0., u;1, 5],g); 
px=abs(px(2,2)-px(1,2));
Nx = 1;
if nargin<5||isempty(h)
  h=px;
  h0=0.;
else
  h=abs(min(h));
  h0=h;
  if h0>0.01*sqrt(L0)
    Nx=2;
    h=[h; px];
  else
    h=px;
    h0=0.;
  end
end

h=reshape(h,length(h),1);
hg=tranproc(utc+sign(defnr)*h,g);





nr = 2;
dt = t(2)-t(1);
R  = spec2cov2(S,nr,Ntime-1,dt);


disp('writing data')
filenames0 = writecov(R);

filenames = {'h.in','reflev.in','dens.out'};
cleanup(filenames{:})

fid = fopen('h.in','wt');
fprintf(fid,'%12.10f\n',hg);
fclose(fid);

fid=fopen('reflev.in','wt');
fprintf(fid,'%12.10E \n',u);
%defnr;

fprintf(fid,'%2.0f \n',defnr);
fprintf(fid,'%2.0f \n',Ntime);
fprintf(fid,'%2.0f \n',N0);
fprintf(fid,'%2.0f \n',nit);
fprintf(fid,'%2.0f \n',speed);
fprintf(fid,'%2.0f \n',SCIS);
fprintf(fid,'%2.0f \n',10^9*rand);  % select a random seed for rind 
fprintf(fid,'%2.0f \n',Nx);
fprintf(fid,'%12.10E \n',dt);
fclose(fid); 
disp('   Starting Fortran executable.')

dos([ wafoexepath, 'cov2acdf.exe']); %compiled cov2acdf.f with rind70.f and intmodule.f

f=createpdf;
switch lower(def)
 case  'tc'
  Htxt = ['Density of Tc with Ac>', num2str(h0), '  u=',num2str(utc)];
  xtxt = 'T [s]';
 case  'tt'
  Htxt = ['Density of Tt with At>', num2str(h0), '  u=',num2str(utc)];
  xtxt = 'T [s]';
 case  'lc'
  Htxt = ['Density of Lc with Ac>', num2str(h0), '  u=',num2str(utc)];
  xtxt = 'L [m]';
 case  'lt'
  Htxt = ['Density of Lt with At>', num2str(h0), '  u=',num2str(utc)];
  xtxt = 'L [m]';
end 
 
f.title=Htxt;
f.labx{1}=xtxt;
ftmp=loaddata('dens.out')/A;
f.x{1}=t*A;

ftmp=reshape(ftmp,Nx,length(t));
if (length(t)>2)
  ftmp(:,2)=0.5*(ftmp(:,3)+ftmp(:,1));
end
if (Nx>1)
  ft_up=ftmp(2,1:Ntime)-ftmp(1,1:Ntime);
else
  ft_up=ftmp(1,1:Ntime);
end
f.f=ft_up;
f.nit=nit;
f.speed=speed;
f.SCIS=SCIS;
f.u=utc;
f.elapsedTime = etime(clock,startTime);

if plotflag 
  pdfplot(f)
end

cleanup(filenames0{:},filenames{:})



