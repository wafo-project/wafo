function [f] = spec2Acdf(spec,utc,def,paramt,h,nit,speed,plotflag)
%SPEC2ACDF  CDF for crests P(Ac<=h) or troughs P(At<=h).
%         
%  CALL:  F   = spec2Acdf(S,u,def,param,h,nit,speed,plotflag);
%
%        F    = cumulative distribution function of crests or trough heights
%        S    = spectral density structure
%        u    = reference level (default the most frequently crossed level).
%       def   = 'Tc',    gives crest height in time, Ac (default).
%               'Tt',    gives trough height in time, At.
%               'Lc' and 'Lt' crest trough in space
%     param  =  [0 tn Nt] defines discretization of crest period (or
%               crest length): tn is the longest period considered while
%               Nt is the number of points, i.e. (Nt-1)/tn is the
%               sampling frequnecy. param = [0 10 51] implies that the
%               crest periods are considered at 51 linearly spaced points
%               in the interval [0,10], i.e. sampling frequency is 5 Hz.
%          h  = vector of  amplitudes; note  h >= 0,
%        nit  =  0,...,9. Dimension of numerical integration  (default 2).
%               -1,-2,-3,... different important sampling type integrations.
%       speed = defines accuraccy of calculations, by choosing different 
%               parameters, possible values: 1,2...,9 (9 fastest, default 4).
%    plotflag = 1, plots cdf  F  and Rayleigh appximation (default)
%               0, no plot 
%       []    = default values are used.
%
% SPEC2ACDF calculates P(Ac<=h) or P(At<=h) in a stationary Gaussian
% transform process X(t) where Y(t) = g(X(t)) (Y zero-mean Gaussian with
% spectrum given in S).
%
% Example: % The upper and lower bounds for F=P(Ac<hi; hi in h ),with
%     S = jonswap;
%     h = 0:0.1:5; paramt=[0 12 51];% is computed by:
%            
%     F = spec2acdf(S,[],'tc',paramt,h);
%
% See also  spec2cov, specnorm, dat2tr, dat2gaus, specdef, wavedef

% Tested on : matlab 5.3
% History: by I. Rychlik 01.10.1998 with name wave_th1.m
% revised by Per A. Brodtkorb 19.09.1999
% revised by I.R. 30.09.1999, bugs removing.
% continued by s.v.i 10.11.1999 by adding crests level in the period densities
% an then calculation of crests distribution.
% changed name and introduced possibility of computation of upper and lower
% bounds by I.R. 17.12.1999.
% Computation of distribution on denser grid I.R. 26 jan. 2000
% Revised by es 000322. Made call with directional spectrum possible.
% revised by IR 6 II 2001 adapted to MATLAB6
% Revised pab Dec 2003
%  replaced code with call to spec2cov2
%  removed tic and toc  
% Revised pab aug2007
% -replaced code with a call to writecov
% GL added plotflag 2017
  
% TODO % Needs further testing  
  
startTime = clock;
if nargin<3||isempty(def)
  def='tc';
end
if strncmpi('l',def,1)
  spec=spec2spec(spec,'k1d');
elseif strncmpi('t',def,1) 
  spec=spec2spec(spec,'freq');
else
  error('Unknown def')
end
switch lower(def)
   case  {'tc','lc'},    defnr = 1; % 'tc' or 'lc'
   case  {'tt','lt'},    defnr =-1; % 'tt' or 'lt'
   otherwise, 
     error('Unknown def')
end			    

[S, xl4, L0, L2]=specnorm(spec);

A = sqrt(L0/L2);
SCIS=0;
if nargin<6||isempty(nit)
  nit=2;
elseif nit<0
 SCIS=min(abs(nit),2);
 nit=1;
end

if isfield(spec,'tr')
   g=spec.tr;
else
   g=[];
end
if isempty(g)
  g=[sqrt(L0)*(-5:0.02:5)', (-5:0.02:5)'];
end

if nargin<2||isempty(utc)
  utc_d = gaus2dat([0, 0],g); % most frequent crossed level 
  utc   = utc_d(1,2);
end

% transform reference level into Gaussian level
uu = dat2gaus([0., utc],g);
u  = uu(2);
disp(['The level u for Gaussian process = ', num2str(u)])

if nargin<7||isempty(speed)
  speed=4;
end			    

if nargin<8, 
   plotflag = 1,
end 
  
ttn    = 1.2*ceil(2*pi*sqrt(L0/L2)*exp(u^2/2)*(0.5+erf(-sign(defnr)*u/sqrt(2))/2));
if nargin<4||isempty(paramt)
  Ntime = 51;
  t0    = ttn;
  tn    = 2*ceil(2*pi*sqrt(L0/L2)*exp(u^2/2)*(0.5+erf(-sign(defnr)*u/sqrt(2))/2));
  paramt=[t0,tn,Ntime];
else
  %t0     = paramt(1);
  tn     = paramt(2);
  Ntime  = paramt(3);
end

%t0=ttn; 
%t0=ttn/A;
tt = tn/A;
t  = levels([0, tt, Ntime]); % normalized times
N0 = 2+ceil(ttn/tn*(Ntime-1)); % the starting point to evaluate
NN0=N0;
if Ntime>101
  disp('nr. of wavelengths >  101., slow')
end
 

nr = 2;

if nargin<5||isempty(h) 
    hg= (0:0.25:5.);
    h=tranproc(u+sign(defnr)*hg',fliplr(g))-utc;
end
% make sure values are positive and in increasing order
%Transform h to Gaussian levels:   
h=reshape(h,length(h),1);
Nx=length(h);
if Nx>101
  error('nr. of amplitudes limited to 101.')
end
hg=tranproc(utc+sign(defnr)*h,g);

dt = t(2)-t(1);
% Calculating covariances
R = spec2cov2(S,nr,Ntime-1,dt);

%NB!!! the spec2thpdf.exe programme is very sensitive to how you interpolate 
%      the covariances, especially where the process is very dependent
%      and the covariance matrix is nearly singular. (i.e. for small t
%      and high levels of u if Tc and low levels of u if Tt)
%     The best is to first interpolate through FFT with a
%     high rate and then interpolate with a spline to obtain the
%     timepoints one want. However for large t
%     it often suffices to interpolate linearly provided that
%     FFT interpolation is good eneough.


disp('writing data')

filenames = {'h.in','reflev.in','dens.out'};
cleanup(filenames{:}) 

fid=fopen('h.in','wt');
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

filenames2 = writecov(R);

disp('   Starting Fortran executable.')


dos([ wafoexepath, 'cov2acdf.exe']);


switch defnr
  case   1, Htxt = 'Distribution of crest height';
  case  -1, Htxt = 'Distribution of trough height';
end 

f = createpdf;
f.title=Htxt;
f.nit=nit;
f.speed=speed;
f.SCIS=SCIS;
f.u=utc;
f.labx1='amplitude [m]';
ftmp=loaddata('dens.out')/A;
f.x{1}=h;




ftmp=reshape(ftmp,Nx,length(t));
if (length(t)>2)
  ftmp(:,2)=0.5*(ftmp(:,3)+ftmp(:,1));
end



%
%% Here we will compute the beginning of the distribution
%% on the denser grid.
%

Ntime = 61;
t0    = 0.;
tn    = ttn;
paramt1=[t0, tn, Ntime];
tt = tn/A;
t  = levels([0, tt, Ntime]); % normalized times
N0 = 1+ceil(t0/tn*(Ntime-1)); % the starting point to evaluate
dt=t(2)-t(1);


cleanup(filenames{2:end}) 

R = spec2cov2(S,nr,Ntime-1,dt);
writecov(R);

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

dos([ wafoexepath, 'cov2acdf.exe']);


ftmp1=loaddata('dens.out')/A;


ftmp1=reshape(ftmp1,Nx,length(t));
if (length(t)>2)
  ftmp1(:,2)=0.5*(ftmp1(:,3)+ftmp1(:,1));
end



X=levels(paramt1);
Nt=paramt(3);
if (NN0<Nt+1)
    t=levels([0, paramt(2), paramt(3)]);
    X=[t(NN0:Nt), X];
    N=length(X);
    [X, indX]=sort(X);
    dx=X;
    dx(1)=0.;
    dx(N)=X(N)-X(N-1);
    dx(2:N-1)=0.5*(X(3:N)-X(1:N-2));
   
    Y=[ftmp(:,NN0:Nt), ftmp1];
    Y=Y(:,indX);
    f.f=min(Y*dx',1);
%       f.f=Y;
%       f.x{1}=X;
   
else
    N=length(X);
    dx=X;
    dx(1)=0.;
    dx(N)=X(N)-X(N-1);
    dx(2:N-1)=0.5*(X(3:N)-X(1:N-2));
   
    Y = ftmp1;
    f.f=min(Y*dx',1);
%       f.f=Y;
%       f.x{1}=X;   
end
f.elapsedTime = etime(clock,startTime);

cleanup(filenames{:},filenames2{:}) 

if plotflag==1
   pdfplot(f)
   hold on
   if abs(u)<0.1
    xu=tranproc(u+f.x{1},g); 
    plot(f.x{1},1-exp(-xu.*xu/2),'k.')
    plot(f.x{1},1-exp(-xu.*xu/2),'k')
   end 
end