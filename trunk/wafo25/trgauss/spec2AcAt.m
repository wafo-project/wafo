function [f] = spec2AcAt(spec,utc,def,paramtc,paramtt,paramt,h1,h2,nit,speed)
%SPEC2ACAT  Survival function for crests and troughs, R(h1,h2)=P(Ac>h1,At>h2).
%         
%  CALL: R = spec2AcAt(S,u,def,paramtc,paramtt,paramtcc,h1,h2,nit,speed,bound);
%
%        R    =  survival function R=P(Ac>h1,At>h2)
%        S    = spectral density structure
%        u    = reference level (default the most frequently crossed level).
%       def   =  'Tcc', P(Ac>h1,At>h2) for waves in time, (default)
%                'Lcc', P(Ac>h1,At>h2) for waves in space.
%    paramtc  = [0 tc Ntc] defines discretization of Tc (or Lc): tc is
%               the longest Tc considered while Ntc is the number of
%               points, i.e. (Ntc-1)/tc is the sampling
%               frequency. paramtc= [0 10 51] implies that the crest
%               periods are considered at 51 linearly spaced points in
%               the interval [0,10], i.e. sampling frequency is 5 Hz.
%   paramtt   = [0 tt Ntt] defines discretization of Tt (or Lt): tt is
%               the longest Tt considered while Ntt is the number of points.
%  paramtcc   = [0 tcc Ntcc] defines discretization of a period Tcc (or
%               Lcc): tcc is the longest period considered while Ntcc is
%               the number of points.
%          h1 = vector of heights of Ac crest; note  h1 > 0,
%          h2 = vector of heights of At trough; note  h2 > 0.
%        nit  =  0,...,9. Dimension of numerical integration  (default 2).
%               -1,-2,-3,... different important sampling type integrations.
%       speed = defines accuraccy of calculations, by choosing different 
%               parameters, possible values: 1,2...,9 (9 fastest, default 4).
%       []    = default values are used.
%
% SPEC2ACAT evaluates P(Ac>h1,At>h2), i.e. probability that crest is
% higher than h1 and trough lower than h2, for waves in time or in space
% in a stationary Gaussian transform process X(t) where Y(t) = g(X(t) (Y
% zero-mean Gaussian with spectrum given in S). 
%
% Example:%  Very accurate approx. of R=P(Ac>h1,At>h2) (waves in
%         %  time with unidirectional JONSWAP spectrum) is computed by:
%   Hm0=7;  Tp=11;
%   S = jonswap(4*pi/Tp,[Hm0 Tp]); 
%   paramtt = [0 8 51]; 
%   paramtc = paramtt; 
%   paramtcc=[0 10 31]; 
% %   paramtt = [0 12 51]; 
% %   paramtc = paramtt; 
% %   paramtcc=[0 18 81]; 
%   h1=0.5:0.5:6; h2=h1;    
%   F = spec2AcAt(S,[],'Tcc',paramtc,paramtc,paramtcc,h1,h2,-2);
%
% See also  spec2cov, specnorm, dat2tr, dat2gaus, datastructures, wavedef

% Tested on : matlab 5.3
% History: by I. Rychlik 01.10.1998 with name wave_t.m
% Revised by I. R.   13.10.1999
% Revised by Sylvie van Iseghem 10.11.1999- adding levels of crests and troughs
% Revised by I. R.   07.12.1999 The upper and lower bounds to  the
%  density included.
% Revised by I. R.   17.12.1999 The case def=-1,-2 are included.
% Revised by es 000322. Made call with directional spectrum possible.
% Revised by IR 29 VI 2000, introducing 'def' to consider also waves in space.
%                           and removing bug from definition of transformation.
% Revised IR 6 II 2001 adapted for MATLAB6
% REvised pab Dec 2003
% replace tic-toc statements with clock and etime  
% Revised pab aug2007
% -replade code with a call to writecov
% -added call to cleanup


% TODO % paramtt is never used remove it???????

startTime = clock;
if nargin<3||isempty(def)
  def='tcc';
end
if strncmpi('l',def,1)
  spec=spec2spec(spec,'k1d');
elseif strncmpi('t',def,1)
  spec=spec2spec(spec,'freq');
else
  error('Unknown def')
end


[S, xl4, L0, L2]=specnorm(spec);
A = sqrt(L0/L2);
SCIS=0;
if nargin<9||isempty(nit)
  nit=2;
elseif nit<0
 SCIS=min(abs(nit),2);
end

if isfield(spec,'tr')
   g=spec.tr;
else
   g=[];
end
if isempty(g)
  g = [sqrt(L0)*(-5:0.02:5)', (-5:0.02:5)'];
end

if nargin<2||isempty(utc)
  utc_d = gaus2dat([0, 0],g); % most frequent crossed level 
  utc   = utc_d(1,2);
end

% transform reference level into Gaussian level
uu = dat2gaus([0., utc],g);
u  = uu(2);
disp(['The level u for Gaussian process = ', num2str(u)])


if nargin<10||isempty(speed)
  speed=4;
end			    
  
if nargin<4||isempty(paramtc)
  defnr=1;
  paramtc=[0., 2*ceil(2*pi*sqrt(L0/L2)*exp(u^2/2)*(0.5+erf(-sign(defnr)*u/sqrt(2))/2)), 51];
end
if nargin<5||isempty(paramtt)
  %defnr=-1;
  %paramtt=[0., 2*ceil(2*pi*sqrt(L0/L2)*exp(u^2/2)*(0.5+erf(-sign(defnr)*u/sqrt(2))/2)), 51];
end
if nargin<6||isempty(paramt)
  Ntime = 51;
  t0    = 0.;
  tn    = 1.5*ceil(2*pi*sqrt(L0/L2)*exp(u^2/2));
  %paramt=[t0,tn,Ntime];
else
  t0     = paramt(1);
  tn     = paramt(2);
  Ntime  = paramt(3);
end
t  = levels([0, tn/A, Ntime]); % normalized times
N0 = 1+ceil(t0/tn*(Ntime-1)); % the starting point to evaluate
if Ntime>101
  disp('Note: nr. of wavelengths (periods) exceeds 101 (slow).')
end

nr = 2; 

px1=gaus2dat([0., u;1, 4],g); 
px1=abs(px1(2,2)-px1(1,2));
px2=gaus2dat([0., u;1, 4],g); 
px2=abs(px2(2,2)-px2(1,2));

if nargin<7||isempty(h1)
  Nx1 = 6;
  h1=levels([0, px1, 6]);
else
  h1=sort(abs(h1(:))); % make sure values are positive and in increasing order
  Nx1=length(h1); 
end
figure(1)
clf
subplot(221)
nit0=nit;
if nit0>=0
  nit0=nit0+3;
end
if strncmpi('l',def,1) 
  FAc=spec2Acdf(spec,utc,'Lc',paramtc,h1,nit0,speed); 
elseif strncmpi('t',def,1)
  FAc=spec2Acdf(spec,utc,'Tc',paramtc,h1,nit0,speed);
else
  error('Unknown def')
end
pdfplot(FAc);
hold on
if nargin<7||isempty(h2)
  Nx2 = 6;
  h2=levels([0, px2, 6]);
else
  h2=sort(abs(h2(:))); % make sure values are positive and  increasing 
  Nx2=length(h2);
end
subplot(222)
if strncmpi('l',def,1)
  FAt=spec2Acdf(spec,utc,'Lt',paramtc,h1,nit0,speed); 
elseif strncmpi('t',def,1)
  FAt=spec2Acdf(spec,utc,'Tt',paramtc,h1,nit0,speed);
else
  error('Unknown def')
end

pdfplot(FAt);
 
 
 
Nx=Nx1*Nx2;

if Nx>101
   disp('Note: nr. of amplitudes exceeds 101, (slow).')
end


dt=t(2)-t(1);
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



%Transform amplitudes to Gaussian levels:   
der1=ones(Nx1,1); % dh/dh=1
h1=reshape(h1,Nx1,1);
hg1=tranproc([utc+h1, der1],g);
%der1=abs(hg1(:,2));
hg1=hg1(:,1); % Gaussian level
der2=ones(Nx2,1); % dh/dh=1
h2=reshape(h2,Nx2,1);   
hg2=tranproc([utc-h2, der2],g);
%der2=abs(hg2(:,2));
hg2=hg2(:,1); % Gaussian level



disp('writing data')

filenames = {'h.in','reflev.in','dens.out'};
cleanup(filenames{:})

fid=fopen('h.in','wt');
fprintf(fid,'%12.10f\n',hg1);
fprintf(fid,'%12.10f\n',hg2);
fclose(fid);


%SCIS=0;


fid=fopen('reflev.in','wt');
fprintf(fid,'%12.10E \n',u);
fprintf(fid,'%2.0f \n',Ntime);
fprintf(fid,'%2.0f \n',N0);
fprintf(fid,'%2.0f \n',nit);
fprintf(fid,'%2.0f \n',speed);
fprintf(fid,'%2.0f \n',SCIS);
fprintf(fid,'%2.0f \n',10^9*rand);  % select a random seed for rind 
fprintf(fid,'%2.0f %2.0f\n',[Nx1, Nx2]);
fprintf(fid,'%12.10E \n',dt);
fclose(fid); 


filenames2 = writecov(R);

disp('   Starting Fortran executable.')

  dos([ wafoexepath, 'cov2tccpdf.exe']); 
  %compiled spec2tccpdf.f with rind60.f and intmodule.f
  %   dos([ wafoexepath, 'sp2tccpdf50.exe']); %compiled spec2tccpdf.f with rind50.f

 
f=createpdf;
f.labx{1}='amplitude [m]';
f.x{1}=h1';
f.x{2}=h2';
f.nit=nit;
f.speed=speed;
f.SCIS=SCIS;
f.u=utc;
tmp=loaddata('dens.out')/A;  
t=t'*A;
dt=t(2)-t(1);
ff1=reshape(tmp(:,1),Nx1*Nx2,length(t));

subplot(212)
fff=zeros(Nx1,Nx2);
for j=1:Nx1
  for i=1:Nx2
    ii1=i+(j-1)*Nx2;
    fff(j,i)=max(dt*sum(ff1(ii1,:))+1.-FAc.f(j,1)-FAt.f(i,1),0.);
    fff(j,i)=min(fff(j,i),1.);
  end
end
f.f=fff';
contour(h1,h2,fff',[1, 0.9, 0.75, 0.5, 0.25, 0.1, 0.05])
axis('square')
%   pdfplot(f);
f.elapsedTime = etime(clock,startTime);

cleanup(filenames{:},filenames2{:})


