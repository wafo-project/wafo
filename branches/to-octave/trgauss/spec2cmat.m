function [f, fmm] = spec2cmat(spec,utc,def,paramt,paramu,nit) 
%SPEC2CMAT Joint intensity matrix for cycles (max,min)-, rainflow- and (crest,trough)
%          
% CALL:  f   = spec2cmat(S,u,def,paramt,paramu,nit); 
% 
%        f    = pdf (density structure) of crests (trough) heights 
%        S    = spectral density structure 
%        u    = reference level (default the most frequently crossed level).
%       def   = 'Mm'  : gives maximum and the following minimum height (default). 
%               'rfc' : gives maximum and the rainflow minimum height.
%               'AcAt': gives (crest,trough) heights  (this option needs
%                       more work). 
%     paramt  = [0 tn Nt] defines discretization of half period: tn is
%               the longest period considered while Nt is the number of
%               points, i.e. (Nt-1)/tn is the sampling frequnecy. 
%               paramt=[0 10 51] implies that the halfperiods are
%               considered at 51 linearly spaced points in the interval
%               [0,10], i.e. sampling frequency is 5 Hz. 
%     paramu  = [u v N] defines discretization of maxima and minima ranges: 
%               u is the lowest minimum considered, v the heighest
%               maximum and N is the number of levles (u,v) included. 
%        nit  =  0,...,9. Dimension of numerical integration (only
%                positive nit are allowed). (default nit=1). 
%       []    = default values are used. 
% 
%      
%  The model for loads is a stationary Gaussian transformed process X(t),
%  where  Y(t) = g(X(t)) is a zero-mean Gaussian with spectrum, S.
%  
%  Note: algorithm uses Markov Chain approximation to the sequence of
%  turning points in Y.   
%
% Example: % The intensity matrix of rainflow cycles is computed by: 
%       S  = jonswap;      
%       L0 = spec2mom(S,1); 
%       paramu = [sqrt(L0)*[-4 4] 41]; 
%       frfc   = spec2cmat(S,[],'rfc',[],paramu); 
% 
% See also  spec2mmtpdf spec2cov, specnorm, dat2tr, dat2gaus, wavedef 
 

% Tested on : matlab 5.3 
% History: by I. Rychlik 01.10.1998 with name minmax.m 
% bounds by I.R. 02.01.2000. 
% Revised by es 000322. Made call with directional spectrum possible.
% revised by ir 000612. Help and plots improved. 
% revised by IR removing error in transformation 29 VI 2000
% revised by I.R. 01.20.2001 Change name minmax to wminmax
% revised by I.R. 6 II 2001 adapted for MATLAB 6
% revised pab 30nov2003

% TODO % AcAt option needs more work  

startTime = clock; 
[S, xl4, L0, L2, L4]=specnorm(spec);

 
A = sqrt(L0/L2); 
%SCIS=0; 
if nargin<6||isempty(nit) 
  nit=1; 
elseif nit<0 
  warning('WAFO:SPEC2CMAT','Only postive nit allowed')
 nit=1; 
end 
 
if isfield(spec,'tr') 
   g = spec.tr; 
else 
   g = []; 
end 
if isempty(g) 
  g = [sqrt(L0)*(-5:0.02:5)', (-5:0.02:5)']; 
end
S.tr = g;

if nargin<2||isempty(utc)
    utc_d = gaus2dat([0, 0],g); % most frequent crossed level 
    utc   = utc_d(1,2);
end

% transform reference level into Gaussian level
uu = dat2gaus([0., utc],g);
u  = uu(2);
disp(['The level u for Gaussian process = ', num2str(u)])


if nargin<4||isempty(paramt) 
  distanceBetweenExtremes = 5*pi*sqrt(L2/L4); %(2.5 * mean distance between extremes)
  paramt = [0 distanceBetweenExtremes,43];   
end
%t0     = paramt(1); 
tn     = paramt(2); 
Ntime  = paramt(3); 
t      = levels([0, tn/A, Ntime]); % normalized times 


if nargin<3||isempty(def) 
    defnr=1;     
else 
 switch lower(def) 
 case  'mm',    defnr = 1;
 case  'rfc',   defnr = 0;
 case  'acat',  defnr =-1; 
 otherwise, error('Unknown def') 
 end 
end			     

 
%nr = 4; 
 
if nargin<5||isempty(paramu)  
   paramu=[-4*sqrt(L0), 4*sqrt(L0), 41]; 
end 
 %Transform amplitudes to Gaussian levels:    
h   = levels(paramu); 
h   = reshape(h,length(h),1); 
Nx  = length(h); 
%der = ones(Nx,1);  
%hg  = tranproc([h, der],g); 
%der = abs(hg(:,2)); 
%hg  = hg(:,1); % Gaussian level 

%if exist('h.in'), delete h.in, end
%fid=fopen('h.in','wt');
%fprintf(fid,'%12.10f\n',hg); 
%fclose(fid);
 
 

%paru=paramu;
% paru(1:2)=paru(1:2)/sqrt(L0);
ftmp = wminmax(S,nit,paramu,t);
 
dh=h(2)-h(1);
ftmp = ftmp *dh*dh;
%for i=1:Nx 
%  ftmp(:,i)=dh*dh*ftmp(:,i);%.*der*der(i);%* sqrt(-R(1,6)/R(1,4))/2/pi;
%end 


f=createpdf;
f.nit=nit; 
f.x{1}=h; 
f.x{2}=h;
if (defnr==0)
  f.f      =fliplr(mctp2rfc(fliplr(ftmp)));%* sqrt(-R(1,6)/R(1,4))/2/pi;
  f.title  ='Joint density of maximum and rainflow minimum';
  f.labx{1}='max [m]';
  f.labx{2}='rainflow min [m]';
elseif (defnr==-1)   
   %CC= normalizing constant= 1/ expected number of u-up-crossings of X
   %CC = 2*pi*sqrt(L0/L2)*exp(0.5D0*u*u/L0); 
   % CC = normalizing constant = 1/ expected number of zero-up-crossings of X' 
   %CC = 2*pi*sqrt(L2/L4); 
  fact = sqrt(L0/L4);
   
  f.f    = fliplr(mctp2tc(fliplr(ftmp*fact),utc,paramu));
  index1 = find(f.x{1}>0);
  index2 = find(f.x{2}<0);
  f.f    = flipud(f.f(index2,index1));
  f.x{1} = f.x{1}(index1);
  f.x{2} = abs(flipud(f.x{2}(index2)));
  f.title  ='Joint density of crest and trough';  
  f.labx{1}='Crest [m]';
  f.labx{2}='Trough [m]';
else  %(defnr==1)
  f.f=ftmp;
  f.title  ='Joint density of maximum and minimum';  
  f.labx{1}='max [m]';
  f.labx{2}='min [m]';
end



try
  [f.cl,f.pl] = qlevels(f.f,[10, 30, 50, 70, 90, 95, 99, 99.9],f.x{1},f.x{2});
end
f.elapsedTime = etime(clock,startTime);


if nargout>1

  fmm=createpdf;
  fmm.title = 'Joint density of maximum and minimum';
  fmm.labx{1}='max [m]';
  fmm.labx{2}='min [m]';
  fmm.nit=nit;
  fmm.x{1}=h;
  fmm.x{2}=h;
  fmm.f=ftmp;
end
 
