function [f] = spec2tccpdf(spec,utc,def,paramt,h1,h2,nit,speed)
%SPEC2TCCPDF Density of crest-to-crest wave -period or -length.  
%         
%  CALL: f   = spec2tccpdf(S,u,def,paramt,h1,h2,nit,speed);
%
%        f    = density structures of wave periods of waves.
%        S    = spectral density structure.
%        u    = reference level (default the most frequently crossed level).
%       def   =  't<', period for waves with Ac<h1, At<h2, (default)
%                't>', period for waves with Ac>h1, At>h2.
%                'L<' and 'L>' ditto for wave length
%     paramt  = [t0 tn Nt] where t0, tn and Nt are the first value, last value 
%               and the number of points, respectively, for which
%               the density will be computed. paramt= [5 5 51] implies
%               that the density is computed only for T=5 and
%               using 51 linearly spaced points in the interval [0,5].
%          h1 = height of Ac crest; note only one positive value h1 > 0,
%          h2 = height of At trough; note only one positive value h2 > 0.
%        nit  =  0,...,9. Dimension of numerical integration  (default 2),
%               -1,-2,... different important sampling type integrations.
%       speed = defines accuraccy of calculations, by choosing different. 
%               parameters, possible values: 1,2,...,9 (9 fastest, default 4).
%       []    = default values are used.
%
% SPEC2TCCPDF calculates pdf of Tcc period  or Lcc length (i.e., crest to crest 
% wave -period or -length) such that Ac>h1 and At>h2, or Ac<h1 and At<h2 in 
% a stationary Gaussian transform process X(t), where Y(t) = g(X(t)).
% (Y zero-mean Gaussian with spectrum given in S). The transformation, g, 
% can be estimated using lc2tr or dat2tr.
%
% Example: %Compute the pdf of Tcc with crest and trough higher then 0.5*Hs:
%    S=jonswap; L=spec2mom(S); NIT=0;
%    f = spec2tccpdf(S,[],'t>',[],[2*sqrt(L(1))],[2*sqrt(L(1))],NIT);
%
% See also  spec2cov2, specnorm, dat2tr, dat2gaus, datastructures, perioddef

% Tested on : matlab 5.3
% History: by I. Rychlik 01.10.1998 with name wave_t.m
% Revised by I. R.   13.10.1999
% Revised by Sylvie van Iseghem 10.11.1999- adding levels of crests and troughs
% Revised by I. R.   07.12.1999 The upper and lower bounds to  the density included.
% Revised by I. R.   27.12.1999 and changed name to spec2tccpdf.
% Revised by es 000322. Made call with directional spectrum possible.
% revised by IR removing error in transformation 29 VI 2000
% revised by IR addapting to MATLAB6.0 8 II 2001
% Revised by pab 
% replaced some code with call to spec2cov2  
% Revised pab Aug2007
% -replaced some code with call to writecov
startTime = clock;

if nargin<3||isempty(def)
  def='t>';    
else
  if strcmpi('l',(def(1)))
    spec=spec2spec(spec,'k1d');
  elseif strcmpi('t',(def(1)))
    spec=spec2spec(spec,'freq');
  else
    error('Unknown def')
  end
  
  switch lower(def)
   case  {'l<','t<'},    defnr = 1; 
   case  {'l>','t>'},    defnr =-1;
   otherwise, error('Unknown def')
  end
end


[S, xl4, L0, L2]=specnorm(spec);
A=sqrt(L0/L2);
SCIS=0;
if nargin<7||isempty(nit)
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
  g=[sqrt(L0)*(-5:0.02:5)', (-5:0.02:5)'];
end


if nargin<2||isempty(utc)
  utc_d = gaus2dat([0, 0],g); % most frequent crossed level 
  utc   = utc_d(1,2);
end

% transform reference level into Gaussian level
uu = dat2gaus([0., utc],g);
u  = uu(2);
disp(['The level u for Gaussian process = ' num2str(u)])


if nargin<8||isempty(speed)
    speed=4;
end			    

if nargin<4||isempty(paramt)
  Ntime = 61;
  t0    = 0.;
  tn    = 2.*ceil(2*pi*sqrt(L0/L2)*exp(u^2/2));
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


 
px1=gaus2dat([0., u;1, 5],g); 
px1=abs(px1(2,2)-px1(1,2));
px2=gaus2dat([0., u;1, 5],g); 
px2=abs(px2(2,2)-px2(1,2));
Nx1 = 1;
Nx2 = 1;
 
if nargin<5||isempty(h1)
  if(defnr>0)
    h1=px1;
  else
    h1=0.;
  end
else
  h1=min(h1);
end  
if nargin<6||isempty(h2)
  if(defnr>0)
    h2=px2;
  else
    h2=0.;
  end
else
  h2=min(h2);
end      
h01=h1;
h02=h2;

if defnr<0 
  if h1<px1
    Nx1=2;
    h1=[h1; px1];
  else
    error('For option ">" the value of h1 is too large')
  end    
  if h2<px2
    Nx2=2;
    h2=[h2; px2];
  else
    error('For option ">" the value of h2 is too large')
  end    
end
%Nx=Nx1*Nx2;
%Transform amplitudes to Gaussian levels:   
der1=ones(Nx1,1); % dh/dh=1
hg1=tranproc([utc+h1, der1],g);
%der1=abs(hg1(:,2));
hg1=hg1(:,1); % Gaussian level
der2=ones(Nx2,1); % dh/dh=1
hg2=tranproc([utc-h2, der2],g);
%der2=abs(hg2(:,2));
hg2=hg2(:,1); % Gaussian level

nr = 2; 
dt = t(2)-t(1);
R = spec2cov2(S,nr,Ntime-1,dt);


disp('writing data')
filenames0 = writecov(R);

filenames = {'h.in','reflev.in'};
cleanup(filenames{:}) 

fid=fopen('h.in','wt');
fprintf(fid,'%12.10f\n',hg1);
fprintf(fid,'%12.10f\n',hg2);
fclose(fid)

fid=fopen('reflev.in','wt');
fprintf(fid,'%12.10E \n',u);
fprintf(fid,'%2.0f \n',Ntime);
fprintf(fid,'%2.0f \n',N0);
fprintf(fid,'%2.0f \n',nit);
fprintf(fid,'%2.0f \n',speed);
fprintf(fid,'%2.0f \n',SCIS);
fprintf(fid,'%2.0f \n',10^9*rand);  % select a random seed for rind 
fprintf(fid,'%2.0f %2.0f\n',[Nx1 Nx2]);
fprintf(fid,'%12.10E \n',dt);
fclose(fid); 



disp('   Starting Fortran executable.')

dos([ wafoexepath 'cov2tccpdf.exe']);



f=createpdf;
 
switch lower(def)
 case  't>'
  Htxt = ['Density of Tcc with Ac>', num2str(h01), ' and At>',  num2str(h02)];
  xtxt = 'T [s]';
 case  'l>'
  Htxt = ['Density of Lcc with Ac>', num2str(h01), ' and At>',  num2str(h02)];
  xtxt = 'L [m]';
 case  't<'
  Htxt = ['Density of Tcc with Ac<', num2str(h01), ' and At<',  num2str(h02)];
  xtxt = 'T [s]';
 case  'l<'
  Htxt = ['Density of Lcc with Ac<', num2str(h01), ' and At<',  num2str(h02)];
  xtxt = 'L [m]';
end 

f.title   = Htxt;
f.labx{1} = xtxt;
t=t'*A;
f.x{1}=t;
f.nit=nit;
f.speed=speed;
f.SCIS=SCIS;
f.u=utc;  
tmp=loaddata('dens.out')/A;
 
if defnr==1    
  f.f=tmp(:,1);
end

if defnr==-1   
  fup=reshape(tmp(:,1),4,length(t));
  f.f=fup(4,:)+fup(1,:)-fup(2,:)-fup(3,:);
end

f.elapsedTime = etime(clock,startTime);

cleanup(filenames0{:},filenames{:}) 