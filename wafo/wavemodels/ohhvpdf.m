function [f,fA,fB] = ohhvpdf(v,h,Hm0,def,norm)
%OHHVPDF Joint (Vcf,Hd) PDF for linear waves with Ochi-Hubble spectra. 
%
%  CALL: f = ohhvpdf(v,h,Hm0,Tp)
% 
%  f   = pdf evaluated at meshgrid(v,h).
%  v,h = vectors of evaluation points.
%  Hm0 = significant wave height [m].
%  def = defines the parametrization of the spectral density (default 1)
%        1 : The most probable spectrum  (default)
%        2,3,...11 : gives 95% Confidence spectra
%
% OHHVPDF approximates the joint distribution of (Vcf, Hd), i.e., crest
% front velocity and wave height, for a Gaussian process with a bimodal
% Ochi-Hubble spectral density. The empirical parameters of the model is
% fitted by least squares to simulated (Vcf,Hd) data for 24 classes of Hm0.
% Between 50000 and 150000 zero-downcrossing waves were simulated for
% each class of Hm0.
% OHHVPDF is restricted to the following range for Hm0: 
%  0.5 < Hm0 [m] < 12
%
% Example:
% Hm0 = 6;Tp = 8;def= 2;
% h = linspace(0,4*Hm0/sqrt(2))'; 
% v = linspace(0,4*2*Hm0/Tp)';
% f = ohhvpdf(v,h,Hm0,def);
% w    = linspace(0,40,5*1024+1).';
% S = ochihubble(w,[Hm0 def]);
% dt = 0.3;
% x = spec2sdat(S,80000,.2); rate = 8;
% [vi,hi] = dat2steep(x,rate,1);
% fk = kdebin([vi,hi],{'L2',.5,'inc',128});
% fk.title = f.title; fk.labx = f.labx; 
% plot(vi,hi,'.'), hold on
% pdfplot(f),pdfplot(fk,'r'),hold off  
%
% See also  ohhpdf, thvpdf

% Reference  
% P. A. Brodtkorb (2004),  
% The Probability of Occurrence of Dangerous Wave Situations at Sea.
% Dr.Ing thesis, Norwegian University of Science and Technolgy, NTNU,
% Trondheim, Norway.   

% Adapted to  cssmooth  by GL Feb 2011    
% History
% By pab 20.12.2000

%error(nargchk(3,5,nargin))
narginchk(3,5)
if nargin<5||isempty(norm), norm=0;end 
if nargin<4||isempty(def), def=1;end 

if any(Hm0>12| Hm0<=0) 
  disp('Warning: Hm0 is outside the valid range')
  disp('The validity of the Hd distribution is questionable')
end

if def>11||def<1 
  warning('DEF is outside the valid range')
  def = mod(def-1,11)+1;
end

global OHHVPAR
if isempty(OHHVPAR)
  OHHVPAR = load('ohhvpar.mat');
end
% Weibull distribution parameters as a function of Hm0 and h2
A11  = OHHVPAR.A11s;
B11  = OHHVPAR.B11s;
Hm00 = OHHVPAR.Hm0;
h2   = OHHVPAR.h2;
[E2 H2] = meshgrid(Hm00,h2);


Hrms = Hm0/sqrt(2);
w    = linspace(0,100,16*1024+1).'; %original spacing
%w    = linspace(0,1000/Tp,8*1024+1).'; % adaptiv spacing
ch   = spec2char(ochihubble(w,[Hm0,def]),{'Tm02','eps2'});

Tm02 = ch(1);
%eps2 = ch(2);
Vrms = 2*Hm0/Tm02;
if nargin<1 ||isempty(v), v=linspace(0,4*Vrms); end
if nargin<2 ||isempty(h), h=linspace(0,4*Hrms); end

if norm==1,
  h = h*Hrms;
  v = v*Vrms;
end
  
%fh = ohhpdf(h(:)/Hrms,Hm0,def,'time',1);
% Fh = fh.f;
[A0 B0 C0] = ohhgparfun(Hm0,def,'time');
Fh = pdfgengam(h(:)/Hrms,A0,B0,C0);

method = '*cubic'; %'spline'
A1 = exp(cssmooth(h2,interp2(E2,H2,log(squeeze(A11(:,def,:))).',Hm0(ones(size(h2))),h2,method),1,h/Hrms,1));
B1 = cssmooth(h2,interp2(E2,H2,squeeze(B11(:,def,:)).',Hm0(ones(size(h2))),h2,method),1,h/Hrms,1);

[V1 A1] = meshgrid(v/Vrms,A1);
[V1 B1] = meshgrid(v/Vrms,B1);

f = createpdf(2);
f.title = 'Joint distribution of (Hd,Vcf) in time';

if norm==0
  f.f = repmat(Fh/Hrms,[1 length(v)]).*pdfweib(V1,A1,B1)/Vrms;
  f.x = {v,h};
  f.norm=0;
  f.labx={'Vcf [m/s]', 'Hd [m]'};
else
  f.f = repmat(Fh,[1 length(v)]).*pdfweib(V1,A1,B1);
  f.x = {v/Vrms,h/Hrms};
  f.labx={'Vcf', 'Hd'};
  f.norm = 1;
end
f.note = ['Ochi-Hubble Hm0=' num2str(Hm0) ' def = ' num2str(def)];
f.f(isinf(f.f)|isnan(f.f))=0;
[f.cl,f.pl] = qlevels(f.f);%,[10:20:90 95 99 99.9 99.99 99.999 99.9999]);
if nargout>1,
  fA      = createpdf(2);
  fA.f    = A11;
  fA.x    = {Hm00,(1:11)' h2};
  fA.labx = {'Hm0','def', 'h'};
  fA.note = ['The conditonal Weibull distribution Parameter A(h)/Hrms' ...
	'for wave heigth as a function of h=Hd/Hrms, def and Hm0 for' ...
	'the ochihubble spectrum'];
    
  ra = range(A11(:));
  st = round(min(A11(:))*100)/100;
  en = max(A11(:));
  fA.cl   = st:ra/20:en;
end
if nargout>2,
  fB      = createpdf(2);
  fB.f    = B11;
  fB.x    = {Hm00,(1:11)',h2};
  fB.labx = {'Hm0','def', 'h'};
  fB.note = ['The conditonal Weibull distribution Parameter B(h)/Hrms' ...
	'for wave heigth as a function of h=Hd/Hrms, def and Hm0 for' ...
	'the ochihubble spectrum'];
  ra = range(B11(:));
  st = round(min(B11(:))*100)/100;
  en = max(B11(:));
  fB.cl   = st:ra/20:en;
end
return






