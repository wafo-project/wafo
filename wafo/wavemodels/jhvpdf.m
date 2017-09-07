function [f,Hrms,Vrms,fA,fB] = jhvpdf(Hd,Vcf,Hm0,Tp,gam,normalizedInput,condon)
%JHVPDF Joint (Vcf,Hd) PDF for linear waves with a JONSWAP spectrum.
%
%  CALL: f = jhvpdf(Hd,Vcf,Hm0,Tp,gamma)
% 
%  f     = pdf evaluated at (Vcf,Hd)
%  Hd    = zero down crossing wave height [m]
%  Vcf   = crest front velocity    [m/s]
%  Hm0   = significant wave height [m]
%  Tp    = Spectral peak period    [s]
%  Gamma = Peakedness parameter of the JONSWAP spectrum
%
% JHVPDF approximates the joint distribution of (Vcf, Hd), i.e., crest
% front velocity (Ac/Tcf) and wave height, for a Gaussian process with a
% JONSWAP spectral density. The empirical parameters of the model is
% fitted by least squares to simulated (Vcf,Hd) data for 13 classes of
% GAMMA between 1 and 7. About 100000 zero-downcrossing waves
% were simulated for each class.
% JHVPDF is restricted to the following range for GAMMA: 
%  1 <= GAMMA <= 7 
%
% NOTE:  The size of f is the common size of the input arguments. 
%
% Example:
% Hm0 = 6;Tp = 9; gam=3.5
% h = linspace(0,4*Hm0/sqrt(2))'; 
% v = linspace(0,4*2*Hm0/Tp)';
% [V,H] = meshgrid(v,h);  
% f = jhvpdf(H,V,Hm0,Tp,gam);
% w = linspace(0,40,5*1024+1).';
% S = jonswap(w,[Hm0, Tp, gam]);
% dt = .3;
% x = spec2sdat(S,80000,dt); rate = 4;
% [vi,hi] = dat2steep(x,rate,1);
% fk = kdebin([vi,hi],{'L2',.5,'inc',128}); 
% plot(vi,hi,'.'), hold on
% contour(v,h,f,fk.cl),
% pdfplot(fk,'r'),hold off
%
% See also  thvpdf

% Reference  
% P. A. Brodtkorb (2004),  
% The Probability of Occurrence of Dangerous Wave Situations at Sea.
% Dr.Ing thesis, Norwegian University of Science and Technolgy, NTNU,
% Trondheim, Norway. 

% Adapted to  cssmooth  by GL Feb 2011  
% History
% revised pab 10 Jan 2004  
% By pab 20.12.2000

%error(nargchk(4,7,nargin))
narginchk(4,7)
if (nargin < 7||isempty(condon)),  condon  = 0; end
if (nargin < 6||isempty(normalizedInput)),  normalizedInput  = 0;end
if (nargin < 5||isempty(gam))
   gam = getjonswappeakedness(Hm0,Tp);
end

multipleSeaStates = any(numel(Hm0)>1|...
			numel(Tp) >1|...
			numel(gam)>1);
if multipleSeaStates
  [icode, Vcf,Hd,Hm0,Tp,gam] = iscomnsize(Vcf,Hd,Hm0,Tp,gam);
else
  [icode, Vcf,Hd] = iscomnsize(Vcf,Hd);
end
if ~icode 
  error('Requires non-scalar arguments to match in size.');
end

displayWarning = 0;
if displayWarning
  if any(any(Tp>5*sqrt(Hm0) | Tp<3.6*sqrt(Hm0)))
    disp('Warning: Hm0,Tp is outside the JONSWAP range')
    disp('The validity of the parameters returned are questionable')
  end
end
k = find(gam<1);
if any(k) 
  if displayWarning,
    warning('Peakedness parameter less than 1. Must be larger than 1.')
  end
  gam(k)=1;
end
k1 = find(7<gam);
if any(k1)
  if displayWarning
  warning('Peakedness parameter larger than 7. The pdf returned is questionable')
  end
  gam(k1) = 7;
end

global JHVPAR
if isempty(JHVPAR)
  JHVPAR = load('jhvpar.mat');
end
% Weibull distribution parameters as a function of e2 and h2
A11 = JHVPAR.A11s;
B11 = JHVPAR.B11s;
e2  = JHVPAR.gam; % gamma
h2  = JHVPAR.h2;  % Hd/Hrms
[E2 H2] = meshgrid(e2,h2);

if 1,
  
  %Tm02 = Tp./(1.30301-0.01698*gam+0.12102./gam);
  %dev = 2e-5;
  c1 =[ 0.16183666835624   1.53691936441548   1.55852759524555];
  c2 =[ 0.15659478203944   1.15736959972513   1];
  Tm02 = Tp.*(polyval(c2,gam)./polyval(c1,gam));
  
   %dev = 2e-4;
  %EPS2cof = [0.00068263671017  -0.01802256231624   0.44176198490431];
  %eps2 = polyval(EPS2cof,gam);
else
  w    = linspace(0,100,16*1024+1).'; % jonswap original spacing
  %Hm0 = 6;
  %gam = linspace(1,7,32);
  Tm02 = zeros(size(gam));
  eps2 = Tm02;
  for ix=1:length(gam(:))
    ch   = spec2char(jonswap(w,[Hm0(ix),Tp(ix) gam(ix)]),{'Tm02','eps2'});
    Tm02(ix) = ch(1);
    eps2(ix) = ch(2);
  end
end


if normalizedInput
  Hrms = 1;
  Vrms = 1;
else
  Hrms = Hm0/sqrt(2);
  Vrms = 2*Hm0./Tm02; % Erms
end

v = Vcf./Vrms;
h = Hd./Hrms;
cSize = size(h); % common size of input



method ='*cubic';
Nh2 = length(h2);
if multipleSeaStates
  h   = h(:);
  v   = v(:);
  Tp  = Tp(:);
  Hm0 = Hm0(:);
  gam = gam(:);
%  eps2 = eps2(:);
  A1 = zeros(length(h),1);
  B1 = A1;
  [gamu,ix,jx] = unique(gam);
  numSeaStates = length(gamu);
  gami = zeros(Nh2,1);
  for iz=1:numSeaStates
    k = find(jx==iz);
    gami(:) = gamu(iz);
    A1(k) = cssmooth(h2,interp2(E2,H2,A11,gami,h2,method),1,h(k),1);
    B1(k) = cssmooth(h2,interp2(E2,H2,B11,gami,h2,method),1,h(k),1);
  end
else
  A1 = cssmooth(h2,interp2(E2,H2,A11,gam(ones(size(h2))),h2,method),1,h,1);
  B1 = cssmooth(h2,interp2(E2,H2,B11,gam(ones(size(h2))),h2,method),1,h,1);
end

% Waveheight distribution in time
% Truncated Weibull  distribution parameters as a function of Tp, Hm0, gam 
[A0, B0, C0] = jhwparfun(Hm0,Tp,gam,'time');

switch condon,
 case 0, % regular pdf is returned 
  f = pdfweibmod(h,A0,B0,C0).*pdfweib(v,A1,B1);
 case 1, %pdf conditioned on x1 ie. p(x2|x1) 
  f = pdfweib(v,A1,B1);
 case 3, % secret option  used by XXstat: returns x2*p(x2|x1) 
  f = v.*pdfweib(v,A1,B1);
 case 4, % secret option  used by XXstat: returns x2.^2*p(x2|x1) 
  f = v.^2.*pdfweib(v,A1,B1);
 case 5, % p(h)*P(V|h) is returned special case used by jhvcdf2
  f = pdfweibmod(h,A0,B0,C0).*cdfweib(v,A1,B1);
 case 6, % P(V|h) is returned special case used by jhvcdf2
  f = cdfweib(v,A1,B1);
 case 7,% p(h)*(1-P(V|h)) is returned special case used by jhvcdf2
  f = pdfweibmod(h,A0,B0,C0).*(1-cdfweib(v,A1,B1));
  otherwise, error('unknown option')
end

if multipleSeaStates
  f = reshape(f,cSize);
end

if condon~=6
  f = f./Hrms./Vrms;
end
f((isnan(f)|isinf(f) ))=0;
if any(size(f)~=cSize)
  disp('Wrong size')
end


if nargout>3,
  fA      = createpdf(2);
  fA.f    = A11;
  fA.x    = {e2,h2};
  fA.labx = {'Gamma', 'h'};
  fA.note = ['The conditonal Weibull distribution Parameter A(h,gamma)/Hrms' ...
	'for Vcf as a function of h=Hd/Hrms and gamma for' ...
	'the Jonswap spectrum'];
    
  ra = range(A11(:));
  st = round(min(A11(:))*100)/100;
  en = max(A11(:));
  fA.cl   = st:ra/20:en;
end
if nargout>4,
  fB      = createpdf(2);
  fB.f    = B11;
  fB.x    = {e2,h2};
  fB.labx = {'Gamma', 'h'};
  fB.note = ['The conditonal Weibull distribution Parameter B(h,gamma)/Hrms' ...
	'for Vcf as a function of h=Hd/Hrms and gamma for' ...
	'the Jonswap spectrum'];
  ra = range(B11(:));
  st = round(min(B11(:))*100)/100;
  en = max(B11(:));
  fB.cl   = st:ra/20:en;
end
return






