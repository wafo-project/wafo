function [f,Hrms,Vrms,fA,fB] = jhvnlpdf(Hd,Vcf,Hm0,Tp,gam,normalizedInput,condon)
%JHVNLPDF Joint (Vcf,Hd) PDF for linear waves with a JONSWAP spectrum.
%
%  CALL: f = jhvnlpdf(Hd,Vcf,Hm0,Tp,gamma)
% 
%  f     = pdf evaluated at (Vcf,Hd)
%  Hd    = zero down crossing wave height [m]
%  Vcf   = crest front velocity    [m/s]
%  Hm0   = significant wave height [m]
%  Tp    = Spectral peak period    [s]
%  Gamma = Peakedness parameter of the JONSWAP spectrum
%
% JHVNLPDF approximates the joint distribution of (Vcf, Hd), i.e., crest
% front velocity (Ac/Tcf) and wave height, for 2nd order Stokes waves with a
% JONSWAP spectral density. The empirical parameters of the model is
% fitted by least squares to simulated (Vcf,Hd) data for 13 classes of
% GAMMA between 1 and 7. Between 140000 and 150000 zero-downcrossing waves
% were simulated for each class.
% JHVNLPDF is restricted to the following range for GAMMA: 
%  1 <= GAMMA <= 7 
%
% Example:
% Hm0 = 6;Tp = 9; gam=3.5
% h = linspace(0,4*Hm0/sqrt(2))'; 
% v = linspace(0,4*2*Hm0/Tp)';
% [V,H] = meshgrid(v,h); 
% f = jhvnlpdf(H,V,Hm0,Tp,gam);
% contourf(v,h,f)    
%
% See also  thvpdf

% Reference  
% P. A. Brodtkorb (2004),  
% The Probability of Occurrence of Dangerous Wave Situations at Sea.
% Dr.Ing thesis, Norwegian University of Science and Technolgy, NTNU,
% Trondheim, Norway. 

% Adapted to  cssmooth  by GL Feb 2011  
% History
% By pab 20.12.2000

error(nargchk(4,7,nargin))
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

global JHVNLPAR
if isempty(JHVNLPAR)
  JHVNLPAR = load('jhvnlpar.mat');
end
% Weibull distribution parameters as a function of e2 and h2
A11 = JHVNLPAR.A11s;
B11 = JHVNLPAR.B11s;




%dev = 2e-5;
c1 =[ 0.16183666835624   1.53691936441548   1.55852759524555];
c2 =[ 0.15659478203944   1.15736959972513   1];
Tm02 = Tp.*(polyval(c2,gam)./polyval(c1,gam));

%w    = linspace(0,10.5,16*1024+1).'; % jonswap original spacing
%ch   = spec2char(jonswap(w,[Hm0,Tp gam]),{'Tm02','eps2'});
%Tm02 = ch(1);
%eps2 = ch(2);


if normalizedInput
  Hrms = 1;
  Vrms = 1;
else
  Hrms = Hm0/sqrt(2);
  Vrms = 2*Hm0./Tm02; % Erms
end

v = Vcf./Vrms;
h = Hd./Hrms;
cSize = size(h); 
 
if gam<1 
  disp('Warning: Peakedness parameter less than 1. Must be larger than 1.')
    gam=1;
elseif gam>7
  disp('Warning: Peakedness parameter larger than 7. The pdf returned is questionable')
end 


h2  = JHVNLPAR.h2(:);  % Hd/Hrms
method ='*cubic';
Nh2 = length(h2);
if 1,
  
  Tpp  = JHVNLPAR.Tp; % gamma
  Hm00 = JHVNLPAR.Hm0;
  [E1, H1, H2] = meshgrid(Tpp,Hm00,h2);
  Nh2 = length(h2);
  if multipleSeaStates
    h   = h(:);
    v   = v(:);
    Tp  = Tp(:);
    Hm0 = Hm0(:);
    gam = gam(:);
    A1 = zeros(length(h),1);
    B1 = A1;
    [TpHm0,ix,jx] = unique([Tp,Hm0],'rows');
    numSeaStates = length(ix);
    Tpi = zeros(Nh2,1);
    Hm0i = zeros(Nh2,1);
    for iz=1:numSeaStates
      k = find(jx==iz);
      Tpi(:)  = TpHm0(iz,1);
      Hm0i(:) = TpHm0(iz,2);
      A1(k) = exp(cssmooth(h2,interp3(E1,H1,H2,log(A11),Tpi,Hm0i,h2,method),...
			 1,h(k),1));
      B1(k) = exp(cssmooth(h2,interp3(E1,H1,H2,log(B11),Tpi,Hm0i,h2,method),...
			 1,h(k),1));
    end
  else
    Tpi  = repmat(Tp,[Nh2,1]);
    Hm0i = repmat(Hm0,[Nh2,1]);
    A1 = exp(cssmooth(h2,interp3(E1,H1,H2,log(A11),Tpi,Hm0i,h2,method),1,h,1));
    B1 = exp(cssmooth(h2,interp3(E1,H1,H2,log(B11),Tpi,Hm0i,h2,method),1,h,1));
  end
else
  e2  = JHVNLPAR.gam; % gamma
  
  [E2 H2] = meshgrid(e2,h2);
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
end
[A0,B0,C0] = jhnlwparfun(Hm0,Tp,gam);

switch condon,
 case 0, % regular pdf is returned 
  f = pdfweibmod(h,A0,B0,C0).*pdfweib(v,A1,B1);
 case 1, %pdf conditioned on x1 ie. p(x2|x1) 
  f = pdfweib(v,A1,B1);
 case 3, % secret option  used by XXstat: returns x2*p(x2|x1) 
  f = v.*pdfweib(v,A1,B1);
 case 4, % secret option  used by XXstat: returns x2.^2*p(x2|x1) 
  f = v.^2.*pdfweib(v,A1,B1);
 case 5, % p(h)*P(V|h) is returned special case used by jhvnlcdf2
  f = pdfweibmod(h,A0,B0,C0).*cdfweib(v,A1,B1);
 case 6, % P(V|h) is returned special case used by jhvnlcdf
  f = cdfweib(v,A1,B1);
 case 7,% p(h)*(1-P(V|h)) is returned special case used by jhvnlcdf
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
	'for wave heigth as a function of h=Hd/Hrms and gamma for' ...
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
	'for wave heigth as a function of h=Hd/Hrms and gamma for' ...
	'the Jonswap spectrum'];
  ra = range(B11(:));
  st = round(min(B11(:))*100)/100;
  en = max(B11(:));
  fB.cl   = st:ra/20:en;
end
return






