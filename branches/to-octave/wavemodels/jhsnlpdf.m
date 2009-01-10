function [f,Hrms,Vrms] = jhsnlpdf(Hd,Scf,Hm0,Tp,gam,normalizedInput,condon)
%JHSNLPDF Joint (Scf,Hd) PDF for nonlinear waves with a JONSWAP spectra.
%
%  CALL: f = jhsnlpdf(Hd,Scf,Hm0,Tp,gam)
% 
%    f     = pdf evaluated at (Vcf,Hd)
%    Hd    = zero down crossing wave height [m]
%    Vcf   = crest front velocity    [m/s]
%    Hm0   = significant wave height [m]
%    Tp    = Spectral peak period    [s]
%    Gamma = Peakedness parameter of the JONSWAP spectrum
%
% JHSNLPDF approximates the joint PDF of (Scf, Hd), i.e., crest front
% steepness (2*pi*Ac/(g*Td*Tcf)) and wave height, for 2nd order nonlinear
% waves with a JONSWAP spectral density. The empirical parameters of the
% model is fitted by least squares to simulated (Scf,Hd) data for 70
% classes of Hm0 and Tp. Between 40000 and 220000 zero-downcrossing
% waves were simulated for each class.
% JHSNLPDF is restricted to the following range: 
% 0.5 < Hm0 < 12 and 3 < Tp < 20
%  
% NOTE: The size of f is the common size of the input arguments.
%  
% Example:
% Hm0 = 6;Tp = 8;
% h = linspace(0,4*Hm0/sqrt(2)); 
% s = linspace(0,6*1.25*Hm0/Tp^2,101);
% [S,H] = meshgrid(s,h); 
% f = jhsnlpdf(H,S,Hm0,Tp);
% contourf(s,h,f)  
%
% See also  jhsnlpdf2

% Reference  
% P. A. Brodtkorb (2004),  
% The Probability of Occurrence of Dangerous Wave Situations at Sea.
% Dr.Ing thesis, Norwegian University of Science and Technolgy, NTNU,
% Trondheim, Norway. 
  
% History
% Validated pab 7 april 2004
% By pab 17.01.2003

error(nargchk(3,7,nargin))


if (nargin < 7||isempty(condon)),  condon  = 0; end
if (nargin < 6||isempty(normalizedInput)),  normalizedInput  = 0;end

if (nargin < 4||isempty(Tp)),  Tp  = 8;end
if (nargin < 3||isempty(Hm0)), Hm0 = 6;end
if (nargin < 5||isempty(gam))
%   gam = getjonswappeakedness(Hm0,Tp);
end

gam = getjonswappeakedness(Hm0,Tp);
multipleSeaStates = any(numel(Hm0)>1||...
			numel(Tp) >1||...
			numel(gam)>1);
if multipleSeaStates
  [csize, Scf,Hd,Hm0,Tp,gam] = comnsize(Scf,Hd,Hm0,Tp,gam);
else
  [csize, Scf,Hd] = comnsize(Scf,Hd);
end
if any(isnan(csize))
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
    warning('WAFO:JHSNLPDF','Peakedness parameter less than 1. Must be larger than 1.')
  end
  gam(k)=1;
end
k1 = find(7<gam);
if any(k1)
  if displayWarning
  warning('WAFO:JHSNLPDF','Peakedness parameter larger than 7. The pdf returned is questionable')
  end
  gam(k1) = 7;
end

global JHSNLPAR
if isempty(JHSNLPAR)
  %JHSNLPAR = load('jhsnlpar.mat');
  JHSNLPAR = load('jhsnlpar21-Jul-2004.mat');
end
% Gamma distribution parameters as a function of e2 and h2
A11 = JHSNLPAR.A11s;
B11 = JHSNLPAR.B11s;
C1 = 1;



if normalizedInput,
  Hrms = 1;
  Vrms = 1;
else
  %Tm02 = Tp./(1.30301-0.01698*gam+0.12102./gam);
  %dev = 2e-5;
  c1 =[ 0.16183666835624   1.53691936441548   1.55852759524555];
  c2 =[ 0.15659478203944   1.15736959972513   1];
  Tm02 = Tp.*(polyval(c2,gam)./polyval(c1,gam));
  Hrms = Hm0/sqrt(2);
  Vrms = 1.25*Hm0./(Tm02.^2); % Erms
end

h = Hd./Hrms;
v = Scf./Vrms;
cSize = size(h); % common size

h2  = JHSNLPAR.h2(:);  % Hd/Hrms
Nh2 = length(h2);
method ='*cubic';
if 1,
  Tpp  = JHSNLPAR.Tp; % gamma
  Hm00 = JHSNLPAR.Hm0;
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
      A1(k) = exp(smooth(h2,interp3(E1,H1,H2,log(A11),Tpi,Hm0i,h2,method),...
			 1,h(k),1));
      B1(k) = exp(smooth(h2,interp3(E1,H1,H2,log(B11),Tpi,Hm0i,h2,method),...
			 1,h(k),1));
    end
  else
    Tpi  = repmat(Tp,[Nh2,1]);
    Hm0i = repmat(Hm0,[Nh2,1]);
    A1 = exp(smooth(h2,interp3(E1,H1,H2,log(A11),Tpi,Hm0i,h2,method),1,h,1));
    B1 = exp(smooth(h2,interp3(E1,H1,H2,log(B11),Tpi,Hm0i,h2,method),1,h,1));
  end
else % old call
  e2  = JHSNLPAR.gam(:); % gamma
  h2  = JHSNLPAR.h2(:);  % Hd/Hrms
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
      A1(k) = exp(smooth(h2,interp2(E2,H2,log(A11),gami,h2,method),1,h(k),1));
      B1(k) = exp(smooth(h2,interp2(E2,H2,log(B11),gami,h2,method),1,h(k),1));
    end
  else
    A1 = exp(smooth(h2,interp2(E2,H2,log(A11),...
			       gam(ones(size(h2))),h2,method),1,h,1));
    B1 = exp(smooth(h2,interp2(E2,H2,log(B11),...
			       gam(ones(size(h2))),h2,method),1,h,1));
  end
end
% Waveheight distribution in time
% Truncated Weibull  distribution parameters as a function of Tp, Hm0, gam 
[A0, B0, C0] = jhnlwparfun(Hm0,Tp,gam,'time');
switch condon,
 case 0, % regular pdf is returned 
  f = pdfweibmod(h,A0,B0,C0).*pdfgengam(v,A1,B1,C1);
 case 1, %pdf conditioned on x1 ie. p(x2|x1) 
  f = pdfgengam(v,A1,B1,C1);
 case 3, % secret option  used by XXstat: returns x2*p(x2|x1) 
  f = v.*pdfgengam(v,A1,B1,C1);
 case 4, % secret option  used by XXstat: returns x2.^2*p(x2|x1) 
  f = v.^2.*pdfgengam(v,A1,B1,C1);
 case 5, % p(h)*P(V|h) is returned special case used by thscdf
  f = pdfweibmod(h,A0,B0,C0).*cdfgengam(v,A1,B1,C1);
 case 6, % P(V|h) is returned special case used by thscdf
  f = cdfgengam(v,A1,B1,C1);
 case 7,% p(h)*(1-P(V|h)) is returned special case used by thscdf
  f = pdfweibmod(h,A0,B0,C0).*(1-cdfgengam(v,A1,B1,C1));
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

return




