function [f,Hrms,Vrms] = jhspdf(Hd,Scf,Hm0,Tp,gam,normalizedInput,condon)
%JHSPDF Joint (Scf,Hd) PDF for linear waves with JONSWAP spectra.
%
%  CALL: f = jhspdf(Hd,Scf,Hm0,Tp,gam)
% 
%    f  = pdf
%   Hd  = zero down crossing wave height
%   Scf = crest front steepness
%   Hm0 = significant wave height
%   Tp  = Spectral peak period 
%
% JHSPDF approximates the joint PDF of (Scf, Hd), i.e., crest front
% steepness (2*pi*Ac/(g*Td*Tcf)) and wave height, for a Gaussian process with a
% JONSWAP spectral density. The empirical parameters of the model is
% fitted by least squares to simulated (Scf,Hd) data for 13 classes of
% GAMMA. Between 47000 and 55000 zero-downcrossing waves were
% simulated for each class of GAMMA.
% JHSPDF is restricted to the following range for GAMMA: 
%  1 <= GAMMA <= 7 
%  
% NOTE: The size of f is the common size of the input arguments.
%  
% Example:
% Hm0 = 6;Tp = 8;
% h = linspace(0,4*Hm0/sqrt(2)); 
% s = linspace(0,6*1.25*Hm0/Tp^2,101);
% [S,H] = meshgrid(s,h); 
% f = jhspdf(H,S,Hm0,Tp);
% contourf(s,h,f)  
%
% See also  jhspdf2

  
% Reference  
% P. A. Brodtkorb (2004),  
% The Probability of Occurrence of Dangerous Wave Situations at Sea.
% Dr.Ing thesis, Norwegian University of Science and Technolgy, NTNU,
% Trondheim, Norway.  

% Adapted to  cssmooth  by GL Feb 2011  
% History
% revised pab March 2007
% By pab 17.01.2003

error(nargchk(3,7,nargin))


if (nargin < 7||isempty(condon)),  condon  = 0; end
if (nargin < 6||isempty(normalizedInput)),  normalizedInput  = 0;end

if (nargin < 4||isempty(Tp)),  Tp  = 8;end
if (nargin < 3||isempty(Hm0)), Hm0 = 6;end
if (nargin < 5||isempty(gam))
   gam = getjonswappeakedness(Hm0,Tp);
end

multipleSeaStates = any(numel(Hm0)>1|...
			numel(Tp) >1|...
			numel(gam)>1);
if multipleSeaStates
  [icode, Scf,Hd,Hm0,Tp,gam] = iscomnsize(Scf,Hd,Hm0,Tp,gam);
else
  [icode, Scf,Hd] = iscomnsize(Scf,Hd);
end
if ~icode
    error('WAFO:JHSPDF','Requires non-scalar arguments to match in size.');
end
displayWarning = 0;
if displayWarning
  if any(any(Tp>5*sqrt(Hm0) | Tp<3.6*sqrt(Hm0)))
    warning('WAFO:JHSPDF','Hm0,Tp is outside the JONSWAP range. \nThe validity of the parameters returned are questionable.')
  end
end
k = find(gam<1);
if any(k) 
  if displayWarning,
    warning('WAFO:JHSPDF','Peakedness parameter less than 1. Must be larger than 1.')
  end
  gam(k)=1;
end
k1 = find(7<gam);
if any(k1)
  if displayWarning
  warning('WAFO:JHSPDF','Peakedness parameter larger than 7. The pdf returned is questionable')
  end
  gam(k1) = 7;
end

global JHSPAR
if isempty(JHSPAR)
  JHSPAR = load('jhspar18-Jan-2006.mat');
 %  JHSPAR = load('jhspar13-Jul-2004.mat');
%  JHSPAR = load('jhspar.mat');
end
%Generalized Gamma distribution parameters as a function of e2 and h2
A11 = JHSPAR.A11s;
B11 = JHSPAR.B11s;
C11 = 1.0;
%C11 = 1.5;
e2  = JHSPAR.gam(:); % gamma
h2  = JHSPAR.h2(:);  % Hd/Hrms
[E2 H2] = meshgrid(e2,h2);


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

Nh2 = length(h2);
method ='*cubic';
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
    A1(k) = exp(cssmooth(h2,interp2(E2,H2,log(A11),gami,h2,method),1,h(k),1));
    B1(k) = exp(cssmooth(h2,interp2(E2,H2,log(B11),gami,h2,method),1,h(k),1));
  end
else
  A1 = exp(cssmooth(h2,interp2(E2,H2,log(A11),...
			     gam(ones(size(h2))),h2,method),1,h,1));
  B1 = exp(cssmooth(h2,interp2(E2,H2,log(B11),...
			     gam(ones(size(h2))),h2,method),1,h,1));
end
% Waveheight distribution in time
% Truncated Weibull  distribution parameters as a function of Tp, Hm0, gam 
[A0, B0, C0] = jhwparfun(Hm0,Tp,gam,'time');
switch condon,
 case 0, % regular pdf is returned 
  f = pdfweibmod(h,A0,B0,C0).*pdfgengam(v,A1,B1,C11);
 case 1, %pdf conditioned on x1 ie. p(x2|x1) 
  f = pdfgengam(v,A1,B1,C11);
 case 3, % secret option  used by XXstat: returns x2*p(x2|x1) 
  f = v.*pdfgengam(v,A1,B1,C11);
 case 4, % secret option  used by XXstat: returns x2.^2*p(x2|x1) 
  f = v.^2.*pdfgengam(v,A1,B1,C11);
 case 5, % p(h)*P(V|h) is returned special case used by thscdf
  f = pdfweibmod(h,A0,B0,C0).*cdfgengam(v,A1,B1,C11);
 case 6, % P(V|h) is returned special case used by thscdf
  f = cdfgengam(v,A1,B1,C11);
 case 7,% p(h)*(1-P(V|h)) is returned special case used by thscdf
  f = pdfweibmod(h,A0,B0,C0).*(1-cdfgengam(v,A1,B1,C11));
  otherwise
    error('WAFO:JHSPDF','unknown option')
end
if multipleSeaStates
  f = reshape(f,cSize);
end

if condon~=6
  f = f./Hrms./Vrms;
end
f((isnan(f)|isinf(f) ))=0;
if any(size(f)~=cSize)
  warning('WAFO:JHSPDF','Wrong size')
end

return



