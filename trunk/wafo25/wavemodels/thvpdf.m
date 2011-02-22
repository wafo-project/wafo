function [f,Hrms,Vrms,fA,fB] = thvpdf(Hd,Vcf,Hm0,Tp,normalizedInput,condon)
%THVPDF Joint (Vcf,Hd) PDF for linear waves with Torsethaugen spectra. 
%
%  CALL: f = thvpdf(Hd,Vcf,Hm0,Tp)
% 
%  f   = pdf evaluated at (Hd,Vcf).
%  Hd  = zero down crossing wave height
%  Vcf = crest front velocity
%  Hm0 = significant wave height [m]. (default Hm0 = 6)
%  Tp  = Spectral peak period    [s]. (default  Tp = 8)
%
% THVPDF approximates the joint PDF of (Vcf, Hd), i.e., crest
% front velocity (Ac/Tcf) and wave height, for a Gaussian process with a
% Torsethaugen spectral density. The empirical parameters of the model is
% fitted by least squares to simulated (Vcf,Hd) data for 600 classes of
% Hm0 and Tp. Between 50000 and 150000 zero-downcrossing waves were
% simulated for each class of Hm0 and Tp.
% THVPDF is restricted to the following range for Hm0 and Tp: 
% 0.5 < Hm0 [m] < 12,  3.5 < Tp [s] < 20,  and  Hm0 < (Tp-2)*12/11.
%  
% NOTE: The size of f is the common size of the input arguments.
%
% Example:
% Hm0 = 6;Tp = 8;
% h = linspace(0,4*Hm0/sqrt(2))'; 
% v = linspace(0,4*2*Hm0/Tp)';
% [V,H] = meshgrid(v,h);  
% f = thvpdf(H,V,Hm0,Tp);
% contourf(v,h,f)  
%
% See also  thspdf, thsspdf

% Reference  
% P. A. Brodtkorb (2004),  
% The Probability of Occurrence of Dangerous Wave Situations at Sea.
% Dr.Ing thesis, Norwegian University of Science and Technolgy, NTNU,
% Trondheim, Norway.   

% Adapted to  cssmooth  by GL Feb 2011  
% History
% revised pab 10.08.2003  
% -revised pab 28.11.2002
%   extended the example
% By pab 20.12.2000

error(nargchk(3,6,nargin))

if (nargin < 6||isempty(condon)),  condon  = 0; end
if (nargin < 5||isempty(normalizedInput)),  normalizedInput  = 0;end
if (nargin < 4||isempty(Tp)),  Tp  = 8;end
if (nargin < 3||isempty(Hm0)), Hm0 = 6;end

multipleSeaStates = any(numel(Hm0)>1|numel(Tp)>1);
if multipleSeaStates
  [icode, Vcf,Hd,Hm0,Tp] = iscomnsize(Vcf,Hd,Hm0,Tp);
else
  [icode, Vcf,Hd] = iscomnsize(Vcf,Hd);
end
if ~icode 
  error('Requires non-scalar arguments to match in size.');
end
displayWarning = 0;
if displayWarning,
  if any(Hm0>11| Hm0>(Tp-2)*12/11) 
    disp('Warning: Hm0 is outside the valid range')
    disp('The validity of the Joint (Hd,Vcf) distribution is questionable')
  end
  if any(Tp>20|Tp<3)
    disp('Warning: Tp is outside the valid range')
    disp('The validity of the Joint (Hd,Vcf) distribution is questionable')
  end
end

global THVPAR
if isempty(THVPAR)
  THVPAR = load('thvpar.mat');
end
% Weibull distribution parameters as a function of e2 and h2
A11 = THVPAR.A11s;
B11 = THVPAR.B11s;
e2  = THVPAR.e2;
h2  = THVPAR.h2;
Tpp   = THVPAR.Tp;
Hm00  = THVPAR.Hm0;
Tm020 = THVPAR.Tm02;
eps20 = THVPAR.eps2;

[Tp1,Hs1] = meshgrid(Tpp,Hm00);
method = '*cubic'; %'spline'

eps2 = interp2(Tp1,Hs1,eps20,Tp,Hm0,method);

if normalizedInput,
  Hrms = 1;
  Vrms = 1;
else  
  Tm02 = interp2(Tp1,Hs1,Tm020,Tp,Hm0,method);
%  w    = linspace(0,100,16*1024+1).'; % torsethaugen original spacing
%  St = torsethaugen(w,[Hm0,Tp]);
%  ch   = spec2char(St,{'Tm02','eps2'});
%  Tm02 = ch(1);
%  eps2 = ch(2);
  Hrms = Hm0/sqrt(2);
  Vrms = 2*Hm0./Tm02; % Erms
end
  
%Fh = thpdf(h(:)/Hrms,Hm0,Tp,eps2,1);
if displayWarning
  if eps2<0.4 
    disp('Warning: eps2 is less than 0.4. The pdf returned is questionable')
  elseif eps2>1.3
    disp('Warning: eps2 is larger than 1.3. The pdf returned is questionable')
  end
end
h = Hd/Hrms;
v = Vcf/Vrms;
cSize = size(h); % common size of input

[E2 H2] = meshgrid(e2,h2);
Nh2 = length(h2);
if multipleSeaStates
  h   = h(:);
  v   = v(:);
  Tp  = Tp(:);
  Hm0 = Hm0(:);
  eps2 = eps2(:);
  A1 = zeros(length(h),1);
  B1 = A1;
  [eps2u,ix,jx] = unique(eps2);
  numSeaStates = length(ix);
  eps2i = zeros(Nh2,1);
  for iz=1:numSeaStates
    k = find(jx==iz);
    eps2i(:)  = eps2u(iz);
    A1(k) = cssmooth(h2,interp2(E2,H2,A11,eps2i,h2,method),1,h(k),1);
    B1(k) = cssmooth(h2,interp2(E2,H2,B11,eps2i,h2,method),1,h(k),1);
  end
else
  eps2i = repmat(eps2,[Nh2,1]);
  A1 = cssmooth(h2,interp2(E2,H2,A11,eps2i,h2,method),1,h,1);
  B1 = cssmooth(h2,interp2(E2,H2,B11,eps2i,h2,method),1,h,1);
end
% Note if eps2<0.4 then B1 is questionable

% Waveheight distribution in time
% Truncated Weibull  distribution parameters as a function of Tp, Hm0 
[A0, B0, C0] = thwparfun(Hm0,Tp,'time');

switch condon,
 case 0, % regular pdf is returned 
  f = pdfweibmod(h,A0,B0,C0).*pdfweib(v,A1,B1);
 case 1, %pdf conditioned on x1 ie. p(x2|x1) 
  f = pdfweib(v,A1,B1);
 case 3, % secret option  used by XXstat: returns x2*p(x2|x1) 
  f = v.*pdfweib(v,A1,B1);
 case 4, % secret option  used by XXstat: returns x2.^2*p(x2|x1) 
  f = v.^2.*pdfweib(v,A1,B1);
 case 5, % p(h)*P(V|h) is returned special case used by thvcdf
  f = pdfweibmod(h,A0,B0,C0).*cdfweib(v,A1,B1);
 case 6, % P(V|h) is returned special case used by thvcdf
  f = cdfweib(v,A1,B1);
 case 7,% p(h)*(1-P(V|h)) is returned special case used by thvcdf
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
  fA.labx = {'eps2', 'h'};
  fA.title = 'pdfweib parameter A';
  fA.note = ['The conditonal Weibull distribution Parameter A(h,eps2)/Hrms' ...
	'for Vcf given h=Hd/Hrms and eps2 for' ...
	'the Torsethaugen spectrum'];
    
  ra = range(A11(:));
  st = round(min(A11(:))*100)/100;
  en = max(A11(:));
  fA.cl   = st:ra/20:en;
end
if nargout>4,
  fB      = createpdf(2);
  fB.f    = B11;
  fB.x    = {e2,h2};
  fB.labx = {'eps2', 'h'};
  fB.title = 'pdfweib parameter B';
  fB.note = ['The conditonal Weibull distribution Parameter B(h,eps2)/Hrms' ...
	'for Vcf given h=Hd/Hrms and eps2 for' ...
	'the Torsethaugen spectrum'];
  ra = range(B11(:));
  st = round(min(B11(:))*100)/100;
  en = max(B11(:));
  fB.cl   = st:ra/20:en;
end
return






