function [f,Hrms,Vrms] = thsnlpdf(Hd,Scf,Hm0,Tp,normalizedInput,condon)
%THSNLPDF Joint (Scf,Hd) PDF for nonlinear waves with Torsethaugen spectra.
%
%  CALL: f = thsnlpdf(Hd,Scf,Hm0,Tp)
% 
%    f  = pdf
%   Hd  = zero down crossing wave height
%   Scf = crest front steepness
%   Hm0 = significant wave height
%   Tp  = Spectral peak period 
%
% THSNLPDF approximates the joint PDF of (Scf, Hd), i.e., crest
% steepness (2*pi*Ac/(g*Td*Tcf)) and wave height, for for 2nd order
% non-linear waves with a Torsethaugen spectral density. The empirical
% parameters of the model is
% fitted by least squares to simulated (Scf,Hd) data for 600 classes of
% Hm0 and Tp. Between 40000 and 200000 zero-downcrossing waves were
% simulated for each class of Hm0 and Tp.
% THSNLPDF is restricted to the following range for Hm0 and Tp: 
%  0.5 < Hm0 [m] < 12,  3.5 < Tp [s] < 20,  and  Hm0 < (Tp-2)*12/11.
%  The size of f is the common size of the input arguments.
%  
% Example:
% Hm0 = 6;Tp = 8;
% h = linspace(0,4*Hm0/sqrt(2)); 
% s = linspace(0,6*1.25*Hm0/Tp^2,101);
% [S,H] = meshgrid(s,h); 
% f = thsnlpdf(H,S,Hm0,Tp);
% contourf(s,h,f)  
%
% See also  thspdf2, thsspdf

% Reference  
% P. A. Brodtkorb (2004),  
% The Probability of Occurrence of Dangerous Wave Situations at Sea.
% Dr.Ing thesis, Norwegian University of Science and Technolgy, NTNU,
% Trondheim, Norway.   
  
% Adapted to  cssmooth  by GL Feb 2011  
% History
% revised pab 09.08.2003
%  changed input and help header  
% validated 20.11.2002
% By pab 20.12.2000

%error(nargchk(3,6,nargin))
narginchk(3,6)

if (nargin < 6||isempty(condon)),  condon  = 0; end
if (nargin < 5||isempty(normalizedInput)),  normalizedInput  = 0;end
if (nargin < 4||isempty(Tp)),  Tp  = 8;end
if (nargin < 3||isempty(Hm0)), Hm0 = 6;end

multipleSeaStates = any(numel(Hm0)>1|numel(Tp)>1);
if multipleSeaStates
  [icode, Scf,Hd,Hm0,Tp] = iscomnsize(Scf,Hd,Hm0,Tp);
else
  [icode, Scf,Hd] = iscomnsize(Scf,Hd);
end
if ~icode 
    error('Requires non-scalar arguments to match in size.');
end
displayWarning = 0;
if displayWarning,
  if any(Hm0>11| Hm0>(Tp-2)*12/11) 
    disp('Warning: Hm0 is outside the valid range')
    disp('The validity of the Joint (Hd,Scf) distribution is questionable')
  end
  if any(Tp>20|Tp<3)
    disp('Warning: Tp is outside the valid range')
    disp('The validity of the Joint (Hd,Scf) distribution is questionable')
  end
end
  
global THSNLPAR
if isempty(THSNLPAR)
  %THSNLPAR = load('thsnlpar.mat');
  THSNLPAR = load('thsnlpar06-Aug-2004.mat') ; 
end

Tpp  = THSNLPAR.Tp;
Hm00 = THSNLPAR.Hm0;
Tm020 = THSNLPAR.Tm02;
h2   = THSNLPAR.h2;
% Interpolation method
method = '*cubic';% Faster interpolation

if normalizedInput,
  Hrms = 1;
  Vrms = 1;
else
  [Tp1,Hs1] = meshgrid(Tpp,Hm00);
 
  if 1,
    Tm02 = interp2(Tp1,Hs1,Tm020,Tp,Hm0,method);
    %Tp/Tm02
  else
    
    Tm02 = Tp;
    for ix= 1:100
      dTp = (Tm02-interp2(Tp1,Hs1,Tm020,Tp,Hm0,method));
      Tp = Tp+dTp;
      if all(abs(dTp)<0.01)
	%dTp
	%ix
	break
      end
    end
    %Tp,    Tp/Tm02
  end
  
%  w    = linspace(0,100,16*1024+1).'; % torsethaugen original spacing
  %w    = linspace(0,10,2*1024+1).'; 
%  St = torsethaugen(w,[Hm0,Tp]);
%  ch   = spec2char(St,{'Tm02','eps2'});
%  Tm02 = ch(1);
%  eps2 = ch(2);
  Hrms = Hm0/sqrt(2);
  Vrms = 1.25*Hm0./(Tm02.^2); % Erms
end

h = Hd./Hrms;
v = Scf./Vrms;
cSize = size(h); % common size

% Gamma distribution parameters as a function of Tp Hm0 and h2
A11 = THSNLPAR.A11s;
B11 = THSNLPAR.B11s;

[E1, H1, H2] = meshgrid(Tpp,Hm00,h2);
Nh2 = length(h2);
if multipleSeaStates
  h   = h(:);
  v   = v(:);
  Tp  = Tp(:);
  Hm0 = Hm0(:);
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
% Waveheight distribution in time
% Truncated Weibull  distribution parameters as a function of Tp, Hm0 
[A0, B0, C0] = thwparfun(Hm0,Tp,'time');
switch condon,
 case 0, % regular pdf is returned 
  f = pdfweibmod(h,A0,B0,C0).*pdfgam(v,A1,B1);
 case 1, %pdf conditioned on x1 ie. p(x2|x1) 
  f = pdfgam(v,A1,B1);
 case 3, % secret option  used by XXstat: returns x2*p(x2|x1) 
  f = v.*pdfgam(v,A1,B1);
 case 4, % secret option  used by XXstat: returns x2.^2*p(x2|x1) 
  f = v.^2.*pdfgam(v,A1,B1);
 case 5, % p(h)*P(V|h) is returned special case used by thscdf
  f = pdfweibmod(h,A0,B0,C0).*cdfgam(v,A1,B1);
 case 6, % P(V|h) is returned special case used by thscdf
  f = cdfgam(v,A1,B1);
 case 7,% p(h)*(1-P(V|h)) is returned special case used by thscdf
  f = pdfweibmod(h,A0,B0,C0).*(1-cdfgam(v,A1,B1));
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




