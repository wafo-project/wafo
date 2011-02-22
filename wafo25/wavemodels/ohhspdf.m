function [f,Hrms,Vrms,fA,fB] = ohhspdf(Hd,Scf,Hm0,def,normalizedInput, ...
				       condon)
%OHHSPDF Joint (Scf,Hd) PDF for linear waves with Ochi-Hubble spectra. 
%
%  CALL: f = ohhspdf(Hd,Scf,Hm0,def)
% 
%  f   = pdf
%  Hd  = zero down crossing wave height
%  Scf = crest front steepness
%  Hm0 = significant wave height [m].
%  def = defines the parametrization of the spectral density (default 1)
%        1 : The most probable spectrum  (default)
%        2,3,...11 : gives 95% Confidence spectra
%
% OHHSPDF approximates the joint distribution of (Scf, Hd) in time, 
% i.e., crest front steepness (2*pi*Ac/(g*Td*Tcf)) and wave height,
%  for a Gaussian process with a bimodal Ochi-Hubble spectral density
% (ochihubble). The empirical
% parameters of the model is fitted by least squares to simulated
% (Scf,Hd) data for 24 classes of Hm0. Between 50000 and 150000
% zero-downcrossing waves were simulated for each class of Hm0.
% OHHSPDF is restricted to the following range for Hm0: 
% 0.5 < Hm0 [m] < 12
% The size of f is the common size of the input arguments, Hd, Scf and
% Hm0.  
%
% Example:
% Hm0 = 6;Tp = 8;def= 2;
% h = linspace(0,4*Hm0/sqrt(2))'; 
% v = linspace(0,6*1.25*Hm0/Tp^2)';
% [V,H]  = meshgrid(v,h);
% f = ohhspdf(H,V,Hm0,def);
% contourf(v,h,f)
%
% See also  ohhpdf, thspdf

% Reference  
% P. A. Brodtkorb (2004),  
% The Probability of Occurrence of Dangerous Wave Situations at Sea.
% Dr.Ing thesis, Norwegian University of Science and Technolgy, NTNU,
% Trondheim, Norway. 

% Adapted to  cssmooth  by GL Feb 2011
% History
% revised pab 09.09.2003
% By pab 06.02.2001

error(nargchk(3,6,nargin))
if (nargin < 6||isempty(condon)),  condon  = 0;end
if (nargin < 5||isempty(normalizedInput)),  normalizedInput  = 0;end
if (nargin < 4||isempty(def)), def=1;end 

multipleSeaStates = any(numel(Hm0)>1);
if multipleSeaStates
  [icode, Scf,Hd,Hm0] = iscomnsize(Scf,Hd,Hm0);
else
  [icode, Scf,Hd] = iscomnsize(Scf,Hd);
end
if ~icode 
    error('Requires non-scalar arguments to match in size.');
end

if any(Hm0>12| Hm0<=0.5) 
  disp('Warning: Hm0 is outside the valid range')
  disp('The validity of the Hd distribution is questionable')
end

if def>11||def<1 
  warning('DEF is outside the valid range')
  def = mod(def-1,11)+1;
end

global OHHSPAR
if isempty(OHHSPAR)
  OHHSPAR = load('ohhspar.mat');
end


%w    = linspace(0,100,16*1024+1).'; %original spacing
%ch   = spec2char(ochihubble(w,[Hm0,def]),{'Tm02','eps2'});
%Tm02 = ch(1);
%Vrms = 2*Hm0/Tm02;
Hm00  = OHHSPAR.Hm0;

if normalizedInput,
  Hrms = 1;
  Vrms = 1;
else
  method = 'cubic'; %'spline'
  Tm020 = OHHSPAR.Tm02;  
  Hrms = Hm0/sqrt(2);
  Tm02 = interp1(Hm00,Tm020(:,def),Hm0,method);
  Vrms = 1.25*Hm0./(Tm02.^2); % Erms
end


h = Hd./Hrms;
v = Scf./Vrms;
cSize = size(h); % common size

% Logarithm of Weibull distribution parameters as a function of Hm0 and h2
A11  = squeeze(OHHSPAR.A11s(:,def,:));
B11  = squeeze(OHHSPAR.B11s(:,def,:));

h2    = OHHSPAR.h2; %Hd/Hrms
[E2 H2] = meshgrid(Hm00,h2);
method = '*cubic'; %'spline'
if multipleSeaStates
  h   = h(:);
  v   = v(:);
  Hm0 = Hm0(:);
  A1 = zeros(length(h),1);
  B1 = A1;
  [Hm0u,ix,jx] = unique(Hm0);
  numSeaStates = length(ix);
  Nh2 = length(h2);
  Hm0i = zeros(Nh2,1);
  for iz=1:numSeaStates
    k = find(jx==iz);
    Hm0i(:) = Hm0u(iz);
    A1(k) = exp(cssmooth(h2,interp2(E2,H2,A11.',Hm0i,h2,method),...
		       1,h(k),1));
    B1(k) = (cssmooth(h2,interp2(E2,H2,exp(B11).',Hm0i,h2,method),...
		       1,h(k),1));
  end
else
  Hm0i = Hm0(ones(size(h2)));
  A1 = exp(cssmooth(h2,interp2(E2,H2,A11.', ...
			     Hm0i,h2,method),1,h,1));
  B1 = (cssmooth(h2,interp2(E2,H2,exp(B11).', ...
			  Hm0i,h2,method),1,h,1));
end
%fh = ohhpdf(h(:)/Hrms,Hm0,def,'time',1);
% Fh = fh.f;
% Waveheight distribution in time
% Generalized gamma distribution parameters as a function of Hm0 
[A0 B0 C0] = ohhgparfun(Hm0,def,'time');


switch condon,
 case 0, % regular pdf is returned 
  f = pdfgengam(h,A0,B0,C0).*pdfweib(v,A1,B1);
 case 1, %pdf conditioned on x1 ie. p(x2|x1) 
  f = pdfweib(v,A1,B1);
 case 3, % secret option  used by XXstat: returns x2*p(x2|x1) 
  f = v.*pdfweib(v,A1,B1);
 case 4, % secret option  used by XXstat: returns x2.^2*p(x2|x1) 
  f = v.^2.*pdfweib(v,A1,B1);
 case 5, % p(h)*P(V|h) is returned special case used by thscdf
  f = pdfgengam(h,A0,B0,C0).*cdfweib(v,A1,B1);
 case 6, % P(V|h) is returned special case used by thscdf
  f = cdfweib(v,A1,B1);
 case 7,% p(h)*(1-P(V|h)) is returned special case used by thscdf
  f = pdfgengam(h,A0,B0,C0).*(1-cdfweib(v,A1,B1));
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
  fA.x    = {Hm00, h2};
  fA.labx = {'Hm0', 'h'};
  fA.note = sprintf('%s \n',...
		    'The conditonal Weibull distribution Parameter',...
		    'A(h)/Hrms for wave heigth as a function of',...
		    'h=Hd/Hrms and Hm0 for the ochihubble spectrum, ',...
		    sprintf('given def = %d.',def));
    
  ra = range(A11(:));
  st = round(min(A11(:))*100)/100;
  en = max(A11(:));
  fA.cl   = st:ra/20:en;
end
if nargout>4,
  fB      = createpdf(2);
  fB.f    = B11;
  fB.x    = {Hm00,h2};
  fB.labx = {'Hm0','h'};
  fB.note =  sprintf('%s \n',...
		    'The conditonal Weibull distribution Parameter',...
		    'B(h)/Hrms for wave heigth as a function of',...
		    'h=Hd/Hrms and Hm0 for the ochihubble spectrum, ',...
		    sprintf('given def = %d.',def));
  ra = range(B11(:));
  st = round(min(B11(:))*100)/100;
  en = max(B11(:));
  fB.cl   = st:ra/20:en;
end
return







% old call
% Fh = pdfgengam(h(:)/Hrms,A0,B0,C0);
% 
% 
% 
% %C1 = exp(smooth(h2,interp2(E2,H2,log(squeeze(C11(:,def,:))).',Hm0(ones(size(h2))),h2,method),1,h/Hrms,1));
% %figure(1), plot(h/Hrms, A1,h2,squeeze(exp(A11(10,def,:))),'r')
% %figure(2), plot(h/Hrms,B1,h2,squeeze(exp(B11(10,def,:))),'r')
% [V1 A1] = meshgrid(v/Vrms,A1);
% [V1 B1] = meshgrid(v/Vrms,B1);
% %[V1 C1] = meshgrid(v/Vrms,C1);
% 
% f = createpdf(2);
% f.title = 'Joint distribution of (Hd,Scf) in time';
% if norm==0
%   f.f = repmat(Fh/Hrms,[1 length(v)]).*pdfweib(V1,A1,B1)/Vrms;
%   f.x = {v,h};
%   f.norm=0;
%   f.labx={'Scf', 'Hd [m]'};
% else
%   f.f = repmat(Fh,[1 length(v)]).*pdfweib(V1,A1,B1);
%   f.x = {v/Vrms,h/Hrms};
%   f.labx={'Scf', 'Hd'};
%   f.norm = 1;
% end
% f.note = ['Ochi-Hubble Hm0=' num2str(Hm0) ' def = ' num2str(def)];
% if 1,
%   f.f(isinf(f.f)|isnan(f.f))=0;
%   % Some places numerical problems occur
%   %f.f(f.f>30)=0; % Needed to calculate the contour levels
%   [f.cl,f.pl] = qlevels(f.f);%,[10:20:90 95 99 99.9 99.99 99.999 99.9999]);
% end
