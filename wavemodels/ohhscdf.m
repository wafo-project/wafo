function f = ohhscdf(Hd,Scf,Hm0,def,tail)
%OHHSCDF Joint (Scf,Hd) CDF for linear waves with Ochi-Hubble spectra.
%
%  CALL: f = ohhscdf(Hd,Scf,Hm0,def)
% 
%  f   = cdf
%  Hd  = zero down crossing wave height
%  Scf = crest front steepness
%  Hm0 = significant wave height [m].
%  def = defines the parametrization of the spectral density (default 1)
%        1 : The most probable spectrum  (default)
%        2,3,...11 : gives 95% Confidence spectra
% tail = 1 if upper tail is calculated 
%        0 if lower tail is calulated  (default)
%
% OHHSCDF approximates the joint CDF of (Scf, Hd) in time, 
% i.e., crest front steepness (2*pi*Ac/(g*Td*Tcf)) and wave height,
%  for a Gaussian process with a bimodal Ochi-Hubble spectral density
% (ochihubble). The empirical
% parameters of the model is fitted by least squares to simulated
% (Scf,Hd) data for 24 classes of Hm0. Between 50000 and 150000
% zero-downcrossing waves were simulated for each class of Hm0.
% OHHSCDF is restricted to the following range for Hm0: 
% 0.5 < Hm0 [m] < 12
% The size of f is the common size of the input arguments, Hd, Scf and
% Hm0.  
%
% Example:
% Hm0 = 6; def= 2;
% Ec = 0.25;
% Hc = 3;
% lowerTail = 0;
% upperTail = ~lowerTail  
% ohhscdf(Hc,Ec,Hm0,def)           % Prob(Hd<Hc,Scf<Ec)
% ohhscdf(Hc,Ec,Hm0,def,upperTail) % Prob(Hd>Hc,Scf>Ec)  
%
% See also  ohhspdf, thspdf

% Reference 
% P. A. Brodtkorb (2004),  
% The Probability of Occurrence of Dangerous Wave Situations at Sea.
% Dr.Ing thesis, Norwegian University of Science and Technolgy, NTNU,
% Trondheim, Norway.    
  
  
% History
% revised pab 09.09.2003
% By pab 06.02.2001

error(nargchk(3,5,nargin))
if (nargin < 5||isempty(tail)),  tail  = 0;end
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
  Warning('DEF is outside the valid range')
  def = mod(def-1,11)+1;
end

persistent OHHSPAR
if isempty(OHHSPAR)
  OHHSPAR = load('ohhspar.mat');
end
method = 'cubic';
Tm020 = OHHSPAR.Tm02;
Hm00  = OHHSPAR.Hm0;
Hrms = Hm0/sqrt(2);
Tm02 = interp1(Hm00,Tm020(:,def),Hm0,method);
Erms = 1.25*Hm0./(Tm02.^2); % Erms



s = Scf./Erms;
hMax = 10;
h = min(Hd./Hrms,hMax);

eps2 = 1e-6;

hlim    = h;

normalizedInput = 1;
lowerTail = 0;

if 0
  % This is a trick to get the html documentation correct.
  k = ohhspdf(1,1,2,3);
end

if (tail == lowerTail)
  %k       = find(h>2.5);
  hlim(h>2.5) = 2.5;
  f = gaussq(@ohhspdf,0,hlim,eps2/2,[],s,Hm0,def,normalizedInput,5)...
      + gaussq(@ohhspdf,hlim,h,eps2/2,[],s,Hm0,def,normalizedInput,5); 
else % upper tail
  %k       = find(h<2.5);
  hlim(h<2.5) = 2.5;
  f = gaussq(@ohhspdf,h,hlim,eps2/2,[],s,Hm0,def,normalizedInput,7)...
      + gaussq(@ohhspdf,hlim,hMax,eps2/2,[],s,Hm0,def,normalizedInput,7); 
end
return

