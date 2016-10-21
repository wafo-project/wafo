function f = thsscdf(Hd,Scf,Hm0,Tp,tail)
%THSSCDF Joint (Scf,Hd) CDF for linear waves in space with Torsethaugen spectra.
%
%  CALL: f = thsscdf(Hd,Scf,Hm0,Tp,tail)
% 
%  f   = CDF evaluated at (Scf,Hd)
%  Hd  = zero down crossing wave height
%  Scf = crest front steepness
%  Hm0 = significant wave height [m]
%  Tp  = Spectral peak period    [s]
% tail = 1 if upper tail is calculated   
%        0 if lower tail is calulated (default)  
%
% THSSCDF approximates the joint CDF of (Scf, Hd), i.e., crest
% front steepness (Ac/Lcf) and wave height in space, for a Gaussian
% process with a Torsethaugen spectral density. The empirical parameters
% of the model is fitted by least squares to simulated (Scf,Hd) data for
% 600 classes of Hm0 and Tp. Between 100000 and 1000000 zero-downcrossing
% waves were simulated for each class of Hm0 and Tp.
% THSSCDF is restricted to the following range for Hm0 and Tp: 
%  0.5 < Hm0 [m] < 12,  3.5 < Tp [s] < 20,  and  Hm0 < (Tp-2)*12/11.
%
% Example:
% Hm0 = 6;Tp = 8;
% Ec = 0.25;
% Hc = 3;
% lowerTail = 0;
% upperTail = ~lowerTail  
% thsscdf(Hc,Ec,Hm0,Tp)           % Prob(Hd<Hc,Scf<Ec)
% thsscdf(Hc,Ec,Hm0,Tp,upperTail) % Prob(Hd>Hc,Scf>Ec)  
%  
% See also  thscdf

% Reference 
% P. A. Brodtkorb (2004),  
% The Probability of Occurrence of Dangerous Wave Situations at Sea.
% Dr.Ing thesis, Norwegian University of Science and Technolgy, NTNU,
% Trondheim, Norway.  
  
% History
% revised pab 14.03.2004  
% revised pab 09.08.2003
% changed input  
% validated 20.11.2002
% By pab 20.12.2000


error(nargchk(3,5,nargin))  
global THSSPARW
global THSSPARG
if (nargin < 5||isempty(tail)),tail = 0; end
if (nargin < 4||isempty(Tp)),Tp = 8; end
if (nargin < 3||isempty(Hm0)), Hm0 = 6; end

multipleSeaStates = any(numel(Hm0)>1|numel(Tp)>1);
if multipleSeaStates
  [icode, Scf,Hd,Hm0,Tp] = iscomnsize(Scf,Hd,Hm0,Tp);
else
  [icode, Scf,Hd] = iscomnsize(Scf,Hd);
end
if ~icode
  error('Requires non-scalar arguments to match in size.');
end
useWeibull = 1;
if useWeibull
 
  if isempty(THSSPARW)
    THSSPARW= load('thsspar27-Jul-2004.mat');
  end
  Tpp  = THSSPARW.Tp;
  Hm00 = THSSPARW.Hm0;
  Tm020 = THSSPARW.Tm02;
else

  if isempty(THSSPARG)
    THSSPARG = load('thsspar.mat');
  end

  Tpp  = THSSPARG.Tp;
  Hm00 = THSSPARG.Hm0;
  Tm020 = THSSPARG.Tm02;
end
  % Interpolation method
method = '*cubic';% Faster interpolation

[Tp1,Hs1] = meshgrid(Tpp,Hm00);
Tm02 = interp2(Tp1,Hs1,Tm020,Tp,Hm0,method);
%  w    = linspace(0,100,16*1024+1).'; % torsethaugen original spacing
  %w    = linspace(0,10,2*1024+1).'; 
%  St = torsethaugen(w,[Hm0,Tp]);
%  ch   = spec2char(St,{'Tm02','eps2'});
%  Tm02 = ch(1);
%  eps2 = ch(2);
Hrms = Hm0/sqrt(2);
Erms = 2*Hm0./(Tm02); % Srms

s = Scf./Erms;
hMax = 10;
h = min(Hd./Hrms,hMax);

eps2 = 1e-6;

f = zeros(size(Hd));

% Only compute within valid range
%k0 = find((2<=Tp) & (Tp<=21) & (Hm0<=(Tp-2)*12/11) & (Hm0<=12));
%upLimit = 6.5;
loLimit = 2.5;
k0 = find((2<=Tp) & (Tp<=21) & (loLimit*sqrt(Hm0)<Tp) & (Hm0<=12));
if any(k0)
  if multipleSeaStates
    h = h(k0);
    s = s(k0);
    Hm0 = Hm0(k0);
    Tp = Tp(k0);
  else
    k0 = 1:numel(Hd);
  end
  hlim    = h;

  normalizedInput = 1;
  lowerTail = 0;
  if tail==lowerTail
    k       = (h>2.5);
    hlim(k) = 2.5;
    f(k0) = gaussq(@thsspdf,0,hlim,eps2/2,[],s,Hm0,Tp,normalizedInput,5)...
	    + gaussq(@thsspdf,hlim,h,eps2/2,[],s,Hm0,Tp,normalizedInput,5); 
  else % upper tail
    k       = (h<2.5);
    hlim(k) = 2.5;
    f(k0) = gaussq(@thsspdf,h,hlim,eps2/2,[],s,Hm0,Tp,normalizedInput,7)...
	    + gaussq(@thsspdf,hlim,hMax,eps2/2,[],s,Hm0,Tp,normalizedInput,7); 
  end
end
return

