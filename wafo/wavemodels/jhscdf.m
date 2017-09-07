function f = jhscdf(Hd,Scf,Hm0,Tp,gam,tail)
%JHSCDF Joint (Scf,Hd) CDF for linear waves with a JONSWAP spectrum.
%
%  CALL: f = jhscdf(Hd,Scf,Hm0,Tp,Gamma,tail)
% 
%   f   = CDF evaluated at (Scf,Hd)
%   Hd  = zero down crossing wave height [m] 
%   Scf = crest front steepness
%   Hm0 = significant wave height [m]
%   Tp  = Spectral peak period    [s]
% Gamma = Peakedness parameter of the JONSWAP spectrum  
%  tail = 1 if upper tail is calculated   
%         0 if lower tail is calulated (default)
%  
% JHSCDF approximates the joint CDF of (Scf, Hd), i.e., crest front
% steepness (2*pi*Ac/(g*Td*Tcf)) and wave height, for a Gaussian process with a
% JONSWAP spectral density. The empirical parameters of the model is
% fitted by least squares to simulated (Scf,Hd) data for 13 classes of
% GAMMA between 1 and 7. Between 47000 and 55000 zero-downcrossing waves were
% simulated for each class of GAMMA.
% JHSCDF is restricted to the following range for GAMMA: 
%  1 <= GAMMA <= 7 
%
% Example:
% Hm0 = 6;Tp = 9; gam = 3.5;
% ec = 0.25;
% hc = 3;
% lowerTail = 0;
% upperTail = ~lowerTail  
% jhscdf(hc,ec,Hm0,Tp,gam)           % Prob(Hd<Hc,Scf<ec)
% jhscdf(hc,ec,Hm0,Tp,gam,upperTail) % Prob(Hd>Hc,Scf>ec)  
%  
%  % Conditional probability of steep and high waves given seastates
%  % i.e., Prob(Hd>hc,Scf>ec|Hs,Tp)  
%  upperTail = 1;
%  Hs = linspace(2.5,11.5,10);
%  Tp = linspace(4.5,19.5,16);
%  [T,H] = meshgrid(Tp,Hs); 
%  p = jhscdf(hc,ec,H,T,gam,upperTail);
%  v = 10.^(-6:-1);  
%  contourf(Tp,Hs,log10(p),log10(v))
%  xlabel('Tp')
%  ylabel('Hs')  
%  fcolorbar(log10(v))  
%  
% See also  jhsnlcdf

% Reference 
% P. A. Brodtkorb (2004),  
% The Probability of Occurrence of Dangerous Wave Situations at Sea.
% Dr.Ing thesis, Norwegian University of Science and Technolgy, NTNU,
% Trondheim, Norway.     
  
% History
% By pab 20.01.2003


%error(nargchk(2,6,nargin))  
narginchk(2,6)
if (nargin < 6||isempty(tail)),tail = 0; end
if (nargin < 4||isempty(Tp)),Tp = 8; end
if (nargin < 3||isempty(Hm0)), Hm0 = 6; end
if (nargin < 5||isempty(gam)),
   gam = getjonswappeakedness(Hm0,Tp);
end

multipleSeaStates = any(numel(Hm0)>1|numel(Tp)>1);
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
%dev = 2e-5;
c1 =[ 0.16183666835624   1.53691936441548   1.55852759524555];
c2 =[ 0.15659478203944   1.15736959972513   1];
Tm02 = Tp.*(polyval(c2,gam)./polyval(c1,gam));

Hrms = Hm0/sqrt(2);
Vrms = 1.25*Hm0./(Tm02.^2); % Erms

v = Scf./Vrms;

hMax = 6;
eps2 = 1e-6;
normalizedInput = 1;
utprb = gaussq('jhspdf',hMax,2*hMax,eps2/2,[],mean(v(:)),mean(Hm0(:)),mean(Tp(:)),mean(gam(:)),normalizedInput,7);
if eps2<utprb
  warning('WAFO:JHSCDF','Check the accuracy of integration!')
end

h = min(Hd./Hrms,hMax);

f = zeros(size(Hd));
% Only compute
if 0, % haver parametrization
  loLimit = 3.6;
  upLimit = 5;
else
  loLimit = 2.5;
  upLimit = 6.5;
end

k0 = find( (loLimit*sqrt(Hm0)<Tp) & (Tp<upLimit*sqrt(Hm0)) );
if any(k0)
  if multipleSeaStates
    h = h(k0);
    v = v(k0);
    Hm0 = Hm0(k0);
    Tp = Tp(k0);
    gam = gam(k0);    
  else
    k0 = 1:numel(Hd);
  end
else
  return
end

if 0
  % This is a trick to get the html documentation correct.
  k = jhspdf(1,1,2,3);
end
hlim  = h;
lowerTail = 0;
if tail==lowerTail,
  %k       = find(h>2.5);%*v);
  hlim(h>2.5) = 2.5;%*v(k);
  f(k0) = gaussq(@jhspdf,0,hlim,eps2/2,[],v,Hm0,Tp,gam,normalizedInput,5)...
      + gaussq(@jhspdf,hlim,h,eps2/2,[],v,Hm0,Tp,gam,normalizedInput,5); 
else % upper tail
  %k       = find(h<2.5);%*v);
  hlim(h<2.5) = 2.5;%*v(k);
  
  f(k0) = gaussq(@jhspdf,h,hlim,eps2/2,[],v,Hm0,Tp,gam,normalizedInput,7)...
      + gaussq(@jhspdf,hlim,hMax,eps2/2,[],v,Hm0,Tp,gam,normalizedInput,7); 
end
return

