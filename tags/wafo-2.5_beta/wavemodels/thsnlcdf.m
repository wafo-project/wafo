function f = thsnlcdf(Hd,Scf,Hm0,Tp,tail)
%THSNLCDF Joint (Scf,Hd) CDF for non-linear waves with Torsethaugen spectra.
%
%  CALL: f = thsnlcdf(Hd,Scf,Hm0,Tp,tail)
% 
%   f   = CDF evaluated at (Scf,Hd)
%   Hd  = zero down crossing wave height (vector)
%   Scf = crest front steepness (vector) 
%   Hm0 = significant wave height [m]
%   Tp  = Spectral peak period    [s]
%  tail = 1 if upper tail is calculated   
%         0 if lower tail is calulated (default)
%  
% THSNLCDF approximates the joint CDF of (Scf, Hd), i.e., crest front
% steepness (2*pi*Ac/(g*Td*Tcf)) and wave height, for 2nd order
% nonlinear waves with a Torsethaugen spectral density.
%  The empirical parameters of the
% model is fitted by least squares to simulated (Scf,Hd) data for 600
% classes of Hm0 and Tp. Between 40000 and 200000 zero-downcrossing waves
% were simulated for each class of Hm0 and Tp.
% THSNLCDF is restricted to the following range for Hm0 and Tp: 
%  0.5 < Hm0 [m] < 12,  3.5 < Tp [s] < 20,  and  Hm0 < (Tp-2)*12/11.
%
% Example:
% Hm0 = 6;Tp = 8;
% sc = 0.25;
% hc = 3;
% lowerTail = 0;
% upperTail = ~lowerTail  
% thsnlcdf(hc,sc,Hm0,Tp)           % Prob(Hd<hc,Scf<sc)
% thsnlcdf(hc,sc,Hm0,Tp,upperTail) % Prob(Hd>hc,Scf>sc)  
%
%  % Conditional probability of steep and high waves given seastates
%  % i.e., Prob(Hd>hc,Scf>sc|Hs,Tz)  
%  upperTail = 1;
%  Hs = linspace(2.5,11.5,10);
%  Tp = linspace(4.5,19.5,16);
%  [T,H] = meshgrid(Tp,Hs); 
%  p = thsnlcdf(hc,sc,H,T,upperTail);
%  v = 10.^(-6:-1);  
%  contourf(Tp,Hs,log10(p),log10(v))
%  xlabel('Tp')
%  ylabel('Hs')  
%  fcolorbar(log10(v))  
%  
% See also  thsnlpdf

% Reference 
% P. A. Brodtkorb (2004),  
% The Probability of Occurrence of Dangerous Wave Situations at Sea.
% Dr.Ing thesis, Norwegian University of Science and Technolgy, NTNU,
% Trondheim, Norway.  
  
% History
% revised pab 09.08.2003
% changed input  
% validated 20.11.2002
% By pab 20.12.2000


error(nargchk(3,5,nargin))  

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

global THSNLPAR
if isempty(THSNLPAR)
  %THSNLPAR = load('thsnlpar.mat');
  THSNLPAR = load('thsnlpar06-Aug-2004.mat');  
end

Tpp  = THSNLPAR.Tp;
Hm00 = THSNLPAR.Hm0;
Tm020 = THSNLPAR.Tm02;
% Interpolation method
method = '*cubic';% Faster interpolation

[Tp1,Hs1] = meshgrid(Tpp,Hm00);
if 1,
  Tm02 = interp2(Tp1,Hs1,Tm020,Tp,Hm0,method);
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
end


%  w    = linspace(0,100,16*1024+1).'; % torsethaugen original spacing
  %w    = linspace(0,10,2*1024+1).'; 
%  St = torsethaugen(w,[Hm0,Tp]);
%  ch   = spec2char(St,{'Tm02','eps2'});
%  Tm02 = ch(1);
%  eps2 = ch(2);
Hrms = Hm0/sqrt(2);
Erms = 1.25*Hm0./(Tm02.^2); % Erms

s = Scf./Erms;
hMax = 5;
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
  
  if 0
    % This is a trick to get the html documentation correct.
    k = thsnlpdf(1,1,2,3);
  end
  
  if tail==lowerTail
   % k       = find(h>2.5);
    hlim(h>2.5) = 2.5;
    f(k0) = gaussq(@thsnlpdf,0,hlim,eps2/2,[],s,Hm0,Tp,normalizedInput,5)...
	+ gaussq('thsnlpdf',hlim,h,eps2/2,[],s,Hm0,Tp,normalizedInput,5); 
  else % upper tail
    %k       = find(h<2.5);
    hlim(h<2.5) = 2.5;
    f(k0) = gaussq(@thsnlpdf,h,hlim,eps2/2,[],s,Hm0,Tp,normalizedInput,7)...
	+ gaussq('thsnlpdf',hlim,hMax,eps2/2,[],s,Hm0,Tp,normalizedInput,7); 
  end
end
return

