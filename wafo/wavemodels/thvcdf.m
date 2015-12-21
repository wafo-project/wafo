function f = thvcdf(Hd,Vcf,Hm0,Tp,tail)
%THVCDF Joint (Vcf,Hd) CDF for linear waves with Torsethaugen spectra. 
%
%  CALL: f = thvcdf(Hd,Vcf,Hm0,Tp,tail)
% 
%   f   = CDF evaluated at (Vcf,Hd)
%   Hd  = zero down crossing wave height 
%   Vcf = crest front velocity
%   Hm0 = significant wave height [m]
%   Tp  = Spectral peak period    [s]
%  tail = 1 if upper tail is calculated   
%         0 if lower tail is calulated (default)
%  
% THVCDF approximates the joint CDF of (Vcf, Hd), i.e., crest front
% velocity (Ac/Tcf) and wave height, for a Gaussian process with a
% Torsethaugen spectral density. The empirical parameters of the model is
% fitted by least squares to simulated (Vcf,Hd) data for 600 classes of
% Hm0 and Tp. Between 50000 and 150000 zero-downcrossing waves were
% simulated for each class of Hm0 and Tp.
% THVCDF is restricted to the following range for Hm0 and Tp: 
%  0.5 < Hm0 [m] < 12,  3.5 < Tp [s] < 20,  and  Hm0 < (Tp-2)*12/11.
%
% Example:
% Hm0 = 6;Tp = 8;
% vc = 3;
% hc = 3;
% lowerTail = 0;
% upperTail = ~lowerTail  
% thvcdf(hc,vc,Hm0,Tp)           % Prob(Hd<Hc,Vcf<Vc)
% thvcdf(hc,vc,Hm0,Tp,upperTail) % Prob(Hd>Hc,Vcf>Vc)  
%  
%  % Conditional probability of steep and high waves given seastates
%  % i.e., Prob(Hd>hc,Vcf>vc|Hs,Tp)  
%  upperTail = 1;
%  Hs = linspace(2.5,11.5,10);
%  Tp = linspace(4.5,19.5,16);
%  [T,H] = meshgrid(Tp,Hs); 
%  p = thvcdf(hc,vc,H,T,upperTail);
%  v = 10.^(-6:-1);  
%  contourf(Tp,Hs,log10(p),log10(v))
%  xlabel('Tp')
%  ylabel('Hs')  
%  fcolorbar(log10(v))  
%  
% See also  thvpdf

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
  [icode, Vcf,Hd,Hm0,Tp] = iscomnsize(Vcf,Hd,Hm0,Tp);
else
  [icode, Vcf,Hd] = iscomnsize(Vcf,Hd);
end
if ~icode 
  error('Requires non-scalar arguments to match in size.');
end

global THVPAR
if isempty(THVPAR)
  THVPAR = load('thvpar.mat');
end

Tpp  = THVPAR.Tp;
Hm00 = THVPAR.Hm0;
Tm020 = THVPAR.Tm02;
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
Vrms = 2*Hm0./Tm02; % Erms

v = Vcf./Vrms;
f = zeros(size(Hd));

% Only compute within valid range
k0 = find((2<=Tp) & (Tp<=21) & (Hm0<=(Tp-2)*12/11) & (Hm0<=12));
if any(k0)
  hMax = 5;
  eps2 = 1e-6;
  h = min(Hd./Hrms,hMax);
  
  if multipleSeaStates
    h = h(k0);
    v = v(k0);
    Hm0 = Hm0(k0);
    Tp = Tp(k0);
  else
    k0 = 1:numel(Hd);
  end
  if 0
    % This is a trick to get the html documentation correct.
    k = thvpdf(1,1,2,3);
  end
  normalizedInput = 1;
  utprb = gaussq(@thvpdf,hMax,2*hMax,eps2/2,[],mean(v(:)),mean(Hm0(:)),mean(Tp(:)),normalizedInput,7);
  if eps2<utprb
    warning('Check the accuracy of integration!')
  end

 

  hlim    = h;


  lowerTail = 0;
  if tail==lowerTail,
    k       = (h>2.5);
    hlim(k) = 2.5;
    k       = (h>1.3*v);
    hlim(k) = 1.3*v(k);
    f(k0) = gaussq(@thvpdf,0,hlim,eps2/2,[],v,Hm0,Tp,normalizedInput,5)...
	+ gaussq(@thvpdf,hlim,h,eps2/2,[],v,Hm0,Tp,normalizedInput,5); 
  else % upper tail
    k       = find(h<1.3*v);
    hlim(k) = 1.3*v(k);
    f(k0) = gaussq(@thvpdf,h,hlim,eps2/2,[],v,Hm0,Tp,normalizedInput,7)...
	+ gaussq(@thvpdf,hlim,hMax,eps2/2,[],v,Hm0,Tp,normalizedInput,7); 
  end
end
return

