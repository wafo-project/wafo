function [p, eps2] = mk87cdf(Hd,Scf,Hs,Tz,tail)
%MK87CDF Myrhaug and Kjeldsen (1987) joint (Scf,Hd) CDF. 
%
% CALL: F = mk87cdf(Hd,Scf,Hs,Tz,tail)
%
%    F  = CDF
%   Hd  = zero down crossing wave height.
%   Scf = crest front steepness.
%   Hs  = significant wave height.
%   Tz  = average zero down crossing period.
%  tail = 1 if upper tail is calculated   
%         0 if lower tail is calulated (default)
%
% MK87CDF returns the joint CDF of (Scf, Hd) given Hs and Tz,
% i.e., crest front steepness (2*pi*Ac/(g*Td*Tcf)) and wave height, given
% the seastate. The root mean square values of Hd and Scf (Hrms,Erms) are
% related to the significant waveheight and the average zero down
% crossing period by:
%             Hrms = 0.715*Hs;
%             Erms = 0.0202+0.826*Hs/(Tz^2);
%  
%   The size of F is the common size of  the input arguments
%
% Examples: 
%  Hs = 5.5;
%  Tz = 8.5;
%  sc = 0.25;
%  hc = 4;
%  p = mk87cdf(hc,sc,Hs,Tz) % Prob(Hd<=hc,Scf<=sc|Hs,Tz) = 0.59
%  
%  % Conditional probability of steep and high waves given seastates
%  % i.e., Prob(Hd>hc,Scf>sc|Hs,Tz)  
%  upperTail = 1;
%  Hs = linspace(2.5,18.5,17);
%  Tz = linspace(4.5,19.5,16);
%  [T,H] = meshgrid(Tz,Hs); 
%  p = mk87cdf(hc,sc,H,T,upperTail);
%  v = 10.^(-6:-1);  
%  contourf(Tz,Hs,log10(p),log10(v))
%  xlabel('Tz')
%  ylabel('Hs')  
%  fcolorbar(log10(v))  
%  
% See also  mk87pdf, gaussq

%   References:
%   Myrhaug, D. and Kjelsen S.P. (1987) 
%  'Prediction of occurences of steep and high waves in deep water'.
%   Journal of waterway, Port, Coastal and Ocean Engineers,
%   Vol. 113, pp 122--138
%
%   Myrhaug & Dahle (1984) Parametric modelling of joint probability
%   density distributions for steepness and asymmetry in deep water 
%   waves
%      

%tested on matlab 5.2
%history:
% revised pab 09.08.2003
% Changed input + updated help header  
% revised pab 04.11.2000
% - no dependency on stats toolbox anymore.
% by  Per A. brodtkorb 1998

  
error(nargchk(3,5,nargin))  

if (nargin < 5||isempty(tail)),tail = 0; end
if (nargin < 4||isempty(Tz)),Tz = 8; end
if (nargin < 3||isempty(Hs)), Hs = 6; end

multipleSeaStates = any(numel(Hs)>1|numel(Tz)>1);
if multipleSeaStates
[icode, Scf,Hd,Tz,Hs] = iscomnsize(Scf,Hd,Tz,Hs);
else
  [icode, Scf,Hd] = iscomnsize(Scf,Hd);
end
if ~icode 
  error('Requires non-scalar arguments to match in size.');
end
  
Hrms = 0.715*Hs;
Erms = 0.0202 + 0.826*Hs./(Tz.^2); 

s = Scf./Erms;
hMax = 20;
h = min(Hd./Hrms,hMax);

eps2 = 1e-5;


p = zeros(size(Hd));
%k0 = find((Hs<=(Tz-4)*13/6+4));
%upLimit = 6.5/1.4;
loLimit = 2.5/1.26;
k0 = find((loLimit*sqrt(Hs)<Tz));
if any(k0)
  if multipleSeaStates
    h = h(k0);
    s = s(k0);
    Hs = Hs(k0);
    Tz = Tz(k0);
  else
    k0 = 1:numel(Hd);
  end
  hlim    = h;

  normalizedInput = 1;
  lowerTail = 0;
  
  
  
  if 0
    % This is a trick to get the html documentation correct.
    k = mk87pdf(1,1,2,3);
  end
  
  if (tail == lowerTail)
    %k       = find(h>2.5);
    hlim(h>2.5) = 2.5;
    p(k0) = gaussq(@mk87pdf,0,hlim,eps2/2,[],s,Hs,Tz,5,normalizedInput)...
	+ gaussq(@mk87pdf,hlim,h,eps2/2,[],s,Hs,Tz,5,normalizedInput); 
  else
    %k       = find(h<2.5);
    hlim(h<2.5) = 2.5;
    p(k0) = gaussq(@mk87pdf,h,hlim,eps2/2,[],s,Hs,Tz,7,normalizedInput)...
	+ gaussq(@mk87pdf,hlim,hMax,eps2/2,[],s,Hs,Tz,7,normalizedInput); 
  end
end
return

