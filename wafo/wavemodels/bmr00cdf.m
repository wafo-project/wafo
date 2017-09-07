function [p, eps2] = bmr00cdf(Hd,Scf,Hs,Tz,tail)
%BMR00CDF Brodtkorb et.al (2000) joint (Scf,Hd) CDF from North Sea.
%
% CALL: F = bmr00cdf(Hd,Scf,Hs,Tz,tail)
%
%    F  = CDF
%   Hd  = zero down crossing wave height.
%   Scf = crest front steepness.
%   Hs  = significant wave height.
%   Tz  = average zero down crossing period.
%  tail = 1 if upper tail is calculated   
%         0 if lower tail is calulated (default)
%
% BMR00CDF returns the joint CDF of (Scf, Hd) given Hs and Tz,
% i.e., crest front steepness (2*pi*Ac/(g*Td*Tcf)) and wave height, given
% the seastate. The root mean square values of Hd and Scf (Hrms,Erms) are
% related to the significant waveheight and the average zero down
% crossing period by:
%             Hrms = Hs/sqrt(2);
%             Erms = 5/4*Hs/(Tz^2);
%  
%   The size of F is the common size of  the input arguments
%
% Examples: 
%  Hs = 5.5;
%  Tz = 8.5;
%  sc = 0.25;
%  hc = 4;
%  p = bmr00cdf(hc,sc,Hs,Tz) % Prob(Hd<=hc,Scf<=sc|Hs,Tz) = 0.68
%  
%  % Conditional probability of steep and high waves given seastates
%  % i.e., Prob(Hd>hc,Scf>sc|Hs,Tz)  
%  upperTail = 1;
%  Hs = linspace(2.5,18.5,17);
%  Tz = linspace(4.5,19.5,16);
%  [T,H] = meshgrid(Tz,Hs); 
%  p = bmr00cdf(hc,sc,H,T,upperTail);
%  v = 10.^(-6:-1);  
%  contourf(Tz,Hs,log10(p),log10(v))
%  xlabel('Tz')
%  ylabel('Hs')  
%  fcolorbar(log10(v))  
%  
% See also  mk87cdf, gaussq

% Reference 
%  Brodtkorb, P.A. and Myrhaug, D. and Rue, H. (2000)
%  "Joint Distributions of Wave Height and Wave Steepness Parameters",
% In Proc. 27'th Int. Conf. on Coastal Eng., ICCE, Sydney, Australia },
%  vol. 1, pp. 545--558, Paper No. 162


%tested on matlab 5.2
%history:
% by  Per A. brodtkorb July 2004

  
%error(nargchk(3,5,nargin))  
narginchk(3,5)
if (nargin < 5||isempty(tail)),tail = 0; end
if (nargin < 4||isempty(Tz)),Tz = 8; end
if (nargin < 3||isempty(Hs)), Hs = 6; end

multipleSeaStates = any(numel(Hs)>1||numel(Tz)>1);
if multipleSeaStates
  [csize, Scf,Hd,Tz,Hs] = comnsize(Scf,Hd,Tz,Hs);
else
  [csize, Scf,Hd] = comnsize(Scf,Hd);
end
if any(isnan(csize))
  error('Requires non-scalar arguments to match in size.');
end
  
Hrms = Hs/sqrt(2);
Erms = 5/4*Hs./(Tz.^2); 
%Erms = (0.0202 + 0.826*Hs./(Tz.^2));

s = Scf./Erms;
hMax = 20;
h = min(Hd./Hrms,hMax);

eps2 = 1e-3;


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

  h0 = 3.79937912772197;
  normalizedInput = 1;
  lowerTail = 0;
  
   
  if 0
    % This is a trick to get the html documentation correct.
    k = bmr00pdf(1,1,2,3);
  end
  
  
  if (tail == lowerTail)
    %k       = find(h>h0);
    hlim(h>h0) = h0;
    p(k0) = gaussq(@bmr00pdf,0,hlim,eps2/2,[],s,Hs,Tz,5,normalizedInput)...
	+ gaussq(@bmr00pdf,hlim,h,eps2/2,[],s,Hs,Tz,5,normalizedInput); 
  else
    %k       = find(h<h0);
    hlim(h<h0) = h0;
    p(k0) = gaussq(@bmr00pdf,h,hlim,eps2/2,[],s,Hs,Tz,7,normalizedInput)...
	+ gaussq(@bmr00pdf,hlim,hMax,eps2/2,[],s,Hs,Tz,7,normalizedInput); 
  end
end
return

