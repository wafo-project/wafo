function gam = getjonswappeakedness(Hm0,Tp)
%GETJONSWAPPEAKEDNESS Peakedness factor Gamma given Hm0 and Tp for JONSWAP
%
%  CALL: Gamma =  getJonswapPeakedness(Hm0,Tp);
%
% Hm0   = significant wave height [m].
% Tp    = peak period [s]
% gamma = Peakedness parameter of the JONSWAP spectrum
%
%  GETJONSWAPPEAKEDNESS relate GAMMA to Hm0 and Tp.
%  A standard value for GAMMA is 3.3. However, a more correct approach is 
%  to relate GAMMA to Hm0 and Tp:
%        D = 0.036-0.0056*Tp/sqrt(Hm0);
%        gamma = exp(3.484*(1-0.1975*D*Tp^4/(Hm0^2)));
%  This parameterization is based on qualitative considerations of deep water
%  wave data from the North Sea, see Torsethaugen et. al. (1984)
%  Here GAMMA is limited to 1..7.
% 
%  NOTE: The size of GAMMA is the common size of Hm0 and Tp.
%  
% Example
% assert(getjonswappeakedness(7,[5, 11, 25]), [7, 2.38529836797459   1], 1e-10)
% Hm0 = linspace(1,20);
% Tp = Hm0;
% [T,H] = meshgrid(Tp,Hm0);
% gam = getjonswappeakedness(H,T);
% contourf(Tp,Hm0,gam,1:7),fcolorbar(1:7)
%
% Hm0 = 1:10;  
% Tp  = linspace(2,16);
% [T,H] = meshgrid(Tp,Hm0);
% gam =  getjonswappeakedness(H,T); 
% plot(Tp,gam)
% xlabel('Tp [s]')  
% ylabel('Peakedness parameter')
%
% close all
%  
% See also jonswap  
  
% Tested on matlab 6.1
%History
% revised pab 13april2004
% -made sure gamma is 1 if Tp/sqrt(Hm0) > 5.1429
% by pab 11Jan2004  
  
%  error(nargchk(2,2,nargin))
  narginchk(2,2)
  [csize,Hm0,Tp] = comnsize(Hm0,Tp);
  if any(isnan(csize))
    error('Requires non-scalar arguments to match in size.');
  end
  x   = Tp./sqrt(Hm0);
  
  gam = ones(size(x));
  
  k1 = find(x<=5.14285714285714);
  if any(k1), %limiting gamma to [1 7]
    D       = 0.036-0.0056*x(k1); % approx 5.061*Hm0^2/Tp^4*(1-0.287*log(gam));
    gam(k1) = min(exp(3.484*( 1-0.1975*D.*x(k1).^4 ) ),7); % gamma 
  end
  return