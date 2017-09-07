function F = cdflomax(x,alpha,m0)
%CDFLOMAX CDF for local maxima for a zero-mean Gaussian process
% 
% CALL:  F = cdflomax(x,alpha,m0)
%
%       F     = distribution function evaluated at x
%       alpha = irregularity factor 
%       m0    = zero-order spectral moment (variance of the process)
%
%
% The cdf is calculated from an explicit expression involving the 
% standard-normal cdf. This relation is sometimes written as a convolution
%
%       M = sqrt(m0)*( sqrt(1-a^2)*Z + a*R )
%
% where  M  denotes local maximum, Z  is a standard normal r.v.,  
% R  is a standard Rayleigh r.v., and "=" means equality in distribution.
%
% Note that all local maxima of the process are considered, not
% only crests of waves. 
% 
% Example: 
%  S     = jonswap(10);
%  xs    = spec2sdat(S,10000);
%  mM    = tp2mm(dat2tp(xs)); 
%  m0    = spec2mom(S,1);  
%  alpha = spec2char(S,'alpha');  
%  x     = linspace(-10,10,200).';
%  plotedf(mM(:,2),[x,cdflomax(x,alpha,m0)])
% 
% See also  spec2mom, spec2bw

% Tested on Matlab 6.0
% History: 
%Revised pab 2008
% -renamed from lomaxcdf to cdflomax
% Revised pab Feb2004  
% -extended example  
% By jr 31.03.2001
  
%error(nargchk(3,3,nargin))
narginchk(3,3)
c1 = 1/(sqrt(1-alpha^2))*x./sqrt(m0);
c2 = alpha*c1;

F = cdfnorm(c1,0,1)-alpha*exp(-x.^2/2/m0).*cdfnorm(c2,0,1);
