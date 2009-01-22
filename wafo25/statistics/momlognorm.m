function [m,v,sk,ku]= momlognorm(m0,v0)
%MOMLOGNORM Mean and variance for the Lognormal distribution.
% 
% CALL:  [m,v,sk,ku] = momlognorm(m0,v0)
%
%   m, v = the mean and variance, respectively 
%  sk,ku = the skewness and kurtosis, respectively. 
% m0, v0 = parameters of the Lognormal distribution.
%
%  Mean (m) and variance (v) for the Lognormal distribution is
%
%  m=exp(m0+v0/2)  and  v=exp(2*m0+v0)*(exp(v0)-1);
%
% Example:
%   par = {-1,1}
%   X = rndlognorm(par{:},1000,1);
%   [mean(X) var(X),skew(X),kurt(X)] % Estimated mean and variance
%   [m,v,sk,ku] = momlognorm(par{:})            % True values
%
% See also pdflognorm, cdflognorm, invlognorm, rndlognorm, fitlognorm


% Reference: Cohen & Whittle, (1988) "Parameter Estimation in Reliability
% and Life Span Models", p. 59 ff, Marcel Dekker.



% Tested on; Matlab 5.3
% History:
% revised pab 24.10.2000
%  - added comnsize, nargchk
%  - added default values
% added ms 10.08.2000

error(nargchk(0,2,nargin))
if nargin<1||isempty(m0),  m0=0;  end
if nargin<2||isempty(v0),  v0=1;  end

v0(v0<0) = nan;
try
  m = exp(m0+v0/2);
  v = exp(2*m0+v0).*(expm1(v0));
  sk = (exp(v0)+2).*sqrt(expm1(v0));
  ku = exp(4*v0)+2*exp(3*v0)+3*exp(2*v0)-3;
catch
  error ('m and v must be of common size or scalar');
end

