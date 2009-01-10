function [m,v,sk,ku]= mominvnorm(m0,l0)
%MOMINVNORM Mean and variance for the Inverse Gaussian distribution.
% 
% CALL:  [m,v,sk,ku] = mominvnorm(m0,l0)
%
%   m, v = the mean and variance, respectively 
%  sk,ku = the skewness and kurtosis, respectively. 
% m0, l0 = parameters of the Inverse Gaussian distribution (see pdfinvnorm)
%
%  Mean (m) and variance (v) for the Inverse Gaussian distribution is
%
%  m=m0  and  v=m0^3/l0;
%
% Example:
%   par = {1,1}
%   X = rndinvnorm(par{:},1000,1);
%   [mean(X) var(X),skew(X),kurt(X)]         % Estimated mean and variance
%   [m,v,sk,ku] = mominvnorm(par{:}) % True mean and variance
%
% See also pdfinvnorm, cdfinvnorm, invinvnorm, rndinvnorm

% Reference: Cohen & Whittle, (1988) "Parameter Estimation in Reliability
% and Life Span Models", p. 259 ff, Marcel Dekker.

% Tested on; Matlab 5.3
% History: 
% revised pab 24.10.2000
%  - added comnsize, nargchk
% added ms 14.08.2000

error(nargchk(2,2,nargin))
m0(m0<=0) = nan;
l0(l0<=0) = nan;
try
  m = m0;
  v = m0.^3./l0;
  sk = 3*sqrt(m0./l0);
  ku = 3+15*m0./l0;
  [csz,m] = comnsize(m,v);
catch
 error ('m and l must be of common size or scalar');
end

