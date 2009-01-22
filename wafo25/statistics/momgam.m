function [m,v,sk,ku]= momgam(varargin)
%MOMGAM Mean and variance for the Gamma distribution.
% 
% CALL:  [m,v,sk,ku] = momgam(a,b)
%
%   m, v = the mean and variance, respectively 
%  sk,ku = the skewness and kurtosis, respectively. 
%      a = parameter, a>0
%      b = parameter, b>0 (default b=1)
%
%  Mean (m) and variance (v) for the Gamma distribution is
%
%  m=ab  and  v=ab^2;
%
% Example:
%   par = {2,1}
%   X = rndgam(par{:},100000,1);
%   [mean(X) var(X),skew(X),kurt(X)]         % Estimated mean and variance
%   [m,v,sk,ku] = momgam(par{:}) % True mean and variance
%
% See also pdfgam, cdfgam, invgam, rndgam, fitgam

% Reference: Johnson, Kotz and Balakrishnan (1994)
% "Continuous Univariate Distributions, vol. 1", p. 494 ff
% Wiley

% Tested on; Matlab 5.3
% History: 
% revised pab Dec2003
% fixed a bug: k1 -> k3
% revised pab 24.10.2000
%  - added comnsize, nargchk
% added ms 26.06.2000
% added b parameter ms 23.08.2000
Np = 2;
error(nargchk(1,Np,nargin))
options = [];
params = parsestatsinput(Np,options,varargin{:});

[a,b] = deal(params{:});
if isempty(b)
  b = 1;
end


a(a<=0) = nan;
b(b<=0) = nan;
try
  m = a.*b;
  v = m.*b;
  sk = 2/sqrt(a);
  ku = 3+6/a;
catch
  error('a and b must be of common size or scalar.');
end
