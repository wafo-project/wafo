function [m,v,sk,ku]= momray(b)
%MOMRAY Mean and variance for the Rayleigh distribution.
% 
% CALL:  [m,v,sk,ku] = momray(b)
%
%   m, v = the mean and variance, respectively 
%  sk,ku = the skewness and kurtosis, respectively. 
%      b = parameter of the Rayleigh distribution (see cdfray)
%
%  Mean (m) and variance (v) for the Rayleigh distribution is
%
%  m=b*(pi/2)^(1/2)  and  v=(2-pi/2)*b^2;
%
% Example:
%   par = {3}
%   X = rndray(par{:},1000,1);
%   [mean(X) var(X),skew(X),kurt(X)]      % Estimated mean and variance
%   [m,v,sk,ku] = momray(par{:}) % True mean and variance
%
% See also pdfray, cdfray, invray, rndray, fitray


% Reference: Cohen & Whittle, (1988) "Parameter Estimation in Reliability
% and Life Span Models", p. 181 ff, Marcel Dekker.

%tested on: matlab 5.3
% history
% revised pab 24.10.2000
% Added checks on b



error(nargchk(1,1,nargin))

b(b<=0) = nan;

m  = b * sqrt(pi/2);
v  = (2 - pi/2) * b .^ 2;
sk = 2*sqrt(pi)*(pi - 3)/(4-pi)^(3/2);
ku = 3-(6*pi^2 - 24*pi +16)/((4-pi)^2);

sk = sk(ones(size(v)));
ku = ku(ones(size(v)));
sk(isnan(b)) = nan;
ku(isnan(b)) = nan;

