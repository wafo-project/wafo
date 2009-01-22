function [m,v,sk,ku]= momweib(a,c)
%MOMWEIB Mean and variance for the Weibull  distribution.
% 
% CALL:  [m,v,sk,ku] = momweib(a,c)
%
%   m, v = the mean and variance, respectively 
%  sk,ku = the skewness and kurtosis, respectively. 
%   a, c = parameters of the Weibull distribution (see cdfweib).
%
%  Mean (m) and variance (v) for the Weibull distribution is
%
%  m=a*gamma(1+1/c)  and  v=a^2*gamma(1+2/c)-m^2;
%
% Example:
%   par = {1,2}
%   X = rndweib(par{:},10000,1);
%   [mean(X) var(X),skew(X),kurt(X)]        % Estimated mean and variance
%   [m,v,sk,ku] = momweib(par{:}) % True mean and variance
%
% See also pdfweib, cdfweib, invweib, rndweib,  fitweib



% Reference: Cohen & Whittle, (1988) "Parameter Estimation in Reliability
% and Life Span Models", p. 25 ff, Marcel Dekker.


% Tested on; Matlab 5.3
% History: 
% revised pab 23.10.2000
%   - added comnsize
% added ms 15.06.2000

error(nargchk(2,2,nargin))

a(a<=0) = nan;
c(c<=0) = nan;
try
  m =  a .* gamma(1 + (1 ./ c));
  v = a.^ 2 .* gamma(1 + (2 ./ c)) - m.^ 2;
  sk = (a.^3.*gamma(1+(3./c))-3.*m.*v-m.^3)./v.^(3/2);
  ku = (a.^4*gamma(1+4./c)-4*sk.*v.^(3/2).*m-6*m.^2.*v-m.^4)/v.^2;
catch
  error('a and c must be of common size or scalar.');
end
