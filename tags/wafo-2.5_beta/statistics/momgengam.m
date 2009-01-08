function [m,v,sk,ku]= momgengam(varargin)
%MOMGENGAM Mean and variance for the Generalized Gamma distribution.
% 
% CALL:  [m,v,sk,ku] = momgengam(a,b,c)
%
%   m, v = the mean and variance, respectively 
%  sk,ku = the skewness and kurtosis, respectively. 
%  a,b,c = parameters of the Generalized Gamma distribution. (see pdfgengam)
%
%  Mean (m) and variance (v) for the Generalized Gamma distribution is
%
%  m=b*gamma(a+1/c)/gamma(a)  and
%  v=b^2*(gamma(a+2/c)/gamma(a)-gamma(a+1/c)^2/gamma(a)^2);
%
% Example:
%   param = {2,2,0.5};N = 1000;  
%   x = rndgengam(param{:},N,1);
%   [mean(x),var(x),skew(x),kurt(x)]  
%   [m,v,sk,ku] = momgengam(param{:})
%
% See also  pdfgengam, cdfgengam, invgengam, rndgengam, fitgengam

% Reference: Cohen & Whittle, (1988) "Parameter Estimation in Reliability
% and Life Span Models", p. 220 ff, Marcel Dekker.

% Tested on; Matlab 5.3
% History: 
% revised pab Jul2004
% fixed a bug   
% reviseed pab Dec2003
% fixed abug: k1 -> k3
% revised pab 24.10.2000
%  - added comnsize, nargchk + default value for b and c
% added ms 09.08.2000

Np = 3;
error(nargchk(1,Np,nargin))
options = [];
params = parsestatsinput(Np,options,varargin{:});

[a,b,c] = deal(params{:});
if isempty(b)
  b = 1;
end
if isempty(c)
  c = 1;
end

a(a<=0) = nan;
b(b<=0) = nan;
c(c<=0) = nan;

try
  gammalna = gammaln(a);
  m = b.*exp(gammaln(a+1./c)-gammalna);
  v = b.^2.*exp(gammaln(a+2./c)-gammalna )-m.^2;
  % E(X^r) = b^r*gamma(1+r./c)/gamma(a)
  sk = (b.^3.*exp(gammaln(a+3./c)-gammalna)-3.*m.*v-m.^3)./v.^(3/2);
  ku = (b.^4*exp(gammaln(a+4./c)-gammalna)-4*sk.*v.^(3/2).*m-6*m.^2.*v-m.^4)/v.^2;
catch
    error('a b and c must be of common size or scalar.');
end



