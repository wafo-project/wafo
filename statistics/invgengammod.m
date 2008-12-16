function [x,xlo,xup] = invgengammod(F,varargin)
%INVGENGAMMOD Inverse of the Generalized Gamma distribution function
%
% CALL:  x = invgengam(F,m,v,L)
%
%        x = inverse cdf for the Generalized Gamma distribution evaluated at F
%        m = location parameter    (default 0)
%        v = scale parameter (v>0) (default 1)
%        L = shape parameter       (default 0)
%
% The modified Generalized Gamma distribution is defined by its cdf
%
%  F(x) = cdfnorm((log(x)-m)/sqrt(v))                    for x>=0, L==0
%         gammainc(exp(L*(log(x)-m)/sqrt(v))/L^2,1./L^2) for x>=0, L~=0.
%
% The parametrization used here is related to the parameters used 
% in PDFGENGAM by the following equations:
%     m = log(b)+1/c*log(a)
%     v = 1/(c^2*a)
%     L = 1/sqrt(a)
%
% Example:
%   x = linspace(0,1,200);
%   p1 = invgengam(x,0.5,1,1); p2 = invgengam(x,1,2,1);
%   p3 = invgengam(x,2,1,2); p4 = invgengam(x,2,2,2);
%   plot(x,p1,x,p2,x,p3,x,p4)
%
% See also  pdfgengammod, cdfgengammod, rndgengammod, fitgengammod,
%           momgengammod

% Reference: Cohen & Whittle, (1988) "Parameter Estimation in Reliability
% and Life Span Models", p. 220 ff, Marcel Dekker.
% 
% http://www.weibull.com/LifeDataWeb/generalized_gamma_distribution.htm


% Tested on; Matlab 5.3
% History: 
% adapted from stixbox ms 10.08.2000
% revised pab 23.10.2000
%  - added comnsize, nargchk 
% revised pab 4nov2005
%  -improved the starting guess for x.


error(nargchk(2,inf,nargin))
options = struct('covariance',[],'alpha',0.05,...
  'lowertail',true,'logp',false,'abseps',1e-90,'releps',sqrt(eps),'max_iter',100); % default options

Np = 3;
[params,options] = parsestatsinput(Np,options,varargin{:});
if numel(options)>1
  error('Multidimensional struct of distribution parameter not allowed!')
end

[m,v,L] = deal(params{:});
if isempty(m), m=0;end
if isempty(v), v=1;end
if isempty(L), L=0;end


[csize] = comnsize(F,m,v,L);
if any(isnan(csize))
  error('F, m, v and L must be of common size or scalar.');
end


a = 1./L.^2;
c = abs(L)./sqrt(v);
b = exp(m-log(a)./c);
epsilon = 2^-13; % defining L to zero
if isscalar(a) && isscalar(b) && isscalar(c)
  if abs(L)<epsilon
    x = invlognorm(F,m,v,options);
  else
    x = invgengam(F,a,b,c,options);
  end
else
  if isscalar(m), m = m(ones(csize));end
  if isscalar(v), v = v(ones(csize));end  
  if isscalar(L), L = L(ones(csize));end  
  if isscalar(a), a = a(ones(csize));end  
  if isscalar(b), b = b(ones(csize));end  
  if isscalar(c), c = c(ones(csize));end  
  if isscalar(F), F = F(ones(csize));end  
  x = zeros(csize);
  k1 = find(abs(L)<=epsilon);
  if any(k1)
    x(k1) = invlognorm(F(k1),m(k1),v(k1),options);
  end
  k1 = find(abs(L)>epsilon);
  if any(k1)
    x(k1) = invgengam(F(k1),a(k1),b(k1),c(k1),options);
  end
end
% x = invgam(F,a,1,options); 
% b(b<=0) = nan;
% c(c<=0) = nan;
% x = x.^(1./c).*b;

if nargout>1
% TODO % Implement  xlo and xup
 warning('WAFO:INVGENGAMMOD','xlo and xup not implemented yet')
 xlo = nan;
 xup = xlo;
end
