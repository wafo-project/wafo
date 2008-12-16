function f = pdfexp(x,varargin)
%PDFEXP Exponential probability density function
%
% CALL:  f = pdfexp(x,m);
%        f = pdfexp(x,phat);
%
%        f = density function evaluated at x
%        m = mean
%     phat = Distribution parameter struct
%            as returned from FITEXP.  
%  options = struct with fieldnames:
%         .logp : if TRUE, density, p, returned as log(p).
%
% The Exponential distribution is defined by its pdf
%
%        f(x)=exp(-x/m)/m, x>=0, m>0.
% 
% Example: 
%   x = linspace(0,6,200);
%   p1 = pdfexp(x,1); p2 = pdfexp(x,2);
%   plot(x,p1,x,p2),shg 
%
% See also  cdfexp, invexp, rndexp, fitexp, momexp

% Reference: Johnson, Kotz and Balakrishnan (1994)
% "Continuous Univariate Distributions, vol. 1", p. 494 ff
% Wiley

% Tested on; Matlab 5.3
% History:
% revised removed dependence on comnsize
%  Revised pab Dec2003
% fixed abug k1->k3
% revised pab 24.10.2000
%  - added comnsize, nargchk
% added ms 15.06.2000

error(nargchk(2,inf,nargin))
Np = 1;
options = struct('logp',false); % default options
[params,options] = parsestatsinput(Np,options,varargin{:});
if numel(options)>1
  error('Multidimensional struct of distribution parameter not allowed!')
end

m = params{1};

try
  m(m<=0) = nan;
  xn = x./m;
  xn(xn<0) = inf;
  if options.logp
    f = -xn-log(m);
  else
    f = exp(-xn)./m;
  end
catch
  error('x and m must be of common size or scalar.');
 end
% 
% % Initialize f to zero.
% f = zeros(size(x));
% 
% k=find(x >= 0 & m>0);
% if any(k),
%   f(k)=exp(-x(k)./m(k))./m(k);
% end
  
% [icode x m] = comnsize(x,m);
% if ~icode > 0
%     error('x and m must be of common size or scalar.');
% end
% 
% % Initialize f to zero.
% f = zeros(size(x));
% 
% k=find(x >= 0 & m>0);
% if any(k),
%   f(k)=exp(-x(k)./m(k))./m(k);
% end
% 
% k3 = find(m<=0);     
% if any(k3)
%   tmp = NaN;
%   f(k3) = tmp(ones(size(k3)));
% end
% 
% 
