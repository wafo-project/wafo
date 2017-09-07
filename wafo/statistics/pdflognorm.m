function f = pdflognorm(x,varargin)
%PDFLOGNORM Lognormal probability density function
%
% CALL:  f = pdflognorm(x,m,v,options);
%        f = pdflognorm(x,phat,options);
%
%        f = density function evaluated at x
%        m = mean of log(x)     (default 0)
%        v = variance of log(x) (default 1)
%     phat = Distribution parameter struct
%            as returned from FITLOGNORM.  
%  options = struct with fieldnames:
%         .logp   : if TRUE, density, p, returned as log(p).
%
% The Lognormal distribution is defined by its pdf
%
%     f(x) = exp(-(log(x)-m)^2/(2*v))/sqrt(v*2*pi*x^2), x>=0.
%
% Example: 
%   x = linspace(0,3,200);
%   p1 = pdflognorm(x,0,1); p2 = pdflognorm(x,.5,0.25);
%   plot(x,p1,x,p2), shg
%
% See also cdflognorm, invlognorm, rndlognorm, fitlognorm, momlognorm

%
%     This program is free software; you can redistribute it and/or modify
%     it under the terms of the GNU Lesser General Public License as published by
%     the Free Software Foundation; either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU Lesser General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.



% Reference: Cohen & Whittle, (1988) "Parameter Estimation in Reliability
% and Life Span Models", p. 59 ff, Marcel Dekker.

% Tested on; Matlab 5.3
% History: 
% revised pab 24.10.2000
%  - added comnsize, nargchk
%  - added default values
%  - fixed a bug in the parameterization
% added ms 10.08.2000

%error(nargchk(1,6,nargin))
narginchk(1,6)
options = struct('logp',false); % default options
Np = 2;
[params,options] = parsestatsinput(Np,options,varargin{:});
if numel(options)>1
  error('Multidimensional struct of distribution parameter not allowed!')
end

[m,v] = deal(params{:});
if isempty(m),  m=0;  end
if isempty(v),  v=1;  end


v(v<=0) = nan;
x(x<=0) = inf; % trick to set pdf to zero
try
  s = sqrt(v);
  xn = (log(x)-m)./s;
  
  if options.logp
    f = -0.5*xn.^2 - log(x)-log(s.*sqrt(2*pi));
  else
    f = exp(-0.5*xn.^2)./(x.*s.*sqrt(2*pi));
  end
catch
  error ('x, m and v must be of common size or scalar');
end


% f=zeros(size(x));
% 
% k = find (x>0&v>0);
% if any(k)    
%   f(k)=1./sqrt(2*pi*v(k)).*exp(-0.5*(log(x(k))-m(k)).^2./v(k))./x(k);
% end
% 
% k1 = find (v<=0);
% if any (k1)
%   tmp=NaN;
%   f(k1) = tmp(ones(size(k1)));
% end




