function f = pdfinvnorm(x,varargin)
%PDFINVNORM Inverse Gaussian probability density function
%
% CALL:  f = pdfinvnorm(x,m,l,options);
%        f = pdfinvnorm(x,phat,options);
%
%        f = density function evaluated at x
%      m,l = parameters
%     phat = Distribution parameter struct
%            as returned from FITINVNORM.  
%  options = struct with fieldnames:
%         .logp   : if TRUE, density, p, returned as log(p).
%
% The Inverse Gaussian distribution is defined by its pdf
%
%        f(x)=(l/(2*pi*x^3))^(1/2)*exp(-l*(x-m)^2/(2*m^2*x)), x>0.
%
% Example: 
%   x = linspace(0,3,200);
%   p1 = pdfinvnorm(x,1,1); p2 = pdfinvnorm(x,1,.25);
%   plot(x,p1,x,p2), shg
%
% See also cdfinvnorm, invinvnorm, rndinvnorm, fitinvnorm, mominvnorm

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
% and Life Span Models", p. 259 ff, Marcel Dekker.


% Tested on; Matlab 5.3
% History:
% revised pab 24.10.2000
%  - added comnsize, nargchk
% added ms 14.08.2000

options = struct('logp',false); % default options
if (nargin==1 && nargout <= 1 && isequal(x,'defaults'))
  f = options; 
  return
end
%error(nargchk(2,inf,nargin))
narginchk(2,inf)
Np = 2;
[params,options] = parsestatsinput(Np,options,varargin{:});
if numel(options)>1
  error('Multidimensional struct of distribution parameter not allowed!')
end

[m,L] = deal(params{:});
m(m<=0) = nan;
L(L<=0) = nan;
x(x<=0) = inf;
try 
  sqx = sqrt(x);
  logf= -0.5*L.*(sqx./m-1./sqx).^2 +0.5*(log(L/(2*pi))-3*log(x));
  
  if options.logp
    f = logf;
  else
    f = exp(logf);
  end
  %f=(L./(2*pi*x.^3)).^(1/2).*exp(-L.*((x-m)./m).^2./(2*x));
catch
  error('x, m and l must be of common size or scalar');
end

return
% 
% [iscmn, x, m, l] = iscomnsize (x,m, L);
% if ~iscmn
%   error ('x, m and l must be of common size or scalar');
% end
% 
% f=zeros(size(x));
% ok=((m>0)&(l>0));
% k = find (x>0&ok);
% if any(k)  
%   f(k)=(l(k)./(2*pi*x(k).^3)).^(1/2).*exp(-l(k).*(x(k)-m(k)).^2./(2*m(k).^2.*x(k)));
% end
% 
% k1 = find (~ok);
% if any (k1)
%   f(k1) = NaN;
% end
% 

