function f = pdfray(x,varargin)
%PDFRAY Rayleigh probability density function
%
% CALL:  f = pdfray(x,b,options);
%        f = pdfray(x,phat,options);
%
%        f = density function evaluated at x
%        b = scale parameter
%     phat = Distribution parameter struct
%            as returned from FITRAY.  
%  options = struct with fieldnames:
%         .logp : if TRUE, density, p, returned as log(p).
% 
% The Rayleigh distribution is defined by its cdf
%
%  F(x;b) = 1 - exp(-x^2/(2b^2)), x>=0, b>0
%
% Example: 
%   x = linspace(0,4,200);
%   p1 = pdfray(x,1); p2 = pdfray(x,0.5);
%   plot(x,p1,x,p2), shg
% 
% See also cdfray, invray, rndray, fitray, momray

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
% and Life Span Models", p. 181 ff, Marcel Dekker.


% Tested on; Matlab 5.3
% History: 
% revised pab 2007
% -removed dependence on comnsize
% revised pab 24.10.2000
%  - added comnsize, nargchk
% added ms 15.06.2000


%error(nargchk(2,4,nargin))
narginchk(2,4)
Np = 1;
options = struct('logp',false); % default options
[params,options] = parsestatsinput(Np,options,varargin{:});
if numel(options)>1
  error('Multidimensional struct of distribution parameter not allowed!')
end
b = params{1};

x(x<0) = 0;
b(b<0) = nan;
try
  xn = min(x./b,realmax);
  if options.logp
    f=-0.5*xn.^2 + log(xn)-log(b);
  else
    f=xn.*exp(-0.5*xn.^2)./b;
  end
catch
    error ('x and b must be of common size or scalar');
end
% 
% [csize, x, b] = comnsize (x,b);
% if any(isnan(csize))
%   error ('x and b must be of common size or scalar');
% end
% 
% f=zeros(size(x));
% 
% k = find ((x>=0)&(b>0));
% if any (k)  
%   f(k)=x(k).*exp(-x(k).^2./(2*b(k).^2))./b(k).^2;
% end
% 
% k1 = find (b<=0);
% if any (k1)
%   tmp=NaN;
%   f(k1) = tmp(ones(size(k1)));
% end
% 
