function f = pdfweibmod(x,varargin)
%PDFWEIBMOD Truncated Weibull probability density function
%
% CALL:  f = pdfweibmod(x,a,b,c);
%
%        f = density function evaluated at x
%        a = scale parameter      
%        b = shape parameter
%        c = truncation parameter (default 0)
%     phat = Distribution parameter struct
%            as returned from FITWEIBMOD.  
%  options = struct with fieldnames:
%         .logp : if TRUE, density, p, returned as log(p).
%
%
% The Truncated Weibull distribution is defined by its cdf
%
%  F(x;a,c) = 1 -  exp(-((x+|c|)/a)^b+(|c|/a)^b), x>=0
%
% Example: 
%   x = linspace(0,6,200);
%   p1 = pdfweibmod(x,1,1,2); p2 = pdfweibmod(x,2,2,3);
%   plot(x,p1,x,p2), shg
%
% See also cdfweibmod, invweibmod, rndweibmod, fitweibmod, momweibmod

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
% and Life Span Models", p. 25 ff, Marcel Dekker.


% Tested on; Matlab 5.3
% History: 
% -revised pab 2007
% revised pab 12.11.2000
%  - added check on f(k): Replace NaN's with zero 
% revised pab 24.10.2000
%  - added comnsize, nargchk
% rewritten ms 15.06.2000


error(nargchk(2,6,nargin))
Np = 3;
options = struct('logp',false); % default options
[params,options] = parsestatsinput(Np,options,varargin{:});
if numel(options)>1
  error('Multidimensional struct of distribution parameter not allowed!')
end

[b,a,c] = deal(params{:});

if isempty(c)
  c =0;
end
if isempty(b)
  b =1;
end

a(a<=0) = nan;
b(b<=0) = nan;
try
  c = abs(c);
  x(x<0) = inf; % trick to set pdf to zero
  xn = (x+c)./b;
  %f = exp(-xn.^a+(c./b).^a).*xn.^(a-1).*a./b;
  limit1 = (c./b).^a + log(a./b);
  logf   = -xn.^a +(c./b).^a +(a-1).*log(min(xn,realmax)) + log(a./b);
  spcase =  (a==1) & (xn==0);
  if any(spcase)
    if ~isscalar(limit1)
      limit1 = limit1(spcase);
    end
    logf(spcase) = limit1;
  end  
  if options.logp
    f = logf;
  else
    f = exp(logf);
  end
catch
  error ('x, a,b, and c must be of common size or scalar');
end

return

% [errorcode, x, a,b, c] = comnsize (x,a,b, abs(c));
% if (errorcode > 0)
%   error ('x, a, b and c must be of common size or scalar');
% end
% 
% f=zeros(size(x));
% 
% ok = ((b > 0)  & (a > 0));
% 
% k = find (x>=0&ok);
% if any (k)  
%   f(k)=((x(k)+c(k))./a(k)).^(b(k)-1).*b(k)./a(k).*exp(-((x(k)+c(k))./a(k)).^b(k)+abs(c(k)./a(k)).^b(k));
%   k0 = find(isnan(f(k)));
%   if any(k0), 
%     f(k(k0))=0; 
%   end
% end
% 
% k1 = find (~ok);
% if any (k1)
%   tmp=NaN;
%   f(k1) = tmp(ones(size(k1)));
% end
