function [f,options] = pdfweib(x,varargin)
%PDFWEIB Weibull probability density function
%
% CALL:  f = pdfweib(x,a,b,c,options);
%        f = pdfweib(x,phat,options);
%
%        f = density function evaluated at x
%        a = scale parameter    
%        b = shape parameter
%        c = location parameter (default 0)
%     phat = Distribution parameter struct
%            as returned from FITWEIB.  
%  options = struct with fieldnames:
%         .logp : if TRUE, density, p, returned as log(p).
%
% The Weibull distribution is defined by its cdf
%
%  F(x;a,b,c) = 1 -  exp(-((x-c)/a)^b), x>=0, a,b>0
%
%   Some references refer to the Weibull distribution with
%   a single parameter, this corresponds to PDFWEIB with a = 1.
%
% Example: 
%   x = linspace(0,6,200);
%   p1 = pdfweib(x,1,1); p2 = pdfweib(x,2,1); p3 = pdfweib(x,5,1);
%   plot(x,p1,x,p2,x,p3), shg
%
% See also cdfweib, invweib, rndweib, fitweib, momweib

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
% -added parameter for location
% -removed dep. on comnsize
% revised pab 24.10.2000
%  - added comnsize, nargchk
% rewritten ms 15.06.2000


%error(nargchk(2,6,nargin))
narginchk(2,6)
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


a(a<=0) = nan;
b(b<=0) = nan;
try
  xn = (x-c)./b;
  xn(xn<0) = inf; % trick to set pdf to zero
  limit1 = log(a./b);
  logf = -xn.^a +(a-1).*log(min(xn,realmax)) + limit1;
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
  error ('x, a, b and c must be of common size or scalar');
end



% [errorcode, x, a, c] = comnsize (x,a, c);
% if (errorcode > 0)
%   error ('x, a and c must be of common size or scalar');
% end
% 
% f=zeros(size(x));
% 
% ok = ((c > 0)  & (a > 0));
% 
% k = find (x>=0&ok);
% if any (k)  
%   f(k)=(x(k)./a(k)).^(c(k)-1).*c(k)./a(k).*exp(-(x(k)./a(k)).^c(k));
% end
% 
% k1 = find (~ok);
% if any (k1)
%   tmp=NaN;
%   f(k1) = tmp(ones(size(k1)));
% end
 
