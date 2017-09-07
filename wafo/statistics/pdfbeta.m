function f = pdfbeta(x,varargin)
%PDFBETA   Beta probability density function
%
% CALL:  f = pdfbeta(x,a,b,options);
%        f = pdfbeta(x,phat,options);
%
%    f = density function evaluated at x
%    x = matrix
% a, b = distribution shape parameters
% phat = Distribution parameter struct
%            as returned from FITBETA.  
%  options = struct with fieldnames:
%         .logp : if TRUE, density, p, returned as log(p).
%
%  The PDF is defined by:
%
%   f = x^(a-1)*(1-x)^(b-1)/H(a,b)    0<= x <= 1, a>0, b>0
% 
% where H(a,b) is a normalization constant.
% Example: 
%   x = linspace(-1,2,200);
%   p1 = pdfbeta(x,1,1); p2 = pdfbeta(x,2,2);
%   plot(x,p1,x,p2)
%
% See also cdfbeta, rndbeta, invbeta, fitbeta, mombeta

% Copyright (C) 1993 Anders Holtsberg, 2000 Per A. Brodtkorb
%
% This program is free software; you can redistribute it and/or modify
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

% tested on matlab 5.3
%History:
% revised pab sep 2007
% -removed dependence on comnsize
% -replaced beta() with betaln()
%revised pab 29.10.2000
% adapted from stixbox
% -added nargchk, comnsize
%       Anders Holtsberg, 18-11-93
%       Copyright (c) Anders Holtsberg


options = struct('logp',false); % default options
if ischar(x) && strcmpi(x,'defaults')
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

[a,b] = deal(params{:});
if isempty(b)
  error('Shape parameter b undefined!')
end

try
  a(a<=0) = nan;
  b(b<=0) = nan;
  %x(x<0 |1< x) = nan;
  %f = x.^(a-1) .* (1-x).^(b-1) ./ beta(a,b)
  alogx = (a-1).*log(x); 
  alogx( (a>0) & (x<0 | 1<x) ) = -inf;
  alogx( (a==1) & (x==0) ) = 0; 
  
  blog1px = (b-1).*log1p(-x);
  blog1px((b>0) & (x<0 | 1<x)) = -inf;
  blog1px((b==1) & (x==1)) = 0;
  logf = alogx+blog1px - betaloge(a,b);
  %logf =  (a-1).*log(x) + (b-1).*log1p(-x) - betaln(a,b);
 if options.logp
   f = logf;
 else
   f = exp(logf);
   %y = lpdfbeta(x,a,b);
   %plot(x,y,x,f,'.'),shg
   %semilogy(x,abs(y-f)./y,'-'),shg
 end
catch
   error('x, a and b must be of common size or scalar');
end


function y = lpdfbeta(x,a,b)
% This subfunction is based on c-code from   Catherine Loader (2000). 

% Reference
% Catherine Loader (2000). 
% "Fast and Accurate Computation of Binomial Probabilities"; 
% http://www.herine.net/stat/software/dbinom.html.
% @misc{ july-fast,
%   author = "Catherine Loader July",
%   title = "Fast and Accurate Computation of Binomial Probabilities",
%   url = "citeseer.ist.psu.edu/312695.html" }


a(a<=0) = nan;
b(b<=0) = nan;

x(x<=0) = realmin;
x(x>=1) = 1;

if (a<1)
  if (b<1)  %/* a<1, b<1 */
    f = a*b./((a+b)*x.*(1-x));
    p = pdfbin(a,a+b,x);
    
  else % /* a<1, b>=1 */
    f = a./x;
    p = pdfbin(a,a+b-1,x);
   end
  else
   if (b<1) %/* a>=1, b<1 */
     f = b./(1-x);
     p = pdfbin(a-1,a+b-1,x);
    
    else % /* a>=1, b>=1 */
     f = a+b-1;
      p = pdfbin(a-1,(a-1)+(b-1),x);
   end
end
y = p.*f;
 



% Old call kept just in case
% [errorcode x,a,b] = comnsize(x,a,b);
% if errorcode>0,
%   error('x, a and b must be of common size or scalar');
% end
% 
% f = zeros(size(x));
% 
% ok = (a>0 & b>0);
% 
% k = find(x>=0&x<=1 & ok);
% if any(k)
%   f(k) = x(k).^(a(k)-1) .* (1-x(k)).^(b(k)-1) ./ beta(a(k),b(k));  
%   % f(k) = exp((a(k)-1).*log(x(k))+ (b(k)-1).*log1p(-x(k)) - betaln(a(k),b(k)));
% end
% 
% k=find(~ok);
% if (any(k)),
%   f(k)=NaN;
% end
