function f = pdff(x,varargin)
%PDFF Snedecor's F probability density function
%
% CALL:  f = pdff(x,df1,df2,options);
%        f = pdff(x,phat,options);
%
%        f = PDF evaluated at x
%        x = matrix
%  df1,df2 = degrees of freedom (1,2,....)
%     phat = Distribution parameter struct
%            as returned from WFFIT.  
%  options = struct with fieldnames:
%         .logp   : if TRUE, density, p, returned as log(p).
%         .disable: if TRUE disable check on integer DF.
%
% The F PDF is defined as
%
%  f(x;df1,df2) = (1+df2/(df1*x))^(-df1/2)*(1+df1/(df2*x))^(-df2/2)/(x*B)
% 
%  where B = beta(df1/2,df2/2)
%
% Example:
%   x  = linspace(0,6,200);
%   p1 = pdff(x,1,1); p2 = pdff(x,2,2);
%   plot(x,p1,x,p2),shg
%
% See also cdff, invf, rndf, wffit, momf

% tested on matlab 5.3
%History:
% revised pab sep 2007
% -removed dependence on comnsize
% -new implementation -> tackles large df1 and df2 better than before
%revised pab 29.10.2000
% adapted from stixbox
% -added nargchk, comnsize + a check that df1,df2 are positive integers
%        Anders Holtsberg, 18-11-93
%        Copyright (c) Anders Holtsberg

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


options = struct('logp',false,'disable',false); % default options
if ischar(x) && strcmpi(x,'defaults')
  f = options;
  return
end

error(nargchk(2,inf,nargin))
Np = 2;

[params,options] = parsestatsinput(Np,options,varargin{:});
if numel(options)>1
  error('Multidimensional struct of distribution parameter not allowed!')
end

a = params{1};
b = params{2};
if options.disable
  a(a<=0) = nan;
  b(b<=0) = nan;
else
  a(floor(a)~=a | a<=0) = nan;
  b(floor(b)~=b | b<=0) = nan;
end
%x(x<0) = inf;
try 
  c = (b./a);
  xx  = sign(x).*exp(-log1p(c./abs(x)));
catch
  error('x, df1 and df2 must be of common size or scalar');
end
if any(abs(a(:)*eps)>abs(b(:)) | abs(a(:))>realmax/2) %any(abs(k(:)*eps)>2) 
  warning('WAFO:PDFF','Result may not be exact.')
end



%   c = b./a;
%   xx = x./(x+c);
%   tmp = pdfbeta(xx,a/2,b/2);
%   f = tmp./(x+c).^2.*c;
f = pdfbeta(xx,a/2,b/2,'logp',true)-log(c)-2*log1p(x./c);
lrg = sqrt(eps)*c>1 & (x>=0) ;
if any(lrg)
  % c-> inf => a*X -> pdfchi2(a) 
  if ~isscalar(x), x = x(lrg);end
  if ~isscalar(a), a = a(lrg);end
  f(lrg) =log(a)+pdfchi2(a.*x,a,'disable',options.disable,'logp',true);
end




% Check
% y = lpdff(x,a,b,true);
% y2 = log(lpdff(1./x,b,a));
% z = log(a)+pdfchi2(a.*x,a,'disable',options.disable,'logp',true);
% z2 = log(b.*pdfchi2(b./x,b,'disable',options.disable));
% clf,plot(x,f,x,y,'r.'),shg
%  semilogy(x,abs(z-f)./abs(f)),shg
% semilogy(x,abs(f-y)./abs(y)),shg

if ~options.logp
   f = exp(f);
  
end
return


function y = lpdff(x,m,n,logp)

% Reference
% Catherine Loader (2000). 
% "Fast and Accurate Computation of Binomial Probabilities"; 
% http://www.herine.net/stat/software/dbinom.html.% @misc{ july-fast,
%   author = "Catherine Loader July",
%   title = "Fast and Accurate Computation of Binomial Probabilities",
%   url = "citeseer.ist.psu.edu/312695.html" }

n(n<=0) = nan;
m(m<=0) = nan;
x(x<0) = 0;

try  
  f = (n+x.*m);
  %q = n./f;
  p = x.*m./f;
catch
  error('Size mismatch')
end
  if   (m>=2)
    %log(q)-log1p(-p)
    f = log(m/2) + log1p(-p);
    dens = pdfbin( (m-2)/2,(m+n-2)/2, p,'logp',true);
  else
   %f = m*m*q ./ (2*p*(m+n));
   f = log(m./(1+n./m)/2)+log1p(-p) -log(p);
    dens = pdfbin(m/2, (m+n)/2, p,'logp',true);
  end
  if logp
    y = dens+f;
  else
    y =  exp(dens+f);
  end