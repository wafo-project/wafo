function f = pdfgev(x,varargin)
%PDFGEV Generalized Extreme Value probability density function
%
% CALL:  f = pdfgev(x,k,s,m, options);
%
%        f = density function evaluated at x
%        k = shape parameter in the GEV 
%        s = scale parameter in the GEV, s>0  (default 1) 
%        m = location parameter in the GEV    (default 0)
%     phat = Distribution parameter struct
%            as returned from FITGEV.  
%  options = struct with fieldnames:
%         .logp   : if TRUE, density, p, returned as log(p).
% 
% The Generalized Extreme Value distribution is defined by its cdf
%
%                exp( - (1 - k(x-m)/s)^1/k) ),  k~=0
%  F(x;k,s,m) =
%                exp( -exp(-(x-m)/s) ),  k==0
%
% for k*(x-m)<s (when k~=0).
%
% Example: 
%  x = linspace(-5,5,200);
%  p1 = pdfgev(x,1.25,1); p2 = pdfgev(x,1,1);
%  p3 = pdfgev(x,0.75,1); p4 = pdfgev(x,0.5,1);
%  p5 = pdfgev(x,0,1);    p6 = pdfgev(x,-0.5,1);
%  p7 = pdfgev(x,-0.75,1);p8 = pdfgev(x,-1,1);
%
% subplot(211), plot(x,p1,x,p2,x,p3,x,p4), title('k>0')
% subplot(212), plot(x,p5,x,p6,x,p7,x,p8), title('k<=0'), shg
%
% See also cdfgev, invgev, rndgev, fitgev, fitgev

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


% References
%  Johnson  N.L., Kotz S. and Balakrishnan, N. (1994)
%  Continuous Univariate Distributions, Volume 1. Wiley. 

% Tested on; Matlab 5.3
% History: 
% revised pab 16 feb 2006
%  - fixed a bug for k==0
% revised jr 14.08.2001
% - a bug in the last if-statement condition fixed
%   (thanks to D Eddelbuettel <edd@debian.org>)
% revised pab 24.10.2000
% - added  nargchk, comnsize and default values for m, s
% added ms 14.06.2000

error(nargchk(2,inf,nargin))
Np = 3;
options = struct('logp',false,'disable',false); % default options
[params,options] = parsestatsinput(Np,options,varargin{:});
if numel(options)>1
  error('Multidimensional struct of distribution parameter not allowed!')
end
[k,s,m] = deal(params{:});
if isempty(m), m=0;end
if isempty(s), s=1;end

s(s<=0) = NaN;
try
  xn = (x-m)./s;
  kxn = k.*xn;
  kxn(kxn>1) = 1;
  kxn(k==0 & xn==inf) = 0;
  log1mkxn = log1pxdx(-kxn);
  logfs = -exp(-xn.*log1mkxn) + (-xn+kxn).*log1mkxn;
  
  % Handle special cases
  logfs(kxn==-inf | kxn==1) = -inf;      %  
  logfs(k==1 & xn==1) = 0; % 0^0 situation
  
  if options.logp
    f = logfs-log(s);
  else
    f = exp(logfs)./s;
  end
catch
  error('x, k, s and m must be of common size or scalar.');
end
  

return


% old code kept just in case
% s(s<=0) = NaN;
% 
% [icode x k s,m] = comnsize(x,k,s,m);
% if ~icode
%   error('x, k, s and m must be of common size or scalar.');
% end
%   
% epsilon=1e-4; % treshold defining k to zero
% 
% 
% f = zeros(size(x));
% k0 = find( abs(k)<=epsilon & s>0);
% if any(k0),
%   tmp=exp(-(x(k0)-m(k0))./s(k0));
%   f(k0) = exp(-tmp).*tmp./s(k0);
% end
% 
% k1=find(((k.*(x-m)<s)&(abs(k)>epsilon)) & s>0);
% if any(k1),
%   tmp = (1-k(k1).*(x(k1)-m(k1))./s(k1));
%   f(k1)=exp(-tmp.^(1./k(k1))).*tmp.^(1./k(k1)-1)./s(k1);
% end
%   
% 
% 
% k3 = find(s<=0 | isnan(s));
% if any(k3),
%    tmp   = NaN;
%    f(k3) = tmp(ones(size(k3)));
% end
% return




