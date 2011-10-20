function [x,xlo,xup] = invgev(F,varargin)
%INVGEV Inverse of the Generalized Extreme Value distribution function
%
% CALL:  x = invgev(F,k,s,m,options)
%         [x,xlo,xup] = invgev(F,phat,options)
%
%        x = inverse cdf for the GEV distribution evaluated at F
%  xlo,xup = 100*(1-alpha) % confidence bounds of x.
%        k = shape parameter in the GEV 
%        s = scale parameter in the GEV, s>0  (default 1)
%        m = location parameter in the GEV    (default 0)
%     phat = Distribution parameter struct
%            as returned from FITGEV.  
%  options = struct with fieldnames:
%         .lowertail: if TRUE (default), F = Prob[X <= x],
%                     otherwise, F = Prob[X > x].
%         .logp     : if TRUE, probability, p, input as log(p).
%         .alpha    : Confidence coefficent    (default 0.05)
%         .proflog  : if TRUE compute  xlo and xup using proflog
%
%
% The Generalized Extreme Value distribution is defined by its cdf
%
%                exp( - (1 - k(x-m)/s)^1/k) ),  k~=0
%  F(x;k,s,m) =
%                exp( -exp(-(x-m)/s) ),  k==0
%
%  for k(x-m)/s<=1 (when k<=0) and x<m+s/k (when k>0).
%
% Example:
%  k=1;s=1;m= 0;    
%  opt = {'lowertail',false,'logp',false}
%  F0 = [logspace(log10(realmin),-1) linspace(0.2,1-1e-3) logspace(log10(1-sqrt(eps)),log1p(-eps)/log10(10))];
%  %F0 = [logspace(-300,-1) linspace(0.11,0.5)];
%  x  = invgev(F0,k,s,m,opt{:});
%  F  = cdfgev(x,k,s,m,opt{:});
%  semilogy(abs(F-F0)./F0+eps), shg % relative error
%
% See also pdfgev, cdfgev, rndgev, fitgev, momgev

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

% Tested on: Matlab 5.3
% History: 
% Revised by jr 22.12.1999
% revised ms 14.06.2000
% - updated header info
% - changed name to invgev (from gevinv)
% revised pab 24.10.2000
% - added  nargchk, comnsize and default values for m, s

error(nargchk(2,9,nargin))
options = struct('proflog',false,'alpha',0.05,...
  'lowertail',true,'logp',false); % default options

Np = 3;
[params,options,tmp,phat] = parsestatsinput(Np,options,varargin{:});
if numel(options)>1
  error('Multidimensional struct of distribution parameter not allowed!')
end

[k,s,m] = deal(params{:});
if isempty(s),  s=1;end
if isempty(m),  m=0;end

s(s<=0) = NaN;
if options.logp
  F(F>0) = NaN;
  if options.lowertail
    logF = F;
  else
    logF = log(-expm1(F));
    lrg = 0.5 < F;
    logF(lrg) = log1p(-exp(F(lrg)));
  end
else
  F(F<0 | 1<F) = NaN;
  if options.lowertail
    logF = log(F);
  else
    logF = log1p(-F);
  end
end


try 
  q = -log(-logF);
  k1 = find(k~=0 & logF<=0);
  if any(k1)
    if ~isscalar(k), k = k(k1); end
   % q(k1) = (1-(-log(F(k1))).^k(k1))./k(k1); 
    q(k1) = -expm1(-k.*q(k1))./k;
  end
   x = m + s.*q;
catch
  error('k s and m must be of common size or scalar.');
end
if nargout>1
  if isempty(phat)
    error('Must have distribution struct!')
  end
  
  alpha = options.alpha;
  if options.proflog
    xlo = x;
    xup = x;
    for ix =1:numel(x)
      [Lp,CI] = proflog(phat,'i',2,'x',x(ix),'link',@lnkgev,'alpha',alpha);
      xlo(ix) = CI(1);
      xup(ix) = CI(2);
    end
  else
% TODO % Implement  xlo and xup
    %pcov = phat.covariance;
    warning('WAFO:INVGEV','xlo and xup not implemented yet')
    xlo = nan;
    xup = xlo;
  end
end

return
% [icode F k,s,m] = iscomnsize(F,k,s,m);
% if ~icode
%     error('k s and m must be of common size or scalar.');
% end
% 
% % Initialize  x to zero.
% x = zeros(size(k));
% 
% epsilon=1e-4;
% 
% ok = (F>=0 & F<=1 & s>0);
% 
% k1 = find(abs(k)< epsilon & ok);
% if any(k1)
%   x(k1) = m(k1) - s(k1).*log(-log(F(k1)));
% end
% 
% k2 = find(abs(k)>= epsilon & ok);
% if any(k2)
%   x(k2) = m(k2) + s(k2).*(1-(-log(F(k2))).^k(k2))./k(k2);
% end
% 
% 
% tmp=Inf;
% k3 = find(abs(k) < epsilon & F==1&ok);
% if any(k3)
%   x(k3)=tmp(ones(size(k3)));
% end
% k4 = find(abs(k) >= epsilon & F==1&ok);
% if any(k4)
%   x(k4)=m(k4)+s(k4)./k(k4);
% end
% 
% k5=find(F==0&ok);
% if any(k5),
%   x(k5)=-tmp(ones(size(k5)));
% end
% 
% k6=find(~ok);
% if any(k6),
%   tmp=NaN;
%   x(k6)=tmp(ones(size(k6)));
% end