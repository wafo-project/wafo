function [x,xlo,xup] = invfrech(F,varargin)
%INVFRECH Inverse of the Frechet distribution function
%
% CALL:  x = invfrech(F,a,c)
%
%        x = inverse cdf for the Frechet distribution evaluated at F
%     a, c = parameters
%
% The Frechet distribution is defined by its cdf
%
%  F(x;a,c) = exp(-(x/a)^(-c)), x>=0, a,c>0
%
% Example:
%   F = linspace(0,1,100);
%   x = invfrech(F,10,5);
%   plot(F,x);
%
%   close all;
%
% See also pdffrech, cdffrech, rndfrech, fitfrech, momfrech 

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

% Reference: 


% Tested on: Matlab 5.3
% History: 
% Added PJ 10-May-2001

error(nargchk(3,11,nargin))

options = struct('alpha',0.05,...
  'lowertail',true,'logp',false); % default options
Np = 2;
[params,options,tmp,phat] = parsestatsinput(Np,options,varargin{:});
if numel(options)>1
  error('Multidimensional struct of distribution parameter not allowed!')
end

[a,b] = deal(params{:});
if isempty(b), b = 0;end
a(a<=0) = nan;
b(b<=0) = nan;
if options.logp
  F(F>0) = nan;
  if options.lowertail
    q = -F;
  else
    q = -log(expm1(F));
    lrg = 0.5 < F;
    q(lrg) = -log1p(-exp(F(lrg)));
  end
else
  F(F<0 | 1 <F) = nan;
  if options.lowertail
    q = -log(F);
  else
    q = -log1p(-F);
  end
end
try
  x=a.*q.^(-1./b);
catch
  error('F, a and c must be of common size or scalar');
end

if nargout>=2
   %Compute confidence bounds on log scale.
	pcov = phat.covariance
   logx = log(x);
   logq = log(q);
   dA = 1./a;
   dB = -1./(b.^2);
   logxvar = pcov(1,1).*dA.^2 + 2*pcov(1,2).*dA.*dB.*logq + pcov(2,2).*(dB.*logq).^2;
%    deriv = [1./A, -1./(B.^2)];
%    pcov = pcov .* (deriv' * deriv);
%    logxvar = pcov(1,1) + 2*pcov(1,2)*logq + pcov(2,2)*logq.^2;
   if any(logxvar<0)
      error('PCOV must be a positive semi-definite matrix.');
   end
   logxcrit = -invnorm(alpha/2)*sqrt(logxvar);
  
   
   % Convert back to Weibull scale
   xlo = exp(logx - logxcrit);
   xup = exp(logx + logxcrit);
end


% [errorcode, F, a, c] = comnsize(F,a, c);
% if (errorcode > 0)
%   error ('F, a and c must be of common size or scalar');
% end
% 
% x=zeros(size(F));
% 
% ok = ((c > 0)  & (a > 0));
%   
% k = find ((F == 1) & ok);
% if any (k),
%   tmp=inf;
%   x(k) = tmp(ones (size(k)));
% end
%   
% k1 = find ((F > 0) & (F < 1) & ok);
% if any (k1),
%   x(k1)=(-log(F(k1))).^(-1./c(k1)).*a(k1);
% end
% 
% k2 = find(F<0 | F>1 | ~ok);
% if any(k2),
%   tmp=NaN;
%   x(k2)=tmp(ones(size(k2)));
% end
