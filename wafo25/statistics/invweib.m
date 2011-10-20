function [x,xlo,xup] = invweib(F,varargin)
%INVWEIB Inverse of the Weibull distribution function
%
% CALL:  x = invweib(F,a,b,c,options)
%        [x,xlo,xup] = invweib(F,phat,options)
%
%          x = inverse cdf for the Weibull distribution evaluated at F
%    xlo,xup = 100*(1-alpha) % confidence bounds of x.
%      a,b,c = parameters
%       phat = Distribution parameter struct
%              as returned from FITWEIB.  
%    options = struct with fieldnames:
%           .lowertail: if TRUE (default), F = Prob[X <= x],
%                       otherwise, F = Prob[X > x].
%           .logp     : if TRUE, probability, p, input as log(p).
%           .alpha    : Confidence coefficent    (default 0.05)
%           .proflog  : if TRUE compute xlo and xup using proflog
%  
%
% The Weibull distribution is defined by its cdf
%
%  F(x;a,b) = 1 -  exp(-((x-c)/a)^b), x>=0, a,b>0
%
% Example:
%   F = linspace(0,1,100);
%   x = invweib(F,10,5);
%   plot(F,x)
%
% See also pdfweib, cdfweib, rndweib, fitweib, momweib

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


% Tested on: Matlab 5.3
% History: 
% revised pab 24.10.2000
% removed comnsize
% -adde pcov and alpha
% - added comnsize, nargchk
% rewritten ms 15.06.2000

error(nargchk(2,inf,nargin))
options = struct('proflog',false,'alpha',0.05,...
  'lowertail',true,'logp',false); % default options
Np = 3;
[params,options,tmp,phat] = parsestatsinput(Np,options,varargin{:});
if numel(options)>1
  error('Multidimensional struct of distribution parameter not allowed!')
end

[a,b,c] = deal(params{:});
if isempty(c), c = 0;end
% Return NaN for out of range parameters.
a(a <= 0) = NaN;
b(b <= 0) = NaN;

if options.logp
  F(F>0) = nan;
  if options.lowertail
    %R = -expm1(F);
    q = -log1p(-exp(F));
  else
    q = -F;  
  end
else
  F(F<0 | 1<F) = NaN;
  if options.lowertail
    q = -log1p(-F);
  else
    q = -log(F);
  end
end

% Return NaN for out of range parameters.
a(a <= 0) = NaN;
b(b <= 0) = NaN;

try
    x = c+a .* q.^(1./b);
catch
   error ('F, a and c must be of common size or scalar!');
end

% 
if nargout>=2
  if isempty(phat)
    error('Must have distribution struct!')
  end
 
  alpha = options.alpha;
  if options.proflog
    xlo = x;
    xup = x;
    for ix =1:numel(x)
      [Lp,CI] = proflog(phat,'i',1,'x',x(ix),'link',@lnkgenpar,'alpha',alpha);
      xlo(ix) = CI(1);
      xup(ix) = CI(2);
    end
  else
    % Compute confidence bounds using the delta method on log scale.
     pcov = phat.covariance;
     %pcov = options.covariance;
   logx = log(x);
   logq = log(q);
   dA = 1./a;
   dB = -1./(b.^2);
   logxvar = pcov(1,1).*dA.^2 + 2*pcov(1,2).*dA.*dB.*logq + pcov(2,2).*(dB.*logq).^2;
%    deriv = [1./A, -1./(B.^2)];
%    pcov = pcov .* (deriv' * deriv);
%    logxvar = pcov(1,1) + 2*pcov(1,2)*logq + pcov(2,2)*logq.^2;
   if any(logxvar<0)
      error('Covariance must be a positive semi-definite matrix.');
   end
   logxcrit = -invnorm(alpha/2)*sqrt(logxvar);
  
   
   % Convert back to Weibull scale
   xlo = exp(logx - logxcrit);
   xup = exp(logx + logxcrit);
  end
end





