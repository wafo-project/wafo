function [x,xlo,xup] = invweibmod(F,varargin)
%INVWEIBMOD Inverse of the modified Weibull distribution function
%
% CALL:  x = invweibmod(F,a,b,c,options)
%        [x,xlo,xup] = invweibmod(F,phat,options)
%
%          x = inverse cdf for the Weibull distribution evaluated at F
%    xlo,xup = 100*(1-alpha) % confidence bounds of x.
%      a,b,c = parameters
%       phat = Distribution parameter struct
%              as returned from FITWEIBMOD.  
%    options = struct with fieldnames:
%           .lowertail: if TRUE (default), F = Prob[X <= x],
%                       otherwise, F = Prob[X > x].
%           .logp     : if TRUE, probability, p, input as log(p).
%           .alpha    : Confidence coefficent    (default 0.05)
%  
%
% The modified Weibull distribution is defined by its cdf
%
%  F(x;a,b) = 1 -  exp(-((x+|c|)/a)^b +|c/a|.^b), x>=0, a,b>0
%
% Example:
%   F = linspace(0,1,100);
%   x = invweibmod(F,10,5);
%   plot(F,x)
%
% See also pdfweibmod, cdfweibmod, rndweibmod, fitweibmod, momweibmod

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

error(nargchk(3,5,nargin))
options = struct('covariance',[],'alpha',0.05,...
  'lowertail',true,'logp',false); % default options
Np = 3;
[params,options] = parsestatsinput(Np,options,varargin{:});
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
c = abs(c);
try
    x = a .* (q+(c./a).^b).^(1./b)-c;
catch
   error ('F, a and c must be of common size or scalar!');
end

% 
if nargout>=2
   %Compute confidence bounds on log scale.
   xlo = nan;
   xup = nan;
   return
    pcov = options.covariance;
    alpha = options.alpha;
% TODO % Implement confidence bounds on x   
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





