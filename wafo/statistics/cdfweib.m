function [F,Flo,Fup] = cdfweib(x,varargin)
%CDFWEIB Weibull cumulative distribution function
%
% CALL:  F = cdfweib(x,a,b,c,options);
%       [F, Flo,Fup] = cdfweib(x,phat,options);
%
%        F = distribution function evaluated at x
%  Flo,Fup = 100*(1-alpha) % confidence bounds of F.
%     a,b,c = parameters
%     phat = Distribution parameter struct
%            as returned from FITWEIB.  
%  options = struct with fieldnames:
%         .lowertail: if TRUE (default), F = Prob[X <= x],
%                     otherwise, F = Prob[X > x].
%         .logp     : if TRUE, probability, p, returned as log(p).
%         .alpha    : Confidence coefficent    (default 0.05)
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
%   p1 = cdfweib(x,1,1); p2 = cdfweib(x,2,2);
%   plot(x,p1,x,p2);
%
%   close all;
%
% See also pdfweib, invweib, rndweib, fitweib, momweib

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
% revised pab Aug 2007
% -added pcov, alpha
% revised pab Feb2005
%  -removed comnsize
%  -added pcov and alpha to input
% revised pab 24.10.2000
%  - added comnsize
% rewritten ms 15.06.2000


error(nargchk(2,inf,nargin))

options = struct('alpha',0.05,...
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

try
  % Force a zero for negative x.
  x(x < c) = c;
  xn = ((x-c)./a).^b;
catch
   error ('x, a, b and c must be of common size or scalar');
end



% TODO % Flo and Fup not correct for c~=0
if nargout>=2
  % Compute confidence bounds on log scale
  pcov = phat.covariance;
  alpha = options.alpha;
  
   logxn = log(xn);
   dA = -b./a;
   dB = logxn./b;
   logxvar = (pcov(1,1).*dA.^2 + 2*pcov(1,2).*dA.*dB + pcov(2,2).*dB.^2);
   
   % Assume comnsize Then:
   %delta = [dA(:), dB(:)];
   %logxvar = reshape(delta*pcov*delta.',size(xn));

   if any(logxvar(:)<0)
      error('Covariance must be a positive semi-definite matrix.');
   end
   xcrit = -invnorm(alpha/2)* sqrt(logxvar);

   xlo = logxn - xcrit;
   xup = logxn + xcrit;

   % Convert back to Weibull scale
   % Convert back to original scale
   if options.lowertail
     Flo = -expm1(-exp(xlo));
     Fup = -expm1(-exp(xup));
   else
     Flo = exp(-exp(xlo));
     Fup = exp(-exp(xup));
   end
   if options.logp
     Flo = log(Flo);
     Fup = log(Fup);
   end
end

if options.logp
  if options.lowertail
    F = log(-expm1(-xn));
    sml = -xn<-1;
    F(sml) = log1p(-exp(-xn(sml)));
  else
    F = -xn;
  end
elseif options.lowertail
    F=-expm1(-xn);
else
    F= exp(-xn);
end






