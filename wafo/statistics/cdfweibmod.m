function [F,Flo,Fup] = cdfweibmod(x,varargin)
%CDFWEIBMOD Truncated Weibull cumulative distribution function
%
% CALL:  F = cdfweibmod(x,a,b,c,options);
%       [F, Flo,Fup] = cdfweibmod(x,phat,options);
%
%        F = distribution function evaluated at x
%    a,b,c = parameters
%     phat = Distribution parameter struct
%            as returned from FITWEIBMOD.  
%  options = struct with fieldnames:
%         .lowertail: if TRUE (default), F = Prob[X <= x],
%                     otherwise, F = Prob[X > x].
%         .logp     : if TRUE, probability, p, returned as log(p).
%         .alpha    : Confidence coefficent    (default 0.05)
%
% The Truncated Weibull distribution is defined by its cdf
%
%  F(x;a,b,c) = 1 -  exp(-((x+|c|)/a)^b+abs(c/a)^b), x>=0, a>0,b>0
%
%   Some references refer to the Weibull distribution with
%   a single parameter, this corresponds to PDFWEIB with a = 1.
%
% Example: 
%   x = linspace(0,6,200);
%   p1 = cdfweibmod(x,1,1,1); p2 = cdfweibmod(x,2,2,3);
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
% revised pab 24.10.2000
%  - added comnsize
% rewritten ms 15.06.2000


error(nargchk(2,inf,nargin))
options = struct('covariance',[],'alpha',0.05,...
  'lowertail',true,'logp',false); % default options
Np = 3;
[params,options] = parsestatsinput(Np,options,varargin{:});
if numel(options)>1
  error('Multidimensional struct of distribution parameter not allowed!')
end

[a,b,c] = deal(params{:});
if isempty(c), c = 0;end



a(a<=0) = nan;
b(b<=0) = nan;
try
  c = abs(c);
  x(x<0) = 0; % trick to set cdf to zero
  xn = ((x+c)./a).^b;
  logR =  -xn+(c./a).^b;
catch
  error ('x, a,b, and c must be of common size or scalar');
end

% TODO % Flo and Fup not correct for c~=0
if nargout>=2
  % Compute confidence bounds on log scale
  pcov = options.covariance;
  alpha = options.alpha;
  
   logxn = log(xn);
   if any(c~=0)
     if isscalar(c)
       logxn(:) = nan;
     else
       logxn(c~=0) = nan;
     end
   end
   dA = -b./a;
   dB = logxn./b;
   logxvar = (pcov(1,1).*dA.^2 + 2*pcov(1,2).*dA.*dB + pcov(2,2).*dB.^2);
   
   % Assume comnsize Then:
   %delta = [dA(:), dB(:)];
   %logxvar = reshape(delta*pcov*delta.',size(xn));

   if any(logxvar(:)<0)
      error('PCOV must be a positive semi-definite matrix.');
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
    F = log(-expm1(logR));
    sml = -xn<-1;
    F(sml) = log1p(-exp(logR(sml)));
  else
    F = logR;
  end
elseif options.lowertail
    F=-expm1(logR);
else
    F= exp(logR);
end

% [icode, x, a,b,c] = iscomnsize (x,a,b, abs(c));
% if ~icode
%   error ('x, a, b and c must be of common size or scalar');
% end
% 
% F=zeros(size(x));
% 
% ok = ((b > 0)  & (a > 0));
% 
% k = find (x>=0&ok);
% if any (k) 
%   F(k)=1-exp(-((x(k)+c(k))./a(k)).^b(k)+abs(c(k)./a(k)).^b(k));
% end
% 
% k1 = find (~ok);
% if any (k1)
%   tmp=NaN;
%   F(k1) = tmp(ones(size(k1)));
% end


