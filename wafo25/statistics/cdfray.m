function [F, Flo,Fup] = cdfray(x,varargin)
%CDFRAY Rayleigh cumulative distribution function
%
% CALL:  F = cdfray(x,b,options);
%        [F,Flo,Fup] = cdfray(x,phat,options);
%
%        F = distribution function evaluated at x
%  Flo,Fup = 100*(1-alpha) % confidence bounds of F.
%        b = distribution parameter
%     phat = Distribution parameter struct
%            as returned from FITRAY.  
%  options = struct with fieldnames:
%         .lowertail: if TRUE (default), F = Prob[X <= x],
%                     otherwise, F = Prob[X > x].
%         .logp     : if TRUE, probability, p, returned as log(p).
%         .alpha    : Confidence coefficent    (default 0.05)
% 
% The Rayleigh distribution is defined by its cdf
%
%  F(x;b) = 1 - exp(-x^2/(2b^2)), x>=0
%
% Example: 
%   x = linspace(0,4,200);
%   p1 = cdfray(x,1); p2 = cdfray(x,0.5);
%   plot(x,p1,x,p2), shg
%
% See also pdfray, invray, rndray, fitray, momray

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
% and Life Span Models", p. 181 ff, Marcel Dekker.


% Tested on; Matlab 5.3
% History:
% revised pab 2007
% - removed dependence on comnsize
% revised pab 24.10.2000
%  - added comnsize, nargchk
% added ms 15.06.2000

error(nargchk(2,8,nargin))
options = struct('proflog',false,'alpha',0.05,...
  'lowertail',true,'logp',false); % default options
Np = 1;
[params,options,tmp,phat] = parsestatsinput(Np,options,varargin{:});
if numel(options)>1
  error('Multidimensional struct of distribution parameter not allowed!')
end
b = params{1};


x(x<0) = 0;
b(b<0) = nan;
try
  xn = 0.5*(x./b).^2;
catch
    error ('x and b must be of common size or scalar');
end

if nargout>=2
  % Compute confidence bounds on log scale.
   logx = log(xn+realmin);
   bVar = phat.covariance;
   alpha = options.alpha;
   if bVar<0
      error('Variance must be non-negative.');
   end
   xcrit = -invnorm(alpha/2).*2.*sqrt(bVar) ./ b;

   xlo = logx - xcrit;
   xup = logx + xcrit;
   
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



% [csize, x, b] = comnsize (x,b);
% if any(isnan(csize))
%   error ('x and b must be of common size or scalar');
% end
% 
% F=zeros(size(x));
% 
% k = find ((x>=0)&(b>0));
% if any (k)  
%   F(k)=1-exp(-x(k).^2./(2*b(k).^2));
% end
% 
% k1 = find (b<=0);
% if any (k1)
%   tmp=NaN;
%   F(k1) = tmp(ones(size(k1)));
% end




