function [x,xlo,xup] = invraymod(F,varargin)
%INVRAYMOD Inverse of the modified Rayleigh distribution function
%
% CALL:  x = invraymod(F,b,c,options)
%        [x,xlo,xup] = invraymod(F,phat,options)
%
%        x = inverse cdf for the Rayleigh distribution evaluated at F
%  xlo,xup = 100*(1-alpha) % confidence bounds of x.
%      b,c = parameters
%     phat = Distribution parameter struct
%            as returned from FITRAYMOD.  
%  options = struct with fieldnames:
%         .lowertail: if TRUE (default), F = Prob[X <= x],
%                     otherwise, F = Prob[X > x].
%         .logp     : if TRUE, probability, p, input as log(p).
%         .alpha    : Confidence coefficent    (default 0.05)
% 
% The modified Rayleigh distribution is defined by its cdf
%
%  F(x;b,c) = 1 - exp(-(x+c)^2/(2b^2)+0.5(c/b)^2), x>=0
%
%
% Example:
%   F = linspace(0,1,100);
%   x = invraymod(F,1);
%   plot(F,x);
%  
%   close all;
%
% See also pdfraymod, cdfraymod, rndraymod, fitraymod, momraymod

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

% Tested on: Matlab 5.3
% History: 
% revised pab 24.10.2000
% - added comnsize, nargchk
% added ms 15.06.2000


error(nargchk(2,inf,nargin))
options = struct('covariance',[],'alpha',0.05,...
  'lowertail',true,'logp',false); % default options
Np = 2;
[params,options] = parsestatsinput(Np,options,varargin{:});
if numel(options)>1
  error('Multidimensional struct of distribution parameter not allowed!')
end
[b,c] = deal(params{:});
if isempty(c)
  c = 0;
end
if options.logp
  F(F>0) = nan;
  if options.lowertail
    %R = -expm1(F);
    logR = log1p(-exp(F));
  else
    logR = F;  
  end
else
  F(F<0 | 1<F) = NaN;
  if options.lowertail
    logR = log1p(-F);
  else
    logR = log(F);
  end
end

b(b<0) = nan;
try
  c = abs(c);
  x = sqrt(-2*logR+(c./b).^2).*b-c;
catch
  error ('F and b must be of common size or scalar');
end

% TODO % confidence bounds not implemented for c~=0;
if nargout>=2
   % Compute confidence bounds  on log scale.
   bvar = options.covariance;
   logx = log(x);
    if any(c~=0)
     if isscalar(c)
       logx(:) = nan;
     else
       logx(c~=0) = nan;
     end
   end
   logxvar = bvar./b.^2;
   if any(logxvar<0)
      error('PCOV must be a positive semi-definite matrix.');
   end
   logxcrit = -invnorm(options.alpha/2).*sqrt(logxvar);
   
   % Convert back to Rayleigh scale
   xlo = exp(logx - logxcrit);
   xup = exp(logx + logxcrit);
end



return
% [icode, F, b] = iscomnsize(F,b);
% if  ~icode 
%   error ('F and b must be of common size or scalar');
% end
% 
% x=zeros(size(F));
% 
%   
% k = find ((F == 1) & (b>0));
% if any (k),
%   tmp=inf;
%   x(k) = tmp(ones (size(k)));
% end
%   
% k1 = find ((F > 0) & (F < 1) & (b>0));
% if any (k1),
%   x(k1)=sqrt(-2*log(1-F(k1))).*b(k1);
% end
% 
% k2 = find(F<0 | F>1 | (b<=0));
% if any(k2),
%   tmp=NaN;
%   x(k2)=tmp(ones(size(k2)));
% end