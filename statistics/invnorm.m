function [x,xlo,xup] = invnorm(F,varargin)
%INVNORM Inverse of the Normal distribution function
%
% CALL:  x = invnorm(F,m,v,options)
%        [x,xlo,xup] = invnorm(F,phat,options)
%
%  xlo,xup = 100*(1-ALPHA)% confidence intervall for X.
%        x = inverse cdf for the Normal distribution evaluated at F
%        F = probability
%        m = mean     (default 0)
%        v = variance (default 1)
%     phat = Distribution parameter struct
%            as returned from FITNORM.  
%  options = struct with fieldnames:
%         .lowertail: if TRUE (default), F = Prob[X <= x],
%                     otherwise, F = Prob[X > x].
%         .logp     : if TRUE, probability, p, given as log(p).
%         .alpha    : Confidence coefficent        (default 0.05)
%
% Example:
%   F = linspace(0,1,100);
%   x = invnorm(F,2.5,0.6);
%   plot(F,x)
% 
% See also  cdfnorm, pdfnorm, rndnorm, momnorm, fitnorm

% Copyright (C) 2000 WAFO-group
%
% This program is free software; you can redistribute it and/or modify
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



% Tested on: Matlab 6.0, 5.3
% History:
% revised pab aug 2007
% - replaced phiinv with erfcinv -> faster
% revised pab nov2005
%  - removed call to comnsize -> faster code
%  - added xlo and xup to output
% revised pab 23.03.2003
% -replace erfinv with PHIINV which perform more accurate
%  inversion in the lower tail. This also results in faster execution.
% revised jr 03.04.2001
%  - fixed a bug in the last if statement
%  - updated information, example
% revised pab 24.10.2000
%  - added comnsize, nargchk
%  - added default values
% added ms 17.06.2000

error(nargchk(1,inf,nargin))
options = struct('proflog',false,'alpha',0.05,...
  'lowertail',true,'logp',false); % default options

Np = 2;
[params,options,tmp,phat] = parsestatsinput(Np,options,varargin{:});
if numel(options)>1
  error('Multidimensional struct of distribution parameter not allowed!')
end

[m,v] = deal(params{:});
if isempty(m),
  m=0;
end
if isempty(v),
  v=1;
end

% Out of range parameter or probabilities:
v(v<=0)      = NaN;


sgn = (~options.lowertail)*2-1;
if options.logp
  F(0<F) = NaN;
  xn = sgn*sqrt(2)*erfcinv(2*exp(F));
  near_zero = (log1p(-eps)<F);
  if any(near_zero(:))
    xn(near_zero) = -sgn*sqrt(2)*erfcinv(-2*expm1(F(near_zero)));
  end
else  
  F(F<0 | 1<F) = NaN;
  xn = sgn*sqrt(2)*erfcinv(2*F);
end

try
  x = m + sqrt(v).*xn;
catch
  error ('F, m and v must be of common size or scalar');
end


if nargout>=2
  if isempty(phat)
    error('Must have distribution struct!')
  end
  alpha = options.alpha;
  if options.proflog
    xlo = x;
    xup = x;
    for ix =1:numel(x)
      [Lp,CI] = proflog(phat,'i',1,'x',x(ix),'link',@lnknorm,'alpha',alpha);
      xlo(ix) = CI(1);
      xup(ix) = CI(2);
    end
  else
  pcov = phat.covariance;
 
  % Compute confidence bounds using the delta method
   var_x = pcov(1,1) + 2*pcov(1,2)*xn + pcov(2,2)*xn.^2;
   if any(var_x(:)<0)
      error('Covariance must be a positive semi-definite matrix.');
   end
   %xcrit = -invnorm(alpha/2)*sqrt(var_x);
   xcrit = erfcinv(alpha)*sqrt(2*var_x);
  
   xlo = x - xcrit;
   xup = x + xcrit;
  end
end
return


