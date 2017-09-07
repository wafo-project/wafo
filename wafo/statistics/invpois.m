% INVPOIS Inverse of the Poisson distribution function
%  
%   CALL:  x = invpois(F,L,options)
%         [x,xlo,xup] = invpois(F,phat,options)
%  
%          x = inverse cdf for the Poisson distribution evaluated at F
%    xlo,xup = 100*(1-alpha) % confidence bounds of x.
%          L = parameter, L>0 
%       phat = Distribution parameter struct
%              as returned from FITPOIS.  
%    options = struct with fieldnames:
%           .lowertail: if TRUE (default), F = Prob[X <= x],
%                       otherwise, F = Prob[X > x].
%           .logp     : if TRUE, probability, p, input as log(p).
%           .alpha    : Confidence coefficent        (default 0.05)  
%  
%   Example:
%     L = 10; N = 40;  
%     opt = {'lowertail',true,'logp',false};
%     F0 = cdfpois(0:N,L,opt{:});
%     x  = invpois(F0,L,opt{:});
%     F  = cdfpois(x,L,opt{:});
%     semilogy(abs(F-F0)./F0+eps); % relative error
%   
%     opt = {'lowertail',true,'logp',true};
%     x0 = 1:N;
%     F  = cdfpois(x0,L,opt{:});
%     x  = invpois(F,L,opt{:});
%     N0 = length(x0);
%     semilogy(1:N0,abs(x-x0)./x0+eps); % relative error
%
%     assert(x, x0)
%     close all;
%  
%   See also  pdfpois, cdfpois, fitpois, rndpois, mompois


% Copyright (C) 1995, 1996, 1997, 2005, 2006, 2007 Kurt Hornik, Per A. Brodtorb
%
% INVPOIS is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or (at
% your option) any later version.
%
% INVPOIS is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with INVPOIS; see the file COPYING.  If not, see
% <http://www.gnu.org/licenses/>.

% Author: KH <Kurt.Hornik@wu-wien.ac.at>
% Revised pab
% Description: Quantile function of the Poisson distribution
function x = invpois (F, varargin)

%error(nargchk(2,inf,nargin))
narginchk(2,inf)
options = struct('proflog',false,'alpha',0.05,...
  'lowertail',true,'logp',false,'maxiter',10000); % default options

Np = 1;
[params,options] = parsestatsinput(Np,options,varargin{:});
if numel(options)>1
  error('Multidimensional struct of distribution parameter not allowed!')
end
L = params{1};

% For large values of L (>1000),
% the normal distribution with mean = L and variance = L is
% an excellent approximation to the Poisson distribution.
x = max(round(invnorm(F,L,L,options)),0);

[csize, F] = comnsize (F, L);
if any(isnan(csize))
  error('WAFO:INVPOIS','F and L must be of common size or scalar');
end

if options.logp
  F(F>0) = NaN;
  if options.lowertail
    F  = exp(F);
  else
    F = -expm1(F);
  end
else
  F(F<0 | 1<F) = NaN;
  if ~options.lowertail
    F = 1-F;
  end
end

L(L<=0) = nan;
ok = ~((F < 0) | (F > 1) | isnan (F) | ~(L > 0));
k = find (~ok);
if (any (k))
  x(k) = NaN;
end %if

  k = find ((F == 1) & (L > 0));
  if (any (k))
      x(k) = Inf;
  end %if

  k = find (F<1 & ok & (0<L & L < 10000));
  if (any (k))
    x(k) = -2;
    cdf = zeros(size(k));
    iter = 0;
    while (iter<=options.maxiter)
      m = find(cdf <= F(k)+eps);
      if (any(m))
        x(k(m)) = x(k(m)) + 1;
        cdf(m) = cdfpois(x(k(m))+1, L);
      else
        break;
      end %if
    end %while
    if iter==options.maxiter
      warning('WAFO:INVPOIS','Some quantiles did not converge')
    end
  end%if
end %function
