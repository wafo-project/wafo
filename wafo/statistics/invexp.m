function [x,xlo,xup] = invexp(F,varargin)
%INVEXP Inverse of the Exponential distribution function
%
% CALL:  x = invexp(F,m,options)
%       [x,xlo,xup] = invexp(F,phat,options)
%
%        x = inverse cdf for the Exponential distribution evaluated at F
%  xlo,xup = 100*(1-alpha) % confidence bounds of x.
%        m = mean, m>0
%     phat = Distribution parameter struct
%            as returned from FITEXP.  
%  options = struct with fieldnames:
%         .lowertail: if TRUE (default), F = Prob[X <= x],
%                     otherwise, F = Prob[X > x].
%         .logp     : if TRUE, probability, p, input as log(p).
%         .alpha    : Confidence coefficent        (default 0.05)
%         .proflog  : if TRUE compute  xlo and xup using proflog
%
% Example:
%   m=1;    
%   opt = {'lowertail',false,'logp',false};
%   F0 = [logspace(log10(realmin),-1) linspace(0.2,1-1e-3) ...
%         logspace(log10(1-sqrt(eps)),log1p(-eps)/log10(10))];
%   %F0 = [logspace(-300,-1) linspace(0.11,0.5)];
%   x  = invexp(F0,m,opt{:});
%   F  = cdfexp(x,m,opt{:});
%   semilogy(abs(F-F0)./F0+eps); % relative error
%
%   close all;
%
% See also  cdfexp, pdfexp, rndexp, fitexp, momexp

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


% Reference: Johnson, Kotz and Balakrishnan (1994)
% "Continuous Univariate Distributions, vol. 1", p. 494 ff
% Wiley


% Tested on; Matlab 5.3
% History: 
% revised pab nov2005
% -added confidence interval for x
% revised pab 24.10.2000
%  - added comnsize, nargchk
% added ms 10.08.2000



error(nargchk(1,9,nargin))
options = struct('proflog',false,'alpha',0.05,...
  'lowertail',true,'logp',false); % default options

Np = 1;
[params,options,tmp,phat] = parsestatsinput(Np,options,varargin{:});
if numel(options)>1
  error('Multidimensional struct of distribution parameter not allowed!')
end

m = params{1};
if isempty(m), m=1;end

m(m<=0)      = NaN;
if options.logp
  F(F>0) = NaN;
  if options.lowertail
    q = -log(expm1(F));
    lrg = 0.5 < F;
    q(lrg) = -log1p(-exp(F(lrg)));
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


try 
  x  = m.*q;
catch
  error ('F and m must be of common size or scalar');
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
      CI = ciproflog(phat,'i',1,'x',x(ix),'link',@(x,logR,phat,ix) -x./logR,'alpha',alpha);
      xlo(ix) = CI(1);
      xup(ix) = CI(2);
    end
  else
    if 1,
      alpha2  = options.alpha/2;
      n = length(phat.data);
      gamcrit = invgam([alpha2, 1-alpha2],n,1);
      xup = n*x/gamcrit(1);
      xlo = n*x/gamcrit(2);
    else
      % Compute confidence bounds on log scale.
      mVar = phat.covariance;
      logx = log(x+realmin);
      if mVar<0
        error('Variance must be non-negative.');
      end
      lgxcrit = -invnorm(options.alpha/2).* sqrt(mVar ./ (m.^2));
      xlo = exp(logx - lgxcrit);
      xup = exp(logx + lgxcrit);
    end
  end
end


