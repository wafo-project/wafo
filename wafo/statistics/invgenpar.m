function [x,xlo,xup] = invgenpar(F,varargin)
%INVGENPAR Inverse of the Generalized Pareto distribution function
%
% CALL:  x = invgenpar(F,k,s,m,options)
%        [x,xlo,xup] = invgenpar(F,phat,options)
%
%           x = inverse cdf for the GPD evaluated at F	
%     xlo,xup = 100*(1-ALPHA)% confidence intervall for X.
%           F = lower or upper tail probability
%           k = shape parameter in the GPD
%           s = scale parameter in the GPD    (default 1)
%           m = location parameter in the GPD (Default 0)
%     phat = Distribution parameter struct
%            as returned from FITGENPAR.  
%  options = struct with fieldnames:
%         .lowertail: if TRUE (default), F = Prob[X <= x],
%                     otherwise, F = Prob[X > x].
%         .logp     : if TRUE, input probability, p, given as log(p).
%         .alpha    : Confidence coefficent    (default 0.05)
%         .proflog  : if TRUE compute xlo and xup using proflog
%
% The Generalized Pareto distribution is defined by its cdf
%
%                1 - (1-k(x-m)/s)^1/k,  k~=0
%  F(x;k,s,m) =
%                1 - exp(-(x-m)/s),  k==0
% 
%  for x>m (when k<=0) and m<x<s/k (when k>0), s>0.
%
% Example:
%   F = linspace(0,1,100);
%   x = invgenpar(F,0.3,2);
%   plot(F,x)
%
% See also  cdfgenpar, cdfgenpar, rndgenpar fitgenpar, momgenpar

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


% References 
%  Johnson  N.L., Kotz S. and Balakrishnan, N. (1994)
%  Continuous Univariate Distributions, Volume 1. Wiley. 

% Tested on: Matlab 5.3
% History: 
% Revised by jr 22.12.1999
% revised ms 14.06.2000
% - updated header info
% - changed name to invgenpar (from gpdinv)
% revised pab 25.10.2000
% - added nargchk + comnsize
% revised pab aug 2007
% - added pcov, alpha, uppertail to input
% -removed threshold defining k to zero
error(nargchk(2,inf,nargin))
options = struct('proflog',false,'alpha',0.05,...
  'lowertail',true,'logp',false); % default options

Np = 3;
[params,options,tmp,phat] = parsestatsinput(Np,options,varargin{:});
if numel(options)>1
  error('Multidimensional struct of distribution parameter not allowed!')
end

[k,s,m] = deal(params{:});
if isempty(s), s=1;end
if isempty(m), m=0;end



s(s<=0) = NaN;
%epsilon = 1e-4; % treshold defining k to zero

if options.logp
  if options.lowertail
    %R = -expm1(F);
    logR = log1p(-exp(F));
  else
    logR = F;  
  end
  R = exp(logR);
else
  F(F<0 | 1<F) = NaN;
  if options.lowertail
    logR = log1p(-F);
  else
    logR = log(F);
  end
  R = exp(logR);
end
% if options.lowertail
%   if options.logp
%     R = -expm1(F);
%   else
%     R = 1-F; %convert to survival probabilities
%     
%     R(R<0 | 1<R) = NaN;
%   end
%   
%   logR = log(R);
% else
%   if options.logp
%     logR = F;
%     R = exp(logR);
%   else
%     R = F; % survival probabilities given in
%     R(R<0 | 1<R) = NaN;
%     logR = log(R);
%   end
% end


msgID = 'MATLAB:log:logOfZero';
state = warning('off',msgID);

try
  [csize, k,R,logR] = comnsize(k,R,logR);
 
  q = -logR;
 
 
 
  k1 = find( k~=0 & (0<=R));
  if any(k1)
    %if ~isscalar(k), k = k(k1); end
    %q(k1) = ( 1 - R(k1).^k(k1) )./k(k1);
    q(k1) =  -expm1(k(k1).*logR(k1))./k(k1);
  end
 
  x = m + s.*q;
catch
  error('x, k, s and m must be of common size or scalar.');
end
warning(state)


if nargout >= 2 
  if isempty(phat)
    error('Must have distribution struct!')
  end
 
  alpha = options.alpha;
  if options.proflog  
    xlo = x;
    xup = x;
    for ix =1:numel(x)
      [Lp,CI] = proflog(phat,'i',2,'x',x(ix),'link',@lnkgenpar,'alpha',alpha);
      xlo(ix) = CI(1);
      xup(ix) = CI(2);
    end
  else
  % Compute confidence bounds using the delta method.
    logq = log(q);
 
    dqdk = zeros(size(k));
    
     pcov = phat.covariance;
     %pcov = options.covariance;
    if any(k~=0)
      k1 = find( k~=0 & (0<R));
      if any(k1)
        %dqdk(k1) = (-(1 - R(k1).^k(k1))./k(k1) +R(k1).^k(k1).*logR(k1))./k(k1);
        dqdk(k1) = (-q(k1)+R(k1).^k(k1).*logR(k1))./k(k1);
      end
      k2 =  find( k > 0 & (R==0));
      if any(k2)
        dqdk(k2) = -1./k(k2).^2;
      end
      k3 =  find( k < 0 & (R==0));
      if any(k3)
        dqdk(k3) = NaN;
      end
    end
    
    
    % Approximate the variance of x=q*s on the log scale.
    %    dlogx/dk = dlogx/dq * dqdk = dqdk/q
    %    dlogx/ds = 1/s
    logx = logq + log(s);
    varlogx = pcov(1,1).*(dqdk./q).^2 + 2.*pcov(1,2).*dqdk./(s.*q) + pcov(2,2)./(s.^2);
    if any(varlogx(:) < 0)
        error('Covariance must be a positive semi-definite matrix.');
    end
    zcrit = -invnorm(alpha/2)* sqrt(varlogx);
    
    % Convert back to original scale  
    xlo = m+exp(logx - zcrit);
    xup = m+exp(logx + zcrit);
  end
end


