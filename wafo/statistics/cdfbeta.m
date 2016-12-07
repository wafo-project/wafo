function [F,Flo,Fup] = cdfbeta(x,varargin)
%CDFBETA  Beta cumulative distribution function
%
% CALL:  F = cdfbeta(x,a,b,options);
%        [F,Flo,Fup] = cdfbeta(x,phat,options);
%
%        F = distribution function evaluated at x
%  Flo,Fup = 100*(1-alpha) % confidence bounds of F.
%        x = matrix
%      a,b = distribution parameters for shape
%     phat = Distribution parameter struct
%            as returned from FITBETA.  
%  options = struct with fieldnames:
%         .lowertail: if TRUE (default), F = Prob[X <= x],
%                     otherwise, F = Prob[X > x].
%         .logp     : if TRUE, probability, p, returned as log(p).
%         .alpha    : Confidence coefficent    (default 0.05)
%
% The BETA PDF is defined as:
%
%   f = x^(a-1)*(1-x)^(b-1)/H(a,b)    0<= x <= 1, a>0, b>0
% 
% where H(a,b) is a normalization constant.
%
% NOTE: To compute accurate uppertail probabilities use the identity
%      1 - CDFBETA(X,A,B) = CDFBETA(1-X,B,A).
%   
% Example: 
%   x = linspace(0,1,200);
%   p1 = cdfbeta(x,1,1); 
%   p2 = cdfbeta(x,2,2);
%   plot(x,p1,x,p2);
%
%   close all;
%
% See also invbeta, pdfbeta, rndbeta, fitbeta, mombeta

% Copyright (C) 1993 Anders Holtsberg, 2000 Per A. Brodtkorb

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

% tested on matlab 5.3
%History:
% revised pab 
% -removed dependence on comnsize
%revised pab 29.10.2000
% adapted from stixbox
% -added nargchk, comnsize
%       Anders Holtsberg, 18-11-93


options = struct('proflog',false,'alpha',0.05,...
  'lowertail',true,'logp',false); % default options

if ischar(x) && strcmpi(x,'defaults')
  F = options;
  return
end
error(nargchk(2,9,nargin))
Np = 2;
[params,options,tmp,phat] = parsestatsinput(Np,options,varargin{:});
if numel(options)>1
  error('Multidimensional struct of distribution parameter not allowed!')
end

[a,b] = deal(params{:});
if isempty(b),
  error(nargchk(3,9,nargin))
end

a(a<=0) = 0;
b(b<=0) = 0;
x(x<0)  = 0;
x(1<x)  = 1;
a = min(a,realmax);
b = min(b,realmax);

try
  notOK = (a<=0 | b<=0 | isnan(x));
  x(isnan(x)) = 0;
  
  if options.lowertail
    F = betainc(x,a,b);
%     k = find(a==1 &  ~notOK);
%     if any(k)
%       if isscalar(b), bk = b; else bk = b(k);end
%       if isscalar(x), xk = x; else xk = x(k); end
%         F(k) = -expm1(bk.*log1p(-xk));
%     end
%      k = find(b==1 &  ~notOK);
%     if any(k)
%       if isscalar(b), bk = b; else bk = b(k);end
%       if isscalar(x), xk = x; else xk =x(k); end
%         F(k) = -expm1(bk.*log1p(-xk));
%     end
  else
    x = 1-x;
    F = betainc(x,b,a);
  end
  F(notOK) = nan;
  
catch
   error('x, a and b must be of common size or scalar');
end

% Compute confidence bounds if requested.
if nargout >= 2
  pcov = phat.covariance;
  alpha1 = options.alpha;
  
   h = 1e-5;
  % Approximate the variance of F on the logit scale
  if options.lowertail
    logitp = log(F./(1-F));
    dBda = (betainc(x,a+h,b)-F)./h;
    dBdb = (betainc(x,a,b+h)-F)./h;
  else
    logitp = -log(F./(1-F));
    dBda = (betainc(x,b,a+h)-F)./h;
    dBdb = (betainc(x,b+h,a)-F)./h;
  end
 
  
  dL   = 1 ./ (F.*(1-F)+realmin); % derivative of logit(F) w.r.t. F
  dLda = dBda.* dL;      % dlogitp/da = dF/da * dlogitp/dF
  dLdb = dBdb.* dL;    % dlogitp/db = dF/db * dlogitp/dF
  varLogitp = pcov(1,1).*dLda.^2 + 2.*pcov(1,2).*dLda.*dLdb + pcov(2,2).*dLdb.^2;
    
  if any(varLogitp(:) < 0)
    error('Covariance must be a positive semi-definite matrix.');
  end
    
  % Use a normal approximation on the logit scale, then transform back to
  % the original CDF scale
  zcrit = -invnorm(alpha1/2) * sqrt(varLogitp);
  zlo = logitp - zcrit;
  zup = logitp + zcrit;
  if options.logp
    Flo = -log1p(exp(-zlo));
    Fup = -log1p(exp(-zup));
  else
    Flo = 1 ./ (1 + exp(-zlo));
    Fup = 1 ./ (1 + exp(-zup));
  end
end
 
if options.logp
  F = log(F);
end





