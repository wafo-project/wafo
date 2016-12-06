function [F,Flo,Fup] = cdfchi2(x,varargin)
%CDFCHI2 Chi squared cumulative distribution function
%
% CALL:  F = cdfchi2(x,p,options);
%        [F,Flo,Fup] = cdfchi2(x,phat,options);
%
%        F = distribution function evaluated at x
%  Flo,Fup = 100*(1-alpha) % confidence bounds of F.
%       df = degrees of freedom
%     phat = Distribution parameter struct
%            as returned from FITCHI2.  
%  options = struct with fieldnames:
%         .lowertail: if TRUE (default), F = Prob[X <= x],
%                     otherwise, F = Prob[X > x].
%         .logp     : if TRUE, probability, p, returned as log(p).
%         .disable  : if TRUE disable check on integer DF.
%         .alpha    : Confidence coefficent    (default 0.05)
%
% The Chi squared distribution is defined by its pdf
%
%   f(x)=x^(df/2-1)*exp(-x/2)/gamma(df/2)/2^(df/2), x>=0, df=1,2,3,...
%
% Example: 
%   x = linspace(0,15,200);
%   p1 = cdfchi2(x,2); 
%   p2 = cdfchi2(x,3);
%   assert(p1(1:50:end), ...
%   [0.0, 0.848083174802328   0.976921278221860   0.996493953857845], 1e-10)
%   assert(p2(1:50:end), ...
%          [0, 0.712469071706393, 0.943402058537635, 0.989821273700460], 1e-10)
%   plot(x,p1,x,p2);
%
%   close all;
%
% See also pdfchi2, rndchi2, fitchi2, momchi2

% Reference: Johnson, Kotz and Balakrishnan (1994)
% "Continuous Univariate Distributions, vol. 1", p. 415 ff
% Wiley

% Tested on; Matlab 5.3
% History:
% revised pab 25.10.2000
%  - added comnsize, nargchk
% added ms 26.06.2000

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


options = struct('covariance',[],'alpha',0.05,...
  'lowertail',true,'logp',false,'disable',false); % default options
if ischar(x) && strcmpi(x,'defaults')
  F = options;
  return
end
error(nargchk(2,8,nargin))

Np = 1;
[params,options] = parsestatsinput(Np,options,varargin{:});
if numel(options)>1
  error('Multidimensional struct of distribution parameter not allowed!')
end

df = params{1};
ok= ((df >0));
if options.disable
  warntxt =  'df should be a positive.';
else
   ok = (ok & df==round(df));
   warntxt =  'df should be a positive integer'; 
end
if any(~ok) 
  warning('WAFO:PRBCHI2',warntxt)
  df(~ok) = nan;
end
tailnames = {'lower','upper'};
tail = tailnames{2-(options.lowertail~=0)};
try
  x(x<0) = 0;
  F = gammainc(x/2,df/2,tail);
  F(x==0 & isnan(df)) = nan; % make sure output is nan
catch
  error('x and df must be of common size or scalar.');
end

% Compute confidence bounds if requested.
if nargout >= 2
    % Approximate the variance of F on the logit scale
    alpha1 = options.alpha;
    pcov = options.covariance;
    if any(pcov(1) < 0)
      error('PCOV must be a positive semi-definite matrix.');
    end
    logitp = log(F./((1-F)));
    dL   = 1 ./ (F.*(1-F)+realmin); % derivative of logit(F) w.r.t. F
    h = 1e-5;
    dGda = (gammainc(x/2,df/2+h,tail)-F)./h;
   
    dLda = dGda.* dL;      % dlogitp/da = dG/da * dlogitp/dF
   
    varLogitp = pcov(1,1).*dLda.^2;
    
    
    
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


