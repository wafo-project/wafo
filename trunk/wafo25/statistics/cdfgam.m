function [F, Flo,Fup] = cdfgam(x,varargin)
%CDFGAM Gamma cumulative distribution function
%
% CALL:  [F, Flo,Fup] = cdfgam(x,a,b,options);
%
%        F = distribution function evaluated at x
%  Flo,Fup = 100*(1-alpha) % confidence bounds of F.
%        a = shape parameter
%        b = scale parameter (default b=1)
%     phat = Distribution parameter struct 
%            as returned from FITGAM.  
%  options = struct with fieldnames:
%         .lowertail: if TRUE (default), F = Prob[X <= x],
%                     otherwise, F = Prob[X > x].
%         .logp     : if TRUE, probability, p, returned as log(p).
%         .alpha    : Confidence coefficent    (default 0.05)
%
% The Gamma distribution is defined by its pdf
%
%        f(x)=x^(a-1)*exp(-x/b)/gamma(a)/b^a, a,b>0, x>=0.
%
% Example: 
%   x = linspace(0,7,200);
%   p1 = cdfgam(x,1); p2 = cdfgam(x,2);
%   plot(x,p1,x,p2), shg
%
% See also  cdfgengam, pdfgam, invgam, rndgam, fitgam, momgam

% Reference: Johnson, Kotz and Balakrishnan (1994)
% "Continuous Univariate Distributions, vol. 1", p. 494 ff
% Wiley

% Tested on; Matlab 5.3
% History:
% revised pab nov2005
% removed call to comnsize -> faster call
% revised pab 24.10.2000
%  - added comnsize, nargchk
% added ms 26.06.2000
% added b parameter ms 23.08.2000

options = struct('covariance',[],'alpha',0.05,...
  'lowertail',true,'logp',false); % default options
if ischar(x) && strcmpi(x,'defaults')
  F = options;
  return
end
error(nargchk(2,inf,nargin))

Np = 2;
[params,options] = parsestatsinput(Np,options,varargin{:});
if numel(options)>1
  error('Multidimensional struct of distribution parameter not allowed!')
end

[a,b] = deal(params{:});

if isempty(b), 
  b=1;
end

a(a<=0) = NaN;
b(b<=0) = NaN;
tailnames = {'lower','upper'};
tail = tailnames{2-(options.lowertail~=0)};
try
  xn        = x./b;
  xn(xn<0)  = 0;
  if options.logp
    logF = gammaincln(xn,a,tail);
    logF(xn==inf) = log(double(options.lowertail~=0));
    logF(xn==0 & isnan(a)) = nan; % make sure output is nan
    F = exp(logF);
  else
    F = gammainc(xn,a,tail);
    F(xn==inf) = (options.lowertail~=0);
    F(xn==0 & isnan(a)) = nan; % make sure output is nan
  end
catch
  error('x, a and b must be of common size or scalar.');
end



if nargout >= 2
  % Compute confidence bounds on the logit scale
    % Variance of F on the logit scale using delta method
    alpha1 = options.alpha;
    pcov = options.covariance;
    logitp = log(F./((1-F)));
    dL   = 1 ./ (F.*(1-F)+realmin); % derivative of logit(F) w.r.t. F
    %h = 1e-5;
    %dGda = (gammainc(xn,a+h,tail)-F)./h;
    dGda = dgammainc(xn,a,tail);
    dGdb = -pdfgam(xn,a,1).*xn./b;
    dLda = dGda.* dL;      % dlogitp/da = dG/da * dlogitp/dF
    dLdb = dGdb.* dL;      % dlogitp/db = dG/db * dlogitp/dF
    varLogitp = pcov(1,1).*dLda.^2 + 2.*pcov(1,2).*dLda.*dLdb + pcov(2,2).*dLdb.^2;
    
    if any(varLogitp(:) < 0)
        error('PCOV must be a positive semi-definite matrix.');
    end
   
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
  F = logF;
end



