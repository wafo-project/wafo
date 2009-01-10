function [F,Flo,Fup] = cdfgengam(x,varargin)
%CDFGENGAM Generalized Gamma cumulative distribution function
%
% CALL:  F = cdfgengam(x,a,b,c);
%        [F,Flo,Fup] = cdfgengam(x,phat,options);
%
%        F = distribution function evaluated at x
%  Flo,Fup = 100*(1-alpha) % confidence bounds of F.
%        a = shape parameter
%        b = scale parameter (default b=1)
%        c = shape parameter (default c=1)
%     phat = Distribution parameter struct 
%            as returned from FITGENGAM.  
%  options = struct with fieldnames:
%         .lowertail: if TRUE (default), F = Prob[X <= x],
%                     otherwise, F = Prob[X > x].
%         .logp     : if TRUE, probability, p, returned as log(p).
%         .alpha    : Confidence coefficent    (default 0.05)
%
% The generalized Gamma distribution is defined by its pdf
%
% f(x;a,b,c)=c*x^(a*c-1)/b^(a*c)*exp(-(x/b)^c)/gamma(a), x>=0, a,b,c>0
% 
% Example: 
%   x = linspace(0,7,200);
%   p1 = cdfgengam(x,1,2,1); p2 = cdfgengam(x,3,1,1);
%   plot(x,p1,x,p2), shg
%
% See also cdfgam, pdfgengam, invgengam, rndgengam, fitgengam, momgengam

% Reference: Cohen & Whittle, (1988) "Parameter Estimation in Reliability
% and Life Span Models", p. 220 ff, Marcel Dekker.

% Tested on; Matlab 5.3
% History: 
% revised pab aug 2007
% -removed call to comnsize
% revised pab 23.10.2000
%  - added comnsize, nargchk+ default values on b c.
% added ms 09.08.2000


error(nargchk(2,10,nargin))
options = struct('alpha',0.05,...
  'lowertail',true,'logp',false); % default options
Np = 3;
[params,options,tmp,phat] = parsestatsinput(Np,options,varargin{:});
if numel(options)>1
  error('Multidimensional struct of distribution parameter not allowed!')
end

[a,b,c] = deal(params{:});

if isempty(b),  b=1; end
if isempty(c),  c=1; end

a(a<=0) = nan;
b(b<=0) = nan;
c(c<=0) = nan;
tailnames = {'lower','upper'};
tail = tailnames{2-(options.lowertail~=0)};

try 
  xn = (x./b).^c;
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
   error('x, a, b and c must be of common size or scalar.');
end

if nargout >= 2
  % Compute confidence bounds on the logit scale
    % Variance of F on the logit scale using delta method
    alpha1 = options.alpha;
    pcov   = phat.covariance;
    logitp = log(F./((1-F)));
    dL   = 1 ./ (F.*(1-F)+realmin); % derivative of logit(F) w.r.t. F
    %h = 1e-5;
    %dGda = (gammainc(xn,a+h,tail)-F)./h;
    dGda = dgammainc(xn,a,tail);
    dGdb = -pdfgam(xn,a,1).*xn.*c./b;
    dGdc = pdfgam(xn,a,1).*xn.*log(x./b);
    dLda = dGda.* dL;      % dlogitp/da = dG/da * dlogitp/dF
    dLdb = dGdb.* dL;      % dlogitp/db = dG/db * dlogitp/dF
    dLdc = dGdc.* dL;      % dlogitp/dc = dG/dc * dlogitp/dF
    

    deriv = [dLda(:),dLdb(:), dLdc(:)];
    varLogitp = reshape(sum((deriv*pcov).* deriv,2),size(F));
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

return



