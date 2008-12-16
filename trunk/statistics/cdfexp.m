function [F,Flo,Fup] = cdfexp(x,varargin)
%CDFEXP Exponential cumulative distribution function
%
% CALL:  F = cdfexp(x,m,options);
%        [F,Flo,Fup] = cdfexp(x,phat,options);
%
%        F = distribution function evaluated at x
%  Flo,Fup = 100*(1-alpha) % confidence bounds of F.
%        x = matrix
%        m = mean (default 1)
%     phat = Distribution parameter struct
%            as returned from FITEXP.  
%  options = struct with fieldnames:
%         .lowertail  : if TRUE (default), F = Prob[X <= x],
%                       otherwise, F = Prob[X > x].
%         .logp       : if TRUE, probability, p, returned as log(p).
%         .alpha      : Confidence coefficent    (default 0.05)
%        
% The Exponential distribution is defined by its cdf
%
%        F(x)=1-exp(-x/m), x>=0, m>0.
% 
% Example: 
%   x = linspace(0,6,200);
%   p1 = cdfexp(x,1); p2 = cdfexp(x,2);
%   plot(x,p1,x,p2), shg
%
% See also pdfexp, invexp, rndexp, fitexp, momexp

% Reference: Johnson, Kotz and Balakrishnan (1994)
% "Continuous Univariate Distributions, vol. 1", p. 494 ff
% Wiley

% Tested on; Matlab 5.3
% History: 
% revised pab nov 2005
% -removed comnsize -> faster
% revised pab Dec2003
% fixed abug: k1 ->k3
% revised pab 24.10.2000
%  - added comnsize, nargchk
% added ms 15.06.2000

error(nargchk(1,10,nargin))
options = struct('proflog',false,'alpha',0.05,...
  'lowertail',true,'logp',false); % default options
Np = 1;
[params,options,tmp,phat] = parsestatsinput(Np,options,varargin{:});
if numel(options)>1
  error('Multidimensional struct of distribution parameter not allowed!')
end
m = params{1};
if isempty(m)
  m = 1;
end

m(m<=0) = NaN;

try 
  xm = x./m;
catch
  error('x and m must be of common size or scalar.');
end

xm(xm<0)  = 0;
if options.lowertail
  if options.logp
    F = log1p(-exp(-xm));
  else
    F = -expm1(-xm);
  end
else
  if options.logp
    F = -xm;
  else
    F = exp(-xm);
  end
end

if nargout>=2
  % Compute confidence bounds on log scale.
   logx = log(xm+realmin);
   mVar = phat.covariance;
   if mVar<0
      error('Variance must be non-negative.');
   end
   xcrit     = -invnorm(options.alpha/2)* sqrt(mVar ./ (m.^2));
  
   xlo = logx - xcrit;
   xup = logx + xcrit;
   
   % Convert back to original scale
   
   if options.lowertail
     if options.logp
       Flo = log1p(-exp(-exp(xlo)));
       Fup = log1p(-exp(-exp(xup)));
     else
       Flo = -expm1(-exp(xlo));
       Fup = -expm1(-exp(xup));
     end
   else
     if options.logp
       Flo = -exp(xlo);
       Fup = -exp(xup);
     else
       Flo = exp(-exp(xlo));
       Fup = exp(-exp(xup));
     end
   end
end

