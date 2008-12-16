function [F,Flo,Fup] = cdffrech(x,varargin)
%CDFFRECH Frechet cumulative distribution function
%
% CALL:  F = cdffrech(x,a,b,options);
%        [F,Flo,Fup] = cdffrech(x,phat,options);
%
%        F = density function evaluated at x
%        x = matrix
%     a, b = parameters for scale and shape, respectively.
%     phat = Distribution parameter struct
%            as returned from WFRECHFIT.  
%  options = struct with fieldnames:
%         .lowertail: if TRUE (default), F = Prob[X <= x],
%                     otherwise, F = Prob[X > x].
%         .logp     : if TRUE, probability, p, returned as log(p).
%         .alpha    : Confidence coefficent    (default 0.05)
%
% The Frechet distribution is defined by its cdf
%
%  F(x;a,c) = exp(-(x/a)^(-b)), x>=0, a,b>0
%
% Example: 
%   x = linspace(0,6,200);
%   F1 = cdffrech(x,1,1); F2 = cdffrech(x,2,2); F3 = cdffrech(x,2,5);
%   plot(x,F1,x,F2,x,F3), shg
%
% See also pdffrech, invfrech, rndfrech, fitfrech, momfrech 

% Reference: 

% Tested on; Matlab 5.3
% History: 
% revised pab aug 2007
% -removed call to comnsize -> faster
% Added PJ 10-May-2001


error(nargchk(2,9,nargin))
options = struct('covariance',[],'alpha',0.05,...
  'lowertail',true,'logp',false); % default options
Np = 2;
[params,options] = parsestatsinput(Np,options,varargin{:});
if numel(options)>1
  error('Multidimensional struct of distribution parameter not allowed!')
end

[a,c] = deal(params{:});
if isempty(c),
  error(nargchk(3,9,nargin))
end

c(c<=0) = nan;
a(a<=0) = nan; 
x(x<0)  = 0;
try
  xn = (x./a).^(-c);
catch
  error ('x, a and c must be of common size or scalar');
end
if options.logp
  if options.lowertail
    F = -xn;
  else
    F = log1p(-exp(-xn));
  end
elseif options.lowertail
  F = exp(-xn);
else
  F = -expm1(-xn);
end

if nargout>=2
  % Compute confidence bounds on log scale.
  pcov = options.covariance;
  alpha = options.alpha;
  
  logx = log(xn+realmin);
  da = c./a;
  dc = logx./c;
  logxvar = (pcov(1,1).*da.^2 + 2*pcov(1,2).*da.*dc + pcov(2,2).*dc.^2);
  if any(logxvar(:)<0)
    error('PCOV must be a positive semi-definite matrix.');
  end
  logxcrit     = -invnorm(alpha/2) * sqrt(logxvar);
  
  xlo = exp(logx - logxcrit);
  xup = exp(logx + logxcrit);
   
  % Convert back to original scale  
  if options.logp
    if options.lowertail
      Flo = -xlo;
      Fup = -xup;
    else
      Flo = log1p(-exp(-xlo));
      Fup = log1p(-exp(-xup));
    end
  elseif options.lowertail
    Flo = exp(-xlo);
    Fup = exp(-xup);
  else
    Flo = -expm1(-xlo);
    Fup = -expm1(-xup);
  end
end

return

