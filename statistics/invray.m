function [x,xlo,xup] = invray(F,varargin)
%INVRAY Inverse of the Rayleigh distribution function
%
% CALL:  x = invray(F,b,options)
%        [x,xlo,xup] = invray(F,phat,options)
%
%        x = inverse cdf for the Rayleigh distribution evaluated at F
%  xlo,xup = 100*(1-alpha) % confidence bounds of x.
%        b = parameter
%     phat = Distribution parameter struct
%            as returned from FITRAY.  
%  options = struct with fieldnames:
%         .lowertail: if TRUE (default), F = Prob[X <= x],
%                     otherwise, F = Prob[X > x].
%         .logp     : if TRUE, probability, p, input as log(p).
%         .alpha    : Confidence coefficent    (default 0.05)
%         .proflog  : if TRUE compute  xlo and xup using proflog
% 
% The Rayleigh distribution is defined by its cdf
%
%  F(x;b) = 1 - exp(-x^2/(2b^2)), x>=0
%
%
% Example:
%   F = linspace(0,1,100);
%   x = invray(F,1);
%   plot(F,x)
%
% See also pdfray, cdfray, rndray, fitray, momray


% Reference: Cohen & Whittle, (1988) "Parameter Estimation in Reliability
% and Life Span Models", p. 181 ff, Marcel Dekker.

% Tested on: Matlab 5.3
% History: 
% revised pab 24.10.2000
% - added comnsize, nargchk
% added ms 15.06.2000



error(nargchk(2,inf,nargin))
options = struct('proflog',[],'alpha',0.05,...
  'lowertail',true,'logp',false); % default options
Np = 1;
[params,options,tmp,phat] = parsestatsinput(Np,options,varargin{:});
if numel(options)>1
  error('Multidimensional struct of distribution parameter not allowed!')
end
[b] = deal(params{:});

if options.logp
  F(F>0) = nan;
  if options.lowertail
    %R = -expm1(F);
    logR = log1p(-exp(F));
  else
    logR = F;  
  end
else
  F(F<0 | 1<F) = NaN;
  if options.lowertail
    logR = log1p(-F);
  else
    logR = log(F);
  end
end

b(b<0) = nan;
try
    x = sqrt(-2*logR).*b;
catch
  error ('F and b must be of common size or scalar');
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
      [Lp,CI] = proflog(phat,'i',1,'x',x(ix),'link',@lnkray,'alpha',alpha);
      xlo(ix) = CI(1);
      xup(ix) = CI(2);
    end
  else
   % Compute confidence bounds  on log scale.
   bvar = phat.covariance;
   logx = log(x);
  
 
   logxvar = bvar./b.^2;
   if any(logxvar<0)
      error('Covariance must be a positive semi-definite matrix.');
   end
   logxcrit = -invnorm(options.alpha/2).*sqrt(logxvar);
   
   % Convert back to Rayleigh scale
   xlo = exp(logx - logxcrit);
   xup = exp(logx + logxcrit);
  end
end



return
% [icode, F, b] = iscomnsize(F,b);
% if  ~icode 
%   error ('F and b must be of common size or scalar');
% end
% 
% x=zeros(size(F));
% 
%   
% k = find ((F == 1) & (b>0));
% if any (k),
%   tmp=inf;
%   x(k) = tmp(ones (size(k)));
% end
%   
% k1 = find ((F > 0) & (F < 1) & (b>0));
% if any (k1),
%   x(k1)=sqrt(-2*log(1-F(k1))).*b(k1);
% end
% 
% k2 = find(F<0 | F>1 | (b<=0));
% if any(k2),
%   tmp=NaN;
%   x(k2)=tmp(ones(size(k2)));
% end