function [F, Flo,Fup] = cdfraymod(x,varargin)
%CDFRAYMOD Modified Rayleigh cumulative distribution function
%
% CALL:  F = cdfraymod(x,b,c,options);
%        [F,Flo,Fup] = cdfraymod(x,phat,options);
%
%        F = distribution function evaluated at x
%        b = scale parameter
%        c = truncation parameter (default 0)  
%     phat = Distribution parameter struct
%            as returned from FITRAYMOD.  
%  options = struct with fieldnames:
%         .lowertail: if TRUE (default), F = Prob[X <= x],
%                     otherwise, F = Prob[X > x].
%         .logp     : if TRUE, probability, p, returned as log(p).
%         .alpha    : Confidence coefficent    (default 0.05)
%
% The modified Rayleigh distribution is defined by its cdf
%
%  F(x;b,c) = 1 - exp(-(x+|c|)^2/(2b^2)+c^2/(2b^2)), x>=0
%
% Example: 
%   x = linspace(0,4,200);
%   p1 = cdfraymod(x,1); p2 = cdfraymod(x,0.5,-2);
%   plot(x,p1,x,p2),shg
%
% See also pdfraymod, invraymod, rndraymod, fitraymod, momraymod

% Reference: Cohen & Whittle, (1988) "Parameter Estimation in Reliability
% and Life Span Models", p. 181 ff, Marcel Dekker.


% Tested on: Matlab 5.3
% History:
% by pab 03.12.2000
% based on pdfray

error(nargchk(2,9,nargin))
options = struct('covariance',[],'alpha',0.05,...
  'lowertail',true,'logp',false); % default options

Np = 2;
[params,options] = parsestatsinput(Np,options,varargin{:});
if numel(options)>1
  error('Multidimensional struct of distribution parameter not allowed!')
end
[b,c] = deal(params{:});
if isempty(c),c=0;end

x(x<0) = 0;
b(b<0) = nan;
try
  c = abs(c);
  xn = 0.5*((x+c)./b).^2;
  logR = -(xn - 0.5*(c./b).^2);
catch
    error ('x and b must be of common size or scalar');
end

% TODO % Flo and Fup is not correct yet for c~=0
if nargout>=2
  % Compute confidence bounds on log scale.
   logx = log(-logR+realmin);
   bVar = options.covariance;
   alpha = options.alpha;
   if bVar<0
      error('Variance must be non-negative.');
   end
   xcrit = -invnorm(alpha/2).*2.*sqrt(bVar) ./ b;

   xlo = logx - xcrit;
   xup = logx + xcrit;
   
   % Convert back to original scale
   if options.lowertail
     Flo = expm1(-exp(xlo));
     Fup = expm1(-exp(xup));
   else
     Flo = exp(-exp(xlo));
     Fup = exp(-exp(xup));
   end
   if options.logp
     Flo = log(Flo);
     Fup = log(Fup);
   end
end

if options.logp
  if options.lowertail
    F = log(-expm1(logR));
    sml = -xn<-1;
    F(sml) = log1p(-exp(logR(sml)));
  else
    F = logR;
  end
elseif options.lowertail
    F=-expm1(logR);
else
    F= exp(logR);
end


% if nargin<3||isempty(c),c=0;end
% if nargin<4||isempty(a),a=2;end
% [icode, x, b,c] = iscomnsize (x,b,c);
% if ~icode 
%   error ('x, b and c must be of common size or scalar');
% end
% 
% F = zeros(size(x));
% 
% k = find ((x>=0)&(b>0));
% 
% if any(k)  
%   F(k)=(1-exp(-(x(k)-c(k)).^a./(2*b(k).^a)+abs(c(k)).^a./(2*b(k).^a)));
% end
% 
% k1 = find (b<=0);
% if any(k1)
%   tmp=NaN;
%   F(k1) = tmp(ones(size(k1)));
% end