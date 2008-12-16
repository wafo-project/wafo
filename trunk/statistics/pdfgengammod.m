function f = pdfgengammod(x,varargin)
%PDFGENGAMMOD Modified Generalized Gamma probability density function (stable)
%
% CALL:  f = pdfgengammod(x,m,v,L);
%
%        f = density function evaluated at x
%        m = location parameter    (default 0)
%        v = scale parameter (v>0) (default 1)
%        L = shape parameter       (default 0)
%
% The modified Generalized Gamma distribution is defined by its cdf
%
%  F(x) = cdfnorm((log(x)-m)/sqrt(v))                    for x>=0, L==0
%         gammainc(exp(L*(log(x)-m)/sqrt(v))/L^2,1./L^2) for x>=0, L~=0.
%
% The parametrization used here is related to the parameters used 
% in PDFGENGAM by the following equations:
%     m = log(b)+1/c*log(a)
%     v = 1/(c^2*a)
%     L = 1/sqrt(a)
%
% Example: 
%   x = linspace(0,3,200);
%   p1 = pdfgengammod(x,0,1); p2 = pdfgengammod(x,.5,0.25);
%   plot(x,p1,x,p2)
%
% See also cdfgengammod, invgengammod, rndgengammod, fitgengammod


% Reference: 
% Lawless, J.F., Statistical Models And Methods for Lifetime Data,
% John Wiley & Sons, Inc., New York, 1982.

% Tested on; Matlab 7.3
% History: 
% by pab sept 2007

options = struct('logp',false); % default options
if ischar(x) && strcmpi(x,'defaults')
  f = options;
  return
end
error(nargchk(1,inf,nargin))
Np = 3;

[params,options] = parsestatsinput(Np,options,varargin{:});
if numel(options)>1
  error('Multidimensional struct of distribution parameter not allowed!')
end

[m,v,L] = deal(params{:});
if nargin<2||isempty(m),  m=0;  end
if nargin<3||isempty(v),  v=1;  end
if nargin<4||isempty(L),  L=0;  end

epsilon = 2^-13; % defining L to zero
v(v<=0) = nan;
x(x<0)  = inf; % trick to set pdf to zero
try
  s = sqrt(v);
  xn = (log(x)-m)./s;
  xs = x.*s;
  %f = exp(-0.5*xn.^2)./(xs.*sqrt(2*pi));
  logf = -0.5*xn.^2-log(xs)-0.5*log(2*pi);
  
  bigL = abs(L)>epsilon;
  if any(bigL(:))
    % Generalized gamma distribution
    if ~isscalar(L), 
      L = L(bigL);   
      if ~isscalar(xn)
        xn = xn(bigL);
        xs = xs(bigL);
      end
    else
      bigL = bigL(ones(size(xn)));
    end  
    %f(bigL) =abs(L).*exp( (L.*xn-log(L.^2) -exp(L.*xn))./L.^2-gammaln(1./L.^2))./xs;
    logf(bigL) =log(abs(L)) + ((min(L.*xn,realmax)-log(L.^2) -exp(L.*xn))./L.^2-gammaln(1./L.^2)) - log(xs);
  end
  % Handle special cases for x == 0 (or xn==-inf):
  spcase = (xn==-inf ) ; %& (abs(L)<=epsilon));
  if any(spcase)
    %f(spcase) = 0;
    logf(spcase) = -inf;
    if any(L.*s==1)
      logb = m+log(L.^2).*s./L;
      limit1 = abs(L)./s.*exp(-gammaln(1./L.^2)-logb);
      if ~isscalar(limit1)
        limit1 = limit1(spcase & L.*s==1);
      end
      %f(spcase & L.*s==1) = limit1;
      logf(spcase & L.*s==1) = log(limit1);
    end
  end
catch
  error ('x, m and v must be of common size or scalar');
end
if options.logp
  f = logf;
else
  f = exp(logf);
end

% f=zeros(size(x));
% 
% k = find (x>0&v>0);
% if any(k)    
%   f(k)=1./sqrt(2*pi*v(k)).*exp(-0.5*(log(x(k))-m(k)).^2./v(k))./x(k);
% end
% 
% k1 = find (v<=0);
% if any (k1)
%   tmp=NaN;
%   f(k1) = tmp(ones(size(k1)));
% end




