function [F,Flo,Fup] = cdfgengammod(x,varargin)
%CDFGENGAMMOD Modified Generalized Gamma cumulative distribution function
%
% CALL:  F = cdfgengammod(x,m,v,L);
%
%        F = distribution function evaluated at x
%        m = location parameter    (default 0)
%        v = scale parameter (v>0) (default 1)
%        L = shape parameter       (default 0)
%
% The modified Generalized Gamma distribution is defined by its cdf
%
%  F(x) = cdfnorm((log(x)-m)/sqrt(v))                    for x>=0, L==0
%         gammainc(exp(L*(log(x)-m)/sqrt(v))/L^2,1./L^2) for x>=0, L>0.
%         1-gammainc(exp(L*(log(x)-m)/sqrt(v))/L^2,1./L^2) for x>=0, L<0.
%
% The parametrization used here is related to the parameters used 
% in PDFGENGAM by the following equations (and only for L>=0):
%     m = log(b)+1/c*log(a)
%     v = 1/(c^2*a)
%     L = 1/sqrt(a)
%
% Example: 
%   x = linspace(0,3,200);
%   p1 = cdfgengammod(x,0,1); p2 = cdfgengammod(x,.5,0.25);
%   plot(x,p1,x,p2)
%
% See also pdfgengammod, invgengammod, rndgengammod, fitgengammod


% Reference: 
% http://www.weibull.com/LifeDataWeb/generalized_gamma_distribution.htm

% Tested on; Matlab 7.3
% History: 
% by pab sept 2007

error(nargchk(1,inf,nargin))
options = struct('covariance',[],'alpha',0.05,...
  'lowertail',true,'logp',false); % default options
if nargin==1 && ischar(x) && strcmpi(x,'defaults')
  F = options;
  return
end
Np = 3;
[params,options] = parsestatsinput(Np,options,varargin{:});
if numel(options)>1
  error('Multidimensional struct of distribution parameter not allowed!')
end

[m,v,L] = deal(params{:});
if isempty(m),  m=0;  end
if isempty(v),  v=1;  end
if isempty(L),  L=0;  end

epsilon = 2^-13; % defining L to zero

v(v<=0) = nan;
x(x<=0) = 0; % trick to set cdf to zero
try
  s = sqrt(v);
  xn = (log(x)-m)./s;
  
  bigL = abs(L)>epsilon;
  if any(~bigL(:))
    F = cdfnorm(xn,options); % lognormal distribution
  else
    F = xn;
  end
 
 
  if any(bigL(:))
    
    tailnames = {'lower','upper'};
    % Generalized gamma distribution
    if ~isscalar(L), 
      L = L(bigL);   
      if ~isscalar(xn)
        xn = xn(bigL);
      end
    else
      %tail = tailnames{2-(lowertail~=0)};
      bigL = bigL(ones(size(xn)));
    end  
    lowertail =  logical(mod(options.lowertail+L>0,2));
    L2 = L.^2;
    if isscalar(L);
      tail = tailnames{2-(lowertail~=0)};
      if options.logp
        F(bigL) = gammaincln(exp(L.*xn)./(L2),1./L2,tail );
        F(bigL(xn==inf)) = log(lowertail~=0);
      else
        F(bigL) = gammainc(exp(L.*xn)./(L2),1./L2,tail );
        F(bigL(xn==inf)) = (lowertail~=0);
      end
    else
      if options.logp
        k1 = find(lowertail);
        if any(k1)
          F(bigL(k1)) = gammaincln(exp(L(k1).*xn(k1))./(L2(k1)),1./L2(k1),'lower');
        end
        k1 = find(~lowertail);
        if any(k1)
          F(bigL(k1)) = gammaincln(exp(L(k1).*xn(k1))./(L2(k1)),1./L2(k1),'upper');
        end      
        F(bigL(xn==inf)) = log(lowertail~=0);
      else
        k1 = find(lowertail);
        if any(k1)
          F(bigL(k1)) = gammainc(exp(L(k1).*xn(k1))./(L2(k1)),1./L2(k1),'lower');
        end
        k1 = find(~lowertail);
        if any(k1)
          F(bigL(k1)) = gammainc(exp(L(k1).*xn(k1))./(L2(k1)),1./L2(k1),'upper');
        end
         F(bigL(xn==inf)) = (lowertail~=0);
      end
    end
% TODO % Call cdfgam here instead    lowertail =
% xor(options.lowertail,L>0);
  end
  
catch
  error ('x, m, v and L must be of common size or scalar');
end

if nargout>1
  Flo = nan;
  Fup = nan;
end
