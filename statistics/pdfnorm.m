function f = pdfnorm(x,varargin)
%PDFNORM Normal probability density function
%
% CALL:  f = pdfnorm(x,m,v,options);
%        f = pdfnorm(x,phat,options);
%
%        f = density function evaluated at x
%        m = mean     (default 0)
%        v = variance (default 1)
%     phat = Distribution parameter struct
%            as returned from FITNORM.  
%  options = struct with fieldnames:
%         .logp   : if TRUE, density, p, returned as log(p).
% 
% Example: 
%   x = linspace(-3,3,200);
%   p1 = pdfnorm(x,0,1); p2 = pdfnorm(x,.5,0.25);
%   plot(x,p1,x,p2),shg
%
% See also cdfnorm, invnorm, rndnorm, fitnorm, momnorm


% Tested on; Matlab 5.3
% History:
% revised pab 9Aug2003
%  fixed a bug: 
% revised pab 24.10.2000
%  - added comnsize, nargchk
%  - added default values
% added ms 15.06.2000

error(nargchk(1,5,nargin))
options = struct('logp',false); % default options
Np = 2;
[params,options] = parsestatsinput(Np,options,varargin{:});
if numel(options)>1
  error('Multidimensional struct of distribution parameter not allowed!')
end

[m,v] = deal(params{:});
if isempty(m),  m=0;  end
if isempty(v),  v=1;  end

try
  v(v<=0) = nan;
  xn2 = -0.5*(x-m).^2./v;
  if options.logp
    f = xn2-0.5*log(2*pi*v);
  else
    f = exp(xn2)./sqrt(2*pi*v);
  end
catch
  error ('x, m and v must be of common size or scalar');
end

% [errorcode, x, m, v] = comnsize (x,m, v);
% if (errorcode > 0)
% catch
%   error ('x, m and v must be of common size or scalar');
% end
% f=zeros(size(x));
% 
% k = find (v>0);
% if any(k)    
%   f(k)=1./sqrt(2*pi*v(k)).*exp(-0.5*(x(k)-m(k)).^2./v(k));
% end
% 
% k1 = find (v<=0);
% if any (k1)
%   tmp=NaN;
%   f(k1) = tmp(ones(size(k1)));
% end

