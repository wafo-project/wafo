function f = pdfraymod(x,varargin)
%PDFRAYMOD Truncated Rayleigh probability density function
%
% CALL:  f = pdfraymod(x,b,c);
%
%        f = density function evaluated at x
%        b = scale parameter
%        c = truncation parameter (default 0)
%     phat = Distribution parameter struct
%            as returned from FITRAYMOD.  
%  options = struct with fieldnames:
%         .logp   : if TRUE, density, p, returned as log(p).
% 
% The Truncated Rayleigh distribution is defined by its cdf
%
%  F(x;b,c) = 1 - exp(-(x+|c|)^2/(2b^2)+c^2/(2*b^2)), x>=0, b>0
%
% Example: 
%   x = linspace(0,4,200);
%   p1 = pdfraymod(x,1); p2 = pdfraymod(x,0.5,-2);
%   plot(x,p1,x,p2), shg
%
% See also cdfraymod, invraymod, rndraymod, fitraymod, momraymod

% Reference: Cohen & Whittle, (1988) "Parameter Estimation in Reliability
% and Life Span Models", p. 181 ff, Marcel Dekker.


% Tested on; Matlab 5.3
% History: 
% revised pab 24.10.2000
%  - added comnsize, nargchk
% added ms 15.06.2000

options = struct('logp',false); % default options
if (nargin==1 && nargout <= 1 && isequal(x,'defaults'))
  f = options; 
  return
end

error(nargchk(2,inf,nargin))

Np = 2;
[params,options] = parsestatsinput(Np,options,varargin{:});
if numel(options)>1
  error('Multidimensional struct of distribution parameter not allowed!')
end
[b,c] = deal(params{:});
if isempty(c),c=0;end
try
  c = abs(c);
  b(b<=0) = nan;
  x(x<0) = inf;
  xn = min((x+c)./b,realmax);
  if options.logp
    f = -0.5*(xn.^2 +(c./b).^2)+log(xn)-log(b);
  else
    f = xn.*exp(-0.5*(xn.^2 +(c./b).^2))./b;
  end
catch
   error ('x, b and c must be of common size or scalar');
end

% return
% if nargin<3||isempty(c),c=0;end
% [icode, x, b,c] = iscomnsize (x,b,abs(c));
% if ~icode
%   error ('x, b and c must be of common size or scalar');
% end
% 
% f=zeros(size(x));
% 
% k = find ((x>=0)&(b>0));
% if any (k)  
%   f(k)=(x(k)+c(k)).*exp(-((x(k)+c(k)).^2 -c(k).^2)./(2*b(k).^2))./b(k).^2;
% end
% 
% k1 = find (b<=0);
% if any (k1)
%   tmp=NaN;
%   f(k1) = tmp(ones(size(k1)));
% end
% 
