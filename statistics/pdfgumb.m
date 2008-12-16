function f = pdfgumb(x,varargin)
%PDFGUMB Gumbel probability density function.
%
% CALL:  f = pdfgumb(x,a,b,options) 
%        f = pdfgumb(x,phat,options) 
%   
%     f    = Gumbel pdf evaluated at x
%    a, b  = parameters of the Gumbel distribution.
%     phat = Distribution parameter struct
%            as returned from FITGUMB.  
%  options = struct with fieldnames:
%         .logp  : if TRUE, density, p, returned as log(p).
%         .trunc : if FALSE regular gumbel distribution (default)
%                  if TRUE  truncated gumbel distribution 
%
%  Gumbel PDF  is given by :                           
%      f(x) = exp(-(x-b)/a)*exp(-exp(-(x-b)/a))/a    -inf < x < inf, a>0
%  or the truncated
%      f(x) = exp(-(x-b)/a)*exp(-exp(-(x-b)/a))/a/(1-exp(-exp(b/a)))
%          0 < x < inf, a>0 
%
% Example: 
%   x = linspace(-4,8,200);
%   p1 = pdfgumb(x,2,0); p2 = pdfgumb(x,1,1);
%   plot(x,p1,x,p2),shg
%
% See also  cdfgumb, invgumb, rndgumb, fitgumb, momgumb

% Reference: 
%  Johnson  N.L., Kotz S. and Balakrishnan, N. (1994)
%  Continuous Univariate Distributions, Volume 2. Wiley. 


%  tested on: matlab 5.2
% history
% revised pab 24.10.2000
% - reimplemented comnsize
% revised pab 8.11.1999
% updated header info
%   Per A. Brodtkorb 17.10.98

error(nargchk(2,7,nargin))
Np = 2;
options = struct('logp',false,'trunc',false); % default options
[params,options] = parsestatsinput(Np,options,varargin{:});
if numel(options)>1
  error('Multidimensional struct of distribution parameter not allowed!')
end

a = params{1};
b = params{2};

a(a<=0)= nan;  
if options.trunc
  x(x<0) = inf;
end

try
  xn = max((x-b)./a,-realmax);
catch
  error('x, a and b must be of common size or scalar.');
end

if options.trunc
  logf = -xn-exp(-xn)-log(a) - log1p(-exp(-exp(b./a)));
else
  logf = -xn-exp(-xn)-log(a);
end
if options.logp
  f = logf;
else
  f = exp(logf);
end


return
% [icode x a b] = iscomnsize(x,a,b);
% if ~icode
%     error('x, a and b must be of common size or scalar.');
% end
% 
% f = zeros(size(x));
% 
% 
% if trunc,
%   k = find(x > 0 & a>0);
% else 
%   k = find(a>0);
% end
% 
% if any(k),
%   tmp=exp(-(x(k) -b(k))./a(k)  );
%   f(k) =  tmp.* exp(-tmp)./a(k);
%   if trunc,
%     f(k)=f(k)./(-expm1(-exp(b(k)./a(k)  )));
%   end
% end
% 
% k1 = find(a <= 0 );
% if any(k1)
%    tmp   = NaN;
%    f(k1) = tmp(ones(size(k1)));
% end
