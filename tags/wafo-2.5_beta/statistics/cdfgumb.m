function [F,Flo,Fup] = cdfgumb(x,varargin)
%CDFGUMB Gumbel cumulative distribution function.
%
% CALL:  F = cdfgumb(x,a,b,options)
%        [F,Flo,Fup] = cdfgumb(x,phat,options)
%   
%     F    = Gumbel cdf evaluated at x
%  Flo,Fup = 100*(1-alpha) % confidence bounds of F.
%    a, b  = parameters of the Gumbel distribution.
%     phat = Distribution parameter struct
%             as returned from FITGUMB.  
%  options = struct with fieldnames:
%           .lowertail: if TRUE (default), F = Prob[X <= x],
%                       otherwise, F = Prob[X > x].
%           .logp     : if TRUE, probability, p, returned as log(p).
%           .alpha    : Confidence coefficent    (default 0.05)
%           .trunc    : if TRUE truncated gumbel distribution
%                       otherwise regular gumbel distribution (default)
%                        
%  Gumbel CDF  is given by :                           
%       F(x) = exp(-exp(-(x-b)/a))    -inf < x < inf, a>0
%  or the truncated 
%       F(x) = [exp(-exp(-(x-b)/a))-exp(-exp(b/a))]/(1-exp(-exp(b/a))) 
%           0 < x < inf,  a>0
%
% Example: 
%   x = linspace(-4,6,200);
%   p1 = cdfgumb(x,2,0); p2 = cdfgumb(x,1,1);
%   plot(x,p1,x,p2), shg
%
% See also  pdfgumb, invgumb, rndgumb, fitgumb, momgumb, plotgumb



%  tested on: matlab 5.2
% history
% revised pab 8.11.1999
% updated header info
%   Per A. Brodtkorb 17.10.98
% rewritten ms 19.06.2000
% revised pab 25.10.2000
% - added nargchk+comnsize

% Reference: 
%  Johnson  N.L., Kotz S. and Balakrishnan, N. (1994)
%  Continuous Univariate Distributions, Volume 2. Wiley. 


error(nargchk(2,7,nargin))
Np = 2;
options = struct('covariance',[],'alpha',0.05,...
  'lowertail',true,'logp',false,'trunc',false); % default options
[params,options] = parsestatsinput(Np,options,varargin{:});
if numel(options)>1
  error('Multidimensional struct of distribution parameter not allowed!')
end

a = params{1};
b = params{2};
a(a<=0) = nan;



try
  if options.trunc
    x(x<0) = inf;
    expba = exp(b./a);
  end
  xn = (x -b)./a;
catch
   error('x, a and b must be of common size or scalar.');
end

if options.trunc
  logF = -expba+log(expm1(-exp(-xn)+expba))-log1p(-exp(-expba));
else
  logF = -exp(-xn);
end
    
   
if options.lowertail
  if options.logp
    F = logF;
  else
    F = exp(logF);
  end
else
  F = -expm1(logF);
      
  if options.logp
    sml = logF<-1;
    F(sml) = log1p(-exp(logF(sml)));
    lrg = -1<=logF;
    F(lrg) = log(F(lrg));
  end  
end
 


if nargout>1
% TODO % Implement  Flo and Fup
 warning('WAFO:CDFGUMB','Flo and Fup not implemented yet')
 Flo = nan;
 Fup = Flo;
end


% [iscmn x a b] = iscomnsize(x,a,b);
% if ~iscmn
%     error('x, a and b must be of common size or scalar.');
% end
% F=zeros(size(x));
% k1 = find(a>0);
% if any(k1),
%   F(k1)=exp(-exp(-(x(k1) -b(k1))./a(k1)) );
%   if trunc,
%     tmp=exp(-exp(b(k1)./a(k1)));
%     F(k1)=(F(k1)-tmp)./(1-tmp).*(x(k1)>0);
%   end
% end
% 
% k2=find(a<=0);
% if any(k2)
%   F(k2)=NaN;
% end
