function [F,Flo,Fup] = cdfgev(x,varargin)
%CDFGEV Generalized Extreme Value cumulative distribution function
%
% CALL:  F = cdfgev(x,k,s,m,options);
%        [F,Flo,Fup] = cdfgenpar(x,phat,options);
%
%        F = distribution function evaluated at x
%  Flo,Fup = 100*(1-alpha) % confidence bounds of F.
%        k = shape parameter in the GEV 
%        s = scale parameter in the GEV, s>0 (default 1)
%        m = location parameter in the GEV   (default 0)
%     phat = Distribution parameter struct
%            as returned from FITGEV.  
%  options = struct with fieldnames:
%         .lowertail: if TRUE (default), F = Prob[X <= x],
%                     otherwise, F = Prob[X > x].
%         .logp     : if TRUE, probability, p, returned as log(p).
%         .alpha    : Confidence coefficent    (default 0.05)
% 
% The Generalized Extreme Value distribution is defined by its cdf
%
%                exp( - (1 - k(x-m)/s))^1/k) ),  k~=0
%  F(x;k,s,m) =
%                exp( -exp(-(x-m)/s) ),  k==0
%
%  for s>0 and  k*(x-m)/s <= 1.
%
% Example: 
%   x = linspace(0,15,200);
%   m = 5;
%   p1 = cdfgev(x,0.8,1,m); p2 = cdfgev(x,0.8,2,m);
%   p3 = cdfgev(x,0.5,1,m); p4 = cdfgev(x,0.5,2,m);
%   p5 = cdfgev(x,-0.8,1,m); p6 = cdfgev(x,-0.8,2,m);
%   p7 = cdfgev(x,-0.5,1,m); p8 = cdfgev(x,-0.5,2,m);
%  
%   subplot(211),plot(x,p1,x,p2,x,p3,x,p4),title('k>0')
%   subplot(212),plot(x,p4,x,p5,x,p6,x,p8),title('k<0'),shg
%
% See also pdfgev, invgev, rndgev, fitgev, momgev

% References
%  Johnson  N.L., Kotz S. and Balakrishnan, N. (1994)
%  Continuous Univariate Distributions, Volume 1. Wiley. 

% Tested on; Matlab 5.3
% History: 
% Revised by jr 22.12.1999
% revised ms 14.06.2000
% - updated header info
% - changed name to cdfgev (from gevcdf)
% revised pab 24.10.2000
% - added  nargchk, comnsize and default values for m, s 
% revised jr 14.08.2001
% - a bug in the last if-statement condition fixed
%   (thanks to D Eddelbuettel <edd@debian.org>)
% revised pab may 2007
% -

error(nargchk(2,10,nargin))
options = struct('proflog',false,'alpha',0.05,...
  'lowertail',true,'logp',false); % default options
Np = 3;
[params,options,tmp,phat] = parsestatsinput(Np,options,varargin{:});
if numel(options)>1
  error('Multidimensional struct of distribution parameter not allowed!')
end

[k,s,m] = deal(params{:});
if isempty(s), s=1;end
if isempty(m), m=0;end


s(s<=0) = NaN;
%epsilon = 1e-4; % treshold defining k to zero

try
  xm         = (x-m)./s;
  kxm        = min(k.*xm,1);
  [csize, k,xm] = comnsize(k,xm);  
catch
  error('x, k, s and m must be of common size or scalar.');
end

logF = -exp(-xm); % k==0
k1  = find((kxm <= 1 ) & (k~=0)); %(abs(k)>epsilon));
if any(k1),
  %F(k1) = exp(-(1-kxm(k1)).^(1./k(k1)));
  logF(k1) = -exp(log1p(-kxm(k1))./k(k1));
end
if options.lowertail
  F = exp(logF);
else
  F = -expm1(logF);
end
if nargout>1
% TODO % Implement  Flo and Fup
 warning('WAFO:CDFGEV','Flo and Fup not implemented yet')
 Flo = nan;
 Fup = Flo(ones(size(F)));
 Flo = Fup;
end


if options.logp
  if options.lowertail
    F = logF;
  else
    sml = logF<-1;
    F(sml) = log1p(-exp(logF(sml)));
    lrg = -1<=logF;
    F(lrg) = log(-expm1(logF(lrg)));
  end
end



% old call

% [errorcode x k s,m] = comnsize(x,k,s,m);
% if errorcode > 0
%   error('x, k, s and m must be of common size or scalar.');
% end
%   
% epsilon=1e-4; % treshold defining k to zero
% 
% F = zeros(size(x));
% %k0 = find(x>=m & abs(k)<=epsilon & s>0);
% k0 = find( abs(k)<=epsilon & s>0);
% if any(k0),
%   F(k0) = exp(-exp(-(x(k0)-m(k0))./s(k0)));
% end
% 
% %k1=find(x>m&(k.*x<s+k.*m)&(abs(k)>epsilon));
% k1=find((k.*x<s+k.*m)&(abs(k)>epsilon));
% if any(k1),
%   tmp = (1-k(k1).*(x(k1)-m(k1))./s(k1));
%   F(k1)=exp(-tmp.^(1./k(k1)));
% end
% 
% k2=find((k.*x>=s+k.*m)&(k>epsilon));
% if any(k2),
%   F(k2)=ones(size(k2));
% end
% 
% k3 = find(s<=0 );
% if any(k3),
%    tmp   = NaN;
%    F(k3) = tmp(ones(size(k3)));
% end
% return
% 
