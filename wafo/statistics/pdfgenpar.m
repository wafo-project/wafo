function f = pdfgenpar(x,varargin)
%PDFGENPAR Generalized Pareto probability density function
%
% CALL:  f = pdfgenpar(x,k,s,m,options);
%        f = pdfgenpar(x,phat,options);
% 
%        f = density function evaluated at x
%        k = shape parameter in the GPD
%        s = scale parameter in the GPD    (default 1)
%        m = location parameter in the GPD (default 0)  
%     phat = Distribution parameter struct
%            as returned from FITGENPAR.  
%  options = struct with fieldnames:
%         .logp : if TRUE, density, p, returned as log(p).
%
% The Generalized Pareto distribution is defined by its cdf
%
%                1 - (1-k(x-m)/s)^1/k,  k~=0
%  F(x;k,s,m) =
%                1 - exp(-(x-m)/s),  k==0
%
% for s > 0, 0 <= x-m and  k*(x-m)/s < 1.
%
% Example: 
%   x = linspace(0,2,200);
%   p1 = pdfgenpar(x,1.25,1); p2 = pdfgenpar(x,1,1);
%   p3 = pdfgenpar(x,0.75,1); p4 = pdfgenpar(x,0.5,1);
%   p5 = pdfgenpar(x,0,1);    p6 = pdfgenpar(x,-0.5,1);
%   p7 = pdfgenpar(x,-0.75,1);p8 = pdfgenpar(x,-1,1);
%   
%   subplot(211); plot(x,p1,x,p2,x,p3,x,p4); title('k>0');
%   subplot(212); plot(x,p5,x,p6,x,p7,x,p8); title('k<=0');
%
%   close all;
%
% See also cdfgenpar, invgenpar, rndgenpar fitgenpar, momgenpar
 
% References
%  Johnson N.L., Kotz S. and Balakrishnan, N. (1994)
%  Continuous Univariate Distributions, Volume 1. Wiley. 

% Tested on: Matlab 5.3
% History: 
% Revised pab Aug 2007
% -New implementation avoiding threshold defining k to zero.
% -added log option
% Revised pab oct 2005
% - limits m<x<s/k replaced with 0 <= x-m < s/k in help header.
% -fixed a bug for k<0 
% revised pab June 2005
% fixed bug for k==0 now fixed
% revised jr 14.08.2001
%  - a bug in the last if-statement condition fixed
%    (thanks to D Eddelbuettel <edd@debian.org>)
% revised pab 24.10.2000
%  - added comnsize, nargchk
%  - added extra check on the scale parameter s
%  - added m + default values on m and s
% added ms 14.06.2000

%error(nargchk(2,inf,nargin))
narginchk(2,inf)
Np = 3;
options = struct('logp',false); % default options
[params,options] = parsestatsinput(Np,options,varargin{:});
if numel(options)>1
  error('Multidimensional struct of distribution parameter not allowed!')
end

[k,s,m] = deal(params{:});
if isempty(s), s=1;end
if isempty(m), m=0;end

s(s<=0) = nan;

try
  xn = (x-m)./s;
  xn(xn<0)   = inf;  % trick to force f to zero.
  kxn        = k.*xn;
  kxn(kxn>1) = 1;
  kxn(k==0 & xn==inf) = 0;
 
  %f = exp(-xn)./s;                   % for  k==0
  %f = (1-k.*xn).^(1./k-1)/s;        % for  k~=0 
  %f = exp((1./k-1).*log1p(-kxn))/s  % for  k~=0
  
  %f = exp((-xn+kxn).*log1p(-kxn)./(-kxn))/s  % for any k kxn~=-inf
  logfs = (-xn+kxn).*log1pxdx(-kxn);
  % Handle special cases
  logfs(kxn==-inf | kxn==1) = -inf;      %  
  logfs(k==1 & xn==1) = 0; % 0^0 situation
  
 if options.logp
   f = logfs-log(s);
 else
   f = exp(logfs)./s;
 end
  
catch
  error('x, k s and m must be of common size or scalar.');
end


% 
% Old calls kept just in case 
% fs = zeros(size(kxn));
% iy = (xn>=0) & (abs(k)==0) & (s>0);
% if any(iy)
%   fs(iy) = exp(-xn(iy));
% end
% ix = (xn>=0) & (kxn < 1) & (k~=0) & (s>0);
% if any(ix)
%   if ~isscalar(k)
%     k = k(ix);
%   end
%   fs(ix) = exp((1./k-1).*log1p(-kxn(ix)));
% end
% % Handle special cases
% fs(kxn==1 & xn==1) = 1; % 0^0 situation
% f = fs./s;
% return
%

% [errorcode x k s m] = comnsize(x,k,s,m);
% if errorcode > 0
%     error('x, k s and m must be of common size or scalar.');
% end
% epsilon=1e-4; % treshold defining k to zero
% 
% xm = x-m;
% 
% 
% f = zeros(size(x));
% 
% k0 = find((xm>=0) & (abs(k)<=epsilon) & (s>0));
% if any(k0), 
%    f(k0) = exp(-xm(k0)./s(k0))./s(k0);
% end
% 
% k1=find((xm>=0) & (k.*xm < s) & (abs(k)>epsilon) & (s>0));
% if any(k1),
%   f(k1) = (1-k(k1).*xm(k1)./s(k1)).^(1./k(k1)-1)./s(k1);
% end
% 
%   
% k3 = find(s<=0);
% if any(k3),
%    tmp   = NaN;
%    f(k3) = tmp(ones(size(k3)));
% end
% 
% return
% 
