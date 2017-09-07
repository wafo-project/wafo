function [F,Flo,Fup] = cdfgenpar(x,varargin)
%CDFGENPAR Generalized Pareto cumulative distribution function
%
% CALL:  F = cdfgenpar(x,k,s,m,options);
%        [F,Flo,Fup] = cdfgenpar(x,phat,options);
% 
%        F = distribution function evaluated at x
%  Flo,Fup = 100*(1-alpha) % confidence bounds of F.
%        k = shape parameter in the GPD
%        s = scale parameter in the GPD    (default 1)
%        m = location parameter in the GPD (default 0)
%     phat = Distribution parameter struct
%            as returned from FITGENPAR.  
%  options = struct with fieldnames:
%         .lowertail: if TRUE (default), F = Prob[X <= x],
%                     otherwise, F = Prob[X > x].
%         .logp     : if TRUE, probability, p, returned as log(p).
%         .alpha    : Confidence coefficent    (default 0.05)
% 
% The Generalized Pareto distribution is defined by its cdf
%
%                1 - (1-k*(x-m)/s)^1/k,  k~=0
%  F(x;k,s,m) =
%                1 - exp(-(x-m)/s),  k==0
%
% for  s>0 and 0 <= x-m and k*(x-m)/s <= 1.
%
% Example: 
%   x = linspace(0,2,200);
%   p1 = cdfgenpar(x,1.25,1); p2 = cdfgenpar(x,1,1);
%   p3 = cdfgenpar(x,0.75,1); p4 = cdfgenpar(x,0.5,1);
%   p5 = cdfgenpar(x,-1.25,1); p6 = cdfgenpar(x,-1,1);
%   p7 = cdfgenpar(x,-0.75,1); p8 = cdfgenpar(x,-0.5,1);
%
%   subplot(211);plot(x,p1,x,p2,x,p3,x,p4),title('k>0');
%   subplot(212);plot(x,p4,x,p5,x,p6,x,p8),title('k<0');
%
%   close all;
%
% See also pdfgenpar, invgenpar, rndgenpar, fitgenpar, momgenpar

%     This program is free software; you can redistribute it and/or modify
%     it under the terms of the GNU Lesser General Public License as published by
%     the Free Software Foundation; either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU Lesser General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.

% References
%  Johnson  N.L., Kotz S. and Balakrishnan, N. (1994)
%  Continuous Univariate Distributions, Volume 1. Wiley. 

% Tested on: Matlab 5.3
% History: 
% Revised pab 2007
% - removed threshold defining k to zero
% - more accurate evaluation
% - New input options struct.
% Revised pab Nov2005
% - added confidence interval Flo,Fhi
% Revised pab oct 2005
% - limits m<x<s/k replaced with 0 <= x-m < s/k in help header.
% revised pab June 2005
% - distribution for k==0 was wrong, now fixed
% Revised by jr 22.12.1999
% Modified by PJ 08-Mar-2000
%   Hjelptext
% revised ms 14.06.2000
% - updated header info
% - changed name to cdfgenpar (from gpdcdf)
% revised pab 24.10.2000
% - added  nargchk, comnsize and default values for m, s 
% revised jr 14.08.2001
% - a bug in the last if-statement condition fixed
%   (thanks to D Eddelbuettel <edd@debian.org>)

%error(nargchk(2,inf,nargin))
narginchk(2,inf)
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
  xm(xm<0)   = 0;
  kxm        = min(k.*xm,1);
  %kxm(kxm>1) = 1;
  [iscmn, k,xm] = comnsize(k,xm);
catch
  error('x, k, s and m must be of common size or scalar.');
end
logR = -xm;
k1  = find((0<xm) & (kxm <= 1 ) & (k~=0)); %abs(k)>epsilon));
if any(k1),
  logR(k1) = -xm(k1).*log1pxdx(-kxm(k1));
  %logR(k1) = log1p(-kxm(k1))./k(k1);
  %logR(k1) = log((1-kxm(k1)).^(1./k(k1)));
end
if options.lowertail
  F = -expm1(logR); 
else
  F = exp(logR); 
end


if nargout >= 2
  % Compute confidence bounds on the logit scale
  pcov = phat.covariance;
  if options.proflog
    Flo = F;
    Fup = F;
    CI = zeros(numel(F),2);
    for ix =1:numel(F)
      [CI(ix,:)] = ciproflog(phat,'i',2,'logR',logR(ix),'link',@link,'check',false);
    end
    if options.logp
      if options.lowertail
         Flo(:) = log1p(-exp(CI(:,1)));
         Fup(:) = log1p(-exp(CI(:,2)));
      else
         Flo(:) = (CI(:,1));
        Fup(:) = (CI(:,2));
      end
    else
      if options.lowertail
        Flo = -expm1(CI(:,1));
        Fup = -expm1(CI(:,2));
      else
        Flo(:) = exp(CI(:,1));
        Fup(:) = exp(CI(:,2));
      end
    end
  else
  alpha1 = options.alpha;
  
  % Approximate the variance of F on the logit scale
  if options.lowertail
    logitp = log(F./(1-F));
  else
    logitp = -log(F./(1-F));
  end
  
  dL   = 1 ./ (F.*(1-F)+realmin); % derivative of logit(F) w.r.t. F
  dLdK = dGPDdK(xm,k) .* dL;      % dlogitp/dk = dF/dk * dlogitp/dF
  dLdS = dGPDdS(xm,k,s) .* dL;    % dlogitp/ds = dF/ds * dlogitp/dF
  varLogitp = pcov(1,1).*dLdK.^2 + 2.*pcov(1,2).*dLdK.*dLdS + pcov(2,2).*dLdS.^2;
    
  if any(varLogitp(:) < 0)
    error('Covariance must be a positive semi-definite matrix.');
  end
    
  % Use a normal approximation on the logit scale, then transform back to
  % the original CDF scale
  zcrit = -invnorm(alpha1/2) * sqrt(varLogitp);
  zlo = logitp - zcrit;
  zup = logitp + zcrit;
  
  if options.logp
    Flo = -log1p(exp(-zlo));
    Fup = -log1p(exp(-zup));
  else
    Flo = 1 ./ (1 + exp(-zlo));
    Fup = 1 ./ (1 + exp(-zup));
  end
  end
end
 
if options.logp
  if options.lowertail
    sml = logR<-1;
    F(sml) = log1p(-exp(logR(sml)));
    lrg = -1<logR;
    F(lrg) = log(-expm1(logR(lrg)));
  else
    F = logR;
  end
end


return


function dk = dGPDdK(xm,k)
%DGPDDK Derivative of GPD cdf w.r.t shape parameter, k.

%epsilon = 1e-4; % treshold defining k to zero

kxm = k.*xm;
tmp = log1p(-kxm+realmin).*xm;
dk  = -exp(-xm).*tmp; % k==0


k1  = find((0<xm) & (kxm < 1 ) & k~=0); %(abs(k)>epsilon));
if any(k1),
  tmp1 = -exp(log1p(-kxm(k1))./k(k1)).*tmp(k1);
  tmp1(kxm(k1)==-inf) = 0;      %  
  dk(k1) = tmp1;
  %dk(k1) = -(1-kxm(k1)).^(1./k(k1)).*tmp(k1);
end
k2 = find((kxm>=1));
if any(k2),
  dk(k2) = 0;
end


function ds = dGPDdS(xm,k,s)
%DGPDDK Derivative of GPD cdf w.r.t scale parameter, s.

%epsilon = 1e-4; % treshold defining k to zero

xms = xm./s;
ds  = -exp(-xm).*xms; % k==0

kxm = k.*xm;
k1  = find((0<xm) & (kxm < 1 ) & k~=0); %(abs(k)>epsilon));
if any(k1),
  tmp = -exp((1./k(k1)-1).*log1p(-kxm(k1))).*xms(k1);
  tmp(kxm(k1)==-inf) = 0;      %  
  ds(k1) = tmp;
  %ds(k1) = -(1-kxm(k1)).^((1-k(k1))./k(k1)).*xms(k1);
end
k2 = find((kxm>=1));
if any(k2),
  ds(k2) = 0;
end




function  s  = link(x,logR,phat,ix)

if numel(phat)<3
  u = 0;
else
  u = phat(3);
end

if phat(1)~=0
  s =  -(x-u).*phat(1)./expm1(phat(1).*logR);
else
  s =  -(x-u)./logR;
end
