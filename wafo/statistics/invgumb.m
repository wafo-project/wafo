function [x,xlo,xup] = invgumb(F,varargin)
%INVGUMB Inverse of the Gumbel distribution function.
%
% CALL:  x = invgumb(F,a,b,options) 
%        [x,xlo,xup] = invgumb(F,phat,options)
%
%     x    = Inverse Gumbel evaluated at F
%  xlo,xup = 100*(1-alpha) % confidence bounds of x.
%    a, b  = parameters of the Gumbel distribution.
%        F = lower or upper tail probability 
%     phat = Distribution parameter struct
%             as returned from FITGUMB.  
%   options = struct with fieldnames:
%           .lowertail: if TRUE (default), F = Prob[X <= x],
%                       otherwise, F = Prob[X > x].
%           .logp     : if TRUE, probability, p, input as log(p).
%           .alpha    : Confidence coefficent    (default 0.05)
%           .proflog  : if TRUE compute  xlo and xup using proflog
%           .trunc    : if TRUE truncated gumbel distribution
%                       otherwise regular gumbel distribution (default)
%
%    Gumbel cdf  is given by :    
%            F(x) = exp(-exp(-(x-b)/a)    -inf < x < inf,  a>0
%    or the truncated
%           F(x) = [exp(-exp(-(x-b)/a)) -exp(-exp(b/a)) ]/(1-exp(-exp(b/a)))
%       0 < x < inf, a>0     
%
% Example: 
%  a=1;b=2;    
%  opt = {'lowertail',false,'logp',false};
%  F0 = [logspace(-300,-1) linspace(0.11,0.5)];
%  x  = invgumb(F0,a,b,opt{:});
%  F  = cdfgumb(x,a,b,opt{:});
%  semilogy(abs(F-F0)./F0+eps); % relative error
%
%   close all;
%
% See also  pdfgumb, cdfgumb, rndgumb, fitgumb, momgumb, plotgumb

%
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


% Reference: 
%  Johnson  N.L., Kotz S. and Balakrishnan, N. (1994)
%  Continuous Univariate Distributions, Volume 2. Wiley. 


%  tested on: matlab 5.2
% history
% revised pab 8.11.1999
% updated header info
%   Per A. Brodtkorb 17.10.98
% rewritten ms 19.06.2000
% revised pab 25.10.2000
% - added nargchk+comnsize

error(nargchk(3,7,nargin))
Np = 2;
options = struct('proflog',false,'alpha',0.05,...
  'lowertail',true,'logp',false,'trunc',false); % default options
[params,options,tmp,phat] = parsestatsinput(Np,options,varargin{:});
if numel(options)>1
  error('Multidimensional struct of distribution parameter not allowed!')
end

a = params{1};
b = params{2};
a(a<=0) = nan;

if options.logp
  F(F>0) = NaN;
  if options.lowertail
    logF = F;
  else
    logF = log(expm1(F));
    lrg = 0.5 < F;
    logF(lrg) = log1p(-exp(F(lrg)));
  end
else
  F(F<0 | 1<F) = NaN;
  if options.lowertail
    logF = log(F);
  else
    logF = log1p(-F);
  end
end

try
  if options.trunc
    expba = exp(b./a);
    %tmp=exp(-exp(b(k1)./a(k1)));
    x =-a.*log(-log(-expm1(-expba).*exp(logF) +exp(-expba)) ) + b; 
  else  
    x =-a.*log(-logF) + b;
  end
catch
   error('F, a and b must be of common size or scalar.');
end

if nargout>1
  if options.proflog
    xlo = x;
    xup = x;
    for ix =1:numel(x)
      [Lp,CI] = proflog(phat,'i',1,'x',x(ix),'link',@lnkgumb,'alpha',alpha);
      xlo(ix) = CI(1);
      xup(ix) = CI(2);
    end
  else
% TODO % Implement  xlo and xup 
    warning('WAFO:INVGUMB','xlo and xup not implemented yet')
    xlo = nan;
    xup = xlo;
  end
end



% [icode F a b] = iscomnsize(F,a,b);
% if ~icode
%     error('F, a and b must be of common size or scalar.');
% end
% x=zeros(size(F));
% 
% ok = (0<=F& F<=1 & a>0);
% 
% k1=find((F>0)&(F<1) & ok);
% 
% if any(k1)
%   if trunc,
%     tmp=exp(-exp(b(k1)./a(k1)));
%     x(k1) =-a(k1).* log(-log((1-tmp).* F(k1) +tmp) ) + b(k1); 
%   else
%     x(k1) =-a(k1).* log(-log( F(k1)) ) + b(k1);
%   end
% end
% tmp=Inf;
% k2=find(F==0 & ok);
% if any(k2)
%   if trunc
%     x(k2)=zeros(size(k2));
%   else
%     x(k2)=-tmp(ones(size(k2)));
%   end
% end
% k3=find(F==1& ok);
% if any(k3)
%   x(k3)=tmp(ones(size(k3)));
% end
% 
% k4= find(~ok);
% if any(k4)
%   x(k4)=NaN;
% end
