function [x, xlo,xup] = invf(F,varargin)
%INVF  Inverse of the Snedecor's F distribution function
%
% CALL:  x = invf(F,df1,df2,options)
%        [x,xlo,xup] = invf(F,phat,options)
%
%   x      = inverse cdf for the F distribution evaluated at F.
%  xlo,xup = 100*(1-ALPHA)% confidence interval for X.
% df1, df2 = degrees of freedom (1,2,....)
%     phat = Distribution parameter struct
%            as returned from WFFIT.  
%  options = struct with fieldnames:
%         .lowertail: if TRUE (default), F = Prob[X <= x],
%                     otherwise, F = Prob[X > x].
%         .logp     : if TRUE, probability, p, input as log(p).
%         .alpha    : Confidence coefficent        (default 0.05)
%         .abseps   : Requested absolute error     (default 1e-90)
%         .releps   : Requested relative error     (default sqrt(eps))
%         .max_iter : Maximum number of iterations (default 100)
%
% Example:
%   df1=1;df2=2;    
%   opt = {'lowertail',true,'logp',false};
%   F0 = [logspace(-300,-1) linspace(0.11,0.5)];
%   x  = invf(F0,df1,df2,opt{:});
%   F  = cdff(x,df1,df2,opt{:});
%   semilogy(abs(F-F0)./F0+eps); % relative error
%
%   close all;
%
% See also pdff, cdff, rndf, wffit, momf


% tested on matlab 5.3
%History:
%revised pab 29.10.2000
% adapted from stixbox
% -added nargchk, comnsize
%        Anders Holtsberg, 18-11-93
%        Copyright (c) Anders Holtsberg

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


%error(nargchk(3,14,nargin))
narginchk(3,14)
options = struct('covariance',[],'alpha',0.05,...
  'lowertail',true,'logp',false,'abseps',1e-90,'releps',sqrt(eps),'max_iter',100); % default options

Np = 2;
[params,options] = parsestatsinput(Np,options,varargin{:});
if numel(options)>1
  error('Multidimensional struct of distribution parameter not allowed!')
end

[a,b] = deal(params{:});


a(a<=0 | floor(a)~=a) = nan;
b(b<=0 | floor(b)~=b) = nan;
try
  if nargout>1
    opt = options;
    opt.covariance = options.covariance/4;
    [q,qlo,qup] = invbeta(F,a/2,b/2,options);
     xlo = qlo.*b./((1-qlo).*a);
     xup = qup.*b./((1-qup).*a);
  else
    q = invbeta(F,a/2,b/2,options);
  end
  x = q.*b./((1-q).*a);
catch
  error('x, df1 and df2 must be of common size or scalar');
end

% return
% [icode F,a,b] = iscomnsize(F,a,b);
% if ~icode
%   error('x, df1 and df2 must be of common size or scalar');
% end
% 
% x = zeros(size(F));
% 
% ok = (a>0 & b>0 & floor(a)==a & floor(b)==b);
% 
% k = find(F>0&F<1 & ok);
% if any(k)
%   tmp = invbeta(F(k),a(k)/2,b(k)/2);
%   x(k) = tmp.*b(k)./((1-tmp).*a(k));
% end
% 
% 
% k2=find(F==1&ok);
% if any(k2)
%   x(k2)=inf;
% end
% 
% 
% k3=find(~ok);
% if any(k3)
%   tmp=NaN;
%   x(k3)=tmp(ones(size(k3)));
% end
