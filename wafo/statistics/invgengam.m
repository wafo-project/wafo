function [x,xlo,xup] = invgengam(F,varargin)
%INVGENGAM Inverse of the Generalized Gamma distribution function
%
% CALL:  x = invgengam(F,a,b,c)
%
%        x = inverse cdf for the Generalized Gamma distribution evaluated at F
%   a,b,c  = parameters (see pdfgengam)
%  xlo,xup = 100*(1-alpha) % confidence bounds of x.
%        a = parameter, a>0
%        b = parameter, b>0 (default b=1)
%        c = 
%     phat = Distribution parameter struct
%            as returned from FITGAM.  
%  options = struct with fieldnames:
%         .lowertail: if TRUE (default), F = Prob[X <= x],
%                     otherwise, F = Prob[X > x].
%         .logp     : if TRUE, probability, p, input as log(p).
%         .alpha    : Confidence coefficent        (default 0.05)
%         .abseps   : Requested absolute error     (default 1e-90)
%         .releps   : Requested relative error     (default sqrt(eps))
%         .max_iter : Maximum number of iterations (default 500)
%
%
% Example:
%   x = linspace(0,1,200);
%   p1 = invgengam(x,0.5,1,1); p2 = invgengam(x,1,2,1);
%   p3 = invgengam(x,2,1,2); p4 = invgengam(x,2,2,2);
%   plot(x,p1,x,p2,x,p3,x,p4);
%
%   close all;
%
% See also  pdfgengam, cdfgengam, rndgengam, fitgengam, momgengam

% Reference: Cohen & Whittle, (1988) "Parameter Estimation in Reliability
% and Life Span Models", p. 220 ff, Marcel Dekker.

% Tested on; Matlab 5.3
% History: 
% adapted from stixbox ms 10.08.2000
% revised pab 23.10.2000
%  - added comnsize, nargchk 
% revised pab 4nov2005
%  -improved the starting guess for x.

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


error(nargchk(2,inf,nargin))
options = struct('proflog',false,'alpha',0.05,...
  'lowertail',true,'logp',false,'abseps',1e-90,'releps',sqrt(eps),'max_iter',100); % default options

Np = 3;
[params,options,tmp,phat] = parsestatsinput(Np,options,varargin{:});
if numel(options)>1
  error('Multidimensional struct of distribution parameter not allowed!')
end

[a,b,c] = deal(params{:});
if isempty(b), b=1;end
if isempty(c), c=1;end
[icode] = iscomnsize(F,a,b,c);
if ~icode 
  error('F, a b and c must be of common size or scalar.');
end
  
x = invgam(F,a,1,options); 
b(b<=0) = nan;
c(c<=0) = nan;
x = x.^(1./c).*b;

if nargout>1
  if isempty(phat)
    error('Must have distribution struct!')
  end
  alpha = options.alpha;
  if options.proflog
    error('Confidence interval usin proflog not implemented!')
    xlo = x;
    xup = x;
    for ix =1:numel(x)
      CI = ciproflog(phat,'i',1,'x',x(ix),'link',@lnkgengam,'alpha',alpha);
      xlo(ix) = CI(1);
      xup(ix) = CI(2);
    end
  else
      pcov = phat.covariance;

% TODO % Implement  xlo and xup
    warning('xlo and xup not implemented yet')
    xlo = nan;
    xup = xlo;
  end
end
