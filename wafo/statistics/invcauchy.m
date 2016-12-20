%INVCAUCHY Inverse of the Cauchy distribution function.
%
% CALL:  x = invcauchy(F,a,c,options);
%       [x, xlo,xup] = invcauchy(F,phat,options);
%
%        x = inverse cdf for the Cauchy distribution evaluated at F
%  xlo,xup = 100*(1-alpha) % confidence bounds of F.
%     a,c  = scale and location parameters (default a = 1, c =0)
%     phat = Distribution parameter struct
%            as returned from FITWEIB.  
%  options = struct with fieldnames:
%         .lowertail: if TRUE (default), F = Prob[X <= x],
%                     otherwise, F = Prob[X > x].
%         .logp     : if TRUE, probability, p, returned as log(p).
%         .alpha    : Confidence coefficent    (default 0.05)
%
% The Cauchy distribution is defined by its cdf
%
%  F(x;a,c) =atan((x-c)/a))/pi+1/2, a>0
%
%
% Example: 
%  a=1;b=2;
%  opt = {'lowertail',false,'logp',false};
%  F0 = [logspace(-300,-1) linspace(0.11,0.5)];
%  x  = invcauchy(F0,a,b,opt{:});
%  F  = cdfcauchy(x,a,b,opt{:});
%  semilogy(abs(F-F0)./F0+eps) % relative error
%  
%  close all;
%
% See also pdfcauchy, cdfcauchy, invcauchy, rndcauchy, fitcauchy, momcauchy


% Copyright (C) 2007 Per A. Brodtkorb
%
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


% tested on matlab 7
%History:
% by pab 2007


function [x,xlo,xup] = invcauchy(F, varargin)
error(nargchk(1,inf,nargin))

options = struct('covariance',[],'alpha',0.05,...
  'lowertail',true,'logp',false); % default options
Np = 2;
[params,options] = parsestatsinput(Np,options,varargin{:});
if numel(options)>1
  error('Multidimensional struct of distribution parameter not allowed!')
end

[scale,loc] = deal(params{:});
if isempty(scale)
  scale = 1;
end
if isempty(loc),
  loc = 0;
end

scale(scale <= 0) = NaN;
if options.logp
  if options.lowertal
    F = exp(F);
  else
    F = -expm1(F);
  end
elseif ~options.lowertail
  F = 1-F;
end

F(F<0 | 1<F) = nan;
try
  x = loc+scale.*tan(pi*(F-0.5));
catch
   error ('F, a and c must be of common size or scalar');
end


% TODO % Implement xlo and xup 
if nargout>=2
  % Compute confidence bounds on log scale
  xlo = nan;
  xup = nan;
end

