%CDFCAUCHY Cauchy cumulative distribution function
%
% CALL:  F = cdfcauchy(x,a,c,options);
%       [F, Flo,Fup] = cdfcauchy(x,phat,options);
%
%        F = distribution function evaluated at x
%  Flo,Fup = 100*(1-alpha) % confidence bounds of F.
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
%   x = linspace(0,6,200);
%   p1 = cdfcauchy(x,1,1); p2 = cdfcauchy(x,2,2);
%   plot(x,p1,x,p2), shg
%
% See also pdfcauchy, invcauchy, rndcauchy, fitcauchy, momcauchy


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


function [F,Flo,Fup] = cdfcauchy(x, varargin)
%error(nargchk(1,inf,nargin))
narginchk(1,inf)
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

try
  xn = ((x-loc)./scale);
catch
   error ('x, a and c must be of common size or scalar');
end


% TODO % Implement Flo and Fup 
if nargout>=2
  % Compute confidence bounds on log scale
  Flo = nan;
  Fup = nan;
end

if options.lowertail
  F  = atan(xn)/pi+0.5;
else
  F  = atan(-xn)/pi+0.5;
end

if options.logp
  F = log(F);
end