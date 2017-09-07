%PDFCAUCHY  Cauchy's probability density function
%
% CALL:  f = pdfcauchy(x,a,c,options)
%
%          f = density function evaluated at x
%          x = matrix
%          a = scale parameter    (default 1)
%          c = location parameter (default 0)
%       phat = Distribution parameter struct
%              as returned from FITWEIB.  
%    options = struct with fieldnames:
%           .logp : if TRUE, density, p, returned as log(p).
%  
%   The Cauchy distribution is defined by its pdf
%
%   f(x;a,c) = 1/(pi*a*(1+((x-c)/a).^2))
%
% Example:
%   x = linspace(-5,5,200);
%   p1 = pdfcauchy(x,1); p2 = pdfcauchy(x,5);
%   plot(x,p1,x,p2,x,pdfnorm(x)),shg
%
% See also cdfcauchy, invcauchy, rndcauchy, fitcauchy, momcauchy

% Copyright (C) 2007 Per A. Brodtkorb
%
% This file is part of WAFO.
%
% WAFO is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2, or (at your option)
% any later version.
%
% WAFO is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with WAFO; see the file COPYING.  If not, write to the Free
% Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
% 02110-1301, USA.


% tested on matlab 7
%History:
% by pab 2007

function f = pdfcauchy(x,varargin)
%error(nargchk(1,inf,nargin))
narginchk(1,inf)
Np = 2;
options = struct('logp',false,'disable',false); % default options
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
  xn2 = ((x-loc)./scale).^2;
catch
   error ('x, a and c must be of common size or scalar');
end

if options.logp
  f = -log1p(xn2) - log(pi*scale);
else
  f = 1./(pi*scale*(1+xn2));
end

