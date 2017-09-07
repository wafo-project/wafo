function f = pdfgam(x,varargin)
%PDFGAM Gamma probability density function
%
% CALL:  f = pdfgam(x,a,b,options);
%        f = pdfgam(x,phat,options);
%
%        f = density function evaluated at x
%        a = parameter
%        b = parameter (default b=1)
%     phat = Distribution parameter struct
%            as returned from FITGAM.  
%  options = struct with fieldnames:
%         .logp   : if TRUE, density, p, returned as log(p).
%
% The Gamma distribution is defined by its pdf
% 
%        f(x)=x^(a-1)*exp(-x/b)/gamma(a)/b^a, a,b>0, x>=0.
%
% Example: 
%   x = linspace(0,7,200);
%   p1 = pdfgam(x,1); p2 = pdfgam(x,2);
%   plot(x,p1,x,p2), shg
%
% See also  pdfgengam, cdfgam, invgam, rndgam, fitgam, momgam

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


% Reference: Johnson, Kotz and Balakrishnan (1994)
% "Continuous Univariate Distributions, vol. 1", p. 494 ff
% Wiley

% Tested on; Matlab 5.3
% History: 
% revised pab 24.10.2000
%  - added comnsize, nargchk
%  - replaced code with a call to pdfgengam -> maintainance easier.
% added ms 26.06.2000
% added b parameter ms 23.08.2000

options = struct('logp',false); % default options
if ischar(x) && strcmpi(x,'defaults')
  f = options;
  return
end
%error(nargchk(2,inf,nargin))
narginchk(2,inf)
Np = 2;
[params,options] = parsestatsinput(Np,options,varargin{:});
if numel(options)>1
  error('Multidimensional struct of distribution parameter not allowed!')
end

[a,b] = deal(params{:});
if isempty(b),  b=1; end
c = 1;
try
  f = pdfgengam(x,a,b,c,options);
catch
  error('x, a and b must be of common size or scalar.');
end


function y = lpdfgam(x,r,b)


% Reference
% Catherine Loader (2000). 
% "Fast and Accurate Computation of Binomial Probabilities"; 
% http://www.herine.net/stat/software/dbinom.html.
% @misc{ july-fast,
%   author = "Catherine Loader July",
%   title = "Fast and Accurate Computation of Binomial Probabilities",
%   url = "citeseer.ist.psu.edu/312695.html" }

if (r<1);
  y = pdfpois(r,x./b,'disable',true).*r./x;
else
  y = pdfpois(r-1,x./b,'disable',true)/b;
end
  



