function f = pdfpois(x,varargin)
%PDFPOIS  Poisson probability mass function
%
%  CALL:  f = pdfpois(x,L,options)
%         f = pdfpois(x,phat,options)
% 
%        f = probability of observing x givel L
%        x = 
%        L = mean and variance of the poisson distribuion
%     phat = Distribution parameter struct
%            as returned from FITPOIS.  
%  options = struct with fieldnames:
%         .logp : if TRUE, density, p, returned as log(p).
%
% The Poisson probability mass function is defined by:
%
%        f(x) = L^x exp(-L)/x!, 0<=L, x=0,1,2,....
%
%  Example
% -log(pdfpois(0:7, 1) .* gamma(1+ (0:7))) % =1
%
% See also cdfpois, invpois, rndpois, fitpois, mompois

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



% Reference
% Catherine Loader (2000). 
% "Fast and Accurate Computation of Binomial Probabilities"; 
% http://www.herine.net/stat/software/dbinom.html.
% @misc{ july-fast,
%   author = "Catherine Loader July",
%   title = "Fast and Accurate Computation of Binomial Probabilities",
%   url = "citeseer.ist.psu.edu/312695.html" }

options = struct('logp',false,'disable',false); % default options
if ischar(x) && strcmpi(x,'defaults')
  f = options;
  return
end
error(nargchk(2,inf,nargin))

Np = 1;
[params,options] = parsestatsinput(Np,options,varargin{:});
if numel(options)>1
  error('Multidimensional struct of distribution parameter not allowed!')
end

[lb] = deal(params{:});
 lb(lb<0) = nan;
  x(x<0) = inf;
 if ~options.disable
   %x = round(x); % force int
   x(x~=floor(x)) = inf;
 end
try
  f = exp(-lb).* (x==0);
  if options.logp
    f = log(f);
  end
catch
  error('x and lb must be of common size or scalar.')
end

px =x>0 & lb>0;
if any(px) 
  PI2 = 6.283185307179586476925286;
  if ~isscalar(x), x = x(px) ; end
  if ~isscalar(lb), lb = lb(px) ; end
  if options.logp
    f(px) = -stirlerr(x)-bd0(x,lb) - 0.5.*log(PI2.*x);
  else
    f(px) = exp(-stirlerr(x)-bd0(x,lb))./sqrt(PI2.*x);
  end
end


