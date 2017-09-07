function  f = pdfbin(x,varargin)
%PDFBIN  Binomial probability mass function
%
%  CALL:  f = pdfbin(x,n,p,options)
%         f = pdfbin(x,phat,options)
% 
%        f = probability of observing x successes in n independent trials,
%            where the probability of success in any given trial is p.
%        x = number of successes          (0<=x<=n)
%        n = number of independent trials
%        p = probability of succes in a given trial. (0<=p<=1)
%     phat = Distribution parameter struct
%            as returned from FITBIN.  
%  options = struct with fieldnames:
%         .logp : if TRUE, density, p, returned as log(p).
%
% The binomial probability function is defined by:
%
%        f(x) = n!/(n-x)!/x!*p^x*(1-p)^(n-x)
%
%  Example
%   n = 50;        % Number of objects to test in a day
%   p = 0.05;      % Defect probability
%   defects = 0:n; % All possible outcomes of a day
%   f = pdfbin(defects,n,p);
%   stairs(defects,f);
%   [x,i] = max(f);
%   assert(defects(i), 2)  % Most likely number of defective objects found in one day.
%
% See also cdfbin, invbin, rndbin, fitbin, mombin

%       Anders Holtsberg, 16-03-95
%       Copyright (c) 1995 Anders Holtsberg, Per A. Brodtkorb 2007
%
% This program is free software; you can redistribute it and/or modify
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

[n,p] = deal(params{:});
if true
  % more accurate for large n
 f = dbinom(x,n,p,options);
else
try
  x(x<0 | x>n) = nan;
  p(p<0 | p>1) = nan;

  %f = binom(n,x) .* p.^x.* (1-p).^(n-x);
  uselog = true;
  logf = binom(n,x,uselog)+x.*log(p) + (n - x).*log1p(-p);
  %logf = gammaln(n + 1) - gammaln(x + 1) - gammaln(n - x + 1) ...
  %  +x.*log(p) + (n - x).*log1p(-p);
catch
  error('x, n and p must be of common size or scalar.');
end

if options.logp
  f = logf;
else
  f = exp(logf);
end
end
function y = dbinom(x,n,p,options)

% Reference
% Catherine Loader (2000). 
% "Fast and Accurate Computation of Binomial Probabilities"; 
% http://www.herine.net/stat/software/dbinom.html.
% @misc{ july-fast,
%   author = "Catherine Loader July",
%   title = "Fast and Accurate Computation of Binomial Probabilities",
%   url = "citeseer.ist.psu.edu/312695.html" }

% Translated from c to matlab pab 2007


try
  p(p<0 | 1<p) = nan;
  x(x<0 | n<x) = nan;
  y = (x==0).*exp(n.*log1p(-p)) +(x==n).*exp(n.*log(p));
  if any(p==0)
    if isscalar(p) && isscalar(x)
        if p==0 && x ==0
          y(:) = 1;
        else
          y(:) = 0;
        end
    else
      y(p==0 & x==0) = 1;
      y(p==0 & x~=0) = 0;
    end
  end
  if any(p==1)
    y(p==1 & x==n) = 1;
    y(p==1 & x~=n) = 0;
  end
  if options.logp
    y = log(y);
  end
catch
  error('x, n and p must be of common size or scalar.')
end


% if (p==0.0) return( (x==0) ? 1.0 : 0.0);
% if (p==1.0) return( (x==n) ? 1.0 : 0.0);
% if (x==0) return(exp(n.*log1p(-p)));
% if (x==n) return(exp(n.*log(p)));
  
inside = (0<p & p<1 & 0<x & x < n);
if any(inside)
  if ~isscalar(n), n = n(inside);end
  if ~isscalar(x), x = x(inside);end
  if ~isscalar(p), p = p(inside);end
  lc = stirlerr(n) - stirlerr(x) - stirlerr(n-x) ...
    - bd0(x,n.*p) - bd0(n-x,n.*(1-p));
  PI2 = 6.283185307179586476925286;
  if options.logp
    y(inside) = lc+0.5.*log(n./(PI2*x.*(n-x)));
  else
    y(inside) = (exp(lc).*sqrt(n./(PI2*x.*(n-x))));
  end
end



