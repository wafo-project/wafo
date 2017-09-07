function R = rndgengammod(varargin)
%RNDGENGAMMOD Random matrices from a Generalized Modified Gamma distribution.
%
% CALL:  R = rndgengammod(m,v,L,sz);
%
%         R = matrix of random numbers
%    a,b,c  = parameters (see pdfgengam)
%      phat = Distribution parameter struct
%             as returned from FITGENGAMMOD.  
%        sz = size(R)    (Default common size of a and b)
%             sz can be a comma separated list or a vector 
%             giving the size of R (see zeros for options). 
% Example:
%   R=rndgengammod(1,2,4,1,100);
%   plotweib(R);
%
%   close all;
%
% See also pdfgengammod, cdfgengammod, invgengammod, fitgengammod, momgengammod 

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
% Cohen & Whittle, (1988) "Parameter Estimation in Reliability
% and Life Span Models", p. 220 ff, Marcel Dekker.
%
% http://www.weibull.com/LifeDataWeb/generalized_gamma_distribution.htm

% Tested on; Matlab 5.3
% History: 
% By pab 2007

%error(nargchk(1,inf,nargin))
narginchk(1,inf)
Np = 3;
options = []; % default options
[params,options,rndsize] = parsestatsinput(Np,options,varargin{:});
if numel(options)>1
  error('Multidimensional struct of distribution parameter not allowed!')
end

[m,v,L] = deal(params{:});

if isempty(m), m=0;end
if isempty(v), v=1;end
if isempty(L), L=0;end



if isempty(rndsize)
  [csize] = comnsize(m,v,L);
else
  [csize] = comnsize(m,v,L,zeros(rndsize{:}));
end
if any(isnan(csize))
  error('m, v and L must be of common size or scalar.');
end

a = 1./L.^2;
c = abs(L)./sqrt(v);
b = exp(m-log(a)./c);
epsilon = 2^-13; % defining L to zero
if isscalar(a) && isscalar(b) && isscalar(c)
  if abs(L)<epsilon
    R = rndlognorm(m,v,rndsize{:});
  else
    R = rndgengam(a,b,c,rndsize{:});
  end
else
  if isscalar(m), m = m(ones(csize));end
  if isscalar(v), v = v(ones(csize));end  
  if isscalar(L), L = L(ones(csize));end  
  if isscalar(a), a = a(ones(csize));end  
  if isscalar(b), b = b(ones(csize));end  
  if isscalar(c), c = c(ones(csize));end  
  R = zeros(csize);
  k1 = find(abs(L)<=epsilon);
  if any(k1)
    R(k1) = rndlognorm(m(k1),v(k1));
  end
  k1 = find(abs(L)>epsilon);
  if any(k1)
    R(k1) = rndgengam(a(k1),b(k1),c(k1));
  end
end


