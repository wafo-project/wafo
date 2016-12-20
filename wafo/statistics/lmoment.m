function [lmom,err] = lmoment(Q,m)
%LMOMENT L-moment based on order statistics.
%
% CALL  Lmom = lmoment(data,m)
%       Lmom = lmoment(Q,m)
% 
%   Lmom = [L1,L2,T3,T4,...Tm] vector of L-moments , size 1, m
%   data = vector of data
%   Q    = quantile function given as function handle, string or inline. 
%   m    = max order of moments
%
%  LMOMENT estimates the L-moments. For m>2 the scaling invariant measures
%  of moments, Ti = Li/L2 are calculated.
%  Analogously to conventional moments L1, L2, T3 and T4 measures the
%  location, scaling, skewness and kurtosis, respectively. Especially, the
%  L-skewness (T3) and L-kurtosis (T4) are found useful in several
%  applications because they are more reliable than the moments-based
%  skewness and kurtosis. L-moments have lower sample variances and are
%  more robust against outliers compared to conventional moments.
%  Another benefit is that the distribution only needs finite mean for the
%  L-moments to exist for all orders.
%  
%
% Example
%  N = 500;
%  R = rndray(1,N,1);
%  lmom = lmoment(R,4);                   % Estimated 
%  tlmom = lmoment(@(x)invray(x,1),4);    % True values
%  assert(lmom, tlmom, 0.3);
%  emom = [mean(R),var(R),skew(R),kurt(R)]; % Estimated
%  [mom{1:4}] = momray(1);                 % True values
%  tmom = [mom{:}]; 
%  assert(lmom, tlmom, 0.3);
%
% See also mean, std, skew, kurt

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
% Juha Karvanen (2006)
% Estimation of quantile mixtures via L-moments and Trimmed L-moments.
% Computational statistics and Data analysis, 51, 947-959.

% History
% by pab 2007

if nargin<2 
  m = 2;
else
  m = max(m,2);
end
lmom = zeros(1,m);

if isnumeric(Q)
  x = sort(Q(:)).';
  n = numel(x);
  
  lmom(1) = mean(x);
  if m>1
    temp = ((1-n):2:(n-1));
    p = temp./n;
    lmom(2) = mean(x.*p);
    if m>2
       p1   = ones(1,n);
       tmp = p;
    end
  end

  for i = 3:m
    p2 = p1;
    p1 = p;
    p = ( (2*i-3)*tmp.*p1 - (i-2) * p2 ) /((i-1));
    %p = ( (2*i-3)*temp.*p1 - (i-2) * (n + i - 2)  * p2 ) /((i-1) * (n - i+1));
    lmom(i) = mean(x.*p)/lmom(2); 
  end
else
  err = lmom;
  [lmom(1) err(1)] = gaussq(Q,0,1);
  p = [2 -1];
  Li =   @(x) feval(Q,x).*polyval(p,x);
  [lmom(2), err(2)] = gaussq(Li,0,1);
  reltol = 1e-5;
  temp = [2 -1 ];
  p1 = 1;
   for i = 3:m
    p2 = p1;
    p1 = p;
    p = polysub((2*i-3)*polymul(temp,p1), (i-2)* p2 ) /(i-1) ;
    Li =   @(x) feval(Q,x).*polyval(p,x);
     [lmom(i), err(i)] = gaussq(Li,0,1,reltol);
   end
  lmom(3:end) = lmom(3:end)./lmom(2);
  err(3:end)  = err(3:end)./lmom(2);
end