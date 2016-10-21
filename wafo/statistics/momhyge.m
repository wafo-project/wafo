function [m,v,sk,ku]= momhyge(n,K,N)
%MOMHYGE Mean and variance for the Hypergeometric distribution.
% 
% CALL:  [m,v,sk,ku] = momhyge(n,K,N)
%
%   m, v = the mean and variance, respectively 
%  sk,ku = the skewness and kurtosis, respectively. 
%      n = sample size without replacement (n<=N)
%      K = total population of item K      (K<=N)
%      N = total population
%
%  Mean (m) and variance (v) for the Binomial distribution is
%
%  m=n*K/N  and  v=n*K/N*(1-K/N)*(N-n)/(N-1));
%
% Example:
%   par = {10,30,100}
%   X = rndhyge(par{:},1000,1);
%   [mean(X) var(X),skew(X),kurt(X)]        % Estimated mean and variance
%   [m,v,sk,ku] = momhyge(par{:}) % True mean and variance
%
% See also pdfhyge, cdfhyge, invhyge, rndhyge, fithyge

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


error(nargchk(3,3,nargin))

try

  pk = K/N;
  m = n*pk;
  v = n*pk*(1-pk)*(N-n)/(N-1);
  sk = ((N-2*K).*sqrt(N-1).*(N-2*n))./(sqrt(n.*K.*(N-K).*(N-n)).*(N-2));
  ku = 3+	(N.^2.*(N-1))./(n.*(N-2).*(N-3).*(N-n)).*...
((N.*(N+1)-6.*N.*(N-n))./(K.*(N-K)) +(3*n.*(N-n).*(N+6))./(N.^2)-6);


catch
  error('n, K and N must be of common size or scalar.');
end
