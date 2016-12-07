function  F = cdfhyge(k,varargin)
%CDFHYGE  The hypergeometric cumulative probability function
%
% CALL  F = cdfhyge(x,n,K,N)
%
%       F = cumulative probability
%       x = Number of items belonging to class K (x = 0,1,2,...min(n,K))
%       n = sample size without replacement (n<=N)
%       K = total population of item K      (K<=N)
%       N = total population
%
%   The Hypergeometric distribution is defined by its pdf
% 
%   f(x;n,K,N) = K!/(K-x)!/x!*(N-K)!/(N-K-n+x)!/(n-x)! (N-n)!*n!/N!
%
% Example
%  n = 10; K = 30; N = 100;
%  x = 0:min(n,K);
%  f = cdfhyge(x,n,K,N);
%  bar(x,f);
%
%  close all;
%
% See also pdfhyge,  invhyge, rndhyge, whypgfit, momhyge

%       Anders Holtsberg, 18-11-93
%       Copyright (c) Anders Holtsberg
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


error(nargchk(2,9,nargin))
options = struct('covariance',[],'alpha',0.05,...
  'lowertail',true,'logp',false); % default options
Np = 3;
[params,options] = parsestatsinput(Np,options,varargin{:});
if numel(options)>1
  error('Multidimensional struct of distribution parameter not allowed!')
end

[n,K,N] = deal(params{:});
try
  if max([numel(n) numel(K) numel(N)]) > 1
   
    [csize,k,n,K,N] = comnsize(k,n,K,N);
    [nKN,ii,jj] = unique([n(:),K(:),N(:)],'rows');
    idx = 1:size(nKN,1);
    grp = idx(jj);
    F = zeros(csize);
    for ix = idx
      F(grp==ix) = cdfhyge1(k(grp==ix),nKN(ix,1),nKN(ix,2),nKN(ix,3),options);
    end
  else
    F = cdfhyge1(k,n,K,N,options);
  end
catch
  error('F, n and p must be of common size or scalar.');
end




function F = cdfhyge1(k,n,K,N,options)
if max([numel(n) numel(K) numel(N)]) > 1
   error('Sorry, this is not implemented');
end

kk = (0:n)';
cdf = max(0,min(1,[0; cumsum(pdfhyge(kk,n,K,N))]));
cdf(n+2) = 1;
F = k;
F(:) = cdf(max(1,min(n+2,floor(k(:))+2)));

if options.logp
  if options.lowertail
    F = log(F);
  else
    F = log1p(-F);
  end
elseif ~options.lowertail
  F = 1-F;
end
