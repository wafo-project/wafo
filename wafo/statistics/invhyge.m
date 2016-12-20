function  k = invhyge(F,varargin)
%INVHYGE Inverse of the Hypergeometric distribution function.
%
%        k = invhyge(F,n,K,N)
%
%        Gives the smallest integer k so that P(X <= k) >= F.
%
% Example
%  p = [-1 0.5, 2];
%  n = 10;K = 30; N = 100;
%  k = invhyge(p,n,K,N);
%  assert(k, [nan, 3, nan]);
%
% See also pdfhyge, cdfhyge, rndhyge, whypgfit, momhyge


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


% The algorithm contains a nice vectorization trick which
% relies on the fact that if two elements in a vector
% are exactely the same then matlab's routine SORT sorts them
% into the order they had. Do not change this, Mathworks!
error(nargchk(4,9,nargin))
options = struct('covariance',[],'alpha',0.05,...
  'lowertail',true,'logp',false); % default options
Np = 3;
[params,options] = parsestatsinput(Np,options,varargin{:});
if numel(options)>1
  error('Multidimensional struct of distribution parameter not allowed!')
end

[n,K,N] = deal(params{:});

if options.logp
  if options.lowertail
    F = exp(F);
  else
    F = -expm1(F);
  end
elseif ~options.lowertail
  F = 1-F;
end

try
   if max([numel(n) numel(K) numel(N)]) > 1
    [csize,F,n,K,N] = comnsize(F,n,K,N);
    [nKN,ii,jj] = unique([n(:),K(:),N(:)],'rows');
    idx = 1:size(nKN,1);
    grp = idx(jj);
    k = zeros(csize);
    for ix = idx
      k(grp==ix) = invhyge1(F(grp==ix),nKN(ix,1),nKN(ix,2),nKN(ix,3));
    end
  else
    k = invhyge1(F,n,K,N);
  end
catch
  error('F, n, K and N must be of common size or scalar.');
end

function k  = invhyge1(p,n,K,N)
if max([numel(n) numel(K) numel(N)]) > 1
   error('Sorry, this is not implemented');
end
% if any(any(abs(2*p-1)>1))
%    error('A probability should be 0<=p<=1, please!')
% end
k = p;
if any(n < 1 | K < 1 | N < 1 );
 k(:) = nan;  
else
  p(p<0 | 1 <p) = nan;
  
  lowerlim = max(0,n-(N-K));
  upperlim = min(n,K);
  kk = (lowerlim:upperlim)';
  %nk = length(kk);
  cdf = max(0,min(1,cumsum(pdfhyge(kk,n,K,N))));
  cdf(length(cdf)) = 1;
  [pp,j] = sort(p(:));
  np = length(pp);
  [S,i] = sort([pp;cdf]);
  i = find(i<=np) - (1:np)' + lowerlim;
  j(j) = (1:np)';
  p(:) = i(j);
  k = p;
  k(k>n) = nan;
end
