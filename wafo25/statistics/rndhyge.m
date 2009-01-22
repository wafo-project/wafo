function  X = rndhyge(varargin)
%RNDHYGE Random numbers from the Hypergeometric distribution
%
% CALL X = rndhyge(n,K,N,sz)
%
%   x =  Number of items belonging to class K (x = 0,1,2,...min(n,K))
%   n = sample size without replacement (n<=N)
%   K = total population of item K      (K<=N)
%   N = total population
%
% Example
%  sz = [100,1];
%  R = rndhyge(20,39,100,sz)
% 
% See also pdfhyge, cdfhyge, invhyge, fithyge, momhyge

%       Anders Holtsberg, 18-11-93
%       Copyright (c) Anders Holtsberg

error(nargchk(3,inf,nargin))
Np = 3;
%options = struct; % default options
options = [];
[params,options,rndsize] = parsestatsinput(Np,options,varargin{:});
% if numel(options)>1
%   error('Multidimensional struct of distribution parameter not allowed!')
% end

[n,K,N] = deal(params{:});

if isempty(rndsize)
  csize = comnsize(n,K,N);
else
  csize = comnsize(n,K,N,zeros(rndsize{:}));
end
if any(isnan(csize))
    error('n, K and N must be of common size or scalar.');
end

X = invhyge(rand(csize),n,K,N);


