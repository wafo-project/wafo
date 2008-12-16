function  x = invbin(F,varargin)
%INVBIN Inverse of the Binomial distribution function.
%
% CALL  x = invbin(F,n,p)
%       x = invbin(F,phat)
%
%        F = probability of observing x or less successes in n independent
%            trials.
%        x = number of successes          (0<=x<=n)
%        n = number of independent trials
%        p = probability of succes in a given trial. (0<=p<=1)
%     phat = Distribution parameter struct
%            as returned from WBINOMFIT.  
%  options = struct with fieldnames:
%         .lowertail: if TRUE (default), F = Prob[X <= x],
%                     otherwise, F = Prob[X > x].
%         .logp     : if TRUE, probability, p, returned as log(p).
%         .alpha    : Confidence coefficent    (default 0.05)
%
% Example
%  n = 10; p = 0.05;
%  xo = -1:n+1
%  F = cdfbin(xo,n,p);
%  x = invbin(F,n,p);
%  plot(abs((x-xo))), shg
%
% See also pdfbin, cdfbin, rndbin, invbin, fitbin, mombin


%       Anders Holtsberg, 16-03-95
%       Copyright (c) Anders Holtsberg

% The algorithm contains a nice vectorization trick which
% relies on the fact that if two elements in a vector
% are exactely the same then matlab's routine SORT sorts them
% into the order they had. Do not change this, Mathworks!

error(nargchk(2,9,nargin))
options = struct('covariance',[],'alpha',0.05,...
  'lowertail',true,'logp',false); % default options
Np = 2;
[params,options] = parsestatsinput(Np,options,varargin{:});
if numel(options)>1
  error('Multidimensional struct of distribution parameter not allowed!')
end

[n,p] = deal(params{:});

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
  if numel(p)>1 || numel(n)>1
    [csize,F,n,p] = comnsize(F,n,p);
    [np,ii,jj] = unique([n(:),p(:)],'rows');
    idx = 1:size(np,1);
    grp = idx(jj);
    x = zeros(csize);
    for ix = idx
      x(grp==ix) = invbin1(F(grp==ix),np(ix,1),np(ix,2));
    end
  else
    x = invbin1(F,n,p);
  end
catch
  error('F, n and p must be of common size or scalar.');
end


function x = invbin1(F,n,p)
if max([numel(n) numel(p)]) > 1
   error('Sorry, this is not implemented');
end

x = F;
if any(n < 0 | p < 0 | p > 1 );
 x(:) = nan;  
else
  F(F<0 | 1 <F) = nan;

  kk = (0:n)';
  cdf = max(0,min(1,cumsum(pdfbin(kk,n,p))));
  cdf(n+1) = 1;
  [pp,J] = sort(F(:));
  np = length(pp);
  [S,I] = sort([pp;cdf]);
  I = find(I<=np) - (1:np)';
  J(J) = (1:np)';
  x(:) = I(J);
  x(x>n) = nan;
end
