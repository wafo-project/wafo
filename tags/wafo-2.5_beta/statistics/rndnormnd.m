function r = rndnormnd(mu,sa,cases,method,cutoff)
%RNDNORMND Random vectors from a multivariate Normal distribution
%
% CALL:  r = rndnormnd(mu,S,n,method)
%
%       r = matrix of random numbers from the multivariate normal
%           distribution with mean mu and covariance matrix S.
%       n = number of sample vectors
%  method = 'svd'  Singular value decomp.  (stable, quite fast) (default)
%           'chol'  Cholesky decomposition (fast, but unstable) 
%           'sqrtm'  sqrtm                 (stable and slow) 
%           'genchol'
%  
%   S must be a symmetric, semi-positive definite matrix with size equal
%   to the length of MU. METHOD used for calculating the square root of S 
%   is either svd, cholesky or sqrtm. (cholesky is fastest but least accurate.)
%   When cholesky is chosen and S is not positive definite, the svd-method 
%   is used instead.
%
% Example
%   mu = [0 5]; S = [1 0.45; 0.45 0.25];
%   r = rndnormnd(mu,S,100)';
%   plot(r(:,1),r(:,2),'.')
%
%   d = 40; rho = 2*rand(1,d)-1;
%   mu = zeros(0,d);
%   S = (rho.'*rho-diag(rho.^2))+eye(d);
%   r = rndnormnd(mu,S,100,'genchol')'; 
%
% See also  chol, svd, sqrtm, genchol

% cutoff = cut off value for eigenvalues

% History
% Revised pab 22.11.2003
% -added option genchol  
% revised jr 19.09.2001
% - Changed *name* in function call. 
% - Fixed a bug: the spurious default option 'swd' -> 'svd'
% revised jr 18.06.2001
% - Changed default from 'cholesky' to 'svd'
% revised pab 17.01.2000 
% - updated documentation
%  by Per A. Brodtkorb 17.09.98,20.08.98

[m n] = size(sa);
if m ~= n
  error('S must be square');
end
if  nargin<1 || isempty(mu);
  mu=zeros(m,1);
end

musiz = size(mu);
rows = max(musiz);
if prod(musiz) ~= rows
   error('Mu must be a vector.');
end
mu=mu(:); % make sure it is a column vector


if m ~= rows
  error('The length of mu must equal the number of rows in S.');
end
if nargin<5||isempty(cutoff), 
  cutoff=0;
else
  cutoff=abs(cutoff); % make sure cutoff is positive
end

if nargin<4||isempty(method), 
  method='svd'; % default method
end

switch lower(method(1:2)), %
  case 'sq', % SQRTM slow but stable
    T=sqrtm(full(sa));%squareroot of S
    
  case 'sv', % SVD stable quite fast
    if issparse(sa) % approximate method
      nz=nnz(sa);
      K=floor(min([10 3*sqrt(n) n/4 4*nz/n]));
      %[U S V]=svds(sa,floor(p/2));
      [U S V]=svds(sa,K,'L');
    else % exact method
      [U S V]=svd(full(sa),0);
    end
    L = n;
    if cutoff>0
      L=sum((diag(S)>=cutoff));
      disp(['cutting off ' num2str(L) ' eigenvalues'])
    end
    
    T=(U(:,1:L)*sqrt(S(1:L,1:L))*V(:,1:L)'); %squareroot of S
    
  case 'ch', % CHOL not stable , but fast.
    [T, p] = chol(sa);
    if p ~= 0
      disp('S is not a positive definite matrix.');
      disp('Cholesky factorization failed; trying svd instead.');
      [U S V]=svd(full(sa),0);
      T=(U*sqrt(S)*V'); %squareroot of S
    end
 case 'ge', % GENCHOL stable
  [T,p,rank1] = genchol(sa,eps);
  r = T*randn(rows,cases)  + mu(p,ones(1,cases));
  [tmp, ix] = sort(p);
  r = r(ix,:);
  return
 otherwise, error('unknown option')
end


r = T*randn(rows,cases)  + mu(:,ones(1,cases));

