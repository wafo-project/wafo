function s = skew(X,dim)
%SKEW Computes sample skewness
%
% CALL:  k = skew(X,dim);
%
%        k = sample skewness (third central moment divided by second^(3/2))
%        X = data vector or matrix
%      dim = dimension to sum across. (default 1'st non-singleton 
%                                              dimension of X)
%
% Example:
%   R=rndgumb(2,2,100,2);
%   skew(R)
%
% See also  mean, var, kurt

% Tested on: Matlab 5.3
% History:
% revised pab 24.10.2000
% - made it more general: accepts any size of X
% - added dim, nargchk
% added ms 16.06.2000

error(nargchk(1,2,nargin))
sz = size(X);
if nargin<2||isempty(dim),
  % Use 1'st non-singleton dimension or dimension 1
  dim = find(sz~=1, 1 ); 
  if isempty(dim), dim = 1; end
end

rsz = ones(size(sz)); rsz(dim)=sz(dim);
mu  = mean(X,dim);
if isscalar(mu)
    Xmu2  = (X-mu).^2;
else
    Xmu2  = (X-repmat(mu,rsz)).^2; 
end
s   = mean(Xmu2.^1.5,dim)./mean(Xmu2,dim).^(1.5);
end


