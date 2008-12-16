function k = kurt(X,dim)
%KURT Computes sample kurtosis
%
% CALL:  k = kurt(X,dim);
%
%        k = sample kurtosis (fourth central moment divided by squared second)
%        X = data vector or matrix
%      dim = dimension to sum across. (default 1'st non-singleton 
%                                              dimension of X)
%
% Example:  
%   R=rndgumb(2,2,100,2);
%   kurt(R)
%
% See also   mean, varskew, skew

% Tested on: Matlab 5.3
% History:
% revised pab 24.10.2000
% - made it more general: accepts any size of X
% - added dim, nargchk
% added ms 16.06.2000

error(nargchk(1,2,nargin))
sz = size(X);
if nargin<2||isempty(dim)
  % Use 1'st non-singleton dimension or dimension 1
  dim = find(sz~=1, 1 ); 
  if isempty(dim), dim = 1; end
end
rsz = ones(size(sz)); rsz(dim)=sz(dim);
mu  = mean(X,dim);
if isscalar(mu)
    Xmu2 = (X-mu).^2;
else
    Xmu2 = (X-repmat(mu,rsz)).^2;
end
k   = mean(Xmu2.^2,dim)./mean(Xmu2,dim).^2;
end






