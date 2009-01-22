function r = iqrange(X,dim)
%IQRANGE Computes the Inter Quartile Range
%
% CALL:  r = iqrange(X,dim);
%
%        r = abs(diff(percentile(X,[0.25 .75]))),
%        X = data vector or matrix
%      dim = dimension to sum across. (default 1'st non-singleton 
%                                              dimension of X)
% IQRANGE is a robust measure of spread.
% The use of interquartile range guards against outliers if 
% the distribution have heavy tails.
%
% Example:
%   sz = [100,2];
%   R=rndgumb(2,2,sz);
%   iqrange(R)
%
% See also  std

% Tested on: Matlab 5.3
% History:
% revised pab 24.10.2000

error(nargchk(1,2,nargin))
sz = size(X);
if nargin<2||isempty(dim),
  % Use 1'st non-singleton dimension or dimension 1
  dim = find(sz~=1, 1 ); 
  if isempty(dim), dim = 1; end
end

if dim~=1, 
  iorder=1:length(sz);
  tmp=iorder(dim);
  iorder(dim)=iorder(1);
  iorder(1)=tmp;
  X = permute(X,iorder);
end
r = abs(diff(percentile(X,[0.25 0.75])));

if dim~=1, 
  iorder=1:length(sz);
  tmp=iorder(dim);
  iorder(dim)=iorder(1);
  iorder(1)=tmp;
  r=ipermute(r,iorder);
end


%sz(dim)=sz(dim);



