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
%   R = rndgumb(2,2,100,2);
%   sk = skew(R);
%   
%   assert(skew(1:10), 0)
%
% See also  mean, var, kurt

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
    Xmu  = (X-mu);
else
    Xmu  = (X-repmat(mu,rsz)); 
end

s   = mean(Xmu.^3,dim)./mean(Xmu.^2,dim).^(1.5);
end


