function z = rndboot(x,B)
%RNDBOOT    Simulate a bootstrap resample from a sample.
%
%	  CALL: Z = rboot(X,B)
%
%	  Give a resample of the same size as X, which is assumed
%	  to have one independent realisation in every row.
%	  RESAMPLE(X,B) gives B columns of resamples. This form works
%	  only for X one-dimensional, ie X column vector.
%
% Example
%
% See also ciboot

%       Anders Holtsberg, 14-12-94
%       Copyright (c) Anders Holtsberg
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



if min(size(x)) == 1
   x = x(:);
end

if nargin > 1 && size(x,2) > 1
   if B > 1, error('X multidimensional and B > 2'), end
elseif nargin < 2
   B = 1;
end

n = size(x,1);
nn = n*B;
I = ceil(n*rand(nn,1));
if B > 1
   z = zeros(n,B);
   z(:) = x(I);
else
   z = x(I,:);
end
