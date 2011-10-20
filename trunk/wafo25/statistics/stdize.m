function x = stdize(x, opt)
%STDIZE   Standardize columns to have mean 0 and standard deviation 1.
%
%	  x = stdize(x)
%	  x = stdize(x, 1)
%
%	  A 1 as second argument gives normalization with N instead 
%	  of N-1 as the one argument version does.

%  Copyright (c) Anders Holtsberg, 1998
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



n = size(x,1);
m = mean(x);
if nargin < 2
   s = std(x);
else
   s = std(x, opt);
end
x = (x - m(ones(n,1),:)) ./ s(ones(n,1),:);
