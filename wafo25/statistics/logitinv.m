function p = logitinv(z,b)
%LOGITINV Inverse logit function.
%
%	  p = logitinv(z)
%
%	  The function is equal to 1/(1+exp(-z)) if one argument is given.
%	  If lodds(X,b) is given the same thing is computed for z = X*b.
%	  If length(b) is one element too much then z = X*b(1:n-1) + b(n)
%	  is assumed.

%       Anders Holtsberg, 14-12-94
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



if nargin>1
   X = z;
   n = length(b);
   if size(X,2) < n
      z = b(n) + X*b(1:n-1);
   else
      z = X*b;
   end
end
p = 1 ./ (1+exp(-z));
