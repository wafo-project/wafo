function r = ranktrf(x)
%RANKTRF  Rank transformation of data material.
%
%	CALL  r = ranktrf(x)
%
%  Replaces the data vector  x  and gives rank numbers.
%	 In case  x  is a matrix every column is rank transformed.
%  
%  Ties in the data are resolved by using average rank value
%
% Example
%  a = rand(5,1);
%  r = ranktrf(a(:));
%  
%  a2 = a;
%  a2(1) = a(2);
%  r2 = ranktrf(a2(:));
%  
% 
% See also spearman, testmean1r, testmean2r

%      Copyright (c) Anders Holtsberg, 1998
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



tr = 0;
if size(x,1) == 1
   x = x';
   tr = 1;
end

%n = size(x,1);
%v = (1:n)';

[dummy, i]  = sort(x);
[dummy, r1] = sort(i);
[dummy, i]  = sort(flipud(x));
[dummy, r2] = sort(i);
r = (r1+flipud(r2))/2;

if tr
   r = r';
end
