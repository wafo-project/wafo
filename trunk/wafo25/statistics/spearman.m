function rho = spearman(x,y)
%SPEARMAN Spearman's rank correlation coefficient.
%
%	  rho = spearman(x,y)
%
%	  This is the correlation coefficient between rank
%	  transformed data. 

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


if nargin < 2
   rho = corrcoef(ranktrf(x));
else
   rho = corrcoef(ranktrf(x),ranktrf(y));
end