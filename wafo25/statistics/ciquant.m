function ci = ciquant(x,p,C)
%CIQUANT  Nonparametric confidence interval for quantile
%
%	  ci = ciquant(x,p,C);
%
%         Input C is confidence level for the interval, with default 
%	  0.95, p is the probability. The interval is of the form
%         [LeftLimit, PointEstimate, RightLimit]'. The interval is
%	  constructed conservatively in both ends, that is 
%	  P(q<LeftLimit) <= C/2 and similarly for the upper limit.
%	  If x is a matrix the procedure is colonwise.
%
%	  See also percentile

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




if size(x,1) == 1, x = x'; end
[n,m] = size(x);
x = [-Inf*ones(1,m); sort(x); Inf*ones(1,m)];
pr = cdfbin(0:n-1,n,p);
a = (1-C)/2;
J = pr >= a & pr <= (1-a);
L = x(find(J==1, 1 ),:);
H = x(find(J==1, 1, 'last' )+2,:);
ci = [L;percentile(x,p,2);H];
