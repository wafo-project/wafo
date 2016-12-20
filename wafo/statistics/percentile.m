function  q = percentile(x,p,method)
%PERCENTILE Empirical quantile (percentile).
%
%  CALL:  q = percentile(x,p,method)
%
%      q = empirical quantile
%      x = input vector or matrix
%      p = probability
% method = 1 Interpolation so that F(X_(k)) == (k-0.5)/n. (default)
%          2 Interpolation so that F(X_(k)) == k/(n+1).
%          3 Based on the empirical distribution.
%
%  If input  x  is a matrix then the quantile is computed for 
%  every column. Input  p  may be vector also. It even 
%  accepts  x  being a vector and  p  a matrix!
%
% Example
%  method = 3;
%  x = rndgumb(2,2,1,100);
%  q  = percentile(x,[0.25 0.5 0.75],method); % 25% 50% and 75% quantile
% 
% See also plotqq

% References: 
%  Holtsberg, Anders (1999)
%  Stixbox. A statistics toolbox for Matlab and Octave. 
%  Lund University
%  http://www.maths.lth.se/matstat/stixbox

% Tested on: Matlab 5.3
% History:
% revised pab 2007
% -added example
% revised pab 
% - added nargchk
% - updated help header to conform to wafo style
% by Anders Holtsberg 1994, 1998

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


error(nargchk(2,3,nargin))
if nargin<3||isempty(method), method=1; end
if min(size(x)) == 1
   x = x(:);
   q = zeros(size(p));
else
   if min(size(p)) > 1 
      error('Not both matrix x and matrix p input')
   end
   q = zeros(length(p),size(x,2));
end
if any(any((p>1|p<0)))
   error('Input p is not probability')
end

x = sort(x); 
p = p(:);
n = size(x,1);
if method == 3
   qq1 = x(ceil(max(1,p*n)),:); 
   qq2 = x(floor(min(p*n+1,n)),:);
   qq = (qq1+qq2)/2;
else                         
   x = [x(1,:); x; x(n,:)];
   if method == 2
      % This method is from Hjort's "Computer
      % intensive statistical methods" page 102
      i = p*(n+1)+1;
   else % Metod 1
      i = p*n+1.5;
   end
   iu = ceil(i);
   il = floor(i);
   d = (i-il)*ones(1,size(x,2));
   qq = x(il,:).*(1-d)+x(iu,:).*d;
end

q(:) = qq;







