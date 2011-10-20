function y = bd0(x,np)
% BD0 Evaluate deviance term x*log(x/np) + np - x
%
% CALL y = bd0(x,np)
%
% 
% See also stirlerr, pdfbin, pdfpois

% Reference
% Catherine Loader (2000). 
% "Fast and Accurate Computation of Binomial Probabilities"; 
% http://www.herine.net/stat/software/dbinom.html.% @misc{ july-fast,
%   author = "Catherine Loader July",
%   title = "Fast and Accurate Computation of Binomial Probabilities",
%   url = "citeseer.ist.psu.edu/312695.html" }

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


try
  y = x.*log(x./np)+np-x;
catch
  error('x and np must be of common size or scalar')
end

r1 = abs(x-np)<0.1*(x+np);
if any(r1)
  if ~isscalar(x),x = x(r1);end
  if ~isscalar(np),np = np(r1);end
 
  v = (x-np)./(x+np);
  
  s1 = (x-np).*v;
  s = s1;
  ej = 2*x.*v;
  %v2 = v.*v;
  v = v.*v;
  j = 0;
  ix = 1:numel(s1);
  while any(ix)
    j = j+1;
    s(ix)  = s1(ix);
    ej(ix) = ej(ix).*v(ix);
    s1(ix) = s(ix)+ej(ix)./(2*j+1);
    ix = find(s1~=s);
  end
  y(r1) = s1;
end
