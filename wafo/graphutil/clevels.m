function [CL,c]=clevels(c)
%CLEVELS  Extract the contour levels from the contour matrix
% 
% CALL: [CL, C1] = clevels(C)
%
% CL = [level1 level2, ... ] vector of contour levels (like the one you
%                            give into CONTOUR when you want manual
%                            control of the contour levels).
% C1 = [NaN x1 x2 x3 ... NaN x1 x2 x3 ...;
%       NaN y1 y2 y3 ... NaN y1 y2 y3 ...]
%       contour matrix with levels and pairs set to NaN's.
% C  = [level1 x1 x2 x3 ... level2 x1 x2 x3 ...;
%       pairs1 y1 y2 y3 ... pairs2 y1 y2 y3 ...]
%      contour matrix as described in CONTOURC
%
% Example:
% 
% c = contour(peaks);
% [cl, c1] = clevels(c);
% plot(c1(1,:),c1(2,:));
% cltext(cl)
%
% See also  contourc, fcolorbar

% History
% By pab 2000,2007
% -based on function findzlevel in the old plotspec function.

% Copyright (C) 2000  Per A. Brodtkorb
% 
%  This file, CLEVELS.M, is part of WAFO.
% 
%     CLEVELS is free software; you can redistribute it and/or modify
%     it under the terms of the GNU Lesser General Public License as published by
%     the Free Software Foundation; either version 3 of the License, or
%     (at your option) any later version.
% 
%     CLEVELS is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU Lesser General Public License for more details.
% 
%     You should have received a copy of the GNU Lesser General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
    


error(nargchk(1,1,nargin));
limit = size(c,2);
if limit>0
  i=1;
  j=1;
  Nmin = max(round(sqrt(limit)),5);
  CL = zeros(1,Nmin);

  while (i <= limit)
    CL(j)    = c(1,i);
    npoints  = c(2,i);
    if nargout >1
      c(:,i) = NaN;
    end
    i      = i+npoints+1;
    j      = j+1;
  end
  CL(j:end) = [];
  % sort and remove duplicate levels
  CL = unique(CL);
else
  CL = [];
end



