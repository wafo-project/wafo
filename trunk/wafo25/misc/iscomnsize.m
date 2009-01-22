function [iscsize,varargout] = iscomnsize(varargin)
%ISCOMNSIZE True if all non-scalar arguments are of common size.
%
% CALL:  [iscsize,y1,y2,...,yN] = comnsize(x1,x2,...,xN);
%
%  iscsize   = TRUE if all non-scalar arguments match in size.
%  y1,...,yN = Same as x1,...,xN, except that scalars are transformed to
%              a constant matrix with same size as the other inputs.
%  x1,...,xN = Input arguments.
%
%  ISCOMNSIZE returns TRUE if all non-scalar arguments match in size and 
%  makes sure that the output arguments Y1,Y2,...,YN, are of common size.
%  This is useful for implementing functions where arguments can either
%  be scalars or of common size. 
%
%  NOTE:  If the iscsize is FALSE, then yi = xi.
%
% Examples: 
%   A = rand(4,5);B = 2;C = rand(4,5);
%   [iscsize,A1,B1,C1] = iscomnsize(A,B,C);
%   iscsize = iscomnsize(A,1:2);
%
% See also comnsize

% Copyright (C) 1999, 2007 Per A. Brodtkorb
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

% Tested on: matlab 5.3,7.4
% History:
% by pab 2007

[csize,varargout{1:nargout-1}] = comnsize(varargin{:});
iscsize = ~any(isnan(csize));
