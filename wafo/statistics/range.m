function [r,xmax,xmin] = range(x)
%RANGE  Computes the range between the maximum and minimum values. 
% 
% CALL: [r,xmax,xmin] = range(x) 
%
%   r =  max(x) - min(x), i.e. a scalar or vector giving the range of x
%   x = vector or matrix
%
% Example:
%   x=rndexp(5,1,10)
%   [r,xmax,xmin] = range(x);
%
% See also  max, min, iqrange

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


% By pab 06.02.2000
xmin = min(x);
xmax = max(x);
r = xmax - xmin;
