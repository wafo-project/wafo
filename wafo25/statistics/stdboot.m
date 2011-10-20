function [s,y] = stdboot(varargin)
%STDBOOT  Bootstrap estimate of the standard deviation of a parameter.
%
%	  s = stdboot(X,fun,B)
%	  
%         The function is equal to sqrt(diag(covboot(X,fun)))
%
%	  See also std, stdjack, covboot, and ciboot.

%       Anders Holtsberg, 02-03-95
%       Copyright (c) Anders Holtsberg
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




[C,y] = feval(@covboot,varargin{:});
s = sqrt(diag(C));
