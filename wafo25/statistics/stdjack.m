function [s,y] = stdjack(varargin)
%STDJACK  Jackknife estimate of the standard deviation of a parameter 
%	  estimate.
%
%	  s = stdjack(X,fun)
%	  
%	  The function is equal to sqrt(diag(covjack(X,'T')))

%       Anders Holtsberg, 28-02-95
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




[C,y] = feval(@covjack,varargin{:});
s = sqrt(diag(C));