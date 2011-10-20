function z = logit(p)
%LOGIT    Logit function.
%
%	  z = lodds(p)
%
%	  The logit function is is log(p/(1-p)) and is an important part of
%	  logistic regression

%       Anders Holtsberg, 14-12-94
%       Copyright (c) Anders Holtsberg
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



if any(any(abs(2*p-1)>=1))
   error('A probability input please')
end
z = log(p./(1-p));
