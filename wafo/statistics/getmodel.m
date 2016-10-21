function [model, ind] =  getmodel(ef,id,treshLow,treshUp)
%GETMODEL Return the model parameters
%
% CALL: model = getmodel(ef,id,treshLow,treshUp)
%
%  model = character array of model parameters
%  ind   = indices to model parameters in id, i.e., model = id(ind,:);  
%     id = identification vector of main and interaction effects.
%     ef = vector of average response, main effects and interaction effects. 
% treshLow = lower treshhold (negative)
% treshUp  = upper treshhold  
%
% Example:
% y = [60 72 54 68 52 83 45 80]; % Responses to design D.
% [ef,id] = yates(y);
% nplot(ef,id) 
% model = getmodel(ef,id,-1); 
%
% See also yates, fitmodel, getmodel, alias, cdr, sudg, ffd ,plotresponse, nplot
 
  
%History 
% revised pab Dec2003
% removed some unused code
% by pab 200?

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


error(nargchk(3,4,nargin))

treshLow = -abs(treshLow);
if nargin<4||isempty(treshUp), treshUp = -treshLow; end

ind = find( ef(2:end)< treshLow  | treshUp<ef(2:end));
model = [];
if any(ind)
  model = id(ind,:);
end
return
