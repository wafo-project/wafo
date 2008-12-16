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

error(nargchk(3,4,nargin))

treshLow = -abs(treshLow);
if nargin<4||isempty(treshUp), treshUp = -treshLow; end

ind = find( ef(2:end)< treshLow  | treshUp<ef(2:end));
model = [];
if any(ind)
  model = id(ind,:);
end
return
