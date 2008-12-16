function [fit,res,sd,dof] = fitmodel(y,model)
%FITMODEL  Fits response by polynomial
%
% CALL:  [fit,res,sd,dof ] = fitmodel(y,model)
%
%  fit   = fitted response
%  res   = residual, i.e., y-fit
%  sd    = standard deviation of residual
%  dof   = degrees of freedom
%  y     = Response in standard order
%  model = character array of model parameters
%
% Example
%   D = ffd(3);                    % complete 2^3 design in standard order.
%   y = [60 72 54 68 52 83 45 80]; % Responses to design D.
%   [ef, id] = yates(y);           % Calculate effects
%   nplot(ef,id)                   % Identify model
%   model = strvcat('A','B','AC'); % model parameters
%   [fit,res,sd,dof] = fitmodel(y,model);
%   plotnorm(res)                 % Diagnostic check on fitted model.
% 
% See also  nplot, yates

error(nargchk(2,2,nargin))
sz = size(y);
n  = length(y); 
if prod(sz) == n, 
  y = y(:);       % Make sure it is a column vector
else   
  n = sz(1);      % Number of runs
end
if isnumeric(model)
  model = cnr2cl(model); % Transform into columnlabels
end
model = sort(model,2, 'descend' );
p = size(model,1);

[ef,id] = yates(y); % Calculate the effects
id  = sort(id,2, 'descend' );
ind = ones(n,1);
for ix=1:p
  k = strmatch(model(ix,:),id,'exact');
  if any(k),
    ind(k+1)=0;
  else
    warning('WAFO:FITMODEL','Something wrong!')
  end
end

ef2      = ef;
ind(1)   = 0;
%ind      = find(ind);
ef2(ind~=0,:) = 0;       % Neglect effects from variables not in the model

fit = ryates(ef2);  % Calculate the fit
res = y-fit;        % Residual
 
sd  = std(res);
dof = n-p-1;




