function [f,k,s] = fitgenparml(x,data)
%FITGENPARML Internal routine for fitgenpar (ML estimates for GPD data) 
%
% CALL:  [f,k,s] = fitgenpar_ml(x,data)
%
% f = function values.
% k = shape parameter of GPD.
% s = scale parameter of GPD.
%
% This function is used by fitgenpar for numerical solution of 
% the ML estimate, i.e. solve f=0 for x.
%   data = rndgenpar(0.3,1,0,200,1);
%   x_ML = fzero('fitgenpar_ml',0,[],data);
%   [f,k_ML,s_ML] = fitgenpar_ml(x_ML,data)  % Estimates k_ML and s_ML
%
% See also  fitgenpar

% References
%
%  Davidson & Smith (1990)
%  Models for Exceedances over high Threholds.
%  Journal of the Royal Statistical Society B,52, pp. 393-442.

% Tested on; Matlab 5.3
% History: 
% Created by PJ 22-Jun-2000
% Revised by PJ 10-Oct-2000
% - Help text added w*

% In order to avoid boundary problems in numerical solution we use a transformation
%   Transformation: x = log(1/max_data - t),   -Inf < t < 1/max_data
%   Inverse Trans.: t = 1/max(data) - exp(x),  -Inf < x < Inf

t = 1/max(data) - exp(x); % Inverse Transformation

N = length(data);

k = zeros(size(x));
f1 = zeros(size(x));
for ix =1:numel(x)
  k(ix) = -sum(log1p(-t(ix)*data))/N; % Shape parameter
  % Evaluate function
  f1(ix) = sum(data./(1-t(ix)*data));
end

s = k./t;                     % Scale parameter
f = (1./k-1).*f1-N./t;

i0 = find(t==0);
if any(i0)
  s(i0) = f1/N;
  f(i0) = 0;
end
  
