function [f,k,s] = fitgenpar_mld(x,data)
%FITGENPAR_MLD Parameter estimates for GPD data
%
% CALL:  [f,k,s] = wgpdfit_mld(x,data)
%
% f = function values.
% k = shape parameter of GPD.
% s = scale parameter of GPD.
%
% This function is used by cmat2extralc for numerical solution of 
% the ML estimate, i.e. solve f=0 for x.
%   data0 = wgpdrnd(0.3,1,0,200,1);
%   x_ML = fzero('wgpdfit_ml',0,[],data0);
%   [f,k_ML,s_ML] = wgpdfit_ml(x_ML,data0) % Estimates k_ML and s_ML
%   data1 = floor(data0*10)/10;
%   x=(0:0.1:(max(data1)+0.1))';
%   N = histc(data1+0.05,x);
%   x_MLD = fzero('wgpdfit_mld',0,[],[x, N(:)]);
%   [f,k_MLD,s_MLD] = wgpdfit_mld(x_MLD,[x N(:)]) % Estimates k_ML and s_ML
%
% See also fitgenpar, cmat2extralc

% References
%
%  Davidson & Smith (1990)
%  Models for Exceedances over high Threholds.
%  Journal of the Royal Statistical Society B,52, pp. 393-442.

% Tested on; Matlab 5.3
% History: 
% Created by PJ 10-Oct-2000
% - created from wgpdfit_ml
% -revised pab Aug 2007:
% -updated example in help header

% In order to avoid boundary problems in numerical solution we use a transformation
%   Transformation: x = log(1/max_data - t),   -Inf < t < 1/max_data
%   Inverse Trans.: t = 1/max(data) - exp(x),  -Inf < x < Inf

M = max(data(:,1)); % Max of data
t = 1/M - exp(x); % Inverse Transformation

N = sum(data(:,2));

if t<1/M
  k = -1/N*sum(data(:,2).*log(1-t*data(:,1))); % Shape parameter
  s = k/t;                     % Scale parameter
  % Evaluate function
  f = (1/k-1)*sum(data(:,2).*data(:,1)./(1-t*data(:,1))) - N/t; 
else
  I = (data(:,1)==M);
  k = -1/N*(sum(data(~I,2).*log(1-t*data(~I,1)))+sum(data(I,2).*(x+log(data(I,1))))); % Shape parameter
  s = k/t;                     % Scale parameter
  % Evaluate function
  f = (1/k-1)*(sum(data(~I,2).*data(~I,1)./(1-t*data(~I,1)))+sum(data(~I,2)*exp(-x)))- N/t; 
end


% Evaluate function
%f = (1/k-1)*sum(data(:,2).*data(:,1)./(1-t*data(:,1))) - N/t; 

if isinf(f)
  f=realmax*sign(f);
end
if isnan(f)   %if x<0, f=realmax; else, f=-realmax; end
  f=realmax;
end
