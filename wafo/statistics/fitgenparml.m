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
%
% Example
%   data = rndgenpar(0.3,1,0,200,1);
%   x_ML = fzero(@(p)fitgenpar_ml(p, data),0);
%   [f,k_ML,s_ML] = fitgenparml(x_ML,data);  % Estimates k_ML and s_ML
%
% See also  fitgenpar

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
  
