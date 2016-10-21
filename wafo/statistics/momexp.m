function [m,v,sk,ku]= momexp(varargin)
%MOMEXP Mean and variance for the Exponential distribution.
% 
% CALL:  [m,v,sk,ku] = momexp(m0)
%
%   m, v = the mean and variance, respectively 
%  sk,ku = the skewness and kurtosis, respectively. 
%     m0 = parameter of the Exponential distribution, m0>0.
%
%  Mean (m) and variance (v) for the Exponential distribution is
%
%  m=m0  and  v=m0^2;
%
% Example:
%   par = {4}
%   X = rndexp(par{:},1000,1);
%   [mean(X) var(X),skew(X),kurt(X)]         % Estimated mean and variance
%   [m,v,sk,ku] = momexp(par{:}) % True mean and variance
%
% See also pdfexp, cdfexp, invexp, rndexp, fitexp

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


% Reference: Johnson, Kotz and Balakrishnan (1994)
% "Continuous Univariate Distributions, vol. 1", p. 494 ff
% Wiley

% Tested on; Matlab 5.3
% History:
% revised pab Dec2003
% fixed a bug k1->k3
% revised pab 24.10.2000
%  - added  nargchk
% added ms 15.06.2000

error(nargchk(1,1,nargin))
Np = 1;
options = [];
params = parsestatsinput(Np,options,varargin{:});
m = params{1};
m(m<=0) = nan;
v = m.^2;
sk = repmat(2,size(v));
ku = repmat(9,size(v));
if any(isnan(m(:)))
  sk(isnan(m)) = nan;
  ku(isnan(m)) = nan;
end


