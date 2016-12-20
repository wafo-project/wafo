function [m,v,sk,ku] = mompois(L)
%MOMPOIS Mean and variance for the Poisson distribution.
% 
% CALL:  [m,v,sk,ku] = mompois(L)
%
%   m, v = the mean and variance, respectively 
%  sk,ku = the skewness and kurtosis, respectively. 
%      L = parameter, L>0
%
%  Mean (m) and variance (v) for the Poisson distribution is
%
%  m=L  and  v=L;
%
% Example:
%   par = {10};
%   X = rndpois(par{:},100000,1);
%   moments = {mean(X) var(X),skew(X),kurt(X)};   % Estimated mean and variance
%   [mom{1:4}] = mompois(par{:}); % True mean and variance
%   assert(moments, mom, -0.25);
%
% See also pdfpois, cdfpois, invpois, rndpois, fitpois

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
% by pab Nov2007

error(nargchk(1,1,nargin))

L(L<=0) = nan;


m = L;
v = L;
sk = 1./sqrt(L);
ku = 3+1./L;
