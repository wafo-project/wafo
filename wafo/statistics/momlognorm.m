function [m,v,sk,ku]= momlognorm(m0,v0)
%MOMLOGNORM Mean and variance for the Lognormal distribution.
% 
% CALL:  [m,v,sk,ku] = momlognorm(m0,v0)
%
%   m, v = the mean and variance, respectively 
%  sk,ku = the skewness and kurtosis, respectively. 
% m0, v0 = parameters of the Lognormal distribution.
%
%  Mean (m) and variance (v) for the Lognormal distribution is
%
%  m=exp(m0+v0/2)  and  v=exp(2*m0+v0)*(exp(v0)-1);
%
% Example:
%   par = {-1, 0.1};
%   X = rndlognorm(par{:},10000,1);
%   moments = {mean(X) var(X),skew(X),kurt(X)};   % Estimated mean and variance
%   [mom{1:4}] = momlognorm(par{:}); % True mean and variance
%   assert(moments, mom, -0.25);
%
% See also pdflognorm, cdflognorm, invlognorm, rndlognorm, fitlognorm

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


% Reference: Cohen & Whittle, (1988) "Parameter Estimation in Reliability
% and Life Span Models", p. 59 ff, Marcel Dekker.



% Tested on; Matlab 5.3
% History:
% revised pab 24.10.2000
%  - added comnsize, nargchk
%  - added default values
% added ms 10.08.2000

error(nargchk(0,2,nargin))
if nargin<1||isempty(m0),  m0=0;  end
if nargin<2||isempty(v0),  v0=1;  end

v0(v0<0) = nan;
try
  m = exp(m0+v0/2);
  v = exp(2*m0+v0).*(expm1(v0));
  sk = (exp(v0)+2).*sqrt(expm1(v0));
  ku = exp(4*v0)+2*exp(3*v0)+3*exp(2*v0)-3;
catch
  error ('m and v must be of common size or scalar');
end

