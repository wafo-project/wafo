function [m,v,sk,ku]= momchi2(varargin)
%MOMCHI2 Mean and variance for the Chi squared distribution.
% 
% CALL:  [m,v,sk,ku] = momchi2(p)
%
%   m, v = the mean and variance, respectively 
%  sk,ku = the skewness and kurtosis, respectively. 
%      p = degrees of freedom, p=1,2,3,....
%
%  Mean (m) and variance (v) for the Chi squared distribution is
%
%  m=p  and  v=2p;
%
% Example: 
%   par = {5};
%   X = rndchi2(par{:},10000,1);
%   moments = {mean(X) var(X),skew(X),kurt(X)};   % Estimated mean and variance
%   [mom{1:4}] = momchi2(par{:}); % True mean and variance
%   assert(moments, mom, 0.2)
%
% See also pdfchi2, cdfchi2, invchi2, rndchi2, fitchi2

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
% "Continuous Univariate Distributions, vol. 1", p. 415 ff
% Wiley

% Tested on; Matlab 5.3
% History: 
% revised pab 25.10.2000
% - found a bug: forgot multiplication
% - added warning+nan's
% added ms 26.06.2000

Np = 1;
%error(nargchk(1,Np,nargin))
narginchk(1,Np)
options = []; %struct; % default options
params = parsestatsinput(Np,options,varargin{:});
p = params{1};
p(p~=round(p) | p<=0) = nan;
m =  p;
v = 2*p;
if any(isnan(p)),
  warning('WAFO:MOMCHI2','p should be a positive integer')
end
sk = 2.*sqrt(2./p);
ku = 3+12/p;


