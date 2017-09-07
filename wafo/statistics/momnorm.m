function [m,v,sk,ku]= momnorm(varargin)
%MOMNORM Mean and variance for the Normal distribution.
% 
% CALL:  [m,v] = momnorm(m0,v0)
%
%   m, v = the mean and variance, respectively 
% m0, v0 = parameters of the Normal distribution.
%
%  Mean (m) and variance (v) for the Normal distribution is
%
%  m=m0  and  v=v0;
%
% Example:
%   par = {-1,1};
%   X = rndnorm(par{:},1000,1);
%   moments = {mean(X) var(X),skew(X),kurt(X)};   % Estimated mean and variance
%   [mom{1:4}] = momnorm(par{:}); % True mean and variance
%   assert(moments, mom, -0.25);
%
% See also pdfnorm, cdfnorm, invnorm, rndnorm, fitnorm

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


Np = 2;
%error(nargchk(1,Np,nargin))
narginchk(1,Np)
options = [];
params = parsestatsinput(Np,options,varargin{:});

[m0,v0] = deal(params{:});
if isempty(v0)
  v0 = 1;
end
if isempty(m0)
  m0=0;
end
v0(v0<=0) = nan;
m =  m0;
v = v0;
sk = 0;
ku = 3;

[iscmn,m,v,sk,ku ] = iscomnsize(m,v,sk,ku);
if ~iscmn
  error ('m and v must be of common size or scalar');
end
