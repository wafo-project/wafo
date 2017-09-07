function [m,v,sk,ku]= momt(df)
%MOMT Mean and variance for the Student's T  distribution.
% 
% CALL:  [m,v,sk,ku] = momt(df)
%
%   m, v = the mean and variance, respectively 
%  sk,ku = the skewness and kurtosis, respectively. 
%   df   = degrees of freedom of the Student's T distribution
%
%  Mean (m) and variance (v) for the T distribution is
%
%  m=0 if df>1  and  v=df/(df-2) if df>2
%
% Example:
%   par = {10};
%   X = rndt(par{:}, 10000,1);
%   moments = {mean(X) var(X),skew(X),kurt(X)};   % Estimated mean and variance
%   [mom{1:4}] = momt(par{:}); % True mean and variance
%   assert(moments, mom, -0.25);
%
% See also pdft, cdft, invt, rndt, fitt

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
% and Life Span Models", Marcel Dekker.


% Tested on; Matlab 5.3
% History: 
% by pab 23.10.2000

%error(nargchk(1,1,nargin))
narginchk(1,1)

m = zeros(size(df));
m(df<=1) = nan;
df(df<=2) = nan;
v = df./(df-2);
df(df<=3) = nan;
sk = m;
sk(isnan(df)) = nan;
df(df<=4) = nan;
ku = 3+6./(df-4);


