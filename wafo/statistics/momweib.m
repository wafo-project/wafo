function [m,v,sk,ku]= momweib(a,b, c)
%MOMWEIB Mean and variance for the Weibull  distribution.
% 
% CALL:  [m,v,sk,ku] = momweib(a, b, c)
%
%   m, v = the mean and variance, respectively 
%  sk,ku = the skewness and kurtosis, respectively. 
%   a, b, c = parameters of the Weibull distribution (see cdfweib).
%
%  Mean (m) and variance (v) for the Weibull distribution is
%
%  m=m0 + c  and  v=a^2*gamma(1+2/b)-m0^2 where m0 = a*gamma(1+1/b).
%
% Example:
%   par = {1, 2, 1};
%   X = rndweib(par{:}, 10000,1);
%   moments = {mean(X) var(X),skew(X),kurt(X)};   % Estimated mean and variance
%   [mom{1:4}] = momweib(par{:}); % True mean and variance
%   assert(moments, mom, -0.25);
%
% See also pdfweib, cdfweib, invweib, rndweib,  fitweib

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
% and Life Span Models", p. 25 ff, Marcel Dekker.


% Tested on; Matlab 5.3
% History: 
% revised pab 23.10.2000
%   - added comnsize
% added ms 15.06.2000

%error(nargchk(2,3,nargin));
narginchk(2,3)
if nargin<3 || isempty(c),
  c=0;
end
a(a<=0) = nan;
b(b<=0) = nan;
try
  m0 =  a .* gamma(1 + (1 ./ b)) + c*0;
  v = a.^ 2 .* gamma(1 + (2 ./ b)) - m0.^ 2;
  sk = (a.^3.*gamma(1+(3./b))-3.*m0.*v-m0.^3)./v.^(3/2);
  ku = (a.^4*gamma(1+4./b)-4*sk.*v.^(3/2).*m0-6*m0.^2.*v-m0.^4)/v.^2;
  m = m0 + c;
catch
  error('a, b and c must be of common size or scalar.');
end
