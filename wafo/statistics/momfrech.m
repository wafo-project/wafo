function [m,v,sk,ku]= momfrech(varargin)
%MOMFRECH Mean and variance for the Frechet distribution.
% 
% CALL:  [m,v,sk,ku] = momfrech(a,c)
%
%   m, v = the mean and variance, respectively 
%  sk,ku = the skewness and kurtosis, respectively. 
%   a, c = parameters of the Frechet distribution (see cdffrech).
%
%  Mean (m) and variance (v) for the Frechet distribution is
%
%  m=a*gamma(1-1/c)  (if c>1) and  v=a^2*(gamma(1-2/c))-m^2  (if c>2)
%
% Example:
%   par = {1,10}
%   X = rndfrech(par{:},10000,1);
%   [mean(X) var(X),skew(X),kurt(X)]        % Estimated mean and variance
%   [m,v,sk,ku] = momfrech(par{:}) % True mean and variance
%
% See also  pdffrech, cdffrech, invfrech, rndfrech, fitfrech 


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

% Reference: 


% Tested on; Matlab 5.3
% History: 
% revised pab 
% - rmoved dependence on comnsize
% Added PJ 10-May-2001

error(nargchk(1,2,nargin))
Np = 2;
options = [];
params = parsestatsinput(Np,options,varargin{:});

[a,c] = deal(params{:});
if isempty(c)
  error('Shape parameter b undefined!')
end


try 
  a(a<=0) = nan;
  c(c<=1) = nan;
  m =  a .* gamma(1 - (1 ./ c));
  c(c<=2) = nan;
  v = a.^ 2 .* gamma(1 - (2 ./ c)) - m.^ 2;
  c(c<=3) = nan;
  sk = (a.^3.*gamma(1-(3./c))-3.*m.*v-m.^3)./v.^(3/2);
  c(c<=4) = nan;
  ku = (a.^4*gamma(1-4./c)-4*sk.*v.^(3/2).*m-6*m.^2.*v-m.^4)/v.^2;
catch
   error('a and c must be of common size or scalar.');
end

