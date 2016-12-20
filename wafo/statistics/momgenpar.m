function [m,v,sk,ku]= momgenpar(varargin)
%MOMGENPAR Mean and variance for the Generalized Pareto distribution.
% 
% CALL:  [m,v,sk,ku] = momgenpar(k,s,m0)
%
%       m, v = the mean and variance, respectively 
%      sk,ku = the skewness and kurtosis, respectively. 
%  k, s, m0  = parameters of the  Generalized Pareto distribution
%              (see cdfgenpar).
%
% Mean (m) and variance (v) for the Generalized Pareto distribution is
%
%  m=m0+s/(1+k)  and
%  v=s^2/((1+k)^2*(1+2k))
%
% The mean does not exist for k<-1, and the variance does not exist for 
% k<-0.5.
%
% Example:
%   par = {0.5,1,1};
%   X = rndgenpar(par{:}, 100000, 1);
%   moments = {mean(X) var(X),skew(X),kurt(X)};   % Estimated mean and variance
%   [mom{1:4}] = momgenpar(par{:}); % True mean and variance
%   assert(moments, mom, -0.2);
%
% See also pdfgenpar, cdfgenpar, invgenpar, rndgenpar, fitgenpar

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

% Reference
% AC Davison (1983)
% "Modelling excesses over high thresholds "
% In Proceedings of the NATO Advanced Study Institute on 
% Statistical Extremes and Applications,
% Vimeiro, Portugal, 31 August-14 September, 1983"--T.p. verso.

% Tested on; Matlab 5.3
% History: 
% revised pab 2007
% -added sk and ku
% -removed dependence on comnsize
% Revised by PJ 02-Apr-2001
%  - Added non-existing mean and var (k<-1 and k<-0.5).
% revised pab 24.10.2000
%  - added comnsize, nargchk + default value for s and m0
% added ms 09.08.2000


Np = 3;
error(nargchk(1,Np,nargin))
options = [];
params = parsestatsinput(Np,options,varargin{:});

[k,s,m0] = deal(params{:});
if isempty(s)
  s = 1;
end
if isempty(m0)
  m0=0;
end


k(k<=-1) = nan;
s(s<=0) = nan;

try
   m = m0+s./(1+k);
   k(k<=-0.5) = nan;
   v =s.^2./((1+k).^2.*(1+2*k));
   k(k<=-1/3) = nan;
   sk = 2.*(1-k).*sqrt(1+2.*k)./(1 +3.*k);
   k(k<=-1/4) = nan;
   
   % E(X^r) = s^r*(-k)^-(r+1)*gamma(1+r)*gamma(-1/k-r)/gamma(1-1/k)
   %  = s^r*gamma(1+r)./( (1+k)*(1+2*k).*....*(1+r*k))
   % E[(1-k(X-m0)/s)^r] = 1/(1+k*r)
   
   %Ex3 = (sk.*sqrt(v)+3*m).*v+m^3
   %Ex3 = 6.*s.^3/((1+k).*(1+2*k).*(1+3*k))
   r = 4;
   Ex4 = s^r*gamma(1+r)./((1+k).*(1+2*k).*(1+3*k).*(1+4.*k));
   m1 = m-m0;
   ku = (Ex4-4*sk.*v.^(3/2).*m1-6*m1.^2.*v-m1.^4)/v.^2;
  
catch
  error('k s and m0 must be of common size or scalar.');
end

