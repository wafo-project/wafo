function [m,v,sk,ku]= momgev(varargin)
%MOMGEV Mean and variance for the GEV distribution.
% 
% CALL:  [m,v,sk,ku] = momgev(k,s,m0)
%
%       m, v = the mean and variance, respectively 
%      sk,ku = the skewness and kurtosis, respectively. 
%   k, s, m0 = parameter of the GEV distribution. (see cdfgev)
%
%  Mean (m) and variance (v) for the GEV distribution is
%
%  m = m0+s*(1-gamma(k+1))/k  and
%  v = (s/k)^2*(gamma(2*k+1)-gamma(k+1))
%
%  Note: mean only exists for k>-1 and variance for k>-0.5
%
% Example:
%   par = {0.5,2,0};  
%   X = rndgev(par{:}, 10000, 1);
%   moments = {mean(X) var(X),skew(X),kurt(X)};   % Estimated mean and variance
%   [mom{1:4}] = momgev(par{:}); % True mean and variance
%   assert(moments, mom, -0.2);
%
% See also pdfgev, cdfgev, invgev, rndgev, fitgev
  
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


% Tested on; Matlab 5.3
% History: 
% revised pab dec 2007
% -fixed a bug for small k
% revised pab 16 Feb 2006
% - variance is now correctly defined for k>-0.5 (not k>0)
%   (This bug is corrected thanks to Gerard Orjubin).
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
  gamk = expm1(gammaln(k+1))./k;
  gamk(abs(k)<=eps) = -0.577215664901532; %psi(1); % k==0
    
   m = m0-s.*gamk;
   
   
   k(k<=-0.5) = nan;
   g = @(n) gamma(n*k+1);
   
   g1 = g(1);
   g2 = g(2);
   g2mg12 = g2-g1.^2;
   sml = abs(k)<=1e-7;
   g2mg12(sml) = k(sml).^2*pi^2/6;

   
   gam2k =  expm1(gammaln(2*k+1)-2*gammaln(k+1))./k.^2;
   gam2k(sml) = pi^2/6; % psi(1,1);
   
   v = (s.*g1).^2.*gam2k;
   
   g3 = g(3);
   g3(k<=-1/3) = nan;
   
  
   
   % skewness
   %sk = sign(k).*(-g3+(3.*g2-2*g1.^2).*g1)./((g2mg12).^(3/2));
   sk = sign(k).*(-g3+(g2+2*g2mg12).*g1)./((g2mg12).^(3/2));
   zeta3 = 1.2020569; % apery's constant
   sml = abs(k)<=eps^0.29;
   sk(sml) = 12*sqrt(6)*zeta3/pi^3;
  
   g4 = g(4);
   g4(k<=-1/4) = nan;
   
   % The kurtosis is:
   ku = (g4+(-4.*g3+3*(g2+g2mg12).*g1).*g1)./((g2mg12).^2);
   %ku = (g4-4*g3.*g1+6*g2.*g1.^2-3*g1.^4)./((g2mg12).^2);
   
   ku(abs(k)<=(eps)^0.23) = 3+12/5;
   [csize,v,sk,ku] = comnsize(v,sk,ku,m);
catch
  error('k s and m0 must be of common size or scalar.');
end


