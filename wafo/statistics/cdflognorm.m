function [F,Flo,Fup] = cdflognorm(x,varargin)
%CDFLOGNORM Lognormal cumulative distribution function
%
% CALL:  F = cdflognorm(x,m,v,options);
%        [F,Flo,Fup] = cdflognorm(x,phat,options);
%
%        F = distribution function evaluated at x
%  Flo,Fup = 100*(1-alpha) % confidence bounds of F.
%      m,v = parameters
%     phat = Distribution parameter struct
%            as returned from FITLOGNORM.  
%  options = struct with fieldnames:
%         .lowertail: if TRUE (default), F = Prob[X <= x],
%                     otherwise, F = Prob[X > x].
%         .logp     : if TRUE, probability, p, returned as log(p).
%         .alpha    : Confidence coefficent    (default 0.05)
%
% The Lognormal distribution is defined by its pdf
%
%        f(x)=(v*2*pi*x^2)^(-1)*exp(-(log(x)-m)^2/(2*v)), x>=0.
%
% Example: 
%   x = linspace(0,6,200);
%   p1 = cdflognorm(x,0,1); p2 = cdflognorm(x,.5,0.25);
%   plot(x,p1,x,p2), shg
%
% See also pdflognorm, invlognorm, rndlognorm, fitlognorm, momlognorm

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



% Reference:
% Cohen & Whittle, (1988) "Parameter Estimation in Reliability
% and Life Span Models", p. 59 ff, Marcel Dekker.



% Tested on; Matlab 5.3
% History:
% revised pab 24.10.2000
%  - added comnsize, nargchk
%  - added default values
% added ms 10.08.2000

%error(nargchk(1,9,nargin))
narginchk(1,9)
options = struct('covariance',[],'alpha',0.05,...
  'lowertail',true,'logp',false); % default options
Np = 2;
[params,options] = parsestatsinput(Np,options,varargin{:});
if numel(options)>1
  error('Multidimensional struct of distribution parameter not allowed!')
end

[m,v] = deal(params{:});

if isempty(m),  m=0;  end
if isempty(v),  v=1;  end

v(v<=0) = nan;
x(x<=0) = 0; % trick to set cdf to zero
try
  s = sqrt(v);
  xn = (log(x)-m)./s;
catch
  error ('x, m and v must be of common size or scalar');
end
F = cdfnorm(xn,0,1,options);
if nargout>1
% TODO % Implement  Flo and Fup
 warning('WAFO:CDFLOGNORM','Flo and Fup not implemented yet')
 Flo = nan;
 Fup = Flo;
end

