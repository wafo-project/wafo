function [x,xlo,xup] = invlognorm(F,varargin)
%INVLOGNORM Inverse of the Lognormal distribution function
%
% CALL:  x = invlognorm(F,m,v,options)
%        [x,xlo,xup] = invlognorm(F,phat,options)
%
%        x = inverse cdf for the Lognormal distribution evaluated at F
%      m,v = parameters     (default 0 and 1, respectively)
%     phat = Distribution parameter struct
%            as returned from FITLOGNORM.  
%  options = struct with fieldnames:
%         .lowertail: if TRUE (default), F = Prob[X <= x],
%                     otherwise, F = Prob[X > x].
%         .logp     : if TRUE, probability, p, returned as log(p).
%         .alpha    : Confidence coefficent    (default 0.05)
%
% Example:
%   F = linspace(0,1,100);
%   x = invlognorm(F,0,1);
%   plot(F,x)
%
% See also cdflognorm, pdflognorm, rndlognorm, fitlognorm, momlognorm

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
%  - fixed a bug: the inversion was not correct 
% added ms 10.08.2000

error(nargchk(1,9,nargin))
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

try
  x=exp(invnorm(F,m,v,options));
catch
  error ('F, m and v must be of common size or scalar');
end

if nargout>1
% TODO % Implement  xlo and xup
 warning('WAFO:INVLOGNORM','xlo and xup not implemented yet')
 xlo = nan;
 xup = xlo;
end


