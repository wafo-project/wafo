function R = rndlognorm(varargin)
%RNDLOGNORM Random matrices from a Lognormal distribution.
%
% CALL:  R = rndlognorm(mu,v,sz,options);
%
%     mu, v = parameters (see pdflognorm) (Default 0 and 1, respectively)
%        sz = size(R)    (Default common size of mu and v)
%             sz can be a comma separated list or a vector 
%             giving the size of R (see zeros for options).
% Example:
%   R = rndlognorm(1,2,100,2);
%   plotnorm(log(R)),shg
%
% See also pdflognorm, cdflognorm, invlognorm, fitlognorm, momlognorm

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
% revised pab 23.10.2000
%  - added comnsize, nargchk
%  - added greater flexibility on the sizing of R
% added ms 10.08.2000


error(nargchk(1,inf,nargin))
Np = 2;
options = struct; % default options
[params,options,rndsize] = parsestatsinput(Np,options,varargin{:});
if numel(options)>1
  error('Multidimensional struct of distribution parameter not allowed!')
end

mu = params{1};
v = params{2};


if isempty(mu),mu=0;end
if isempty(v), v=1;end

if isempty(rndsize)
  csize = comnsize(mu,v);
else
  csize = comnsize(mu,v,zeros(rndsize{:}));
end
if any(isnan(csize))
    error('mu and v must be of common size or scalar.');
end

R = exp(rndnorm(mu,v,csize));

