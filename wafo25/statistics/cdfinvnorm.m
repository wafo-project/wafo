function [F,Flo,Fup] = cdfinvnorm(x,varargin)
%CDFINVNORM Inverse Gaussian cumulative distribution function
%
% CALL:  F = cdfinvnorm(x,m,l,options);
%        [F,Flo,Fup] = cdfinvnorm(x,phat,options);
%
%        F = distribution function evaluated at x
%  Flo,Fup = 100*(1-alpha) % confidence bounds of F.
%      m,l = parameters
%     phat = Distribution parameter struct
%            as returned from FITINVNORM.  
%  options = struct with fieldnames:
%         .lowertail: if TRUE (default), F = Prob[X <= x],
%                     otherwise, F = Prob[X > x].
%         .logp     : if TRUE, probability, p, returned as log(p).
%         .alpha    : Confidence coefficent    (default 0.05)
%
% The Inverse Gaussian distribution is defined by its pdf
%
%        f(x)=(l/(2*pi*x^3))^(1/2)*exp(-l*(x-m)^2/(2*m^2*x)), x>0.
%
% Example: 
%   x = linspace(0,3,200);
%   p1 = cdfinvnorm(x,1,1); p2 = cdfinvnorm(x,1,.25);
%   plot(x,p1,x,p2),shg
%
% See also pdfinvnorm, invinvnorm, rndinvnorm, fitinvnorm, mominvnorm

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
% and Life Span Models", p. 259 ff, Marcel Dekker.



% Tested on; Matlab 5.3
% History: 
% revised pab 24.10.2000
%  - added comnsize, nargchk
% added ms 14.08.2000

options = struct('covariance',[],'alpha',0.05,...
  'lowertail',true,'logp',false); % default options

if (nargin==1 && nargout <= 1 && isequal(x,'defaults'))
  F = options; 
  return
end
error(nargchk(2,inf,nargin))
Np = 2;
[params,options] = parsestatsinput(Np,options,varargin{:});
if numel(options)>1
  error('Multidimensional struct of distribution parameter not allowed!')
end

[m,L] = deal(params{:});
m(m<=0) = nan;
L(L<=0) = nan;
x(x<=0) = 0;
try
   sqx = sqrt(x);
   sqL = sqrt(L);
   y1 = sqL.*(sqx./m-1./sqx);
   y2 = sqL.*(sqx./m+1./sqx);
catch
   error('x, m and l must be of common size or scalar');
end
% TODO % Refine lowertail and logp
F=(cdfnorm(y1,0,1)+exp(2*L./m).*cdfnorm(-y2,0,1));

if options.logp
  if options.lowertail
    F = log(F);
  else
    F = log1p(-F);
  end
elseif ~options.lowertail
  F = 1-F;
end

 if nargout>1
% TODO % Implement  Flo and Fup
   warning('WAFO:CDFINVNORM','Flo and Fup not implemented yet')
   Flo = nan;
   Fup = Flo;
 end
 
return

% [iscmn, x, m, l] = iscomnsize (x,m, l);
% if ~iscmn
%   error ('x, m and l must be of common size or scalar');
% end
% 
% F=zeros(size(x));
% ok=((m>0)&(l>0));
% 
% k=find(x>0& ok);
% if any(k)
%   F(k)=(cdfnorm((l(k)./x(k)).^(1/2).*(x(k)./m(k)-1),0,1)+...
%       exp(2*l(k)./m(k)).*cdfnorm(-(l(k)./x(k)).^(1/2).*(x(k)./m(k)+1),0,1));
% end
% 
% k1 = find (~ok);
% if any (k1)
%   tmp=NaN;
%   F(k1) = tmp(ones(size(k1)));
% end