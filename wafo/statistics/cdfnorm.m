function [F,Flo,Fup] = cdfnorm(x,varargin)
%CDFNORM Normal cumulative distribution function
%
% CALL:  F = cdfnorm(x,m,v,options);
%        [F,Flo,Fup] = cdfnorm(x,phat,options);
%
%        F = distribution function evaluated at x
%  Flo,Fup = 100*(1-alpha) % confidence bounds of F.
%        x = matrix
%        m = mean     (default 0)
%        v = variance (default 1)
%     phat = Distribution parameter struct
%            as returned from FITNORM.  
%  options = struct with fieldnames:
%         .lowertail: if TRUE (default), F = Prob[X <= x],
%                     otherwise, F = Prob[X > x].
%         .logp     : if TRUE, probability, p, returned as log(p).
%         .alpha    : Confidence coefficent    (default 0.05)
% 
% Example: 
%   x = linspace(-3,3,200);
%   p1 = cdfnorm(x,0,1); p2 = cdfnorm(x,.5,0.25);
%   plot(x,p1,x,p2), shg
% 
% See also pdfnorm, momnorm, rndnorm, fitnorm

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
% revised pab 2007
% revised pab nov2005
% -removed comnsize -> faster execution
% revised pab 23.03.2003
% -changed call from erf to erfc in order 
%  to get more accurate lower probabilities
% -added a fix up for a bug in erfcore   
% revised pab 24.10.2000
%  - added comnsize, nargchk
%  - added default values
% added ms 15.06.2000
options = struct('proflog',false,'alpha',0.05,...
  'lowertail',true,'logp',false); % default options
if ischar(x) && strcmpi(x,'defaults')
  F = options;
  return
end
%error(nargchk(1,inf,nargin))
narginchk(1,inf)
Np = 2;
[params,options,tmp,phat] = parsestatsinput(Np,options,varargin{:});
if numel(options)>1
  error('Multidimensional struct of distribution parameter not allowed!')
end

[m,v] = deal(params{:});

if isempty(m),  m=0;  end
if isempty(v),  v=1;  end

sgn = (~options.lowertail)*2-1;
v(v<=0) = nan;
try
  z = sgn*(x-m)./sqrt(2*v);
catch
  error ('x, m and v must be of common size or scalar');
end


if options.logp
  F = log(0.5.*erfc(z));
  small_z = (z<-5);
  if any(small_z)
    F(small_z) = log1p(-0.5*erfc(-z(small_z)));
  end
else
  F = 0.5.*erfc(z);
end
% fix up for a bug in erfcore (Matlab R11 and earlier)
F(isnan(z)) = NaN;


if nargout >= 2
  if isempty(phat)
    error('Must have distribution struct!')
  end
  pcov = phat.covariance;
  alpha = options.alpha;
  % Compute confidence bounds using the delta method
   var_z = 0.5*(pcov(1,1) + 2*pcov(1,2)*z + pcov(2,2)*z.^2)./v;
   if any(var_z(:)<0)
      error('Covariance must be a positive semi-definite matrix.');
   end
   %xcrit = -invnorm(alpha/2)*sqrt(var_x);
   zcrit = erfcinv(alpha)*sqrt(2*var_z);
  
   zlo = z + zcrit;
   zup = z - zcrit;
   if options.logp
     Fup = log(0.5.*erfc(zup));
     Flo = log(0.5.*erfc(zlo));
     %small_z = (z<-5);
     if any(small_z)
       Flo(small_z) = log1p(-0.5*erfc(-zlo(small_z)));
       Fup(small_z) = log1p(-0.5*erfc(-zup(small_z)));
     end
   else
     Fup = 0.5.*erfc(zup);
     Flo = 0.5.*erfc(zlo);
   end
   % fix up for a bug in erfcore (Matlab R11 and earlier)
   Fup(isnan(z)) = NaN;
   Flo(isnan(z)) = NaN;
end
