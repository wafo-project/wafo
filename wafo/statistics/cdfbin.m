function  [F,Flo,Fup] = cdfbin(x,varargin)
%CDFBIN  Binomial cumulative probability function
%
%  CALL  F = cdfbin(x,n,p,options)
%        F = cdfbin(x,phat,options)
%
%        F = probability of observing x or less successes in n independent
%            trials.
%        x = number of successes          (0<=x<=n)
%        n = number of independent trials
%        p = probability of succes in a given trial. (0<=p<=1)
%     phat = Distribution parameter struct
%            as returned from FITBIN.  
%  options = struct with fieldnames:
%         .lowertail: if TRUE (default), F = Prob[X <= x],
%                     otherwise, F = Prob[X > x].
%         .logp     : if TRUE, probability, p, returned as log(p).
%         .alpha    : Confidence coefficent    (default 0.05)
%
% Example
%  n = 10; p = 0.05;
%  x = -1:n+1;
%  F = cdfbin(x,n,p);
%  stairs(x,F);
%
%  close all;
%
% See also cdfbin, invbin, rndbin, fitbin, mombin

%       Anders Holtsberg, 27-07-95
%       Copyright (c) Anders Holtsberg
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


options = struct('covariance',[],'alpha',0.05,...
  'lowertail',true,'logp',false); % default options

if ischar(x) && strcmpi(x,'defaults')
  F = options;
  return
end

error(nargchk(2,inf,nargin))

Np = 2;
[params,options] = parsestatsinput(Np,options,varargin{:});
if numel(options)>1
  error('Multidimensional struct of distribution parameter not allowed!')
end

[n,p] = deal(params{:});
try
  %F = betainc(1-p,n-x,x+1);
  nmx = min(max(n-x,realmin),n);
  xp1 = min(max(x+1,realmin),n);
  
  F = cdfbeta(1-p,nmx,xp1,options);
  
  prbZero = (xp1<1 & nmx>-1 & p>-1);
  prbOne =  (nmx<1 & xp1>-1 & p>-1);
  if options.lowertail
    if options.logp
      F(prbZero) = -inf;
      F(prbOne)  = 0;
    else
      F(prbZero) = 0;
      F(prbOne)  = 1;
    end
  else
    if options.logp
      F(prbOne)  = -inf;
      F(prbZero) = 0;     
    else
      F(prbOne)  = 0;
      F(prbZero) = 1;
    end
  end
  if nargout>1
% TODO % Confidence bands for F, i.e., Flo and Fup.
% TODO % Implement  Flo and Fup
    warning('WAFO:PRBBIN','Flo and Fup not implemented yet')
    Flo = nan;
    Fup=Flo;   
  end
  
catch
  error('x, n and p must be of common size or scalar.');
end
return
% old call
%       Anders Holtsberg, 27-07-95
%       Copyright (c) Anders Holtsberg
% try
%   [m,v] = mombin(n,p);
%   F = cdfnorm(x+0.5,m,v,options); % normal approximation with cont. correction
%   
%   % Test if n is large enough so that the normal approximation is OK
%   useexact = (m-3*sqrt(v)<0 | n < m+3*sqrt(v));
%   
%   if any(useexact(:))
%     if numel(useexact)>1
%       if numel(p)>1, p = p(useexact); end
%       if numel(n)>1, n = n(useexact); end
%       if numel(x)>1, x = x(useexact); end
%     
%       [csize,x,n,p] = comnsize(x,n,p);
%     
%       [np,ii,jj] = unique([n(:),p(:)],'rows');
%       idx = 1:size(np,1);
%       grp = idx(jj);
%       F0 = zeros(csize);
%       for ix = idx
%         F0(grp==ix) = cdfbin1(x(grp==ix),np(ix,1),np(ix,2),options);
%       end
%       F(useexact) = F0;
%     else
%       F = cdfbin1(x,n,p,options);
%     end
%   end
% catch
%   error('x, n and p must be of common size or scalar.');
% end
% 
% 
% function F = cdfbin1(x,n,p,options)
% 
% if max([numel(n) numel(p)]) > 1
%    error('Sorry, this is not implemented');
% end
% F = x;
% if any(n < 0 | p < 0 | p > 1 );
%  F(:) = nan;  
% else
%   kk = (0:n)';
%   cdf = max(0,min(1,[0; cumsum(pdfbin(kk,n,p))]));
%   cdf(n+2) = 1;
% 
%   F(:) = cdf(max(1,min(n+2,floor(x(:))+2)));
% end
% if options.logp
%   if options.lowertail
%     F = log(F);
%   else
%     F = log1p(-F);
%   end
% elseif ~options.lowertail
%   F = 1-F;
% end
