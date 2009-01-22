function [T, pvalue, Tn,H]=logps(phat,data,varargin)
%LOGPS Moran's negative log Product Spacings statistic
% 
% CALL:  [T,pvalue,Tn] = logps(phat,data,opt1,opt2,...,optN,dist);
%
%       T    = -sum(log(D_i(phat))), i.e., the negative log-product spacing
%     pvalue = Prob(X>Tn), probability of obtaining a result at least as 
%              extreme as Tn, assuming it was the result of chance alone.
%              Large p-value -> selected model gives a good fit to data.
%              Small p-value -> selected model gives a lesser fit to data.
%       Tn   = Normalized test statistic, which is approximately Chi2
%              distibuted with Nd degrees of freedom.
%       phat = [A1,A2,...,Am], vector of distribution parameters.
%       data = x1, or cellarray of vectors, ie., {x1,x2,..xN} 
%              where each xi are column vectors of data points (of length Nd).
%  opt1,opt2,
%  ...,optN  = additional input to distribution
%       dist = function handle or string defining the CDF to use.
%           
% LOGPS is a utility for parameter estimation as well as goodness of fit test.
% This function works on any continous CDF and PDF having the following 
% calling syntax:
%
%       F = prbtest(x1,x2,..,xN,A1,A2,..,Am,opt1,opt2,...,optN)
%       f = pdftest(x1,x2,..,xN,A1,A2,..,Am,opt1,opt2,...,optN)
%
% where x1,x2,...,xN contain the points to be evaluated and A1,A2... are
% the distribution parameters. 
%
%  The benefit of MPS estimation are: 
% 1) 	MPS estimators are consistent under much more general conditions than
%     ML-estimators. MPS method gives consistent estimators even in 
%     situations where ML fails typically for  distributions where 
%     the endpoints are unknown and the density J-shaped.
% 2) 	MPS estimators are asymptotically normal and asymptotically as 
%     efficient as ML estimators when these exist. 
%
% Example: %MPS and asymptotic covariance of phat:
%   R = rndweib(1,3,0,100,1);                      
%   phat0 = [1.3 2.5];                             %initial guess
%   phat = fminsearch('logps',phat0,[],R,'cdfweib')
%   [L, cov] = loglike(phat,R,'pdfweib')
%   [T,pval,Tn] = logps(phat,R,'cdfweib')
%   Tncrit = invchi2(1-0.05,length(R));
%
% See also loglike, proflog


% Copyright (C) 2000 Per A. Brodtkorb
%
% This program is free software; you can redistribute it and/or modify
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
%
% R. C. H. Cheng; N. A. K. Amin (1983)
% "Estimating Parameters in Continuous Univariate Distributions with a
% Shifted Origin.", 
% Journal of the Royal Statistical Society. Series B (Methodological),
%  Vol. 45, No. 3. (1983), pp. 394-403. 
%
% R. C. H. Cheng; M. A. Stephens (1989)
% "A Goodness-Of-Fit Test Using Moran's Statistic with Estimated
% Parameters", Biometrika, 76, 2, pp 385-392 
%
% Wong, T.S.T. and Li, W.K. (2006)
% "A note on the estimation of extreme value distributions using maximum
% product of spacings.",
% IMS Lecture Notes Monograph Series 2006, Vol. 52, pp. 272-283


%Tested on: matlab 5.3
% History:
% revised pab May 2007
% -added xmin
% -made sure LL is a number if x is NaN or infinite
% -dist can now also be a function handle.
% by pab 31.10.2000


error(nargchk(3,inf,nargin))

params = num2cell(phat(:).',1);



options  = varargin(1:end-1); % cell array of vectors with data points
dist   = varargin{end};
switch class(dist)
  case {'char','function_handle'} % OK
  otherwise
  error('Distribution is unspecified')
end
if isnumeric(data)
  if issorted(data(:))
    data= {data(:)};
  else
    data = {sort(data(:))};
  end
  dx = diff(data{1});
  tie = find(dx<=0);
else
  d = size(data{1});
  if any(d(2:end)>1)
    for ix=1:length(data)
      data{ix}  = data{ix}(:); %% make sure it is a vector.
    end
  end
  if ~issorted([data{:}],'rows')
    data = num2cell(sortrows([data{:}]),1);
  end
  dx = diff([data{:}]);
  tie = find(all(dx==0,2));
end

%xmin = eps;

prb = feval(dist,data{:},params{:},options{:});
if numel(data)>1
  [prb,ix1] = sort(prb(:));
end
dprb = diff(prb);
lowertail = any(dprb>0);
if lowertail
  logD = log([prb(1);max(dprb(:),0);1-prb(end)]);
else
  logD = log(-[prb(1)-1;min(dprb(:),0);-prb(end)]);
end
if any(tie)
  if false %nargout>1
% TODO % implement this method for treating ties in data:
% Assume measuring error is delta. Then compute
% yL = F(xi-delta,theta)
% yU = F(xi+delta,theta)
% and replace 
% logDj = log((yU-yL)/(r-1)) for j = i+1,i+2,...i+r-1

  else
    % This is OK when only minimization of T is wanted
    model = getdistname(dist);
    tiedata = cell(1,numel(data));
    for ix = 1:numel(data)
      tiedata{ix} = data{ix}(tie);
    end
    if numel(data)>1
      ix2 = 1:length(prb);
      ix2  = ix2(ix1);
      logD(ix2(tie)+1)= feval(['pdf',model],tiedata{:},params{:},options{:},'logp',true);
    else
      logD(tie+1)= feval(['pdf',model],tiedata{:},params{:},options{:},'logp',true);
    end
  end
end
if any(~isfinite(logD))
  T =  -sum(logD(isfinite(logD))) + 100*log(realmax)*sum(~isfinite(logD));
else
  T = - sum(logD); %Moran's negative log product spacing statistic 
end
if nargout>1
  if any(tie)
    warning('WAFO:LOGPS','pvalue is on the conservative side (i.e. too large) due to ties in the data!')
  end
  np1 = size(logD,1);
  n = np1-1;
  k = length(params);
  isParUnKnown = true;
  m = (np1)*(log(np1)+0.57722)-0.5-1/(12*(np1));
  v = (np1)*(pi^2/6-1)-0.5-1/(6*(np1));
  C1 = m-sqrt(0.5*n*v);
  C2 = sqrt(v/(2*n));
  Tn = (T + 0.5*k*isParUnKnown-C1)./C2; % chi2 with n degrees of freedom
  pvalue = cdfchi2(Tn,n,'lowertail',false);
  
  if nargout>3
    mylogps = @(phat1)logps(phat1,data,options{:},dist);
    %[H,err] = hessian(mylogps,phat,'RombergTerms',4);
     H = hessian(mylogps,phat,'RombergTerms',4);
  end
 
end

end
