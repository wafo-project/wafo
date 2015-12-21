function [Lp,CI,nphat]=proflog(phat1,varargin)
%PROFLOG Profile Log- likelihood or Product Spacing-function.
% 
% CALL:  [Lp,CI,nphats] = proflog(phats,options);
%
%         Lp = Profile log-likelihood function with parameters phat given
%              the data, phat(i) and quantile x (if given), i.e., 
%                Lp = max(log(f(phat|data,phat(i)))),
%              or
%                Lp = max(log(f(phat|data,phat(i),x))) 
%         CI = Confidence interval for either phat(i) or quantile x.
%     nphats = struct with updated distribution parameters.
%      phats = struct with ML or MPS estimated distribution parameters.
%    options = struct with fieldnames:
%         .i           : integer defining which distribution parameter to 
%                        profile, i.e. which parameter to keep fixed (default 1)
%         .pmin, .pmax : Interval for either the parameter, phat(i) or x, used in the
%                        optimization of the profile function (default is
%                        based on the 100*(1-alpha)% confidence interval
%                        computed using the delta method.)
%         .N           : Number of points used in Lp (default 100)
%         .x           : quantile (return value)
%         .logR        : log survival probability,i.e., R = Prob(X>x;phat)
%         .link        : function connecting the quantile (x) and the 
%                        survival probability (R) with the fixed distribution
%                        parameter, i.e.: phat(i) = link(x,logR,phat,i), 
%                        where logR = log(Prob(X>x;phat)).
%                        This means that if link:
%                         1) and x is not empty then x is profiled 
%                         2) and logR is not empty then logR is profiled
%                         3) is empty then phat(i) is profiled (default)
%         .alpha       : confidence coefficent (default 0.05)
%         .plotflag    :     
%         .optimset    : optimset structure defining performance of the
%                        optimization routine (see optimset for details)
%             
%  PROFLOG is a utility function for making inferences either on a particular
%  component of the vector phat or the quantile, x, or the probability, R.
%  This is usually more accurate than using the delta method assuming 
%  asymptotic normality of the ML estimator or the MPS estimator.
%  This works on any PDF and/or CDF having the following calling syntax:
%
%       f = pdftest(x,A1,A2,..,Am,opt1)
%       F = cdftest(x,A1,A2,..,Am,opt2)
%       x = invtest(F,A1,A2,..,Am,opt2)
%
% where x is the quantile and A1,A2... are the distribution parameters. 
% opt1 and opt2 are options structs which must have the fields 'logp' and 
% 'lowertail', respectively
%
% Examples % MLE and better CI for phat.params(1)
%   R = rndweib(1,3,0,100,1);                      
%   phat = fitweib(R)         
%   opt = proflog('defaults');
%   opt.i = 2;opt.plotflag=1;
%   [Lp,CI] = proflog(phat,opt)
%   [x,xlo,xup] = invweib(1/990,phat,'lowertail',false);
%   [Lp0,CI0] = proflog(phat,'i',2,'x',x,'link',@lnkweib,'plot',1);
%
% % MPS and better confidence interval for the quantile, x.
%   R2 = rndgenpar(0.3,1,0,200,1);
%   m = 1;
%   [phat2] = fitgenpar(R2(R2>m)-m,'method','mps');
%   [x,xlo,xup] = invgenpar(1/9900,phat2,'lowertail',false);
%   [Lp2,CI2] = proflog(phat2,'i',2,'x',x,'link',@lnkgenpar,'plot',1);
%   [Lp3,CI3] = proflog(phat2,'i',2,'logR',-log(9900),'link',@lnkgenpar,'plot',1);
%
% See also loglike, logps, optimset


% Copyright (C) 2007 Per A. Brodtkorb
%
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

% Reference 
% Stuart Coles (2004)
% "An introduction to statistical modelling of extreme values".
% Springer series in statistics

%Tested on: matlab 7.4
% History:
% by pab Sep 2007
if nargin==1 && ischar(phat1)
  Lp = ciproflog(phat1);
  return
end
switch nargout
  case {0,1,2}
    [CI,Lp]=ciproflog(phat1,varargin{:});
  otherwise
    [Lp,CI,nphat]=ciproflog(phat1,varargin{:});
end

return

% error(nargchk(1,inf,nargin))
% options = struct('i',1,'pmin',[],'pmax',[],'N',100,'x',[],'logR',[],'link',[],'lowertail',false,'alpha',0.05,...
%   'plotflag',0,'check',true,'warning',true,'optimset', optimset('disp','off','TolX',1e-5,'TolFun',1e-5,'MaxIter',500));
% if ischar(phat1) && strcmpi(phat1,'defaults')
%   Lp = options;
%   return
% end
% if isnumeric(phat1)
% %% Secret call kept just in case
%   % CALL    [Lp,CI,nphats] = proflog(phatv,data,opt1,opt2,...,optN,dist,options);
%   %
%   %      phatv = [A1,A2,...,Am], vector of distribution parameters.
%   %       data = x1, or cellarray of vectors, ie., {x1,x2,..xN}
%   %              where each xi are column vectors of data points.
%   %  opt1,opt2,
%   %  ...,optN  = additional input to distribution
%   %       dist = function handle or string defining the PDF to use.
%   error(nargchk(3,inf,nargin))
% 
%   phatv = phat1;
%   ni = nargin-1;
%   if isstruct(varargin{end})
%     options = parseoptions(options,varargin{ni});
%     ni = ni-1;
%   end
%   dist  = varargin{ni};
%   data1  = varargin{1}; % cell array of vectors with data points
% 
%   if isnumeric(data1)
%     data1 = {data1};
%   end
%   for ix=1:length(data1),
%     data1{ix}  = data1{ix}(:); %% make sure it is a vector.
%   end
%  
%   phat1 = mlest(dist,phatv,data1,{'method','ML','search',false},varargin{2:ni-1});
% else
%   options = parseoptions(options,varargin{:});
% end
% 
% 
% %% initializing
% logfun = @loglike;
% likestr  ='likelihood';
% if isa(phat1,'struct') || isa(phat1,'fdata')
%  
%   best_phatv  = phat1.params;
%   dist  = phat1.distribution;
%   distname = getdistname(dist);
%   if isempty(phat1.pdfoptions)
%     pdfoptions = {};
%   else
%     pdfoptions = {phat1.pdfoptions};
%   end
%   if iscell(phat1.data)
%     data1 = phat1.data;
%   else
%     data1 = {sort(phat1.data(:))};
%   end
%   if strcmpi(phat1.method,'ML')
%     dist = ['pdf' distname ];
%     LLmax = phat1.loglikemax;   
%   else
%     likestr  ='product spacing';
%     dist = ['cdf' distname ];
%     LLmax = phat1.logpsmax;
%     logfun = @logps;
%     if ~strcmpi(phat1.method,'MPS')
%       warning('WAFO:PROFLOG','PHAT is apparently not a ML- or MPS- estimator!\n CI might be inaccurate.')
%     end
%   end
% else
%  error('Something strange happened')
% end
% if isempty(phat1.fixpar)
%   isnotfixed = repmat(true,size(phat1.params));
% else
%   isnotfixed = ~isfinite(phat1.fixpar);
% end
% 
% i_notfixed = find(isnotfixed);
% if 1
%   i_fixed = options.i;
%   if ~isnotfixed(i_fixed)
%     freetxt = sprintf(' %d,',i_notfixed);
%     freetxt(end) = '';
%     error('options.i must be one equal to an index to one of the free parameters, i.e., one of%s',freetxt)
%   end
% else
%   numfree = sum(isnotfixed);  % number of free parameters.
%   if numfree<options.i
%     error('options.i must be less or equal to %d',numfree)
%   end
%   i_fixed = i_notfixed((1:length(i_notfixed))==options.i);
% end
% isfree     = isnotfixed;
% isfree(i_fixed) = false;
% 
% phatv = best_phatv;
% 
% 
% LLrange     = 0.5 * invchi2(options.alpha, 1,'lowertail',false);
% 
% cross_level = LLmax - LLrange;
% lowLevel    = cross_level -LLrange/7;
% 
% link = options.link;  
% 
% 
% %% Check that phatv are actually at the optimum
% foundNewphat = false;
% phatfree = phatv(isfree);
% if options.check
%   if ~isempty(phatfree)
%     phatfree = fminsearch(@mylogfun,phatfree,options.optimset,phatv(i_fixed),[]);
%   end
%   LLt = -mylogfun(phatfree,phatv(i_fixed),[]);
%   if LLt>LLmax
%     foundNewphat = true;
%     best_phatv(isfree) = phatfree;
%     phatv = best_phatv;
%     if options.warning
%       warning('WAFO:proflog','newphat = %s', sprintf('%g ',best_phatv))
%     end
%     dL = LLmax-LLt;
%     lowLevel= lowLevel-dL;
%   
%     cross_level = cross_level-dL;
%     LLmax = LLt;
%   end
% end
% 
% %% Set up variables to profile
% if ~isempty(link) 
%   prbfun = str2func(['cdf',distname]);
%   invfun = str2func(['inv' distname]);
%   optsfun = {pdfoptions{:},'lowertail',options.lowertail,'logp',true};
%   if isempty(options.x)
%     if isempty(options.logR)
%       error('You must supply a non-empty quantile (x) or probability (logR) in order to profile it!')
%     else
%       isProfileX = false;
%        p_opt = options.logR;
%        qntl = feval(invfun,p_opt,phat1,optsfun{:});
%        xlabtxt = 'log(R)';
%     end
%   else
%      isProfileX = true;
%      p_opt = options.x;
%      prb = feval(prbfun,p_opt,phat1,optsfun{:});
%      xlabtxt = 'x';
%   end
% else
%   xlabtxt = sprintf('phat(%d)',i_fixed);
%   p_opt = phatv(i_fixed);
% end
% %% Find proper interval for the variable to profile.
% if isempty(options.pmin) || isempty(options.pmax)   
%   if ~isempty(link)
%     if isProfileX
%       drl = gradest(@myinvfun,phatv(isnotfixed));
%     else
%       drl = gradest(@myprbfun,phatv(isnotfixed));
%     end
%     pcov = phat1.covariance(isnotfixed,isnotfixed);
%     pvar = sum((drl*pcov).*drl,2);
%   else
%     pvar = phat1.covariance(i_fixed,i_fixed);
%   end
% %   if pvar<=0
% %     disp('PROFLOG: variance zero')
% %   end
%   p_crit = -invnorm(options.alpha/2)* sqrt(pvar)*1.5;
%   if isempty(options.pmin)
%     options.pmin = p_opt - 5 * p_crit;
%   end
%   if isempty(options.pmax)
%     options.pmax = p_opt + 5 * p_crit;
%   end
%   N4 = floor(options.N/4);
%   pvec = unique([ linspace(options.pmin,p_opt - p_crit,N4+1),...
%       linspace(p_opt - p_crit,p_opt + p_crit,options.N-2*N4),...
%       linspace(p_opt + p_crit,options.pmax,N4+1)]);
% 
% else
%   pvec = linspace(options.pmin,options.pmax,options.N);
% end
% 
% 
% %% Find Profile log likelihood or Product spacing function LL
% LL = pvec;
% LL(:) = lowLevel;
% 
%   
% k1 = find(pvec>p_opt, 1,'first');
% for ix = k1:-1:1
%   if ~isempty(phatfree)
%     phatfree = fminsearch(@mylogfun,phatfree,options.optimset,pvec(ix),link);
%   end
%   LL(ix) = -mylogfun(phatfree,pvec(ix),link);
%   if LL(ix)<cross_level
%     break
%   end
% end
% phatfree = phatv(isfree);
% for ix = k1+1:options.N
%   if ~isempty(phatfree)
%     phatfree = fminsearch(@mylogfun,phatfree,options.optimset,pvec(ix),link);
%   end
%   LL(ix) = -mylogfun(phatfree,pvec(ix),link);
%   if LL(ix)<cross_level
%     break
%   end
% end
% 
% 
% %% Check that phatv are actually at the optimum
% if options.check
%   [LLt,ix] = max(LL);
%   iz = max(ix-1,1):min(ix+2,options.N-1);
%   try
%     dLL = diff(LL([iz iz(end)+1]));
%     ind = findcross(dLL,0);
%   catch
%     ind = [];
%   end
%   if any(ind)
%     p_opt = ecross((pvec(iz)+pvec(iz+1))/2,dLL,ind,0);
%   else
%     p_opt = pvec(ix);
%   end
%   
%   if isempty(link)
%     phatv(i_fixed) = p_opt;
%   elseif isProfileX
%     phatv(i_fixed) = feval(link,p_opt,prb,phatv,i_fixed);
%   else
%     phatv(i_fixed) = feval(link,qntl,p_opt,phatv,i_fixed);
%   end
%   phatfree = phatv(isfree);
%   if ~isempty(phatfree)
%     phatfree = fminsearch(@mylogfun,phatfree,options.optimset,phatv(i_fixed),[]);
%   end
%   LLt = -mylogfun(phatfree,phatv(i_fixed),[]);
% 
%   if LLt>LLmax
%     foundNewphat = true;
%     best_phatv(i_fixed) = phatv(i_fixed);
%     best_phatv(isfree) = phatfree;
%     if options.warning
%       warning('WAFO:proflog','newphat = %s', sprintf('%g ',best_phatv))
%     end
%     dL = LLmax-LLt;
%     lowLevel        = lowLevel-dL;
%     LL(LL<lowLevel) = lowLevel;
%     cross_level     = cross_level-dL;
%     LLmax           = LLt;
%   end
% end
% if nargout>2
%   if foundNewphat
%     fixpar = phat1.fixpar;
%     opts1 = {'method',phat1.method,'search',false,'fixpar',fixpar};
%     nphat = mlest(dist,best_phatv(isnotfixed),data1,opts1);
%     nphat.dataname = phat1.dataname;
%   else
%     nphat = phat1;
%   end
% end
% 
% %% Confidence interval (CI)     
% ind = findcross(LL,cross_level);
% switch length(ind)
%   case 0
%     CI  = [options.pmin, options.pmax];
%     warning('WAFO:PROFLOG','Upper bound for %s is larger!',xlabtxt)
%     warning('WAFO:PROFLOG','Lower bound for %s is smaller!',xlabtxt)
%   case 1
%     x0 = ecross(pvec,LL,ind,cross_level);
%     isUpcrossing = LL(ind)<LL(ind+1);
%     if isUpcrossing
%       CI = [x0 options.pmax];
%       warning('WAFO:PROFLOG','Upper bound for %s is larger!',xlabtxt)
%     else
%       CI = [options.pmin x0];
%       warning('WAFO:PROFLOG','Lower bound for %s is smaller!',xlabtxt)
%     end
%   case 2
%     CI = ecross(pvec,LL,ind,cross_level);
%   otherwise
%     error('Number of crossings too large')
% end
% 
% %% Prettify result
% ix = find(LL==lowLevel);
% LL(ix) = [];
% pvec(ix) = [];
% ind = findcross(LL,lowLevel-LLrange/10);
% if any(ind)
%   % Strip off any very small values because it will only obscure the profile
%   % loglikelihood plot
%   p0 = ecross(pvec,LL,ind,lowLevel-LLrange/10);
%   isDown = LL(ind)>LL(ind+1);
%   ind = ind(:)+isDown(:);
%   LL(ind) = lowLevel-LLrange/10;
%   pvec(ind) = p0;
% end
% %% Put result into one object
% Lp = createwdata('data',LL,'args',pvec,...
%   'labels',{sprintf('%s (%s)',xlabtxt,distname),...
%   sprintf('%g %s CI',100*(1-options.alpha),'%')},...
%   'title',sprintf('Profile Log-%s ',likestr),...
%   'dataCI',repmat([LLmax cross_level],length(pvec),1),...
%   'workspace',options);
% 
% if ~isoctave
%   Lp = wdata(Lp);
%   if (options.plotflag)
%     plot(Lp,options.plotflag)
%     % hline([LLmax cross_level],'g--')
%     %vline(CI, 'r')
%     hold on
%     plot(repmat(CI,2,1),repmat([lowLevel;LLmax],1,2),'r'),hold off
%   end
% end
% %% Nested functions
% %   function L = mylogfun(isfree,phatfix)
% %     %  phatv = [ki,sc];
% %     mphat =  phatv;
% %     mphat(i_fixed) = phatfix; 
% %     mphat(isfree) = phatfree;
% %     L = logfun(mphat,data1,dist);
% %   end % function mylogfun
%   function L = mylogfun(phatfree1,fix_p,linkfun1)
%     %Loglikelihood given a fixed distribution parameter or return value (quantile)
%     %
%     % L = myloglikRL(thetafree,fix_p,linkfun)
%     % 
%     % L        = negative loglike or logps function
%     % phatfree = distribution parameters that are not fixed.
%     % fix_p    = fixed parameter, i.e. either quantile (return level),
%     %            probability (return period) or
%     %            distribution parameter
%     % linkfun  = linkfunction connecting the fixed return quantile (fix_p) 
%     %            with the fixed distribution paramter, i.e.: 
%     %                phatv(i_fixed) = linkfun(fix_p,phatv,prb,i_fixed)
%     %            if empty linkfun, it is assumed that fix_p is equal to
%     %            phatv(i_fixed)
%     
%     mphat =  phatv;
%     mphat(isfree) = phatfree1;
%     if isempty(linkfun1) % fix_p = distribution parameter
%       mphat(i_fixed) = fix_p;
%     elseif isProfileX    % fix_p = quantile 
%       mphat(i_fixed) = feval(linkfun1,fix_p,prb,mphat,i_fixed);
%     else                 % fix_p = probability
%       mphat(i_fixed) = feval(linkfun1,qntl,fix_p,mphat,i_fixed);
%     end
%     L = logfun(mphat,data1,pdfoptions{:},dist);
%   end % function mylogfun
%   function x1 = myinvfun(phatnotfixed)
%     mphat = phatv;
%     mphat(isnotfixed) = phatnotfixed;
%     cphat = num2cell(mphat,1);
%     x1 = feval(invfun,prb,cphat{:} ,optsfun{:});
%   end %myinvfun
%   function prb1 = myprbfun(phatnotfixed)
%     mphat = phatv;
%     mphat(isnotfixed) = phatnotfixed;
%     cphat = num2cell(mphat,1);
%     prb1 = feval(prbfun,qntl,cphat{:} ,optsfun{:});
%   end %myinvfun
end % main
