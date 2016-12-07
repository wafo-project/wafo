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
%   phat = fitweib(R);
%   opt = proflog('defaults');
%   opt.i = 2;opt.plotflag=1;
%   [Lp,CI] = proflog(phat,opt)
%   [x,xlo,xup] = invweib(1/990,phat,'lowertail',false);
%   [Lp0,CI0] = proflog(phat,'i',2,'x',x,'link',@lnkweib,'plotflag',1);
%
% % MPS and better confidence interval for the quantile, x.
%   R2 = rndgenpar(0.3,1,0,200,1);
%   m = 1;
%   [phat2] = fitgenpar(R2(R2>m)-m,'method','mps');
%   [x,xlo,xup] = invgenpar(1/9900,phat2,'lowertail',false);
%   [Lp2,CI2] = proflog(phat2,'i',2,'x',x,'link',@lnkgenpar,'plotflag',1);
%   [Lp3,CI3] = proflog(phat2,'i',2,'logR',-log(9900),'link',@lnkgenpar,'plotflag',1);
%
%   close all;
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
end % main
