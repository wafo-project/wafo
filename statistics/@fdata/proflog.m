function [varargout]=proflog(self,varargin)
%PROFLOG Profile Log- likelihood or product spacing function.
% 
% CALL:  [Lp,CI,nphat] = proflog(self,options);
%
%     Lp = Profile log-likelihood function with parameters phat given
%          the data, phat(i) and quantile x (if given), i.e., 
%            Lp = max(log(f(phat|data,phat(i)))),
%          or
%            Lp = max(log(f(phat|data,phat(i),x))) 
%     CI = Confidence interval for either phat(i) or quantile x.
% nphats = struct with updated distribution parameters.
%   self = FDATA object with fields
%     .data         = data vector
%     .distribution = string containing the name of the PDF or function handle.
%     .params       = distribution parameters (default estimated from data)
%     .dist         = function handle or string defining the PDF to use.
%    options = struct with fieldnames:
%         .i           : integer defining which distribution parameter to 
%                        profile, i.e. which parameter to keep fixed (default 1)
%         .pmin, .pmax : Interval for either the parameter, phat(i) or x, used in the
%                        optimization of the profile function (default is
%                        based on the 100*(1-alpha)% confidence interval
%                        computed using the delta method.)
%         .N           : Number of points used in Lp (default 50)
%         .x           : quantile (return value)
%         .link        : function connecting the quantile (x) with the
%                        fixed distribution paramter, i.e.:
%                        phat(i) = linkfun(x,phat,R,i), where R = Prob(X>x;phat) is 
%                        the survival probability. This means that if
%                        linkfun:
%                         1) is not empty, x is profiled.
%                         2) is empty, phat(i) is profiled (default)
%         .alpha       : confidence coefficent (default 0.05)
%         .plotflag    :     
%         .optimset    : optimset structure defining performance of the
%                        optimization routine (see optimset for details)
%             
%  PROFLOG is a utility function for making inferences on a particular
%  component phat(i) of the vector phat or the quantile, x. This is usually 
%  more accurate than using the delta method and the asymptotic normality 
%  of the ML estimator or the MPS estimator.
%  This works on any PDF and/or CDF having the following calling syntax:
%
%       f = pdftest(x,A1,A2,..,Am,opt1)
%       F = prbtest(x,A1,A2,..,Am,opt2)
%       x = invtest(F,A1,A2,..,Am,opt2)
%
% where x is the quantile and A1,A2... are the distribution parameters. 
% opt1 and opt2 are options structs which must have the fields 'logp' and 'lowertail',
% respectively
%
% Examples % MLE and better CI for phat.params(1)
%   R = rndweib(1,3,0,100,1);                      
%   phat = fitweib(R)         
%   opt = proflog('defaults');
%   opt.i = 2;opt.plotflag=1;
%   [Lp,CI] = proflog(phat,opt)
%
% % MPS and better confidence interval for the quantile, x.
%   R2 = rndgenpar(0.3,1,0,200,1);
%   m = 1;
%   [phat2] = fitgenpar(R2(R2>m)-m,'method','mps');
%   x = invgenpar(1/9900,phat2,'lowertail',false);
% % Reorganizing w.r.t. phat(2) = scale, Eq. 4.13 and 4.14, pp 81 in Coles (2004) gives
%   link = @(x,phat,prb,ix) -(x-phat(3)).*phat(1)./expm1(phat(1).*log(prb)); 
%   [Lp2,CI2] = proflog(phat2,'i',2,'x',x,'link',link,'plot',1);
%
% See also loglike, logps, optimset



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
% Stuart Coles (2004)
% "An introduction to statistical modelling of extreme values".
% Springer series in statistics

%Tested on: matlab 7.3
% History:
% by pab Sep 2007

[varargout{1:nargout}] = proflog(struct(self),varargin{:});
return
