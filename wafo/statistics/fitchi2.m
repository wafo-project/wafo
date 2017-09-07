function [phat]=fitchi2(data,varargin)
%FITCHI2 Parameter estimates for Chi squared data.
%
% CALL: phat = fitchi2(data, options)
%
%     phat = Struct with estimated parameters 
%     data = one-dimensional data set
%  options = struct with fieldnames
%     .method : 'ml'  : Maximum Likelihood method (default)
%               'mps' : Maximum Product of Spacing method
%               'mom  : Moments method
%     .plotflag : 1, plot the empiricial distribution
%                   function and the estimated cdf 
%                 0, do not plot
%     .alpha    : Confidence coefficent             (default 0.05)
%     .optimset : optimset structure defining performance of the
%                 optimization routine (see optimset for details)
%          
% Example:
%   R=rndchi2(2,1,100);
%   phat = fitchi2(R);
%   plotfitsumry(phat);
%
%   close all;
%
% See also  pdfchi2, cdfchi2, invchi2, rndchi2, fitchi2, momchi2

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

% No Reference
%  

% tested on: matlab 5.3
% History:
% By  pab 24.10.2000
% revised PJ 03-Apr-2001
%  - fmins changed name to fminsearch for version >= 5.3
% revised pab sept 2007
% -fixed a bug: disable was not passed correctly to pdfchi2

global WAFO_WSTATS_DEFAULT_PLOTFLAG
%error(nargchk(1,inf,nargin))
narginchk(1,inf)
options = struct('method','ML','alpha',0.05,...
  'plotflag', WAFO_WSTATS_DEFAULT_PLOTFLAG,...
  'optimset',optimset('disp','off')); % default options
if (nargin==1 && nargout <= 1 && isequal(data,'defaults'))
  phat = options; 
  return
end



data = data(:); % make sure it is a vector

mu   = mean(data);
sa   = std(data)^2;

% Supply a starting guess with method of moments:
phat0 = max(round((mu+sa/2)/2),1);
ismom = strcmpi(options.method,'ml');
if  strcmpi(options.method,'ml') || strcmpi(options.method,'mps') || ismom,  % Maximum Likelihood
  if  ~ismom
    plotflag = options.plotflag;
    options.plotflag = 0;
    phat =  mlest(@pdfchi2,phat0,data,options,'disable',true);
  
    phat0 = max(round(phat.params),1);
    options.plotflag = plotflag;
  end
  phat = mlest(@pdfchi2,phat0,data,{options,'search',false},'disable',true);
  phat.dataname = inputname(1);
  phat.pdfoptions = struct('disable',true);
else
  error(['Unknown method ' options.method '.']);
end 
return
% disable = 1;
% phat1 = round(fminsearch(@loglike,phat0,options.optimset,data,'disable',disable,@pdfchi2));
% 
% 
% 
% phat1 = [phat1-1 phat1 phat1+1] +(phat1==1);
% Np = length(phat1);
% LL = zeros(1,Np);
% for ix=1:Np
%   LL(ix) = loglike(phat1(ix),data,@pdfchi2);
% end
% [Y,ind] = min(LL);
% phat = phat(ind);
% 
% 
% if nargout>1,
%   [LL, pcov] = loglike(phat,data,@pdfchi2);
% end
% if nargout>2
%   alpha2 = ones(1,1)*0.05/2;
%   var = pcov;
%   pci = invnorm([alpha2;1-alpha2], [phat;phat],[var;var]);
% end
% 
% if plotflag 
%   fitsummaryplot(data,'pdfchi2',phat,plotflag,'MLE')
% end
