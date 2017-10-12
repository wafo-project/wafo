function [phat]=fitgam(data,varargin)
%FITGAM Parameter estimates for Gamma data.
%
% CALL: phat = fitgam(data, options)
%
%     phat = Struct with estimated parameters 
%     data = one-dimensional data set
%  options = struct with fieldnames
%     .method : a string, describing the method of estimation
%                'ml'  = Maximum Likelihood method (default)
%                'mps' = Maximum Product Spacing
%     .plotflag : 1, plot the empiricial distribution
%                   function and the estimated cdf 
%                 0, do not plot
%     .alpha    : Confidence coefficent             (default 0.05)
%     .optimset : optimset structure defining performance of the
%                 optimization routine (see optimset for details)
%          
% Example:
%   R = rndgam(5,1,1,100);
%   phat = fitgam(R);
%   plotfitsumry(phat);
%
%   close all;
%
% See also  pdfgam, cdfgam, invgam, rndgam, momgam

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


% Reference:
% Bowman & Shenton, "Properties of estimators for the Gamma distribution"
% Marcel Dekker

% tested on: matlab 5.3
% History:
% revised pab Aug 2007
% - added wgamafit as an subfunction
% revised pab 21.01.2004
% revised pab 24.10.2000
% - added check on fzero in order to run on matlab 5.2
% - added pci
% added ms 23.08.2000

global WAFO_WSTATS_DEFAULT_PLOTFLAG

%error(nargchk(1,inf,nargin))
narginchk(1,inf)
% Add these options?: 'shape',nan,'scale',nan,'location',0, 
options = struct('method','ML','fixpar',[],'alpha',0.05,...
  'plotflag', WAFO_WSTATS_DEFAULT_PLOTFLAG,...
  'optimset',optimset('disp','off')); % default options
if (nargin==1 && nargout <= 1 && isequal(data,'defaults'))
  phat = options; 
  return
end
options        = parseoptions(options,varargin{:});
options.method = upper(options.method);


if any(data<=0)
  error('data must be strictly positive!')
end

data=data(:); % make sure it is a vector

meanData    = mean(data);
logData     = log(data);
meanLogData = mean(logData);



G = log(meanData)-meanLogData;
% Error of ahat0 within 1.5% from the true value
ahat0  = (3-G+sqrt((G-3).^2+24*G))./(12*G);
%ahat1 = -(meanLogData-mean(data.*logData)/meanData)^(-1);
anyparfixed = ~isempty(options.fixpar) && any(isfinite(options.fixpar));

if strcmpi(options.method,'ml'),  % Maximum Likelihood

  dosearch = ~anyparfixed;
  %start = 1./(2*G);
  start = ahat0;
  %start = ahat1;
  
  ahat = fzero(@(x)localfitgam(x,G),start,options.optimset);
  
elseif  strcmpi(options.method,'mps'),  % Maximum product spacing
  dosearch = true;
  ahat = ahat0;
else
  error(['Unknown method ' options.method '.']);
  %ahat = ahat0;  
end

bhat = meanData/ahat;
phat0=[ahat, bhat];
if anyparfixed
  phat0 = phat0(~isfinite(options.fixpar));
end

phat = mlest(@pdfgam,phat0,data,{options,'search',dosearch});
phat.dataname = inputname(1);
return
% 
% if nargout>1,
%   [LL, cov] = loglike(phat,data,@pdfgam);
% end
% if nargout>2
%   alpha2 = ones(1,2)*0.05/2;
%   var = diag(cov).';
%   pci = invnorm([alpha2;1-alpha2], [phat;phat],[var;var]);
% end
% 
% if plotflag 
%   fitsummaryplot(data,'pdfgam',phat,plotflag)
% end

function L=localfitgam(a,G)
%FITGAM Is an internal routine for fitgam
%


%h=10^(-5);
%L=log(a)-(gammaln(a+h)-gammaln(a))/h-G;
% wpsi = digamma function
L=log(a)-wpsi(a)-G;

