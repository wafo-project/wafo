function [phat] = fitgumb(data,varargin)
%FITGUMB Parameter estimates for Gumbel data.
%
% CALL: [phat] = fitgumb(data,plotflag) 
%
%     phat = the maximum likelihood estimates of the  
%            parameters of the Gumbel distribution given the data.
%     cov  = asymptotic covariance matrix of estimates
%     data = data vector
% plotflag = 0, do not plot
%          > 0, plot the empiricial distribution function and the
%               estimated cdf (see plotfitsumry for options)(default)
%
% Example:
%   R = rndgumb(1,2,100,1,'trunc',1);
%   phat = fitgumb(R);
%   plotfitsumry(phat);
%
%   close all;
%
% See also  pdfgumb, cdfgumb, invgumb, rndgumb, momgumb plotgumb

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
%  Johnson  N.L., Kotz S. and Balakrishnan, N. (1994)
%  Continuous Univariate Distributions, Volume 2. Wiley. 


% tested on: matlab 5.3
% rewritten ms 05.07.2000
% revised jr 01.09.2000
% - ahat and bhat were reversed in the covariance matrix.
% -revised pab added nargchk, pci + safer call to fzero  
% - made sure data is  vector
% revised PJ 02-Apr-2001
%   Fixed problem with no or empty plotflag

global WAFO_WSTATS_DEFAULT_PLOTFLAG

%error(nargchk(1,inf,nargin))
narginchk(1,inf)
% Add these options?: 'shape',nan,'scale',nan,'location',0, 
options = struct('method','ML','alpha',0.05,...
  'plotflag', WAFO_WSTATS_DEFAULT_PLOTFLAG,...
  'trunc',false, 'start',[],... 
  'optimset',optimset('disp','off','TolX',1e-5,'TolFun',1e-5,'MaxIter',500)); % default options
if (nargin==1 && nargout <= 1 && isequal(data,'defaults'))
  phat = options; 
  return
end
options        = parseoptions(options,varargin{:});
options.method = upper(options.method);

data=sort(data(:)); % make sure it is a vector

if ~isempty(options.start)
  start = options.start;
else
  start=6^(1/2)/pi*std(data); % Moment estimate of scale parameter a
end

ahat=fzero(@(p)myfitgumb(p,data),start,options.optimset);

bhat=-ahat*log(mean(exp(-data/ahat)));
phat0=[ahat, bhat];

if strcmpi(options.method,'ml') || ...
    strcmpi(options.method,'mps')  % Maximum Likelihood
  if options.trunc || ...
    strcmpi(options.method,'mps')
    phat = mlest('pdfgumb',phat0,data,options,options);
    phat.dataname = inputname(1);
    phat.pdfoptions = struct('trunc',options.trunc);
  else
    pcov=[0.60793,0.25696;0.25696,1.10867]*ahat^2/length(data);
    
     zcrit = -invnorm(options.alpha/2);
     pvar = diag(pcov).';
     ciL = phat0-sqrt(pvar).*zcrit;
     ciU = phat0+sqrt(pvar).*zcrit;
     
     [LPS,pvalue] = logps(phat0,data,'cdfgumb');
     phat = createfdata(options,'dist','pdfgumb',...
    'params',phat0,'lower',ciL,'upper',ciU,...
    'covar',pcov,'var',diag(pcov).',...
    'dataname',inputname(1),'data',data,...
    'loglikemax', -loglike(phat0,data,'pdfgumb'),...
    'logpsmax',-LPS,'pvalue',pvalue,...
    'pdfoptions',struct('trunc',options.trunc), 'note',sprintf('Moran''s statistic on fit: pvalue = %g',pvalue));

    phat = fdata(phat);

    if options.plotflag
      plotfitsumry(phat,options.plotflag)
    end
  end
else
  error(['Unknown method ' options.method '.']);
  %ahat = ahat0;  
end



function L=myfitgumb(a,data)
%MYFITGUMB Is an internal routine for fitgumb
%

L=a-mean(data)+mean(data.*exp(-data/a))/mean(exp(-data/a));


