function [phat2]=fitt(data,varargin)
%FITT Parameter estimates for Student's T data.
%
% CALL:  phat = fitt(data, options)
%
%      phat = Struct with estimated parameters 
%   data    = Onedimensional dataset
%   options = struct with fieldnames
%     .method : a string, describing the method of estimation
%               'ml'  = Maximum Likelihood method (default)
%               'mps' = Maximum Product Spacing method 
%     .plotflag : 1, plot the empiricial distribution
%                   function and the estimated cdf 
%                 0, do not plot
%     .alpha    : Confidence coefficent             (default 0.05)
%     .search   : true (default) 
%     .optimset : optimset structure defining performance of the
%                 optimization routine (see optimset for details)
%          
% Example:
%   R=rndt(2,1,100);
%   phat = fitt(R) 
%   plotfitsumry(phat)
%
% See also  cdft, plotfitsumry

% No Reference
%  

% tested on: matlab 5.3
% History:
% By  pab 24.10.2000
% revised PJ 03-Apr-2001
%  - fmins changed name to fminsearch for version >= 5.3

global WAFO_WSTATS_DEFAULT_PLOTFLAG
error(nargchk(1,inf,nargin))
% Add these options?: 'shape',nan,'scale',nan,'location',0, 
options = struct('method','ML','alpha',0.05,...
 'plotflag', WAFO_WSTATS_DEFAULT_PLOTFLAG,'optimset',optimset); % default options
if (nargin==1 && nargout <= 1 && isequal(data,'defaults'))
  phat2 = options; 
  return
end
options        = parseoptions(options,varargin{:});
options.method = upper(options.method);

data=data(:); % make sure it is a vector

sa = std(data)^2;

% Supply a starting guess with method of moments:
phat0 = max(round(2*sa/(sa-1)),1);

pdfoptions = struct('disable',true);
if strcmpi(options.method,'ml'),  % Maximum Likelihood
  dist = @pdft;
  likfun = @loglike;
elseif strcmpi(options.method,'mps')  %Maximum product spacing
 dist  = @cdft;
 likfun = @logps;
else
  error(['Unknown method ' options.method '.']);
end



phat = round(fminsearch(likfun,phat0,options.optimset,data,pdfoptions,dist));

%  phat = round(fmins('loglike',phat0,[],[],data,1,'pdft'));


phat = [phat-1 phat phat+1] +(phat==1);
Np = length(phat);
LL = zeros(1,Np);

for ix=1:Np
  LL(ix) = likfun(phat(ix),data,dist);
end
[Y,ind] = min(LL);
phat = phat(ind);
% 


[LL, pcov] = loglike(phat,data,'pdft');
[LPS, pvalue] = logps(phat,data,'cdft');

zcrit = -invnorm(options.alpha/2);
ciL   = phat-zcrit*sqrt(pcov);
ciU   = phat+zcrit*sqrt(pcov);


phat2 = createfdata(options,'dist',dist,...
  'params',phat,'lower',ciL,'upper',ciU,...
  'covar',pcov,'var',pcov,...
  'dataname',inputname(1),'data',data,...
  'pdfoptions',pdfoptions,...
  'loglikemax', -LL,'logps',-LPS,'pvalue',pvalue);


if options.plotflag 
  plotfitsumry(phat,options.plotflag)
end
