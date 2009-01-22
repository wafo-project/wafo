function phat = fitraymod(data,varargin)
%FITRAYMOD Parameter estimates for Truncated Rayleigh data.
%
% CALL:  phat = fitraymod(data, options)
%
%     data = one-dimensional data set
%  options = struct with fieldnames
%     .method : a string, describing the method of estimation
%                'ml' Maximum Likelihood method (default)
%               'mps' Maximum product of spacings method
%     .plotflag : 1, plot the empiricial distribution
%                   function and the estimated cdf 
%                 0, do not plot
%     .alpha    : Confidence coefficent             (default 0.05)
%     .optimset : optimset structure defining performance of the
%                 optimization routine (see optimset for details)
%          
% Example:
%   R=rndraymod(2,2,1,100);
%   phat = fitraymod(R)
%   plotfitsumry(phat)
%
% See also  pdfraymod, cdfraymod, rndraymod, invraymod, momraymod


% Reference: Cohen & Whittle, (1988) "Parameter Estimation in Reliability
% and Life Span Models", p. 181 ff, Marcel Dekker.

%tested on: matlab 5.x
% History:
%  by Per A. Brodtkorb 17.10.98
% revised ms 15.06.2000
% - updated header info
% - changed name to fitray (from raylfit)
% revised ms 11.08.2000
% - changed to standard *fit form
% -revised pab 24.10.2000
%  - replaced gamma with gammaln -> more robust
%  - added nargchk
% revised PJ 03-Apr-2001
%  - fmins changed name to fminsearch for version >= 5.3
% Revised pab Dec2003

global WAFO_WSTATS_DEFAULT_PLOTFLAG
error(nargchk(1,inf,nargin))
% Add these options?: 'shape',nan,'scale',nan,'location',0, 
options = struct('method','ML','fixpar',[],'alpha',0.05,...
 'plotflag', WAFO_WSTATS_DEFAULT_PLOTFLAG,'optimset',optimset); % default options
if (nargin==1 && nargout <= 1 && isequal(data,'defaults'))
  phat = options; 
  return
end
options        = parseoptions(options,varargin{:});
options.method = upper(options.method);

data = data(:);
n = length(data);
phat0 = sqrt(sum(data.^2)/n/2); % Initial guess (MLE with c=0)

phat0 = [phat0 0];
if ~isempty(options.fixpar)
  phat0 = phat0(~isfinite(options.fixpar));
end

if strcmpi(options.method,'ml') || strcmpi(options.method,'mps'),  % Maximum Likelihood
  pdf = @pdfraymod;
  phat = mlest(pdf,phat0,data,options);
  phat.dataname = inputname(1);
else
  error(['Unknown method ' options.method '.']);
end




