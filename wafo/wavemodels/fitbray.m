function [phat]=fitbray(data1,varargin)
%FITBRAY Parameter estimates for Beta-Rayleigh data. 
%
% CALL: phat = fitbray(data,alpha);
%
%   phat  = [a, b, c] = maximum likelihood estimates of the
%           parameters of the Beta-Rayleigh distribution (see
%           braylpdf) given the data.
%   cov   = asymptotic covariance matrix of estimates
%   pci   = 100*(1-alpha) percent confidense intervals
%   data  = data matrix
%   alpha = confidence level (default 0.05 corresponding to 95% CI)
%
% Example:
%  a = .9; b = 105; sz = [100,1]
%  R = sort(rndbeta(a,b,sz));
%  phat = fitbray(R,'fixpar',[nan 105,nan]);
%  plotedf(R,[R pdfbray(R,phat)]);
%
%  close all;
%
% See also  pdfbray, cdfbray, fitbeta

% tested on: matlab 5.2
%History:

% revised pabnov 2004
% -replaced fmins with fminsearch  
% by Per A. Brodtkorb 14.02.99
%   Reference:

global WAFO_WSTATS_DEFAULT_PLOTFLAG
%error(nargchk(1,inf,nargin))
narginchk(1,inf)
% Add these options?: 'shape',nan,'scale',nan,'location',0, 
options = struct('method','ML','fixpar',[],'alpha',0.05,...
 'plotflag', WAFO_WSTATS_DEFAULT_PLOTFLAG,'optimset',optimset); % default options
if (nargin==1 && nargout <= 1 && isequal(data1,'defaults'))
  phat = options; 
  return
end
options        = parseoptions(options,varargin{:});
options.method = upper(options.method);
%method         = options.method;

error(nargchk(1,inf,nargin))

data1=data1(:);

c=sqrt(2)*max(data1);
phat0 = fitbeta((data1./c).^2);
pinit=[phat0.params c];
if ~isempty(options.fixpar)
  pinit(isfinite(options.fixpar)) = [];
end

phat = mlest(@pdfbray,pinit,data1,options);
phat.dataname = inputname(1);

