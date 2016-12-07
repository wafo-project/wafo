function [phat2] = fitray(data,varargin)
%FITRAY Parameter estimates for Rayleigh data.
%
% CALL:  [phat] = fitray(data, plotflag)
%
%     phat = Struct with estimated parameters 
%     data = one- or multi-dimensional data set
%  options = struct with fieldnames
%     .method : 'ml'  Maximum Likelihood method (default)
%               'mps' Maximum product of spacings method
%     .plotflag : 1, plot the empiricial distribution
%                   function and the estimated cdf 
%                 0, do not plot
%     .alpha    : Confidence coefficent             (default 0.05)
%     .optimset : optimset structure defining performance of the
%                 optimization routine (see optimset for details)
%
% Example:
%   R=rndray(2,100,2);
%   phat=fitray(R);
%   plotfitsumry(phat);
%
%   close all;
%
% See also  pdfray, cdfray, invray, rndray, momray

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
%method         = options.method;

sz = size(data);
Nsz=length(sz);
dim = find(sz~=1, 1 );  %1st non-singleton dimension
% make sure dim=1 is the first non-singleton dimension
if isempty(dim) || dim ~= 1, 
  order = [dim 1:dim-1 dim+1:Nsz];
  data  = permute(data,order);
  sz    = size(data);
end
m = prod(sz(2:end));
n =sz(1);

phat=sqrt(sum(data.^2)/n/2).'; % MLE

if strcmpi(options.method,'ml'),  % Maximum Likelihood
 
elseif strcmpi(options.method,'mps')  %Maximum product spacing
  for ix = 1:m
    phat(ix,:) = fminsearch(@logps,phat(ix,:),options.optimset,data(:,ix),@cdfray);
  end
else
  error(['Unknown method ' options.method '.']);
end



if (n<200),
  pvar=phat.^2.*(n-(exp(gammaln(1/2+n)-gammaln(n))).^2)/n;
else    %n too large for evaluating gamma
  pvar=phat.^2.*0.25/n;
end


zcrit = -invnorm(options.alpha/2);

ciL = phat-sqrt(pvar).*zcrit;
ciU = phat+sqrt(pvar).*zcrit;

phat2 = createfdata;
phat2.distribution = 'pdfray';
phat2.alpha = options.alpha;
phat2.method = upper(options.method);
phat2.dataname = inputname(1);
if m>1 % expand
  [phat2(1:m) ] = deal(phat2);
end
phat = num2cell(phat,2);
ciL = num2cell(ciL,2);
ciU = num2cell(ciU,2);
pvar = num2cell(pvar,2);
[phat2.params] = deal(phat{:});
[phat2.lowerbound] = deal(ciL{:});
[phat2.upperbound] = deal(ciU{:});
[phat2.variance] = deal(pvar{:});  
  

for ix = 1:m
  [phat2(ix).logpsmax,phat2(ix).pvalue] = logps(phat2(ix),data(:,ix),'cdfray');
  phat2(ix).covariance = diag(phat2(ix).variance);
  phat2(ix).data = data(:,ix);
  phat2(ix).loglikemax = -loglike(phat2(ix).params,data(:,ix),@pdfray);
end

if options.plotflag 
  plotfitsumry(phat2,options.plotflag)
end
if ~isoctave
  phat2 = fdata(phat2);
end


