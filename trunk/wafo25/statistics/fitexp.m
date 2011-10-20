function [phat2] = fitexp(data,varargin)
%FITEXP Parameter estimates for Exponential data.
%
% CALL: phat = fitexp(data, options)
%
%     phat = Struct with estimated parameters 
%     data = one- or multi-dimensional data set
%  options = struct with fieldnames
%     .method   : 'ml'  Maximum Likelihood method (default)
%               : 'mps' Maximum Product of Spacing method
%     .plotflag : 1, plot the empiricial distribution
%                   function and the estimated cdf 
%                 0, do not plot
%     .alpha    : Confidence coefficent             (default 0.05)
%
% Example:
%   R=rndexp(2,100,1);
%   [mhat]=fitexp(R,'plotflag',1)
%   R=rndexp(2,100,3);
%   [mhat2]=fitexp(R)
%   plotfitsumry(mhat2,3)
%
% See also  pdfexp, cdfexp, invexp, rndexp, momexp

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



% Reference: Johnson, Kotz and Balakrishnan (1994)
% "Continuous Univariate Distributions, vol. 1", p. 494 ff
% Wiley


%tested on: matlab 5.x
% History:
% revised pab 
% -added alpha to input
% -fixed a bug for CI for phat.
% revised pab 24.10.2000
% - added  nargchk + 95% CI for phat
% - fixed some bugs when data is a matrix 
% added ms 16.08.2000

global WAFO_WSTATS_DEFAULT_PLOTFLAG
error(nargchk(1,inf,nargin))
% Add these options?: 'shape',nan,'scale',nan,'location',0, 
options = struct('method','ML','alpha',0.05,...
 'plotflag', WAFO_WSTATS_DEFAULT_PLOTFLAG,...
 'optimset',optimset('disp','off','TolX',1e-5,'TolFun',1e-5,'MaxIter',500)); % default options
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


if ~isscalar(data) && (isempty(dim) || dim ~= 1), 
  order = [dim 1:dim-1 dim+1:Nsz];
  data  = permute(data,order);
  sz    = size(data);
end
m = prod(sz(2:end));
n =sz(1);
phat=mean(data);

if  strcmpi(options.method,'ml'),  % Maximum Likelihood
  

elseif  strcmpi(options.method,'mps')  %Maximum product spacing
  for ix = 1:m
    phat(ix,:) = fminsearch(@logps,phat(ix,:),options.optimset,data(:,ix),@cdfexp);
  end
else
   error(['Unknown method ' options.method '.']);
end

pvar=phat.^2/n;
% phat ~ 1/gamma(n,1/(phat*n))
alpha2  = options.alpha/2;
gamcrit = invgam([alpha2 1-alpha2],n,1);
ciU = n*phat/gamcrit(1);
ciL = n*phat/gamcrit(2);

phat2 = createfdata;
phat2.distribution = 'pdfexp';
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
  [phat2(ix).logpsmax,phat2(ix).pvalue] = logps(phat2(ix),data(:,ix),'cdfexp');
  phat2(ix).logpsmax = -phat2(ix).logpsmax;
  phat2(ix).covariance = diag(phat2(ix).variance);
  phat2(ix).data = data(:,ix);
  phat2(ix).loglikemax = -loglike(phat2(ix).params,data(:,ix),@pdfexp);
end


if options.plotflag 
  plotfitsumry(phat2,options.plotflag)
end
if ~isoctave
  phat2 = fdata(phat2);
end