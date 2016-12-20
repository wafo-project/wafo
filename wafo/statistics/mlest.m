function phat = mlest(dist,phat0,data,options1,varargin)
% MLEST Maximum Likelihood or Maximum Product Spacing estimator
%
% CALL phat = mlest(dist,phat0,data,options,pdfoptions)
%
%      phat = Struct with estimated parameters 
%      dist = pdf function (string or function handle)
%     phat0 = vector of initial phat-values for non-fixed parameters.
%      data = one-dimensional data set
%pdfoptions = list of options passed directly to pdf-function
%   options = struct with fieldnames
%     .method   : a string, describing the method of estimation
%               'ml'  = Maximum Likelihood method (default)
%               'mps' = Maximum Product Spacing method 
%     .fixpar   : vector giving the fixed parameters. (Must be empty or 
%                 have the same length as the number of parameters. 
%                 Non fixed parameters must be given as NaN's)
%     .plotflag : 1, plot the empiricial distribution
%                   function and the estimated cdf 
%                 0, do not plot
%     .alpha    : Confidence coefficent             (default 0.05)
%     .search   : true (default) 
%     .optimset : optimset structure defining performance of the
%                 optimization routine (see optimset for details)
% Example
%   R = rndbeta(2,0.5,1,100);
%   phat0 = [2,2];
%   phat =  mlest(@pdfbeta,phat0,R);
%   plotfitsumry(phat);
%   fixpar = [2,nan];
%   phat2 = mlest(@pdfbeta,phat0(2),R,{'fixpar',fixpar});
%   plotfitsumry(phat2);
%
%   close all;
%
% See also loglike, logps, fminsearch, optimset

% Copyright (C) 2007 Per A. Brodtkorb
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

global WAFO_WSTATS_DEFAULT_PLOTFLAG

options = struct('method','ML','fixpar',[],'search',true,'alpha',0.05,...
  'plotflag', WAFO_WSTATS_DEFAULT_PLOTFLAG,...
  'optimset',optimset('disp','off')); % default options
if (nargin==1 && nargout <= 1 && isequal(dist,'defaults'))
  phat = options; 
  return
end
error(nargchk(3,inf,nargin))
if nargin>3 && ~isempty(options1)
  if iscell(options1)
    options        = parseoptions(options,options1{:});
  else
    options        = parseoptions(options,options1);
  end
end

msgId = 'WAFO:PARSEOPTIONS';
if ~isoctave,
  warnstate = warning('query',msgId);
  warning('off',msgId);
end
if isnumeric(data)
  data = {data(:)};
else
  d = size(data{1});
  if any(d(2:end)>1)
    for ix=1:length(data)
      data{ix}  = data{ix}(:); %% make sure it is a vector.
    end
  end
end
allfixed  = isempty(phat0);
somefixed = ~(isempty(options.fixpar) || all(isnan(options.fixpar)));
if somefixed
  phat2 = options.fixpar;
  notfixed = find(~isfinite(phat2));
  i_fixed  = find(isfinite(phat2));
  if length(phat0)~=length(notfixed)
    error('WAFO:MLEST','The length of PHAT0 must equal the number of non-fixed parameters given in options.fixpar!')
  end
end

model = getdistname(dist);
pdf = ['pdf' model ];
cdf = ['cdf' model ];
if options.search && ~allfixed
  if strcmpi(options.method,'ML')
    pdfOrCdf = pdf;
    if somefixed
      logfun = @(x)fxloglike(x, phat2, notfixed, data,varargin{:},pdfOrCdf);
    else
      logfun = @(x)loglike(x, data,varargin{:},pdfOrCdf);
    end
    
  elseif strcmpi(options.method,'MPS')
    data = num2cell(sortrows([data{:}]),1);
    pdfOrCdf = cdf;
    if somefixed
      logfun = @(x)fxlogps(x, phat2, notfixed, data,varargin{:},pdfOrCdf);
    else
      logfun = @(x)logps(x,data,varargin{:},pdfOrCdf);
    end
    
  else
    error('Unknown method: %s',options.method)
  end
  if isoctave,
    phat1 = fminsearch(logfun,phat0,options.optimset);
  else
    [phat1 , dummy,Converged]= fminsearch(logfun,phat0,options.optimset);
    
    if ~Converged 
       warning('WAFO:MLEST','Minimization-routine did not converge!'); 
    end
  end
  if somefixed
      phat2(notfixed) = phat1;
      phat1 = phat2;
   end
elseif ~allfixed && somefixed
  phat2(notfixed) = phat0;
  phat1 = phat2;
else
  phat1 = phat0;
end
[LL, pcov,H]  = loglike(phat1,data,varargin{:},pdf);
try
    if strcmpi(options.method,'MPS')
        [LPS, pvalue, Tn, H] = logps(phat1,data,varargin{:},cdf);
        pcov = -pinv(H);
    else
        [LPS, pvalue] = logps(phat1,data,varargin{:},cdf);
    end
catch
  LPS = nan;
  pvalue = nan;
end
if somefixed
  pcov(:,i_fixed) = 0;
  pcov(i_fixed,:) = 0;
  pcov(notfixed,notfixed) = -pinv(H(notfixed,notfixed));
end


pvar = diag(pcov).';
zcrit = -invnorm(options.alpha/2);
ciL   = phat1-zcrit*sqrt(pvar);
ciU   = phat1+zcrit*sqrt(pvar);

if numel(data)==1
  data = data{1};
end
phat = createfdata(options,'dist',dist,...
  'params',phat1,'lower',ciL,'upper',ciU,...
  'covar',pcov,'var',pvar,...
  'dataname',inputname(3),'data',data,...
  'loglikemax', -LL,'logps',-LPS,'pvalue',pvalue);


if options.plotflag 
  plotfitsumry(phat,options.plotflag)
end
phat = fdata(phat);
if ~isoctave
  warning(warnstate);
end

end % mlest


function L = fxloglike(phat10, phat2, notfixed, varargin)
  %notfixed = find(isnan(options.fixpar));
  phat2(notfixed) = phat10;
  L = loglike(phat2,varargin{:});
end
function L = fxlogps(phat10, phat2, notfixed, varargin)
  %notfixed = find(isnan(options.fixpar));
  phat2(notfixed) = phat10;
  L = logps(phat2,varargin{:});
end