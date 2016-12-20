function [phat] = fitgev(data,varargin) 
%FITGEV Parameter estimates for GEV data
%
% CALL:  [phat] = fitgev(data,options)
%
%        phat  =  Struct with estimated parameters 
%                [k s m] = [shape scale location]   (see cdfgev)
%        data  = a one-dimensional data set
%   options = struct with fieldnames
%    .method   : a string, describing the method of estimation
%                'PWM' = Probability Weighted Moments (default)
%                'ML'  = Maximum Likelihood estimates
%                'MPS' = Maximum Product of Spacing
%     .fixpar   : vector giving the fixed parameters. (Must be empty or 
%                 have the same length as the number of parameters. 
%                 Non-fixed parameters must then be given as NaN's)
%                 (default [nan nan nan])
%    .start    : starting values for 'ML' method
%                (default 'PWM' estimate) 
%    .plotflag = 0, do not plot
%              > 0, plot the empiricial distribution function and the
%                   estimated cdf (see plotfitsumry for options)(default)
%                  
% The variances of the ML estimates are usually smaller than those of
% the PWM estimates. However, it is recommended to first use the PWM
% method since it works for a wider range of parameter values. 
% If method = 'ml' and  start  is missing, then  pwm  is used to 
% give the starting values  start.
%
% Example:  
%   R = rndgev(0.2,2,7.5,200,1);
%   phat= fitgev(R,'method','pwm');
%   phat2 = fitgev(R,'method','ml');
%
% See also  cdfgev

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

% References
%
%  Prescott, P. and Walden, A.T. (1980)
%  Maximum likelihood estimation of the parameters of the generalized
%  extreme-value distribution
%  Biometrika (67), pp. 723-724
%
%  Hosking, J.R.M, Wallis, J.R. and Wood E.F. (1985)
%  Estimation of the generalized extreme-value distribution by the
%  method of probability-weighted moments
%  Technometrics (27), pp. 251-261
%
%  Borg, S. (1992)
%  XS - a statistical program package in Splus for extreme-value
%  analysis. 
%  Dept. of Mathematical Statistics, Lund University, 1992:E2

% Tested on: Matlab 5.3
% History: 
% revised pab 2007
% -replaced fminsearch with a call to mlest
% Revised by pab 13.06.2001
% - Replaced 
% [phat,dummy,Converged] = fminsearch('likgev',mlstart,....);
% with
%  [phat,dummy,Converged] = feval('fminsearch','likgev',mlstart,...);
% to avoid "Unknown function referenced: fminsearch" error on matlab 5.2
% and below. 
%
% Revised by PJ 02-Apr-2001
%   Method 'ml' now works with new format of fminsearch.
% Revised by jr 30.09.1999
% Modified by PJ 08-Mar-2000
%   Added 'hold off' in plotting.
%   Added routine 'gevll'. Now method 'ml' works.
% revised ms 14.06.2000
% - updated header info
% - changed name to fitgev (from gevfit)
% - added w* to used WAFO-files
% - enabled consistent use of 3 and 4 arguments
% - enabled use of empty start in ML
% revised pab 29.10.2000
%  - added nargchk

global WAFO_WSTATS_DEFAULT_PLOTFLAG
error(nargchk(1,inf,nargin))
% Add these options?: 'shape',nan,'scale',nan,'location',0, 
options = struct('method','PWM','fixpar',[],'start',[],'alpha',0.05,...
  'plotflag', WAFO_WSTATS_DEFAULT_PLOTFLAG,...
  'optimset',optimset('disp','off','TolX',1e-5,'TolFun',1e-5,'MaxIter',500)); % default options
if (nargin==1 && nargout <= 1 && isequal(data,'defaults'))
  phat = options; 
  return
end
options        = parseoptions(options,varargin{:});
options.method = upper(options.method);
method         = options.method;


data = data(:); % make sure it is a column vector
somefixed = ~isempty(options.fixpar);

switch lower(method),
  case 'pwm',
    if somefixed
      error('Can not fix parameters when using PWM method')
    end
    [phat] = pwmfit(data,options);
  case {'ml','mps'},
    % ml.gev <- function(data, shape=NA, scale=NA, location=NA) {
    % J. Statist. Comput. Simul. 16, 241-250
    if isempty(options.start)
       % compute starting values by  pwm
      mlstart=pwmfit2(data);% Added ms
      if somefixed
        isnotfixed = (~isfinite(options.fixpar));
        mlstart = mlstart(isnotfixed);
      end
    else
      mlstart=options.start;
    end
    % Solve the ML-equation
    phat = mlest('pdfgev',mlstart,data,options);
  otherwise
    error('Uknown method: %s',method)
end
phat.dataname = inputname(1); 
phat = fdata(phat);

if options.plotflag 
  plotfitsumry(phat,options.plotflag)
end

function [phat] = pwmfit(data,options)

[phat1,pcov] = pwmfit2(data);
[LL]    = loglike(phat1,data,options,@pdfgev);
[LPS, pvalue] = logps(phat1,data,options,@cdfgev);

pvar = diag(pcov).';
zcrit = -invnorm(options.alpha/2);
ciL   = phat1-zcrit*sqrt(pvar);
ciU   = phat1+zcrit*sqrt(pvar);

phat = createfdata(options,'dist','pdfgev',...
  'params',phat1,'lower',ciL,'upper',ciU,...
  'covar',pcov,'var',pvar,...
  'dataname',inputname(1),'data',data,...
  'loglikemax', -LL,'logps',-LPS,'pvalue',pvalue);



function [phat,pcov] = pwmfit2(data)
%PWMFIT
w=[ -0.4, 1.6637,  1.3355, 1.1405, 1.8461, 1.1628, 2.9092;
      -0.3, 1.4153,  0.8912, 0.5640, 1.2574, 0.4442, 1.4090;
      -0.2, 1.3322,  0.6727, 0.3926, 1.0013, 0.2697, 0.9139;
      -0.1, 1.2915,  0.5104, 0.3245, 0.8440, 0.2240, 0.6815;
      0.0, 1.2686,  0.3704, 0.2992, 0.7390, 0.2247, 0.5633;
      0.1, 1.2551,  0.2411, 0.2966, 0.6708, 0.2447, 0.5103;
      0.2, 1.2474,  0.1177, 0.3081, 0.6330, 0.2728, 0.5021;
      0.3, 1.2438, -0.0023, 0.3297, 0.6223, 0.3033, 0.5294;
      0.4, 1.2433, -0.1205, 0.3592, 0.6368, 0.3329, 0.5880];
    
    n = length(data);
    sortdata = sort(data).';
    koeff1 = ((1:n)-1)/(n-1);
    koeff2 = koeff1.*((1:n)-2)/(n-2);
    b2 = (koeff2*sortdata')/n;
    b1 = (koeff1*sortdata')/n;
    b0 = mean(sortdata);
    z = (2*b1-b0)/(3*b2-b0)-log(2)/log(3);
    shape = 7.8590*z+2.9554*z^2;
    scale = (2*b1-b0)*shape/(gamma(1+shape)*(1-2^(-shape)));
    location = b0+scale*(gamma(1+shape)-1)/shape;
    
    phat=[shape scale location];
    % if (abs(shape)>=0.5) changed by ms to deal with shapes in (.4,.5)
    if (abs(shape)>=0.4)
      warning('WAFO:FITGEV',' The estimate of shape is not within the range where\n the PWM estimator is valid.')
      % elseif (abs(shape)<=0.5) changed by ms to deal with shapes in (.4,.5)
      pcov = repmat(nan,3,3);
    elseif (abs(shape)<=0.4) && nargout>1
      % calculate covariance matrix of estimates with
      % linear interpolation between rows
      if nargout>1
      i1 = sum((w(:,1)<=shape));
      i2 = 10-sum((w(:,1)>=shape));
      w_s = w(i1,:)+(shape-w(i1,1))*(w(i2,:)-w(i1,:))/(w(i2,1)-w(i1,1));
      W=[w_s(7),w_s(6),w_s(4);w_s(6),w_s(5),w_s(3);w_s(4),w_s(3),w_s(2)];
      pCov1=[1,scale,scale;scale,scale^2,scale^2;scale,scale^2,scale^2];
      pcov=W.*pCov1/n;
      end
    end

