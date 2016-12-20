function [phat] = fitgenpar(data,varargin) 
%FITGENPAR Parameter estimates for Generalized Pareto data
%
% CALL:  phat = fitgenpar(data,options)
%
%     phat = Struct with estimated parameters 
%     data = one-dimensional data set
%  options = struct with fieldnames
%     .method   : a string, describing the method of estimation
%                'ls'  = Least squares on log scale (robust if n large)
%                'mom' = Moment method
%                'ml'  = Maximum Likelihood method
%                'mps' = Maximum Product of Spacings method. (default)
%                'pkd' = Pickands' estimator 
%                'pwm' = Probability Weighted Moments 
%     .fixpar   : vector giving the fixed parameters. (Must be empty or 
%                 have the same length as the number of parameters. 
%                 Non-fixed parameters must then be given as NaN's)
%                 (default [nan nan 0])
%     .ksign    : Restriction on region for shape, k (LS method only): 
%                 1  k>0                  
%                 0  k can be any value (default)
%                -1  k<0
%     .plotflag : 1, plot the empiricial distribution
%                   function and the estimated cdf 
%                 0, do not plot
%     .alpha    : Confidence coefficent             (default 0.05)
%     .optimset : optimset structure defining performance of the
%                 optimization routine (see optimset for details)
%                  
% FITGENPAR estimates the shape and scale parameters in the Generalized
% Pareto Distribution (GPD). The location parameter is assumed to be zero.
% The default method is MPS since it works for all values of the shape 
% parameter and have the same asymptotic properties as the ML method 
% (when it is valid). PKD and LS also works for any value of the shape.
%  The ML is only valid when shape<=1, the PWM when shape>-0.5, 
% the MOM when shape>-0.25. The variances of the ML estimates are usually 
% smaller than those of the other estimators. However, for small sample 
% sizes it is recommended to use the PWM, MOM or MPS if they are valid.
%
% Example:  
%   R = rndgenpar(0.3,1,0,40,1);
%   
%   phat = fitgenpar(R,'method','ls');
%   figure(1); plotfitsumry(phat);
%   phat2 = fitgenpar(R,'method','ml');
%   figure(2); plotfitsumry(phat2);
%   phat3 = fitgenpar(R,'method','mps');
%   figure(3); plotfitsumry(phat3);
%
%   close all;
% 
% See also plotfitsumry, pdfgenpar, cdfgenpar, invgenpar, rndgenpar, momgenpar

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
%  Borg, S. (1992)
%  XS - a statistical program package in Splus for extreme-value
%  analysis. 
%  Dept. of Mathematical Statistics, Lund University, 1992:E2
%
%  Davidson & Smith (1990)
%  Models for Exceedances over high Threholds.
%  Journal of the Royal Statistical Society B,52, pp. 393-442.
%
%  Grimshaw, S. D. (1993)
%  Computing the Maximum Likelihood Estimates for the Generalized Pareto Distribution.
%  Technometrics, 35, pp. 185-191.
%
% Wong, T.S.T. and Li, W.K. (2006)
% "A note on the estimation of extreme value distributions using maximum
% product of spacings.",
% IMS Lecture Notes Monograph Series 2006, Vol. 52, 272-283

% Tested on; Matlab 5.3
% History: 
% Revised by jr 22.12.1999
% Modified by PJ 08-Mar-2000
%   Changed 'pgpd' to 'gpdcdf' for method 'pkd'.
%   Hjalptext
%   Added 'hold off' in plotting.
%   Added line: 'data = data(:)';'
% revised ms 14.06.2000
% - updated header info
% - changed name to fitgenpar (from gpdfit)
% - added w* to used WAFO-files
% Updated by PJ 22-Jun-2000
%   Added new method 'ml', maximum likelihood estimation of 
%   parameters in GPD.
% Correction by PJ 19-Jul-2000
%   Found error in the covariance matrix for method 'ml'. Now correct!
%   Some other small changes in method 'ml'
% Correction by PJ 30-Nov-2000
%   Updated method 'ml'. New method for finding start values for numerical solving.
%   Updated help text, and some other small changes.
% Correction by PJ 11-Dec-2000
%   Now method 'ml' works with Matlab <5.3. ('fzero' changed format)
% Revised pab Nov2005
% -added method = pkdlog method: Slightly more robust for estimating the upper tail than pkd
% -added method = ls: More robust for estimating the upper tail than pkd
%   and pwm if n is large n>100;
% revised pab aug 2007
%- moved wgpdfun to here as a subfunction
%- made subfunctions for of each of the methods
%- output is now a struct

global WAFO_WSTATS_DEFAULT_PLOTFLAG


error(nargchk(1,inf,nargin))
% Add these options?: 'shape',nan,'scale',nan,'location',0, 
options = struct('method','MPS','fixpar',[nan, nan, 0],'ksign',0,...
  'alpha',0.05,'plotflag', WAFO_WSTATS_DEFAULT_PLOTFLAG,...
  'monitor',false,'optimset',optimset('disp','off','TolX',1e-6,'TolFun',1e-6,'MaxIter',600)); % default options
if (nargin==1 && nargout <= 1 && isequal(data,'defaults'))
  phat = options; 
  return
end
options        = parseoptions(options,varargin{:});
options.method = upper(options.method);
method         = options.method;

data = sort(data(:)).';  % Make sure data is a row vector, PJ 08-Mar-2000
if isempty(options.fixpar)
    options.fixpar = [nan,nan,nan];
end
somefixed = any(isfinite(options.fixpar));
if somefixed && (all(isfinite(options.fixpar) == [0 0 1]))
  m0 = options.fixpar(3);
  if (strcmpi(method,'mps'))
    [phat0,pcov] = mpsfit(data-m0,options);
  elseif strcmpi(method,'ml'),  % Maximum Likelihood
    [phat0,pcov] = mlfit(data-m0,options);
  elseif (strncmpi(method,'pkd',3))
    [phat0,pcov] = pkdfit(data-m0,options);
  elseif strcmpi(method,'mom'),
    [phat0,pcov] = momfit(data-m0);
  elseif (strcmpi(method,'pwm') || strcmpi(method,'ls')),
    [phat0,pcov] = pwmfit(data-m0,options);
  else
    error(['Unknown method ' method '.']);
  end
  %if m0~=0
    phat0 = [phat0, m0];
    pcov = [pcov,zeros(2,1);zeros(1,3)];
  %end
else
  if somefixed
    isnotfixed = (~isfinite(options.fixpar));
    if isnotfixed(3)
      c0 = options.fixpar(3);
    else
      c0 =  min(data)-0.01*std(data);
    end
  else
    c0 = min(data)-0.01*std(data);
  end
  [phat1] = pkdfit(data-c0,options);
  phat0 =  [phat1,c0];
  if somefixed
    phat0 = phat0(isnotfixed);
  end
  phat = mlest(@pdfgenpar,phat0,data,options);
  phat.dataname= inputname(1);
  return
end

n = numel(data);
tcrit    = -invt(options.alpha/2 ,n-1);
% if 1,
pvar  = diag(pcov).';
ciL = phat0-tcrit.*sqrt(pvar);
ciU = phat0+tcrit.*sqrt(pvar);

%shape = phat0(1);
%scale = phat0(2);
%ciL = [shape - tcrit*sqrt(pcov(1)), scale - tcrit*sqrt(pcov(2,2))];
%ciU = [shape + tcrit*sqrt(pcov(1)), scale + tcrit*sqrt(pcov(2,2))];
%   elseif 1
%     chi2crit = invchi2([alpha2 1-alpha2],n-1);
%     ciL = [shape - tcrit*sqrt(pcov(1)), (scale*(n-1)./chi2crit(2))];
%     ciU = [shape + tcrit*sqrt(pcov(1)), (scale*(n-1)./chi2crit(1))];
%   else
%     halfwidth = tcrit*sqrt((1./scale).^2*pcov(2,2));
%   
%     ciL = [shape - tcrit*sqrt(pcov(1)), scale*exp(-halfwidth)];
%     ciU = [shape + tcrit*sqrt(pcov(1)), scale*exp(halfwidth)];
%   end  

[LPS,pvalue] = logps(phat0,data,'cdfgenpar');
  phat = createfdata(options,'dist','pdfgenpar',...
    'params',phat0,'lower',ciL,'upper',ciU,...
    'covar',pcov,'var',pvar,...
    'dataname',inputname(1),'data',data,...
    'loglikemax', -loglike(phat0,data,'pdfgenpar'),'logpsmax',-LPS,...
    'pvalue',pvalue,'note',sprintf('Moran''s statistic on fit: pvalue = %g',pvalue));
  
if options.plotflag 
  plotfitsumry(phat,options.plotflag)
end

phat = fdata(phat);

%   function H = myhessian(phat1)
%     [L,pcov,H] = loglike(phat1,phat.data,phat.)
%   end
% 
% end

function [phat,pcov] = mpsfit(data,options)
%MPSFIT Maximum Product of spacings estimator
 x = sort(data);
 
% %x = unique(data);
 
 phat0 = [0,mean(data)];
% [phat0] = pwmfit(x,options);
% if phat0(1)<0 || any(isnan(phat0))
%   %[phat0] = mlfit(data,options);
%  % if any(isnan(phat0))
%     [phat0] = pkdfit(x,options);
%   %end
% end
if ~isoctave
  msgId = 'WAFO:PARSEOPTIONS';
  warnstate = warning('query',msgId);
  warning('off',msgId);
 end
phat = fminsearch(@(p)logps(p, x,'lowertail',false,@cdfgenpar),phat0,options.optimset); 
%,x,'lowertail',false,@cdfgenpar);
%phat = fminsearch(@logps,phat0,options.optimset,x,@cdfgenpar);
if ~isoctave
  warning(warnstate);
 end
%fun1 = @(phat1) logps(phat1,x,'lowertail',false,@cdfgenpar);
%Hfun = hessian(fun1,phat0 );
%Gfun = @(phat1) gradest(fun1,phat1 );
%phat = mmfminunc(fun1,phat0,'Hessian',Hfun);
%phat = phat(:).';
 %[tmp, pcov2] = loglike(phat,data,@pdfgenpar);
 [LL,pcov] = likgenpar(phat,data);
if  any(det(pcov)<0) || any(~isfinite(pcov(:)))
 
  %
  %[LL,pcov] = likgenpar(phat,data);
  if phat(1)>-0.5
    n = numel(data);
    shape = phat(1);
    scale = phat(2);
    f = 1/(1+2*shape)/(3+2*shape);
    Vscale = scale^2*(7+18*shape+11*shape^2+2*shape^3);
    Covari = scale*(2+shape)*(2+6*shape+7*shape^2+2*shape^3);
    Vshape = (1+shape)*(2+shape)^2*(1+shape+2*shape^2);
    pcov = f/n*[Vshape Covari; Covari, Vscale];
  else
     n = numel(data);
     shape = phat(1);
     scale = phat(2);
     Vshape = 1-shape;
     Vscale = 2*scale^2;
     Covari = scale;
     pcov = (1-shape)/n*[Vshape Covari; Covari Vscale];
  end
end


function [phat,pcov] = mlfit(data,options)
%MLFIT ML estimator
 % See Davidson & Smith (1990) and Grimshaw (1993)
 
  % Calculate start values
  % The start value can not be less than  1/max_data ,
  %   since if it happens then the upper limit of the GPD is 
  %   lower than the highest data value
  
  % Change variables, in order to avoid boundary problems in numerical solution 
  %   Transformation: x = log(1/max_data - t),   -Inf < t < 1/max_data
  %   Inverse Trans.: t = 1/max(data) - exp(x),  -Inf < x < Inf
  
  
  % New version
  % Find an interval where the function changes sign.
  % Gives interval for start value for the zero-search.
  % Upper and lower limits are given in Grimshaw (1993)
  
  %n = numel(data);
  max_data = max(data);
  X1 = min(data); 
  Xn = max_data; 
  Xmean = mean(data);
  Eps = 1e-6/Xmean;
  t_L = 2*(X1-Xmean)/X1^2;     % Lower limit
  t_U = 1/Xn-Eps;              % Upper limit
  x_L = log(1/max_data - t_L); % Lower limit
  x_U = log(1/max_data - t_U); % Upper limit
  Nx = 10;
  x =linspace(x_U,x_L,Nx);
  [f,k] = fitgenparml(x,data);
  i0 = find(k<=1,1,'first');
  if i0>1
    x_U = x(i0-1);
    x(1:i0-1) = [];
    f(1:i0-1) = [];
  end
  I = findcross(f,0);
  if isempty(I) || f(1)<0
    [x,f]=fplot(@(x)fitgenparml(x,data),[x_U,x_L]);
    i0 = find(f>0,1,'first');
    if i0>1
      x(1:i0-1) = [];
      f(1:i0-1) = [];
    end
    I = findcross(f,0);
  end
  if isempty(I)
      x_start = [];
  else
    i = I(1);
    x_start = [x(i) x(i+1)];
  end
  
  
  if isempty(x_start)
    shape=NaN; scale=NaN; 
    pcov=ones(2)*NaN;
    warning('WAFO:FITGENPAR','Can not find an estimate. Probably the ML-estimator does not exist for this data. Try other methods.')
  else
    % Solve the ML-equation
    
    xx = fzero(@(x)fitgenparml(x,data),x_start,options.optimset);
    %xx = fzero('fitgenparml',x_start,options.optimset,data);
    %xx = fzero('fitgenparml',x_start,optimset('disp','iter'),data);
    
    
    % Extract estimates
    [f,shape,scale] = fitgenparml(xx,data);
    if nargout>1
    % Calculate the covariance matrix by inverting the observed information. 
    if (shape >= 1.0),
      pcov=ones(2)*NaN;
      warning('WAFO:FITGENPAR',[' The estimate of shape (' num2str(shape) ') is not within the range where the ML estimator is valid (shape<1).'])
    else
      if (shape >= 0.5),
        warning('WAFO:FITGENPAR',[' The estimate of shape (' num2str(shape) ') is not within the range where the ML estimator is asymptotically normal (shape<1/2).'])
      end
      %[tmp, pcov2] = loglike([shape,scale],data,@pdfgenpar);
      [LL,pcov] = likgenpar([shape,scale],data);
      % This is wrong:
      %      Vshape = 1-shape;
      %      Vscale = 2*scale^2;
      %      Covari = scale;
      %      pcov = (1-shape)/n*[Vshape Covari; Covari Vscale];
    end
    end
  end
  phat = [shape,scale];
% End MLFIT

function [phat,pcov]=pwmfit(data,options)
  n = numel(data);
  a0 = mean(data);
  xio = sort(data);
  %R1 = n+0.35-(1:n);
  %a1 = sum(p1.*xio)/n/n;
  R1 = (n+0.35-(1:n))/n;
  a1 = mean(R1.*xio);
  %a2 = mean(p1.^2.*xio);
  
  shape = a0/(a0-2*a1)-2;
  scale = 2*a0*a1/(a0-2*a1);
  
  phat = [shape,scale];
  LS = strcmpi(options.method,'ls');
  if LS  
    %R1 = 1-(0.5:(n - 0.5)).'./n;
    if n>10000
      Ri = linspace(R1(1),R1(end-5),10000).';
      xio =interp1(log(R1),xio,log(Ri),'linear');
      %xio =interp1(R1,xio,Ri,'linear');
      R1  = Ri;
    end
   
    %ksign = 0; % k can be any value
    %ksign = 1; % k>0
    %ksign = -1; % k<0
    shape = shape*(options.ksign==0 || sign(options.ksign)==sign(shape));

    phat = fminsearch(@(p)gpdfun(p,xio(:),R1(:),options),[shape,scale],options.optimset);

    shape = phat(1);
    scale = phat(2);
  end
  
  if (shape <= -0.5) 
    %pcov=ones(2)*NaN;
    if ~LS
      warning('WAFO:FITGENPAR',[' The estimate of shape (' num2str(shape) ') is not within the range where the PWM estimator is valid (shape>-0.5).'])
    end
   % [tmp, pcov] = loglike([shape,scale],data,'pdfgenpar');
    [tmp, pcov] = likgenpar([shape,scale],data);
  else 
    f = 1/(1+2*shape)/(3+2*shape);
    Vscale = scale^2*(7+18*shape+11*shape^2+2*shape^3);
    Covari = scale*(2+shape)*(2+6*shape+7*shape^2+2*shape^3);
    Vshape = (1+shape)*(2+shape)^2*(1+shape+2*shape^2);
    pcov = f/n*[Vshape Covari; Covari, Vscale];
  end
 
function [phat,pcov]=momfit(data,options)
%MOMFIT moments estimator
n = numel(data);
m = mean(data); 
s = std(data);
%sk = skew(data);

shape = ((m/s)^2 - 1)/2;
scale = m*((m/s)^2+1)/2;
  
if (shape<=-0.25),
  pcov=ones(2)*NaN;
  warning('WAFO:FITGENPAR',[' The estimate of shape (' num2str(shape) ') is not within the range where the Moment estimator is valid (shape>-0.25).'])
else
  Vshape = (1+2*shape)^2*(1+shape+6*shape^2);
  Covari = scale*(1+2*shape)*(1+4*shape+12*shape^2);
  Vscale = 2*scale^2*(1+6*shape+12*shape^2);
  cov_f = (1+shape)^2/(1+2*shape)/(1+3*shape)/(1+4*shape)/n;
  pcov=cov_f*[Vshape Covari; Covari, Vscale];
end
phat = [shape,scale];


function [phat,pcov] = pkdfit(data,options)
% PKDFIT Pickands' estimator
  epsilon = 1.e-4;
  n  = numel(data);
  xi = -sort(-data);  %  data in descending order
  x  = flipud(xi(:));
  R  = ((n:-1:1).'-0.5)/n;
  
  
  
  % Find the M that minimizes diff. between EDF and G, the estimated cdf.
  n4 = floor(n/4);
  dmax = realmax;
  d  = dmax(ones(1,n4));
  p = inf;
  plog = 1;
  useLogMetric = strncmpi(options.method,'pkdlog',6);
  
  for M=1:n4,
    shape = - log((xi(M)-xi(2*M))/(xi(2*M)-xi(4*M)))/log(2);
    scale = (xi(2*M)-xi(4*M));
    
    if (abs(shape) < epsilon),
      scale = scale/log(2);
    else
      scale = scale*shape/(1-0.5^shape);
    end
    logR1 = max(cdfgenpar(x,shape,scale,0,'lowertail',false,'logp',true),-realmax);
    if useLogMetric
      d(M) = norm((logR1-log(R))/n,plog);
      %d(M) = norm(log((exp(logR1)-R))/n,plog);
    else
      d(M) = norm((exp(logR1)-R)/n,p);
    end
    
  end    % end M-loop
  
  [tmp,M] = min(d);
  %least = find(d(:)==min(d));
  %M = least(1); % possibly more than one M; in that case use the smallest.
  
  shape = - log((xi(M)-xi(2*M))/(xi(2*M)-xi(4*M)))/log(2);
  scale = (xi(2*M)-xi(4*M));
  
  if (abs(shape) < epsilon),
    scale = 1000*scale/log(2);
  else
    scale = scale*shape/(1-1/(2^shape));
  end
phat = [shape,scale];
pcov=ones(2,2)*nan;
%[tmp, pcov] = loglike([shape,scale],data,'pdfgenpar');


function y=gpdfun(phat,x,R1,options)
% GPDFUN Is an internal routine for fitgenpar
%

s =1; m = 0;

k = phat(1);
Np = length(phat);
if Np>1, s = phat(2);end
if Np>2, m = phat(3); end

mlogR1 = -log(R1);
mlogR = min(-cdfgenpar(x,k,s,m,'lowertail',false,'logp',true),realmax);

p = 1;
y = norm(mlogR1-mlogR,p)/numel(x);
if (options.ksign~=0)
  %y = y + 1000*abs(k)*(((k<0)&(options.ksign==1)) | ((k>0)&(options.ksign==-1)));
  y = y + 1000*abs(k)*(sign(options.ksign)~=sign(k));
end

if options.monitor
  semilogx(x,[ mlogR1 mlogR]); drawnow,shg
  disp(['err = ' num2str(y,10)   ' k s m = ' num2str([k,s,m],4) ])
end




