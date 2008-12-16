function [lcEst,Est] = extralc(lc,u,method,plotflag)
%EXTRALC  Extrapolate level crossing spectrum
%
% CALL: [lcEst,Est] = extralc(lc,u,method,plotflag)
%
%   lc     = Level crossing spectrum.                       [n,2]
%   u      = Levels [u_min u_max], extrapolate below u_min and above u_max. 
%   method = A string, describing the method of estimation.
%            Generalized Pareto distribution (GPD):
%             'gpd,ml'   = Maximum likelihood estimator. (default)
%             'gpd,mom'  = Moment method.
%             'gpd,pwm'  = Probability Weighted Moments.
%             'gpd,pkd'  = Pickands' estimator.
%            Exponential distribution (GPD with k=0)
%             'exp,ml'   = Maximum likelihood estimator.
%            Rayleigh distribution 
%             'ray,ml'   = Maximum likelihood estimator.
%   plotflag = 1: Diagnostic plots. (default)
%              0: Don't plot diagnostic plots.
%
%   lcEst    = The estimated LC.     [struct array]
%   Est      = Estimated parameters. [struct array]
%
% Extrapolates the level crossing spectrum for high and for low levels. 
% The tails of the LC is fited to a survival function of a GPD. 
%   H(x) = (1-k*x/s)^(1/k)               (GPD)
% The use of GPD is motivated by POT methods in extreme value theory. 
% For k=0 the GPD is the exponential distribution
%   H(x) = exp(-x/s),  k=0               (Exp)
% The tails with the survival function of a Rayleigh distribution.
%   H(x) = exp(-((x+x0).^2-x0^2)/s^2)    (Ray)
% where x0 is the value where the LC has its maximum.
% The methods 'gpd,...' uses the GPD. We recommend the use of 'gpd,ml'.
% The method 'exp,ml' uses the Exp. 
% The method 'ray,ml' uses Ray, and should be used if the load is a Gaussian process.
%
% Example: 
%   S = jonswap;
%   x = spec2sdat(S,100000,0.1,[],'random');
%   lc = dat2lc(x); s = std(x(:,2));
%   [lcEst,Est] = extralc(lc,s*[-2 2]);
%   [lcEst,Est] = extralc(lc,s*[-2 2],'exp,ml');
%   [lcEst,Est] = extralc(lc,s*[-2 2],'ray,ml');
%
% See also  cmat2extralc, rfmextrapolate, lc2rfmextreme, extralc, fitgenpar

% References:
%
%   Johannesson, P., and Thomas, J-.J. (2000): 
%   Extrapolation of Rainflow Matrices. 
%   Preprint 2000:82, Mathematical statistics, Chalmers, pp. 18. 

% Tested  on Matlab  5.3
%
% History:
% Created by PJ (Pär Johannesson) 10-Mar-2000
% Changed by PJ 14-Mar-2000
%   Added sub-function 'make_increasing'
% Changed by PJ 15-Mar-2000
%   Added sub-function 'make_increasing'
% Updated by PJ 07-Dec-2000
%   Added 'exp,ml' and 'ray,ml'
%   Help text
% Corrected by PJ 17-Feb-2004

% Check input

ni = nargin;
%no = nargout;
error(nargchk(1,4,ni));

if ni < 3, method = []; end
if ni < 4, plotflag = []; end

% Default values
if isempty(method), method = 'gpd,ml'; end
if isempty(plotflag), plotflag = 1; end

% Maximum of lc
[M,I] = max(lc(:,2));  % Corrected by PJ 17-Feb-2004
lc_max = lc(I(1),1);

% Extrapolate LC for high levels
[lcEst.High,Est.High] = extrapolate(lc,u(2),method,u(2)-lc_max);

% Extrapolate LC for low levels

lc1 = [-flipud(lc(:,1)) flipud(lc(:,2))];

[lcEst1,Est1] = extrapolate(lc1,-u(1),method,lc_max-u(1));
lcEst.Low = [-flipud(lcEst1(:,1)) flipud(lcEst1(:,2:end))];
Est.Low = Est1;

if plotflag
  semilogx(lc(:,2),lc(:,1),lcEst.High(:,2),lcEst.High(:,1),lcEst.Low(:,2),lcEst.Low(:,1))
end

%%%
function [lcEst,Est] = extrapolate(lc,u,method,offset)
  % Extrapolate the level crossing specrta for high levels
  
  method = lower(method);
  methodShape = method(1:3);
  methodEst = method(5:end);
  
  % Excedences over level u
  Iu = lc(:,1)>u;
  lc1 = lc(Iu,:);
  lc2 = flipud(lc1);
  lc3 = make_increasing(lc2);
  
  % Corrected by PJ 17-Feb-2004
  lc3=[0 0; lc3]; x=[];
  for k=2:length(lc3(:,1))
      nk = lc3(k,2)-lc3(k-1,2);
      x = [x; ones(nk,1)*lc3(k,1)];
  end
  x = x-u;   
  
  % Estimate tail
  switch methodShape
    
  case 'gpd'
    
    [phat] = fitgenpar(x,'method',methodEst);
    
    Est.k = phat.params(1);
    Est.s = phat.params(2);
    if isnan(phat.fixpar(3))
      Est.cov = phat.covariance;
    else
      Est.cov = phat.covariance(1:2,1:2);
    end
    if Est.k>0,
      Est.UpperLimit = u+Est.s/Est.k;
    end
    
    xF = (0.0:0.01:4)';
    F = cdfgenpar(xF,Est.k,Est.s);
    Est.lcu = interp1(lc(:,1),lc(:,2),u)+1;
    
    % Calculate 90 % confidence region, an ellipse, for (k,s)
    [B,D] =eig(Est.cov);
    b = [Est.k; Est.s];
    
    r = sqrt(-2*log(1-90/100)); % 90 % confidence sphere
    Nc = 16+1;
    ang = linspace(0,2*pi,Nc);
    c0 = [r*sqrt(D(1,1))*sin(ang); r*sqrt(D(2,2))*cos(ang)]; % 90% Circle
%    plot(c0(1,:),c0(2,:))
    
    c1 = B*c0+b*ones(1,length(c0)); % Transform to ellipse for (k,s)
%    plot(c1(1,:),c1(2,:)), hold on
    
    % Calculate conf.int for lcu
    % Assumtion: lcu is Poisson distributed
    % Poissin distr. approximated by normal when calculating conf. int.
    dXX = 1.64*sqrt(Est.lcu); % 90 % quantile for lcu
    
    lcEstCu = zeros(length(xF),Nc); lcEstCl = lcEstCu;
    for i=1:Nc
      k=c1(1,i);
      s=c1(2,i);
      F2 = cdfgenpar(xF,k,s);
      lcEstCu(:,i) = (Est.lcu+dXX)*(1-F2);
      lcEstCl(:,i) = (Est.lcu-dXX)*(1-F2);
    end
    
    lcEstCI = [min(lcEstCl')' max(lcEstCu')'];
    
    lcEst = [xF+u Est.lcu*(1-F) lcEstCI];
    
  case 'exp'
    
    n = length(x);
    s = mean(x);
    cov = s/n;
    Est.s = s;
    Est.cov = cov;
    
    xF = (0.0:0.01:4)';
    F = 1-exp(-xF/s);
    Est.lcu = interp1(lc(:,1),lc(:,2),u)+1;
    
    lcEst = [xF+u Est.lcu*(1-F)];
    
  case 'ray'
    
    n = length(x);
    Sx = sum((x+offset).^2-offset^2);
    s=sqrt(Sx/n);  % Shape parameter
    
    Est.s = s;
    Est.cov = NaN;
    
    xF = (0.0:0.01:4)';
    F = 1 - exp(-((xF+offset).^2-offset^2)/s^2);
    Est.lcu = interp1(lc(:,1),lc(:,2),u)+1;
    
    lcEst = [xF+u Est.lcu*(1-F)];
    
  otherwise
    
    error(['Unknown method: ' method]);
    
  end
  
    
%% End extrapolate
  
function y=make_increasing(x)
  % Makes the signal x strictly increasing. 
  %
  % x = two column matrix with times in first column and values in the second.
  
  n = length(x);
  
  i=1;
  y = x;
  y(1,:) = x(1,:); j=1;
  
  while i<n
    while x(i,2)<=y(j,2) && i<n
      i = i+1;
    end
    if x(i,2)>y(j,2)
      j = j+1;
      y(j,:) = x(i,:);
    end
  end
  
  y = y(1:j,:);
  
%% End make_increasing

