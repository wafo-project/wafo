function [lcEst,Est,R,MSE] = cmat2extralc(param,F,u,method,plotflag)
%CMAT2EXTRALC  Extrapolate level crossing spectrum
%
% CALL: lcEst = cmat2extralc(param,F,u)
%       [lcEst,Est,R,MSE] = cmat2extralc(param,F,u,method,plotflag)
%
%   param  = Parameter vector, [a b n], defines the grid.
%   F      = Cycle matrix. (e.g. Rainflow matrix)          [nxn]
%   u      = Levels [u_min u_max], extrapolate below u_min and above u_max. 
%   method = A string, describing the method of estimation.
%            Generalized Pareto distribution (GPD):
%             'gpd,ml'   = Maximum likelihood estimator. (default)
%             'gpd,mom'  = Moment method.
%             'gpd,pwm'  = Probability Weighted Moments.
%            Exponential distribution (GPD with k=0)
%             'exp,ml'   = Maximum likelihood estimator.
%             'exp,mld'  = Maximum likelihood estimator, discrete.
%             'exp,ls1'  = Least Squares estimator, linear, 1 parameter.
%             'exp,ls4'  = Least Squares estimator, quadratic, 2 parameters. 
%             'exp,wls1' = Weighted Least Squares estimator.
%             'exp,wls4' = Weighted Least Squares estimator.
%            Rayleigh distribution 
%             'ray,ml'   = Maximum likelihood estimator.
%   plotflag = 1: Diagnostic plots. (default)
%              0: Don't plot diagnostic plots.
%
%   lcEst    = The estimated LC.     [struct array]
%   Est      = Estimated parameters. [struct array]
%   R        = Results from internal calculations. [struct array]
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
% The methods 'exp,...' uses the Exp. We recommend the use of 'exp,mld'.
% Also the methods 'exp,ls1', 'exp,wls1', 'exp,ls4', 'exp,wls4' work good.
% The method 'ray,ml' uses Ray, and should be used if the load is a Gaussian process.
%
% For further help use 'type cmat2extralc.txt'.
%
% Example: 
%   [G, Gh] = mktestmat([-1 1 64],[-0.2 0.2], 0.15,1);
%   xD = mctpsim({G Gh}, 2000);
%   Frfc = dtp2rfm(xD, 64, 'CS');
%   [lcEst, Est] = cmat2extralc([-1 1 64], Frfc, [-0.4 0.4]);
%   lcG = cmat2lc([-1 1 64], G/sum(sum(G)));
%   lcF = cmat2lc([-1 1 64], Frfc);
%   clf;
%   semilogx(1000*lcG(:,2), lcG(:,1),...
%            lcF(:,2), lcF(:,1),...
%            lcEst.lc(:,2), lcEst.lc(:,1));
%
%   close all;
%
% See also  rfmextrapolate, lc2rfmextreme, extralc, fitgenpar

% References:
%
%   Johannesson, P., and Thomas, J-.J. (2000): 
%   Extrapolation of Rainflow Matrices. 
%   Preprint 2000:82, Mathematical statistics, Chalmers, pp. 18. 

% Tested  on Matlab  5.3
%
% History:
% Created by PJ (Par Johannesson) 10-Mar-2000
% Changed by PJ 14-Mar-2000
%   Added sub-function 'make_increasing'
% Changed by PJ 30-Mar-2000
%   Added method 'exp' (GPD with k=0)
% Changed by PJ 05-Apr-2000
%   Modification of 'extrapolation'
% Changed by PJ 11-Apr-2000
%   Added method 'exp,ls'
% Changed by PJ 17-Apr-2000
%   Added methods 'exp,ls1','exp,ls2','exp,ls3','exp,wls1'
% Changed by PJ Apr-2000
%   Added methods 'exp,ls4','exp,wls2','exp,wls3','exp,wls4'
% Changed by PJ 19-Jul-2000
%   Diagnostic plots.
%   Tidy the code.
%   Updated helptext. 
% Updated by PJ 10-Oct-2000
%   More efficient implemantation of 'gpd,ml'.
%   Old version is now called 'gpd,ml0'
%   Some additional changes also made.
% Updated by PJ 03-Nov-2000
%   New method 'ray,ml'
% Updated by PJ 15-Nov-2000
%   Updated method 'gpd,ml'.
%   Now more robust. Does not get into infinite loops.
%   If ML-estimator does not exist, then it uses MOM instead.
% Updated by PJ 6-Dec-2000
%   Updated help text.
% Correction by PJ 11-Dec-2000
%   Now method 'ml' works with Matlab <5.3. ('fzero' changed format)

% Check input
ni = nargin;
no = nargout;
error(nargchk(3,5,ni));

if ni < 4, method = []; end
if ni < 5, plotflag = []; end

% Default values
if isempty(method), method = 'gpd,ml'; end
if isempty(plotflag), plotflag = 1; end

% Extract level crossings
n = param(3);
paramD = [1 n n];
lc0 = cmat2lc(paramD,F);
uu = levels(param)';
uuD = levels(paramD)';

% Find the index for the high level
u_high = u(2);
i_high = min(uuD(uu>u_high));

% Find the index for the low level
u_low = u(1);
i_low = max(uuD(uu<u_low));

% Remove cycles with  min>i_high  and  max<i_low
F1 = F;
for i = i_high:n
  F1(i,:) = 0;
end
for j = 1:i_low
  F1(:,j) = 0;
end

% Get level crossings
lc1 = cmat2lc(paramD,F1);

lc=lc1;
%lc=lc0;
[dummy,I]=max(lc(:,2));
i_LCmax = lc(I(1),1);

% Extrapolate LC for high levels
if plotflag, subplot(2,1,1), end

lcHigh = lc(i_high:end,:);
[lcEst.D.High,Est.D.High,MSE.High] = extrapolate(lcHigh,method,plotflag,i_high-i_LCmax);

if plotflag, title([method ': Extrapolation for high levels']),  end

% Extrapolate LC for low levels
if plotflag, subplot(2,1,2), end

lcLow = lc(1:i_low,:);
lcLow1 = [n+1-flipud(lcLow(:,1)) flipud(lcLow(:,2))];

[lcEst1,Est1,MSE.Low] = extrapolate(lcLow1,method,plotflag,i_LCmax-i_low);
lcEst.D.Low = [n+1-flipud(lcEst1(:,1)) flipud(lcEst1(:,2))];
Est.D.Low = Est1;

if isfield(Est.D.Low,'UpperLimit');
  Est.D.Low.LowerLimit = n+1-Est.D.Low.UpperLimit;
  rmfield(Est.D.Low,'UpperLimit');
end

if plotflag, title([method ': Extrapolation for low levels']), end

% Total Level crossing spectrum
lcEst.D.lc = lc;
lcEst.D.lc(1:i_low,:) = lcEst.D.Low;
lcEst.D.lc(i_high:end,:) = lcEst.D.High;

% Convert to correct scales
lcEst.High = [uu(lcEst.D.High(:,1)) lcEst.D.High(:,2)];
lcEst.Low = [uu(lcEst.D.Low(:,1)) lcEst.D.Low(:,2)];
lcEst.lc = [uu(lcEst.D.lc(:,1)) lcEst.D.lc(:,2)];

if no > 1
  Est.method = method;
end

% Supplementary results from internal calculations
if no > 2
  R.lc0 = cmat2lc(param,F);
  R.lc1 = cmat2lc(param,F1);
  R.F1 = F1;
end

% END cmat2extralc
%%%%%

function [lcEst,Est,MSE] = extrapolate(lc,method,plotflag,offset)
%
% Extrapolate the level crossing specrta for high levels
%
  
  iu = lc(1,1);
  
  % Excedences over level u
  lc1 = [lc(:,1)-iu lc(:,2)];
  lc2 = flipud(lc1);
  lc3 = lc2;
  
  method = lower(method);
  methodShape = method(1:3);
  methodEst = method(5:end);
  
  if strcmp(methodShape,'gpd') && ~strcmp(method,'gpd,ml')
  %if strcmp(methodShape,'gpd') | strcmp(method,'exp,ml')| strcmp(method,'exp,mld')
    x = [];
    for k = 2:length(lc3)
      nn = lc3(k,2)-lc3(k-1,2);
      if nn > 0
        x = [x; lc3(k,1)*ones(nn,1)];
      end
    end
    x = x+.5;   
  end
  
  if strcmp(methodShape,'exp') || strcmp(method,'gpd,ml') || strcmp(method,'ray,ml')
    dN = [];
    for k = 2:length(lc3)
      nn = lc3(k,2)-lc3(k-1,2);
      if nn ~= 0
        dN = [dN; lc3(k,1) nn]; 
      end
    end
    NN = [dN(:,1) cumsum(dN(:,2))];
    NN = flipud(NN);
    dN = flipud(dN);
  end
  
  % Estimate tail
  if strcmp(methodShape,'gpd') % GPD
    
    if strcmp(methodEst,'ml')
      
      dNN = [dN(:,1)+0.5 dN(:,2)];
      n = NN(1,2);
      
      % Calculate start values
      X1=dNN(1,1); Xn=dNN(end,1); Xmean = sum(dNN(:,1).*dNN(:,2))/sum(dNN(:,2));
      Eps = 1e-6/Xmean;
      t_L = 2*(X1-Xmean)/X1^2;     % Lower limit
      t_U = 1/Xn-Eps;              % Upper limit
      x_L = log(1/Xn - t_L); % Lower limit
      x_U = log(1/Xn - t_U); % Upper limit
      
      x=linspace(x_L,x_U,10);
      for i=1:length(x)
        f(i)=fitgenpar_mld(x(i),dNN);
      end
      
      I = find(f(1:end-1).*f(2:end) < 0); % Find sign-shift
      
      % Check if any sign-shift was found
      if ~isempty(I) % If ML-estimator exists, then use it! 
        i = I(1);
        x_start = [x(i) x(i+1)];
        if exist('optimset') >= 2 % Function 'optimset' exists ?
          % For Matlab 5.3 and higher ???
          x_ML = fzero(@(x)fitgenpar_mld(x,dNN),x_start,optimset('disp','off'));
        else 
          % For Matlab 5.2 and lower ???
          x_ML = fzero(@(x)fitgenpar_mld(x,dNN),x_start);
        end
        [f,shape,scale] = fitgenpar_mld(x_ML,dNN); % Estimates k_ML and s_ML
        
        % Calculate the (asymptotic) covariance matrix (obtained by inverting the observed information). 
        if (shape >= 1.0),
          cov=[];
          warning([' The estimate of shape (' num2str(shape) ...
              ') is not within the range where the ML estimator is valid (shape<1).'])
        elseif (shape >= 0.5),
          cov=[];
          warning([' The estimate of shape (' num2str(shape) ...
              ') is not within the range where the ML estimator is asymptotically normal (shape<1/2).'])
        else % Calculate the covariance matrix
          Vshape = 1-shape;
          Vscale = 2*scale^2;
          Covari = scale;
          cov = (1-shape)/n*[Vshape Covari; Covari Vscale];
        end
        
      else % If ML-estimator does not exists, then use MOM!
        warning([' ML-estimator does not exist for this data. Using MOM instead.'])

        Xvar = sum((dNN(:,1)-Xmean).^2.*dNN(:,2))/(sum(dNN(:,2))-1);
        shape = (Xmean^2/Xvar-1)/2;
        scale = Xmean*(Xmean^2/Xvar+1)/2;
        cov = [];
      end
      parms = [shape scale];
      
    else
      if strcmp(methodEst,'ml0')
        [parms,cov] = fitgenpar(x,'ml',0);
      else
        [parms,cov] = fitgenpar(x,methodEst,0);
      end
    end
    
    Est.k = parms(1);
    Est.s = parms(2);
    if Est.k>0,
      Est.UpperLimit = iu+Est.s/Est.k;
    end
    Est.cov = cov;
    xF = (0:length(lc1)-1)';
    F = cdfgenpar(xF,Est.k,Est.s);
    
  elseif strcmp(methodShape,'exp') % GPD with k=0
    
    switch methodEst
      
    case 'ml' % Maximum Likelihood
      
      % Old version
      %n = length(x);
      %parms(1) = 0;        % k
      %parms(2) = sum(x)/n; % s
      
      % New version
      n = NN(1,2);
      Sx = sum(dN(:,2).*(dN(:,1)+0.5));
      parms(1) = 0;        % k
      parms(2) = Sx/n;     % s
      
      cov = zeros(2);          % Covariance matrix for estimates.
      cov(2,2) = 1/n*parms(2); % Variance of estimate  s
      
    case 'mld' % Maximum Likelihood, discrete obs.
      
      n = NN(1,2);
      Sx = sum(dN(:,2).*(dN(:,1)+0.0));
      s1 = 1/log(1+n/Sx);
      
      parms(1) = 0;    % k
      parms(2) = s1;   % s
    
      cov = zeros(2);          % Covariance matrix for estimates.
      cov(2,2) = 1/n*parms(2); % Variance of estimate  s
      
      
    case {'ls1','wls1'} % One parameter: Least Squares / Weighted LS
      
      X = NN(:,1);
      Y = log(NN(:,2))-log(NN(1,2));
      n = length(Y);
      
      if strcmp(methodEst,'ls1') % Least Squares
        W = eye(n);
      else % Weighted Least Squares
        W = diag(dN(:,2));
      end
      
      % Compute Estimate
      th = (X'*W*X)\(X'*W*Y);

      parms(1) = 0;          % k
      parms(2) = -1/th(1);   % s
      cov = [];
            
    case {'ls2','wls2'} % Two parameters: Least Squares / Weighted LS
          
      X = [ones(size(NN(:,1))) NN(:,1)];
      Y = log(NN(:,2));
      n = length(Y);
      
      if strcmp(methodEst,'ls2') % Least Squares
        W = eye(n);
      else % Weighted Least Squares
        W = diag(dN(:,2));
      end
      
      % Compute Estimate
      th = (X'*W*X)\(X'*W*Y);

      parms(1) = 0;          % k
      parms(2) = -1/th(2);   % s
      parms(3) = exp(th(1)); % mu0
      cov = [];
                  
    case {'ls3','wls3'} % Three parameters: Least Squares / Weighted LS
          
      X = [ones(size(NN(:,1))) NN(:,1) NN(:,1).^2];
      Y = log(NN(:,2));
      n = length(Y);
      
      if strcmp(methodEst,'ls3') % Least Squares
        W = eye(n);
      else % Weighted Least Squares
        W = diag(dN(:,2));
      end
      
      % Compute Estimate
      th = (X'*W*X)\(X'*W*Y);
      
      parms(1) = 0;          % k
      parms(2) = -1/th(2);   % s
      parms(3) = exp(th(1)); % mu0
      parms(4) = -1/th(3);   % s2
      Est.s2 = parms(4);
      cov = [];
                  
    case {'ls4','wls4'} % Three parameters: Least Squares / Weighted LS
          
      X = [NN(:,1) NN(:,1).^2];
      Y = log(NN(:,2))-log(NN(1,2));
      n = length(Y);
      
      if strcmp(methodEst,'ls4') % Least Squares
        W = eye(n);
      else % Weighted Least Squares
        W = diag(dN(:,2));
      end
      
      % Compute Estimate
      th = (X'*W*X)\(X'*W*Y);
      
      parms(1) = 0;          % k
      parms(2) = -1/th(1);   % s
      parms(4) = -1/th(2);   % s2
      Est.s2 = parms(4);
      cov = [];
            
    end % switch methodEst
    
    Est.k = parms(1);
    Est.s = parms(2);
    Est.cov = cov;
    xF = (0:length(lc1)-1)';
    switch methodEst
    case {'ml', 'mld', 'ls1', 'wls1', 'ls2', 'wls2'}
      F = 1 - exp(-xF/Est.s);
    case {'ls3', 'wls3', 'ls4', 'wls4'}
      F = 1 - exp(-xF/Est.s-xF.^2/Est.s2);
    end
    
  elseif strcmp(methodShape,'ray') % Rayleigh
    
    n = NN(1,2);
    Sx = sum(dN(:,2).*((dN(:,1)+offset+0.5).^2-(offset)^2));
    a=sqrt(Sx/n);  % Shape parameter
    
    Est.s = a;
    
    xF = (0:length(lc1)-1)';
    F = 1 - exp(-((xF+offset).^2-offset^2)/a^2);
    
  else % Unknown method
    
    error(['Unknown method: ' method]);
    
  end
  
  % Estimate of mu0
  switch methodEst
  case {'ls2','wls2','ls3','wls3'}
    Est.lcu = parms(3);
  otherwise
    Est.lcu = lc(1,2);
  end
  
  lcEst = [xF+iu Est.lcu*(1-F)];
  lcEst(end,2) = 0;  % No crossings of the highest level
  
  % Compute MSE (Mean Square Error)
  if strcmp(method,'gpd,ml') || strcmp(method,'exp,ml') || strcmp(method,'exp,mld') || strcmp(method,'ray,mld')
    MSE = 0;
    % Not impemented
  elseif strcmp(methodShape,'gpd')
    n=length(x);
    q1 = ((1:n)-1/2)/n;
    xq1 = invgenpar(q1,Est.k,Est.s)';
    MSE = sum((x-xq1).^2)/n;
  elseif strcmp(method,'ray,ml')
    MSE = 0;
    % Not impemented
  else % (Weighted) Least Squares estimate
    Lam = sqrt(W);
    MSE = (Lam*X*th - Lam*Y)'*(Lam*X*th - Lam*Y)/sum(diag(W));
  end
  
  % Plot Diagnostics 
  if plotflag
    if strcmp(method,'gpd,ml') || strcmp(method,'exp,ml') || strcmp(method,'exp,mld')
      Csum = cumsum(dN(:,2));
      q1 = ([1; Csum(1:end-1)+1; Csum(1:end)]-1/2)/n;
      xq1 = invgenpar(q1,Est.k,Est.s)';
      x = [dN(:,1); dN(:,1)]+0.5;
      plotqq(x,xq1); hold on
      plot([0 max(x)],[0 max(x)]), hold off
      xlabel('data'), ylabel('Quantile')
    elseif strcmp(methodShape,'gpd')
      plotqq(x,xq1); hold on
      plot([0 max(x)],[0 max(x)]), hold off
      xlabel('data'), ylabel('Quantile')
    elseif strcmp(method,'ray,ml')
      Csum = cumsum(dN(:,2));
      q1 = ([1; Csum(1:end-1)+1; Csum(1:end)]-1/2)/n;
      xq1 = sqrt(offset^2-Est.s^2*log(1-q1))-offset;
      x = [dN(:,1); dN(:,1)]+0.5;
      plotqq(x,xq1); hold on
      plot([0 max(x)],[0 max(x)]), hold off
      xlabel('data'), ylabel('Quantile')
    else
      plot(X(:,1),Y,'+'), hold on
      plot(X(:,1),X*th,'r'), hold off
      xlabel('data'), ylabel('log(lc)')
    end
  end
    
% END extrapolate
%%%%%
