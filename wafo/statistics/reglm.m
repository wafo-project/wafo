%REGLM Fit multiple Linear Regression Model.
%
% CALL model = reglm(y, X,options)
%
% model = fitted LM object with methods
%      .compare() : Compare small LM object versus large one
%      .predict() : Predict from a fitted LM object
%      .summary() : Display summary of fitted LM object.
%
% y      = column vector of observed values
% X      = matrix of regressors, with the first column filled with the constant value 1
% options= struct defining performance of REGLM with fieldnames: 
% .constant: if TRUE include constant term in the model otherwise omit it.
% .offset  : vector or offset values (default 0)
% .weights : vector/scalar of weigths, i.e. inverse variance at each Y(i) (default 1)
% .alpha   : Confidence coefficent             (default 0.05)
% .deletecolinear : If true delete colinear covarites (default)
%
% REGLM performs multiple Linear Regression using weighted Least Squares 
% Fit of y on X. The underlying model for the regression is
%
%  y = X * beta + E
%
% where beta is a column vector of regression parameters and e is a column
% vector of random errors 
%
% Plot r and rint to visualize the residual intervals and identify outliers.
%
% NaN values in y and X are removed before calculation begins.
%
% Example 
% 
%  x=(1:10)';  % Covariate
%  y= x+randn(10,1);
%  b = reglm(y,x);
%  b.display() % members and methods
%  b.get()     % return members
%  b.summary()
%  [y1,ylo,yup] = b.predict(x); 
%  plot(x, y,'o',x,y1,x,[ylo,yup],'r','LineWidth',2)
%
%  b2 = reglm(y,[x,x.^2]);
%  b2.compare(b)
%
%  close all;
% 
% See also regsteplm, reglm>predict regglm, reglogit, regnonlm





% Copyright (C) 2005, 2006 William Poetra Yoga Hadisoeseno
%
% REGLM is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3, or (at your option)
% any later version.
%
% REGLM is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Octave; see the file COPYING.  If not, write to the Free
% Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
% 02110-1301, USA.



% References:
% - Matlab 7.0 documentation (pdf)
% - ����ѧ��ѧʵ�顷 ����Դ �� (textbook)
% - http://www.netnam.vn/unescocourse/statistics/12_5.htm
% - wsolve.m in octave-forge
% - http://www.stanford.edu/class/ee263/ls_ln_matlab.pdf

% History
% 
% revised pab dec 2007
% translated from octave to matlab.
% renamed from regress to reglm

function methods_ = reglm (y, X, varargin)
%[b, bint, r, rint, stats] = regress (y, X, alpha)

options = struct('constant',true,'offset',0,'weights',1,'alpha',0.05,'deletecolinear',true);
if nargin==1 && strcmpi(y,'defaults')
  methods_ = options;
  return
end
%error(nargchk(1,inf,nargin))
narginchk(1,inf)
options = parseoptions(options,varargin{:});

if nargin<2 ||isempty(X)
  X = zeros(size(y,1),0);
end

members = mkmemberstruct();
members.options = options;
members.dispersn = 1;
members.date = datestr(now);


check_input()
compute_lm()

methods_.display = @display_;
methods_.fieldnames = @fieldnames_;
methods_.get = @get_;
%methods_.set = @set_;
methods_.predict = @predict;
methods_.summary = @summary;
methods_.compare = @compare;

%% nested functions
  function check_input()
    
    if (~isnumeric(options.weights) || any(options.weights<0))
      error('Negative weights not allowed.')
    end
    
    if (~isnumeric(X))
      error('WAFO:REGLM', 'X must be a numeric matrix');
    else
      [n1,p]=size(X);
      if p>n1
        X=X';
      end
    end



    y = y(:); % make sure y is a column vector
    n = numel(y);
   
    if (n ~= size(X,1))
      error ('WAFO:REGLM','y and X must contain the same number of rows');
    end
    if options.constant==1
      X=[ones(n,1),X];
    end

    offset = options.offset;
    weights = options.weights;

    notnans = ~(any(~isfinite([y, X]), 2) | ~isfinite(offset) | ~isfinite(weights) | weights ==0);
    y = y(notnans,:);
    X = X(notnans,:);
    if numel(offset)==n
      offset = offset(notnans);
      options.offset = offset;
    end
    if numel(weights)==n
      weights = weights(notnans);
      options.weights = weights;
    end

    %% Make sure X is full rank
    s = svd(X);
    tol = max(size(X)) * eps(max(s));
    ix = find(s<=tol);
    if any(ix) && options.deletecolinear
      X(:,ix) = [];
      %p = size(X,2);
      txt = sprintf(' %d,',ix);
      txt(end) = '';
      warning('WAFO:REGLM','Covariate matrix is singular. Removing column(s):%s',txt)
    end
  end % check_XY

  function compute_lm()
    [n,p]=size(X);
    weights = options.weights;
    offset = options.offset;
    if numel(weights)==n
      weights = weights(notnans);
      options.weights = weights;
    else
      weights = weights(ones(n,1));
    end
    W = spdiags(sqrt(weights),0,n,n);



    if ~options.deletecolinear
      pinv_X = pinv(W*X);
    else
      [Xq Xr] = qr(W*X, 0);
      pinv_X = Xr \ Xq';
    end
    b = pinv_X * W*(y-offset);

 

dof = max(n - p,0);
t_alpha_2 = invt(options.alpha / 2, dof);
H = W*X * pinv_X;

r = (eye (n) - H) * W*(y-offset);
SSE = sum (r .^ 2); % Residual Sum of Squares
v = SSE / dof; %  variance of error term also called dispersion

% c = diag(inv (X' * X)) using (economy) QR decomposition
% which means that we only have to use Xr
bcov = (pinv_X*pinv_X')*v;
%bcov = pinv (Xr' * Xr)*v;
   
bstd =  sqrt(diag(bcov));
bcrit = t_alpha_2 * bstd;

bint = [b + bcrit, b - bcrit];

  

  

dof1 = n - p - 1;
h = diag (H);

% From Matlab's documentation on Multiple Linear Regression,
%   sigmaihat2 = norm (r) ^ 2 / dof1 - r .^ 2 / (dof1 * (1 - h));
%   dr = -tinv (1 - alpha / 2, dof) * sqrt (sigmaihat2 .* (1 - h));
% Substitute
%   norm (r) ^ 2 == sum (r .^ 2) == SSE
%   -tinv (1 - alpha / 2, dof) == tinv (alpha / 2, dof) == t_alpha_2
% We get
%   sigmaihat2 = (SSE - r .^ 2 / (1 - h)) / dof1;
%   dr = t_alpha_2 * sqrt (sigmaihat2 .* (1 - h));
% Combine, we get
%   dr = t_alpha_2 * sqrt ((SSE * (1 - h) - (r .^ 2)) / dof1);

dr = t_alpha_2 * sqrt ((SSE * (1 - h) - (r .^ 2)) / dof1);
rint = [r + dr, r - dr];


% s = sqrt(SSE/dof);
% pupper = 1-(options.alpha)/2;
% plower = (options.alpha)/2;
% Is = sqrt(dof ./ invchi2([pupper plower],dof))*s;


 
SST = sum ((y - mean (y)) .^ 2); % Total sum of squares
  
R2 = 1 - SSE / SST ; % coefficient of determination
R2adj =  max(1 - SSE / SST * (n-1)/(n-p-1),0); % adjusted coefficient of determination

pm1 = max(p-1,1);
%    F = (R2 / (p - 1)) / ((1 - R2) / dof);

MSm = (SST-SSE)/pm1; % Mean square model
MSe = SSE/dof;       % Mean square error
Fstat = MSm/MSe;
%F = dof / (pm1) / (1 / R2 - 1);


 eta = X*b+offset;
 
  members.family = 'normal';
  members.link = 'identity';
  members.numvar    = p;
  members.numobs    = n;
  members.df        = dof;
  members.params    = b.';
  members.params_ci  = bint.';
  
  members.dispersnfit = v;
  
  members.params_cov = bcov;
  members.params_std = bstd.';
  members.params_corr      = bcov./(bstd*bstd.');
  members.params_tstat = (b./sqrt(diag(bcov))).';
  members.params_pvalue = 2*cdft(-abs(members.params_tstat), dof);
  members.mu = eta.';
  members.eta = eta.';
  members.X = X;
  members.W = W;
  members.residual = r.'; 
  members.residual_ci = rint.';
  members.residualT = (r./sqrt(v.*(1-h))).';
  members.deviance = SSE; %sum(r.^2);
  members.SSE = SSE;
  members.SST = SST;
  
  members.R2 = R2;
  members.R2adj = R2adj;
  members.model_fstat = Fstat;
  members.model_pvalue = 1 - cdff(Fstat, pm1, dof);

  
  members.options = options;
  members.note = '';
  members.date = datestr(now);
  
  
 % [b, bint, r, rint, stats]
end

 
 
 %   function set_(varargin)
%     members = parseoptions(members,varargin{:});
%   end
  function [value,varargout] = get_(varargin)
    n = numel(members);
    varargout = cell(1,nargout-1);

    switch nargin
      case {0}
        value = members;
      case {1}
        no = max(min(n,nargout),1);
        [value, varargout{1:no-1}] = deal(members(1:no).(varargin{1}));
      otherwise
        for ix = 1:nargin
          [value(1:n).(varargin{ix})] = deal(members.(varargin{ix}));
        end
    end
  end
  function display_()
    display(members)
    display(methods_)
  end
  function t = fieldnames_()
    t = fieldnames(members);
  end

 
  function  localpvalue = compare(object2)
    %REGLM/Compare  small LM versus large one
%
%   CALL     [pvalue] = compare(object2)
%
%	  The standard hypothesis test of a larger linear regression 
%	  model against a smaller one. The standard F-test is used.
%	  The output is the p-value, the residuals from the smaller 
%	  model, and the residuals from the larger model.
%
% Example
%  x=(1:10)';  % Covariate
%  y= x+randn(10,1);
%  b = reglm(y,x);
%  b2 = reglm(y,[x,x.^2]);
%  b2.compare(b)
%
%	See also reglm
    %error(nargchk(1,1,nargin))
    narginchk(1,1)
    try
    fn = object2.fieldnames();
    if any(~ismember({'family','link','options','dispersnfit','deviance','X','df','numvar'},fn));
      error('Apparently not a valid regression object: %s',inputname(1))
    end
    members2 = object2.get();
    catch
       error('Apparently not a valid regression object: %s',inputname(1))
    end
    if members.numvar>members2.numvar
      devL = members.deviance;
      nL   = members.numvar;
      dfL = members.df;
      Al = members.X;
      disprsn = members.dispersnfit;
      devs = members2.deviance;
      ns   = members2.numvar;
      dfs = members2.df;
      As = members2.X;
    else
      devL = members2.deviance;
      nL   = members2.numvar;
      dfL = members2.df;
      Al = members2.X;
      disprsn = members2.dispersnfit;
      devs = members.deviance;
      ns   = members.numvar;
      dfs = members.df;
      As = members.X;
    end
    if any(any((As-Al*(Al\As))>500*eps)) || ~strcmpi(members2.family,members.family) || ~strcmpi(members2.link,members.link)
      warning('WAFO:REGGLM','Small model not included in large model, result is rubbish!')
    end
    
    pmq = abs(nL-ns);
    disp(' ')
    disp('                       Analysis of Deviance')
    if true %options.estdispersn   
      localstat = abs(devL-devs)/disprsn/pmq;
      localpvalue = 1-cdff(localstat,pmq,dfL);
      disp('Model    DF      Residual deviance      F-stat        Pr(>F)')
    else
      localstat = abs(devL-devs)/disprsn;
      localpvalue = 1-cdfchi2(localstat,pmq);
       disp('Model    DF      Residual deviance      Chi2-stat        Pr(>Chi2)')
    end
    
   
    fprintf('Small    %d       %12.4f       %12.4f    %12.4f \n',dfs,devs,localstat,localpvalue)
    fprintf('Full     %d       %12.4f \n',dfL,devL)
    disp(' ')
     if  ( strcmpi(members.family,'binomial') || strcmpi(members.family,'poisson'))
        warning('WAFO:REGGLM','using F-test with a %s distribution is inappropriate',family)
%       elseif ~options.estdispersn
%         warning('WAFO:REGGLM','using F-test with a fixed dispersion is inappropriate')
      end
    
    
  end
 
  function anova()
  
disp(' ')
disp('                       Analysis of Variance')
disp('Source    DF      Sum of Sqrs           Mean Sqr             F-stat            Pr(>F)')
fprintf('Model     %d       %12.4f        %12.4f    %12.4f    %12.4f \n',members.numvar-1,members.SST-members.SSE,(members.SST-members.SSE)/(members.numvar-1),members.model_fstat,members.model_pvalue)
fprintf('Error     %d       %12.4f        %12.4f \n',members.df,members.SSE,members.SSE/(members.df))
fprintf('Total     %d       %12.4f \n',members.numobs-1,members.SST)
disp(' ')
fprintf(' R2 =  %2.4f,     R2adj = %2.4f\n',members.R2,members.R2adj)
disp(' ')
  end
  function summary()
  
disp('Call:')
fprintf('lm(formula = y ~ x, family = %s)\n',members.family)
disp(' ')
disp('Deviance Residuals:')
disp('    Min       1Q         Median       3Q        Max  ')
fprintf('%2.4f     ',percentile(members.residual,[0 0.25 0.5 0.75 1]))
disp(' ')
disp(' Coefficients:')

disp('            Estimate      Std. Error     t value       Pr(>|t|)')
almat = [members.params(:) members.params_std(:),members.params_tstat(:),members.params_pvalue(:)];
if options.constant
  fprintf('(Intercept)  %2.4f        %2.4f        %2.4f        %2.4f\n',almat(1,:))
end
for iz1 = 1:members.numvar-options.constant
  fprintf('x_%d          %2.4f        %2.4f        %2.4f        %2.4f\n',iz1,almat(iz1+options.constant,:))
end
disp(' ')
fprintf('(Dispersion parameter for %s family taken to be %2.2f)\n',members.family,members.dispersnfit)
disp(' ')
% if options.constant
%   fprintf('    Null deviance: %2.4f  on %d  degrees of freedom\n',members.devianceN,members.dfN)
% end
fprintf('Residual deviance: %2.4f  on %d  degrees of freedom\n',members.deviance,members.df)

anova()

end % summary
 
 function [y,ylo,yup] = predict(Xnew,varargin)
%REGLM/PREDICT Predict from a fitted LM object
%
%  CALL [y,ylo,yup] = predict(Xnew,options)
%
%  y        = predicted value
%  ylo,yup  = 100(1-alpha)% confidence interval for y
%  
%  Xnew     =  new covariate
%  options  = options struct defining the calculation
%         .alpha : confidence coefficient (default 0.05) 
%         .cimean : if true return confidence interval for the mean (default true) 
%                   otherwise return confidence interval for a new value
%
% Example
%  x=(1:10)';  % Covariate
%  y= x+rand(10,1)
%  b = reglm(y,x);
%  [y1,ylo,yup] = b.predict(x); 
%  plot(x, y,'o',x,y1,x,[ylo,yup],'r','LineWidth',2)
% 
% See also reglm

opts = struct('alpha',0.05,'cimean',true);
if nargin==2 && strcmpi(Xnew,'defaults')
  y = opts;
  return
end
%error(nargchk(0,inf,nargin))
narginchk(0,inf)
opts = parseoptions(opts,varargin{:});

if nargin<2 || isempty(Xnew)
  Xnew = X;
else
  n = size(Xnew,1);
  if options.constant==1
    Xnew=[ones(n,1),Xnew];
  end
  notnans = ~(any(~isfinite(Xnew), 2) );
  Xnew = Xnew(notnans,:);
end
[n,p] = size(Xnew);

  

if p ~=members.numvar
  error('Number of covariates must match the number of regression coefficients')
end
y = (Xnew*members.params(:)+members.options.offset);


if nargout>1
  pcov =members.params_cov;

  [U S V]=svd(pcov,0);
  R=(U*sqrt(S)*V'); %squareroot of pcov
  %[R,P] = genchol(pcov);
  %R = chol(pcov)
   if opts.cimean
     varxb = sum((Xnew*R).^2,2);
   else
     varxb = sum((Xnew*R).^2,2) + members.dispersnfit;
   end

   crit = -invt(opts.alpha/2,members.df);
 
  ycrit = crit * sqrt(varxb(:));
  ylo = y-ycrit;
  yup = y+ycrit;
end


end %function predict
end % main

 function  compare()
  %REGLM/Compare  small LM versus large one
%
%   CALL     [pvalue] = compare(object2)
%
%	  The standard hypothesis test of a larger linear regression 
%	  model against a smaller one. The standard F-test is used.
%	  The output is the p-value, the residuals from the smaller 
%	  model, and the residuals from the larger model.
%
% Example
%  x=(1:10)';  % Covariate
%  y= x+randn(10,1)
%  b = reglm(y,x);
%  b2 = reglm(y,[x,x.^2]);
%  b2.compare(b)
%
%	See also reglm

 end

 function predict()
%REGLM/PREDICT Predict from a fitted LM object
%
%  CALL [y,ylo,yup] = predict(Xnew,options)
%
%  y        = predicted value
%  ylo,yup  = 100(1-alpha)% confidence interval for y
%  
%  Xnew     =  new covariate
%  options  = options struct defining the calculation
%         .alpha : confidence coefficient (default 0.05) 
%         .predictmean : default 
%
% Example
%  x=(1:10)';  % Covariate
%  y= x+rand(10,1)
%  b = reglm(y,x);
%  [y1,ylo,yup] = b.predict(x); 
%  plot(x, y,'o',x,y1,x,[ylo,yup],'r','LineWidth',2)
% 
% See also reglm

% Do not remove it is a trick in order to get help reglm>predict to work
 end
 
 function members1 = mkmemberstruct()

 
  
fn = {'family','link','options','numvar','numobs','df','params',...
  'params_ci','params_cov','params_std','params_corr','params_tstat',...
  'params_pvalue','mu','eta','X','Y','W','residual','residual_ci',...
  'residualT','deviance','SSE','SST',...
  'dispersnfit','dispersn','R2','R2adj','model_fstat','model_pvalue','note','date'};
fn(2,:) = {[]};
members1 = struct(fn{:});
end

 
 % UNITTEST % 
 % %Longley data from the NIST Statistical Reference Dataset
%  Z = [  60323    83.0   234289   2356     1590    107608  1947
%    61122    88.5   259426   2325     1456    108632  1948
%    60171    88.2   258054   3682     1616    109773  1949
%    61187    89.5   284599   3351     1650    110929  1950
%    63221    96.2   328975   2099     3099    112075  1951
%    63639    98.1   346999   1932     3594    113270  1952
%    64989    99.0   365385   1870     3547    115094  1953
%    63761   100.0   363112   3578     3350    116219  1954
%    66019   101.2   397469   2904     3048    117388  1955
%    67857   104.6   419180   2822     2857    118734  1956
%    68169   108.4   442769   2936     2798    120445  1957
%    66513   110.8   444546   4681     2637    121950  1958
%    68655   112.6   482704   3813     2552    123366  1959
%    69564   114.2   502601   3931     2514    125368  1960
%    69331   115.7   518173   4806     2572    127852  1961
%    70551   116.9   554894   4007     2827    130081  1962 ];
%  % Results certified by NIST using 500 digit arithmetic
%  % b and standard error in b
%  V = [  -3482258.63459582         890420.383607373
%    15.0618722713733         84.9149257747669
%    -0.358191792925910E-01    0.334910077722432E-01
%    -2.02022980381683         0.488399681651699
%    -1.03322686717359         0.214274163161675
%    -0.511041056535807E-01    0.226073200069370
%    1829.15146461355         455.478499142212 ];
%  Rsq = 0.995479004577296;
%  F = 330.285339234588;
%  y = Z(:,1); X = (Z(:,2:end));
%  b = reglm(y, X);
%  b.summary()
%  res = b.get();
%  plot(res.residual), hold on, plot(res.residualci') % residual plot
%  
% assert(all(abs(res.params'-V(:,1))<3e-6));
% assert(abs(res.R2-Rsq)<1e-12);
% assert(abs(res.Fstat-F)<3e-8);
% assert(abs(res.params_std-V(:,2))./V(:,2)<1.e-10);


%      .df           : degrees of freedom for error.
%      .params       : estimated model parameters.
%      .params_ci    : 100(1-alpha)% confidence interval for model  parameters
%      .params_tstat : t statistics for model's estimated parameters.
%      .params_pvalue: p value for model's estimated parameters.
%      .params_std   : standard errors for estimated parameters
%      .params_corr  : correlation matrix for estimated parameters.
%      .mu           : fitted values for the model.
%      .eta          : linear predictor for the model.
%      .residual     : residual for the model (Y-phat.mu).
%      .residualT    : standardized residual
%      .recidual_ci  : 100(1-alpha)% confidence interval for residual.
%      .deviance     : deviance for the model.
%      .dispersnfit  : The estimated error variance
%      .R2           : R^2 statistic
%      .model_fstat  : F statistic
%      .model_pvalue : p value for the full model
