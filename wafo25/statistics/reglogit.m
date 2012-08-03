%REGLOGIT Fit ordinal logistic regression model.
%
%  CALL model = reglogit (y, x,options)
%
%    model = fitted model object with methods
%      .compare() : Compare small LOGIT object versus large one
%      .predict() : Predict from a fitted LOGIT object
%      .summary() : Display summary of fitted LOGIT object.
%
%       y = vector of K ordered categories
%       x = column vectors of covariates 
% options = struct defining performance of REGLOGIT
%      .maxiter    : maximum number of iterations. 
%      .accuracy   : accuracy in convergence.
%      .betastart  : Start value for BETA           (default 0)
%      .thetastart : Start value for THETA          (default depends on Y)
%      .alpha      : Confidence coefficent          (default 0.05)
%      .print      : 1 display summary info about fitted model
%                    2 display convergence info in each iteration
%                      otherwise no action
%      .deletecolinear : If true delete colinear covarites (default)
%
% Suppose Y takes values in K ordered categories, and let
% gamma_i (x) be the cumulative probability that Y
% falls in one of the first i categories given the covariate
% X.  The ordinal logistic regression model is
%
% logit (mu_i (x)) = theta_i + beta' * x,   i = 1...k-1
%
% The number of ordinal categories, K, is taken to be the number
% of distinct values of round (Y).  If K equals 2,
% Y is binary and the model is ordinary logistic regression.  The
% matrix X is assumed to have full column rank.
%
% Given Y only, theta = REGLOGIT(Y) fits the model with baseline logit odds
% only. 
% 
% Example
% X = sort(5*rand(40, 1))
% y = 2*(sin(X)>rand(40))-1
%  b = reglogit(y,X)
%  b.display() % members and methods
%  b.get()     % return members
%  b.summary()
%  [mu,plo,pup] = b.predict();
%  plot(x,mu,'g',x,plo,'r:',x,pup,'r:')
% 
%  y=[1 1 2 1 3 2 3 2 3 3]'
%  x = (1:10)'
%  b = reglogit(y,x)
%  b.display() % members and methods
%  b.get()     % return members
%  b.summary()
%  [mu,plo,pup] = b.predict();
%  plot(x,mu,'g',x,plo,'r:',x,pup,'r:')
%
%  y2 = [zeros(5,1);ones(5,1)];
%  x1 = [29,30,31,31,32,29,30,31,32,33];
%  x2 = [62,83,74,88,68,41,44,21,50,33];
%  X = [x1;x2].';
%  b2 = reglogit(y2,X);
%  b2.summary();
%  b21 = reglogit(y2,X(:,1));
%  b21.compare(b2)
%
% See also regglm, reglm, regnonlm


% Copyright (C) 1995, 1996, 1997, 1998, 1999, 2000, 2002, 2005, 2007
%               Kurt Hornik
%
% Reglogit is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or (at
% your option) any later version.
%
% Reglogit is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Reglogit; see the file COPYING.  If not, see
% <http://www.gnu.org/licenses/>.



% Original for MATLAB written by Gordon K Smyth <gks@maths.uq.oz.au>,
% U of Queensland, Australia, on Nov 19, 1990.  Last revision Aug 3,
% 1992.

% Author: Gordon K Smyth <gks@maths.uq.oz.au>,
% Adapted-By: KH <Kurt.Hornik@wu-wien.ac.at>
% Revised by: pab
%  -renamed from  logistic_regression to reglogit
%  -added predict, summary and compare
% Description: Ordinal logistic regression

% Uses the auxiliary functions logistic_regression_derivatives and
% logistic_regression_likelihood.

function methods_ = reglogit(y, X, varargin)



options = struct('thetastart',[],'betastart',[],'maxiter',500,'accuracy',1e-6,...
  'alpha',0.05,'deletecolinear',true,'print',false);
if nargin==1 && strcmpi(y,'defaults')
  methods_ = options;
  return
end
error(nargchk(1,inf,nargin))
if (nargin < 2) 
  X = [];
end
options = parseoptions(options,varargin{:});

members = mkmemberstruct();
members.options = options;
members.dispersn = 1;
members.date = datestr(now);

logit = @(p) log(p./(1-p));
logitinv = @(x) 1./(exp(-x)+1);

check_options()
check_XY()
compute_logit()

methods_.display = @display_;
methods_.fieldnames = @fieldnames_;
methods_.get = @get_;
%methods_.set = @set_;
methods_.predict = @predict;
methods_.summary = @summary;
methods_.compare = @compare;

  %% nested functions
   function check_options()
    if ~isnumeric(options.accuracy) || options.accuracy<=0
      error('value of accuracy must be > 0')
    end
    if (~isnumeric(options.maxiter) || options.maxiter <= 0)
      error('maximum number of iterations must be > 0')
    end
        
  end % check_options
  function check_XY()
     % check input
  y = round (y);
  [my, ny] = size (y);
  if isempty(X)
    X = zeros (my, 0);
  else
    %% Make sure X is full rank
    s = svd(X);
    tol = max(size(X)) * eps(max(s));
    ix = find(s<=tol);
    if any(ix) && options.deletecolinear
      X(:,ix) = [];
      p = size(X,2);
      txt = sprintf(' %d,',ix);
      txt(end) = '';
      warning('WAFO:REGLOGIT','Covariate matrix is singular. Removing column(s):%s',txt)
    end
  end
  [mx, nx] = size (X);
  if (mx ~= my)
    error ('x and y must have the same number of observations');
  end
    
  end

  function compute_logit()
%family = 'multinomial';
%link = 'logit';


 
  % initial calculations
  %X = -X;
  tol = options.accuracy; 
  incr = 10; decr = 2;
  ymin = min (y); ymax = max (y); yrange = ymax - ymin;
  z  = (y * ones (1, yrange)) == ((y * 0 + 1) * (ymin : (ymax - 1)));
  z1 = (y * ones (1, yrange)) == ((y * 0 + 1) * ((ymin + 1) : ymax));
  z  = z(:, any (z));
  z1 = z1 (:, any(z1));
  [mz, nz] = size (z);
  [mx, nx] = size (X);
  [my, ny] = size (y);
   
   
  g = cumsum (sum (z))' ./ my;
  theta0 = log (g ./ (1 - g));
  beta0  = zeros (nx, 1);
  
  % starting values
  if isempty(options.betastart) 
    beta = beta0;
  else
    beta = options.betastart;
  end
  
  if isempty(options.thetastart)
    theta = theta0;
  else
     theta = options.thetastart;
  end;
 
 
  tb = [theta; beta];

  % likelihood and derivatives at starting values
  [dev,p,g, g1] = logistic_regression_likelihood (y, X, tb, z, z1);
  [dl, d2l] = logistic_regression_derivatives (X, z, z1, g, g1, p);
  epsilon = std (d2l) / 1000;
  if any(beta) || any(theta~=theta0)
    nulldev = logistic_regression_likelihood (y, X, [theta0;beta0], z, z1);
  else
    nulldev = dev;
  end
  
 
  
  % maximize likelihood using Levenberg modified Newton's method
  iter = 0;
  while (abs (dl' * (d2l \ dl) / length (dl)) > tol && iter<=options.maxiter)
    iter = iter + 1;
    tbold = tb;
    devold = dev;
    tb = tbold - d2l \ dl;
    [dev,p,g, g1] = logistic_regression_likelihood (y, X, tb, z, z1);
    if ((dev - devold) / (dl' * (tb - tbold)) < 0)
      epsilon = epsilon / decr;
    else
      while ((dev - devold) / (dl' * (tb - tbold)) > 0)
        epsilon = epsilon * incr;
         if (epsilon > 1e+15)
           error ('epsilon too large');
         end
         tb = tbold - (d2l - epsilon * eye (size (d2l))) \ dl;
         [dev,p,g, g1] = logistic_regression_likelihood (y, X, tb, z, z1);
         disp ('epsilon'); disp (epsilon);
      end %while
    end
    [dl, d2l] = logistic_regression_derivatives (X, z, z1, g, g1, p);
    if (options.print >= 2)
      if options.print==2
        fprintf('Iter: %d,  Deviance: %8.6f \n',iter,dev)
      else
        disp ('Iteration'); disp (iter);
        disp ('Deviance'); disp (dev);
        disp ('First derivative'); disp (dl');
        disp ('Eigenvalues of second derivative'); disp (eig (d2l)');
      end
    end
  end %while

  % tidy up output

  theta = tb (1 : nz, 1);
  beta  = tb ((nz + 1) : (nz + nx), 1);
  pcov = pinv(-d2l);
  se = sqrt (diag (pcov));
  
  if (options.print >= 1)
    summary()
  end

  
    if (nx > 0)
      eta = ((X * beta) * ones (1, nz)) + ((y * 0 + 1) * theta');
    else
      eta = (y * 0 + 1) * theta';
    end
    gammai = diff ([(y * 0), logitinv(eta), (y * 0 + 1)],1,2);
    k0 = min(y);
    mu = (k0-1)+gammai*(1:nz+1)'; % E(Y|X)
    r  = corrcoef([y,mu]);
    R2 = r(1,2).^2; % % coefficient of determination
    R2adj =  max(1 - (1-R2)* (my-1)/(my-nx-nz-1),0); % adjusted coefficient of determination

    res = y-mu;
    
   if nz==1
     members.family = 'binomial';
   else
     members.family = 'multinomial';
   end
   members.link   = 'logit';
   members.options = options;
   members.numvar = nx+nz;
   members.numobs = my;
   members.numk   = nz+1;
   members.df     = max(my-nx-nz,0);
   members.df_null    = my-nz; %nulldf;  nulldf =  n - nz;
   members.params = tb(1:(nz + nx),1).';
   members.params_ci = 1;
   members.params_std = se.';  
   members.params_cov = pcov;
   members.params_tstat=(members.params./members.params_std);
   if false % options.estdispersn %dispersion_parameter=='mean_deviance'
     members.params_pvalue=2.*cdft(-abs(members.params_tstat),members.df);
     bcrit = -se.'*invt(options.alpha/2,members.df);
   else
     members.params_pvalue=2.*cdfnorm(-abs(members.params_tstat));
     bcrit = -se.'*invnorm(options.alpha/2);
   end
   members.params_ci = [members.params+bcrit;members.params-bcrit];
   
   
   members.mu = gammai;
   members.eta = logit(gammai);
   members.X = X;
   
   members.theta = theta;
   members.beta  = beta;
   members.gamma = gammai;
   members.residual  = res.'; 
   members.residualD = sign(members.residual).*sqrt(-2*log(p)).';
   members.deviance  = dev;
   members.deviance_null = nulldev;
   members.d2L = d2l;
   members.dL = dl.';
   members.dispersnfit=1;
   members.dispersn = 1;
   members.R2 = R2;
   members.R2adj = R2adj;
   members.numiter = iter;
   members.convergence = iter<options.maxiter;
   members.note = '';
   members.date = datestr(now);
   
  end % compute_logit()

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
  %Compare  small LOGIT versus large one
%
%   CALL     [pvalue] = compare(object2)
%
%	  The standard hypothesis test of a larger linear regression 
%	  model against a smaller one. The standard Chi2-test is used.
%	  The output is the p-value, the residuals from the smaller 
%	  model, and the residuals from the larger model.
%
%	  See also fitls  
    error(nargchk(1,1,nargin))
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
    if false %options.estdispersn   
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
%      if options.estdispersn && ( strcmpi(family,'binomial') || strcmpi(family,'poisson'))
%         warning('WAFO:REGGLM','using F-test with a %s distribution is inappropriate',family)
% %       elseif ~options.estdispersn
% %         warning('WAFO:REGGLM','using F-test with a fixed dispersion is inappropriate')
%       end
    
    
    end

function anode()
    disp(' ')
    disp('                       Analysis of Deviance')
    if false %options.estdispersn
      localstat = abs(members.deviance_null-members.deviance)/members.dispersnfit/(members.numvar-1);
      localpvalue = 1-cdff(localstat,members.numvar-1,members.df);
      disp('Model    DF      Residual deviance      F-stat        Pr(>F)')
    else
      localstat = abs(members.deviance_null-members.deviance)/members.dispersnfit;
      localpvalue = 1-cdfchi2(localstat,members.numvar-1);
       disp('Model    DF      Residual deviance      Chi2-stat        Pr(>Chi2)')
    end
    
   
    fprintf('Null     %d       %12.4f       %12.4f    %12.4f \n',members.df_null,members.deviance_null,localstat,localpvalue)
    fprintf('Full     %d       %12.4f \n',members.df,members.deviance)
    disp(' ')

    fprintf(' R2 =  %2.4f,     R2adj = %2.4f\n',members.R2,members.R2adj)
    disp(' ')
%      if options.estdispersn && ( strcmpi(family,'binomial') || strcmpi(family,'poisson'))
%         warning('WAFO:REGGLM','using F-test with a %s distribution is inappropriate',family)
% %       elseif ~options.estdispersn
% %         warning('WAFO:REGGLM','using F-test with a fixed dispersion is inappropriate')
%       end
  end
  function summary()
    txtlink = members.link;
    if ~ischar(txtlink)
      txtlink = func2str(members.link{1});
    end
disp('Call:')
fprintf('reglogit(formula = %s(Pr(grp(y)<=i)) ~ theta_i+beta*x, family = %s)\n',txtlink,members.family)
disp(' ')
disp('Deviance Residuals:')
disp('    Min       1Q         Median       3Q        Max  ')
fprintf('%2.4f     ',percentile(members.residualD,[0 0.25 0.5 0.75 1]))
disp(' ')
disp(' Coefficients:')
if false %options.estdispersn
  disp('            Estimate      Std. Error     t value       Pr(>|t|)')
else
  disp('            Estimate      Std. Error     z value       Pr(>|z|)')
end
almat = [members.params(:) members.params_std(:),members.params_tstat(:),members.params_pvalue(:)];
for ix1 = 1:members.numk-1
  fprintf('theta_%d         %2.4f        %2.4f        %2.4f        %2.4f\n',ix1,almat(ix1,:))
end
for ix1 = members.numk:members.numvar
  fprintf(' beta_%d         %2.4f        %2.4f        %2.4f        %2.4f\n',ix1-members.numk+1,almat(ix1,:))
end
disp(' ')
fprintf('(Dispersion parameter for %s family taken to be %2.2f)\n',members.family,members.dispersnfit)
disp(' ')
if true %options.constant
  fprintf('    Null deviance: %2.4f  on %d  degrees of freedom\n',members.deviance_null,members.df_null)
end
fprintf('Residual deviance: %2.4f  on %d  degrees of freedom\n',members.deviance,members.df)

anode()

end % summary
  
%    function summary()
%      fprintf('\n');
%     fprintf('Logistic Regression Results:\n');
%     fprintf('\n');
%     fprintf('Number of Iterations: %d\n', iter);
%     fprintf('Deviance:             %f\n', dev);
%     fprintf('Parameter Estimates:\n');
%     fprintf ('     Theta         S.E.\n');
%     for i = 1 : nz
%       fprintf('   %8.4f     %8.4f\n', tb (i), se (i));
%     end %for
%     if (nx > 0)
%       fprintf('      Beta         S.E.\n');
%       for i = (nz + 1) : (nz + nx)
%         fprintf('   %8.4f     %8.4f\n', tb (i), se (i));
%       end %for
%     end
%    end
 
 
 function [y,ylo,yup] = predict(Xnew,varargin)
%LOGIT/PREDICT Predict from a fitted LOGIT object
%
%  CALL [y,ylo,yup] = predict(Xnew,options)
%
%  y        = predicted value
%  ylo,yup  = 100(1-alpha)% confidence interval for y
%  
%  Xnew     =  new covariate
%  options  = options struct defining the calculation
%         .alpha : confidence coefficient (default 0.05)
%         .size  : size if binomial family (default 1).    
% 
% See also reglogit

opts = struct('size',1,'alpha',0.05,'type','eta','cimean',true);
if nargin==2 && strcmpi(Xnew,'defaults')
  y = opts;
  return
end
error(nargchk(0,inf,nargin))

opts = parseoptions(opts,varargin{:});


 [mx, nx] = size(members.X);
if nargin<1 || isempty(Xnew)
  Xnew = members.X;
else
  notnans = ~(any(~isfinite(Xnew), 2) );
  Xnew = Xnew(notnans,:);
end
[n,p] = size(Xnew);

  
if p ~= nx
  error('Number of covariates must match the number of regression coefficients')
end

nz = members.numk-1;
one = ones(n,1);
 if (nx > 0)
   eta = ((Xnew * members.beta) * ones (1, nz)) + ( one*members.theta');
 else
   eta = one * members.theta';
 end
 y = diff ([zeros(n,1), logitinv(eta), one],1,2);
if nargout>1
  
  pcov = members.params_cov;
  if (nx > 0)
    np = size(pcov,1);

    [U S V]=svd(pcov,0);
    R=(U*sqrt(S)*V'); %squareroot of pcov
    ixb = nz+1:np;

    % Var(eta_i) = var(theta_i+Xnew*b)
    vareta = zeros(n,nz);
    for ix1 = 1:nz
      vareta(:,ix1) = max(sum(([one Xnew]*R([ix1,ixb],[ix1,ixb])).^2,2),eps);
    end
  else
    vareta = diag(pcov);
  end
  crit = -invnorm(opts.alpha/2);

  
  ecrit = crit * sqrt(vareta);
  mulo = logitinv(eta-ecrit);
  muup = logitinv(eta+ecrit);
  ylo1 = diff ([zeros(n,1), mulo , one],1,2);
  yup1 = diff ([zeros(n,1), muup , one],1,2);
  
  ylo = min(ylo1,yup1);
  yup = max(ylo1,yup1);

  for ix1 = 2:members.numk-1
    yup(:,ix1)  = max( [yup(:,ix1),muup(:,ix1)-mulo(:,ix1-1)],[],2);
  end
end
end % predict
end %main function


function model()
%REGLOGIT>MODEL 
%
% Member variables
%      .options      : input options as explained below.
%      .df           : degrees of freedom for error.
%      .params       : estimated model parameters
%      .params_ci    : 100(1-alpha)% confidence interval for model parameters
%      .params_tstat : t statistics for model's estimated parameters.
%      .params_pvalue: p value for model's estimated parameters.
%      .params_std   : standard errors for estimated parameters
%      .params_corr  : correlation matrix for estimated parameters.
%      .mu           : fitted values for the model.
%      .eta          : linear predictor for the model.
%      .residual     : residual for the model (Y-E(Y|X)).
%      .dispersnfit  : The estimated error variance
%      .deviance     : deviance for the model equal minus twice the log-likelihood.
%      .d2L          : Hessian matrix (double derivative of log-likelihood)
%      .dL           : First derivative of loglikelihood w.r.t. THETA and BETA.
%
% Methods
%      .predict   : Predict from a fitted LOGIT object
%      .summary   : Display summary of fitted LOGIT object.
%      .compare   : Compare small LOGIT versus large one
%      .get       : Return member variables
%
% See also reglogit

% Do not remove 
% This is a trick in order to get help regglm>model to work
end



function [dl, d2l] = logistic_regression_derivatives (x, z, z1, g, g1, p)
% [dl, d2l] = logistic_regression_derivatives (@var{x}, @var{z}, @var{z1}, @var{g}, @var{g1}, @var{p})
% Called by logistic_regression.  Calculates derivates of the
% log-likelihood for ordinal logistic regression model.


% Copyright (C) 1995, 1996, 1997, 1998, 1999, 2000, 2002, 2005, 2007
%               Kurt Hornik
%
%
% Reglogit is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or (at
% your option) any later version.
%
% Reglogit is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Reglogit; see the file COPYING.  If not, see
% <http://www.gnu.org/licenses/>.

% Author: Gordon K. Smyth <gks@maths.uq.oz.au>
% Adapted-By: KH <Kurt.Hornik@wu-wien.ac.at>
% Description: Derivates of log-likelihood in logistic regression


  % first derivative
  v = g .* (1 - g) ./ p; 
  v1 = g1 .* (1 - g1) ./ p;
  dlogp = [(dmult (v, z) - dmult (v1, z1)), (dmult (v - v1, x))];
  dl = sum (dlogp)';

  % second derivative
  w = v .* (1 - 2 * g); w1 = v1 .* (1 - 2 * g1);
  d2l = [z, x]' * dmult (w, [z, x]) - [z1, x]' * dmult (w1, [z1, x]) ...
      - dlogp' * dlogp;

end %function

function [dev,p,g, g1] = logistic_regression_likelihood (y, x, beta, z, z1)
% [g, g1, p, dev] = logistic_regression_likelihood (y ,x,beta,z,z1)
% Calculates likelihood for the ordinal logistic regression model.
% Called by logistic_regression.



% Copyright (C) 1995, 1996, 1997, 1998, 1999, 2000, 2002, 2005, 2007
%               Kurt Hornik
%
%
% Reglogit is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or (at
% your option) any later version.
%
% Reglogit is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Reglogit; see the file COPYING.  If not, see
% <http://www.gnu.org/licenses/>.

% Author: Gordon K. Smyth <gks@maths.uq.oz.au>
% Adapted-By: KH <Kurt.Hornik@wu-wien.ac.at>
% Description: Likelihood in logistic regression
  e = exp ([z, x] * beta); 
  e1 = exp ([z1, x] * beta);
  g = e ./ (1 + e); g1 = e1 ./ (1 + e1);
  g = max (y == max (y), g);
  g1 = min (y > min(y), g1);

  p = g - g1;
  dev = -2 * sum (log (p));

end %function

function M = dmult(A,B)
% PURPOSE: computes the product of diag(A) and B
% -----------------------------------------------------
% USAGE:     m = dmult(a,b)
%  where:    a = a matrix
%            b = a matrix
% -----------------------------------------------------
% RETURNS:  m = diag(A) times B
% -----------------------------------------------------             
% NOTE: a Gauss compatability function
% -----------------------------------------------------

% written by:
%  Gordon K Smyth, U of Queensland, Australia, gks@maths.uq.oz.au
% Nov 19, 1990.  Last revision Aug 29, 1995.

% documentation modifications made by:
% James P. LeSage, Dept of Economics
% University of Toledo
% 2801 W. Bancroft St,
% Toledo, OH 43606
% jpl@jpl.econ.utoledo.edu


[mb,nb] = size(B);
M=A(:,ones(1,nb)).*B;
end


function members1 = mkmemberstruct()

fn = {'family','link','options','numvar','numobs','numk','df','df_null','params',...
  'params_ci','params_cov','params_std','params_corr','params_tstat',...
  'params_pvalue','mu','eta','X','Y','theta','beta','residual',...
  'residualD','deviance','deviance_null','d2L','dL',...
  'dispersnfit','dispersn','R2','R2adj','numiter','convergence','note','date'};
fn(2,:) = {[]};
members1 = struct(fn{:});
end
