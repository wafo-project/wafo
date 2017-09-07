function methods_ = regnonlm(x,y,pin,F,varargin)
%REGNONLM Non-Linear Model Regression   
%
% CALL model = regnonlm(x,y,pin,F,options)
%
% model = fitted NONLM object with methods
%      .predict()   : Predict from a fitted NONLM object
%      .summary()   : Display summary of fitted NONLM object.
%      .compare()   : Compare small NONLM object versus large one
%
%       x = column vector or matrix of independent variables, 1 row per
%           observation: x = [x0 x1....xm].
%       y = vector of observed values, same number of rows as x.
%     pin = vector of initial parameters to be adjusted by REGNONLM.
%           All Zero guesses not acceptable.
%       F = name of function or function handle of the form y=f(x,p),
%           with y, x, p of the form y, x, pin as described above.
% options = options structure defining performance of REGNONLM with fieldnames:
%    .weights : vector of statistical weights proportional to 1/var(y). The
%               constant of  proportionality will be estimated. (default = 1)
%    .dp      : fractional increment of p for numerical partial derivatives;
%               default = .001*ones(size(pin))
%              dp(j) > 0 means central differences on j-th parameter p(j).
%              dp(j) < 0 means one-sided differences on j-th parameter p(j).
%              dp(j) = 0 holds p(j) fixed i.e. the initial guess: pin(j)
%    .jacobian: partial derivative function of the form prt=dfdp(x,f,p,dp,F)
%               (default is 'dfdp' a slow but general partial derivatives
%               function; (see regnonlm>dfdp details)
%    .accuracy: accuracy in convergence for the improvement in scalar sum of
%               squares = sum((wt.*(y-f))^2); (default 0.0001)
%    .maxiter : scalar maximum number of iterations. (default 100)
%    .pprec   : desired precision in parameter estimates. Iterations are
%               terminated if  all(abs(P-Pprevious) < pprec)  on two
%               consecutive iterations.], (default zeros()). 
%    .maxstep : maximum step change in parameter vector between successive
%               iterations, ie. abs(chg(i))=abs(min([chg(i) maxstep.*current param estimate])).],
%               (default  inf*ones())
%
% REGNONLM provide a framework for fitting non-linear models (NONLM) to data.
% It is assumed that each outcome of the dependent variables, Y, is
% generated from a Normal-distribution, where the mean, MU, of the distribution 
% depends on the independent variables, X,  through:   
%  
%      E(Y) = MU = F(X,P) 
%  
% where F and P is the inverse link-function and the unknown parameters,
% respectively.  The parameters are fitted by Levenberg-Marquardt
% algorithm (LMA).  
%  
% Example
%  x=(1:10)';  % Covariate
%  y= x+randn(10,1)
%  b = reglm(y,x);
%  fun1 = @(x,p) p(1) + p(2)*x;
%  fun2 = @(x,p) p(1) + p(2)*x+p(3)*x.^2;
%  pin = [0.5 1]'; % initial guess
%  b1 = regnonlm(x,y,pin,fun1);
%  b2 = regnonlm(x,y,[pin;0],fun1);
%  [y0,ylo,yup] = b.predict(x); 
%  [y1,ylo1,yup1] = b1.predict(x); 
%  plot(x, y,'o',x,y0,x,[ylo,yup],'r','LineWidth',2), hold on
%  plot(x, y,'o',x,y1,'.',x,[ylo1,yup1],'r.','LineWidth',2), hold off
%
%  b.summary()
%  b1.summary()  % Compare output
%
%  b1.compare(b2) % compare models
%  
% See also regglm, reglm, reglogit, regnonlm>dfdp

% Copyright (C) 1992-2007 Richard Shrager, Arthur Jutan, Ray Muzic, Per A. Brodtkorb
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301
% USA



% A modified version of Levenberg-Marquardt
% Non-Linear Regression program previously submitted by R.Schrager.
% This version corrects an error in that version and also provides
% an easier to use version with automatic numerical calculation of
% the Jacobian Matrix. In addition, this version calculates statistics
% such as correlation, etc....
%
% Version 3 Notes
% Errors in the original version submitted by Shrager (now called version 1)
% and the improved version of Jutan (now called version 2) have been corrected.
% Additional features, statistical tests, and documentation have also been
% included along with an example of usage.  BEWARE: Some the the input and
% output arguments were changed from the previous version.
%
%     Ray Muzic     <rfm2@ds2.uh.cwru.edu>
%     Arthur Jutan  <jutan@charon.engga.uwo.ca>

% Richard I. Shrager (301)-496-1122
% Modified by A.Jutan (519)-679-2111
% Modified by Ray Muzic 14-Jul-1992
%       1) add maxstep feature for limiting changes in parameter estimates
%          at each step.
%       2) remove forced columnization of x (x=x(:)) at beginning. x could be
%          a matrix with the ith row of containing values of the 
%          independent variables at the ith observation.
%       3) add verbose option
%       4) add optional return arguments covp, stdresid, chi2
%       5) revise estimates of corp, stdev
% Modified by Ray Muzic 11-Oct-1992
%	1) revise estimate of Vy.  remove chi2, add Z as return values
% Modified by Ray Muzic 7-Jan-1994
%       1) Replace ones(x) with a construct that is compatible with versions
%          newer and older than v 4.1.
%       2) Added global declaration of verbose (needed for newer than v4.x)
%       3) Replace return value var, the variance of the residuals with covr,
%          the covariance matrix of the residuals.
%       4) Introduce options as 10th input argument.  Include
%          convergence criteria and maxstep in it.
%       5) Correct calculation of xtx which affects coveraince estimate.
%       6) Eliminate stdev (estimate of standard deviation of parameter
%          estimates) from the return values.  The covp is a much more
%          meaningful expression of precision because it specifies a confidence
%          region in contrast to a confidence interval..  If needed, however,
%          stdev may be calculated as stdev=sqrt(diag(covp)).
%       7) Change the order of the return values to a more logical order.
%       8) Change to more efficent algorithm of Bard for selecting epsL.
%       9) Tighten up memory usage by making use of sparse matrices (if 
%          MATLAB version >= 4.0) in computation of covp, corp, stdresid.
% Modified by Francesco Potortì
%       for use in Octave
% Revised pab 2007
% changed name from leasqr to regnonlm
% -replaced input options with a struct.
% -Added predict, summary and compare subfunctions
% -output is now one struct object


% References:
% Bard, Nonlinear Parameter Estimation, Academic Press, 1974.
% Draper and Smith, Applied Regression Analysis, John Wiley and Sons, 1981.
%
%set default args


% stol,niter,wt,dp,dFdp,options
% [f,p,kvg,iter,corp,covp,covr,stdresid,Z,R2]
% [f,p,kvg,iter,corp,covp,covr,stdresid,Z,R2]= ...
%
% covr = diag(covariance matrix of the residuals).
% Z = matrix that defines confidence region (see comments in the source).


options = struct('weights',1,'jacobian',@dfdp,'dp',0.001,'maxiter',100,'accuracy',0.0001,'pprec',0,'maxstep',inf,'alpha',0.05,'verbose',false);
 
% argument processing
if nargin==1 && strcmpi(x,'defaults')
  methods_ = options;
  return
end
%error(nargchk(4,inf,nargin))
narginchk(4,inf)
options = parseoptions(options,varargin{:});

members = mkmemberstruct();
members.options = options;
members.dispersn = 1;
members.date = datestr(now);


check_input()
compute_nonlm()

methods_.display = @display_;
methods_.fieldnames = @fieldnames_;
methods_.get = @get_;
%methods_.set = @set_;
methods_.predict = @predict;
methods_.summary = @summary;
methods_.compare = @compare;

%% Nested functions
  function check_input()

y=y(:); pin=pin(:); %change all vectors to columns
% check data vectors- same length?
m=length(y); n=length(pin);
[m1,m2]=size(x);
if m1~=m ,
  error('Data x and y must have same number of rows!') 
end;
defval = [0.001 0 inf];
fnames = {'dp','pprec','maxstep'};
for ix= 1:3
  fn = fnames{ix};
  ndp = numel(options.(fn));
  if ndp==0
    options.(fn) = defval(ix);
  elseif ndp~=n && ndp>1
    error('Options: %s must have same number of rows as the parameter matrix, PIN!',fn),
  end
end

  ndp = numel(options.weights);
  if ndp==0
    options.weigths = 1;
  elseif ndp~=m && ndp>1
    error('Options: %s must have same number of rows as the response matrix, Y!','weights'),
  end
  end % check_input


  function compute_nonlm()
    
% set up for iterations
%
m = numel(y);
n=length(pin);
wt = zeros(m,1);
pprec = zeros(n,1);
dp    = zeros(n,1);
maxstep = zeros(n,1);
wt(:) = sqrt(options.weights(:));
pprec(:) = options.pprec(:);
dp(:)    = options.dp(:);
maxstep(:) = options.maxstep(:);
p=pin;
 
f=feval(F,x,p);
fbest=f; 
pbest=p;
r  = wt.*(y-f);
ss = r'*r; % sum of squares
sbest=ss;
nrm=zeros(n,1);
chgprev=Inf*ones(n,1);
kvg=0;
epsLlast=1;
epstab=[.1, 1, 1e2, 1e4, 1e6];

dFdp = options.jacobian;
% do iterations
%
for iter=1:options.maxiter,
  pprev = pbest;
  prt   = feval(dFdp,x,fbest,pprev,dp,F);
  r     = wt.*(y-fbest);
  sprev=sbest;
  sgoal=(1-options.accuracy)*sprev;
  for j=1:n,
    if dp(j)==0,
      nrm(j)=0;
    else
      prt(:,j)=wt.*prt(:,j);
      nrm(j)=prt(:,j)'*prt(:,j);
      if nrm(j)>0,
        nrm(j)=1/sqrt(nrm(j));
      end;
    end
    prt(:,j)=nrm(j)*prt(:,j);
  end;
% above loop could ? be replaced by:
% prt=prt.*wt(:,ones(1,n)); 
% nrm=dp./sqrt(diag(prt'*prt)); 
% prt=prt.*nrm(:,ones(1,m))';
  [prt,s,v]=svd(prt,0);
  s=diag(s);
  g=prt'*r;
  for jjj=1:length(epstab),
    epsL = max(epsLlast*epstab(jjj),1e-7);
    se=sqrt((s.*s)+epsL);
    gse=g./se;
    chg=((v*gse).*nrm);
%   check the change constraints and apply as necessary
    ochg=chg;
    idx = ~isinf(maxstep);
    limit = abs(maxstep(idx).*pprev(idx));
    chg(idx) = min(max(chg(idx),-limit),limit);
    if (options.verbose && any(ochg ~= chg)),
      disp(['Change in parameter(s): ', ...
         sprintf('%d ',find(ochg ~= chg)), 'were constrained']);
    end;
    aprec=abs(pprec.*pbest);       %---
% ss=scalar sum of squares=sum((wt.*(y-f))^2).
    if (any(abs(chg) > 0.1*aprec)),%---  % only worth evaluating function if
      p=chg+pprev;                       % there is some non-miniscule change
      f=feval(F,x,p);
      r=wt.*(y-f);
      ss=r'*r;
      if ss<sbest,
        pbest=p;
        fbest=f;
        sbest=ss;
      end;
      if ss<=sgoal,
        break;
      end;
    end;                          %---
  end;
  epsLlast = epsL;
  if (options.verbose),
    plot(x(:,1),y,''+'',x(:,1),f); shg
    %eval(plotcmd);
  end;
  if ss<eps,
    break;
  end
  aprec=abs(pprec.*pbest);
%  [aprec, chg, chgprev]
  if (all(abs(chg) < aprec) && all(abs(chgprev) < aprec)),
    kvg=1;
    if (options.verbose),
      fprintf('Parameter changes converged to specified precision\n');
    end;
    break;
  else
    chgprev=chg;
  end;
  if ss>sgoal,
    break;
  end
end

% set return values
%
p=pbest;
f=fbest;
%ss=sbest;
kvg=((sbest>sgoal)|(sbest<=eps)|kvg);
if kvg ~= 1 , disp(' CONVERGENCE NOT ACHIEVED! '), end;

% CALC VARIANCE COV MATRIX AND CORRELATION MATRIX OF PARAMETERS
% re-evaluate the Jacobian at optimal values
jac=feval(dFdp,x,f,p,dp,F);
msk = dp ~= 0;
n = sum(msk);           % reduce n to equal number of estimated parameters
jac = jac(:, msk);	% use only fitted parameters

%% following section is Ray Muzic's estimate for covariance and correlation
%% assuming covariance of data is a diagonal matrix proportional to
%% diag(1/wt.^2).  
%% cov matrix of data est. from Bard Eq. 7-5-13, and Row 1 Table 5.1 

if exist('sparse','file')  % save memory
  Q=sparse(1:m,1:m,1./wt.^2);
  Qinv=sparse(1:m,1:m,wt.^2);
else
  Q=diag((0*wt+1)./(wt.^2));
  Qinv=diag(wt.*wt);
end
resid=y-f;                                    %un-weighted residuals
covr=resid'*Qinv*resid*Q/(m-n);               %covariance of residuals
Vy=covr/(1-n/m);  % Eq. 7-13-22, Bard         %covariance of the data 

jtgjinv=pinv(jac'*Qinv*jac);			%argument of inv may be singular
covp=jtgjinv*jac'*Qinv*Vy*Qinv*jac*jtgjinv; % Eq. 7-5-13, Bard %cov of parm est
stdp=sqrt(abs(diag(covp)));
corp=covp./(stdp*stdp');

% if exist('sparse','file')
%   covr=spdiags(covr,0);
%   stdresid=resid./sqrt(spdiags(Vy,0));
%else
  covr=diag(covr);                 % convert returned values to compact storage
  stdresid=resid./sqrt(diag(Vy));  % compute then convert for compact storage
%end
Z=((m-n)*jac'*Qinv*jac)/(n*resid'*Qinv*resid);

%%% alt. est. of cov. mat. of parm.:(Delforge, Circulation, 82:1494-1504, 1990
%%disp('Alternate estimate of cov. of param. est.')
%%acovp=resid'*Qinv*resid/(m-n)*jtgjinv

%Calculate R^2 (Ref Draper & Smith p.46)
%
r=corrcoef([y(:),f(:)]);
R2=r(1,2).^2;
R2adj =  max(1 - (1-R2)* (m-1)/(m-n-1),0); % adjusted coefficient of determination

% if someone has asked for it, let them have it
%
if (options.verbose), 
  plot(x(:,1),y,''+'',x(:,1),f); shg
  %eval(plotcmd);
  disp(' Least Squares Estimates of Parameters')
  disp(p')
  disp(' Correlation matrix of parameters estimated')
  disp(corp)
  disp(' Covariance matrix of Residuals' )
  disp(covr)
  disp(' Correlation Coefficient R^2')
  disp(R2)
  sprintf(' 95%% conf region: F(0.05)(%.0f,%.0f)>= delta_pvec''*Z*delta_pvec',n,m-n)
  Z
%   runs test according to Bard. p 201.
  n1 = sum((f-y) < 0);
  n2 = sum((f-y) > 0);
  nrun=sum(abs(diff((f-y)<0)))+1;
  if ((n1>10)&&(n2>10)), % sufficent data for test?
    zed=(nrun-(2*n1*n2/(n1+n2)+1)+0.5)/(2*n1*n2*(2*n1*n2-n1-n2)...
      /((n1+n2)^2*(n1+n2-1)));
    if (zed < 0),
      prob = erfc(-zed/sqrt(2))/2*100;
      disp([num2str(prob),'% chance of fewer than ',num2str(nrun),' runs.']);
    else
      prob = erfc(zed/sqrt(2))/2*100;
      disp([num2str(prob),'% chance of greater than ',num2str(nrun),' runs.']);
    end
  end
end
dev = sbest;
df = max(m-n,0);
%if options.estdispersn %dispersion_parameter=='mean_deviance'
psi_value = dev/df; % better estimate for psi_value?

members.family    = 'normal';
members.link      = 'nonlinear';
members.numvar    = n;
members.numobs    = m;
members.df        = df;
members.params    = p.';
members.params_ci  = 1;
members.params_cov       = covp;
members.params_std = stdp.';
members.params_corr      = corp;

members.params_tstat=(members.params./members.params_std);
if true %options.estdispersn %dispersion_parameter=='mean_deviance'
  members.params_pvalue=2.*cdft(-abs(members.params_tstat),members.df);
  bcrit = -stdp.'*invt(options.alpha/2,members.df);
else
  members.params_pvalue=2.*cdfnorm(-abs(members.params_tstat));
  bcrit = -stdp.'*invnorm(options.alpha/2);
end
members.params_ci = [members.params+bcrit;members.params-bcrit];

%eta = X*members.params+offset;
%mu = clink{3}(eta);
members.mu    = f.';
members.eta   = [];
members.X     = x;
%members.W     = W;
members.residual= resid.';
members.residualT = stdresid.';
members.residualP=resid.';
members.residualD=resid.';
members.residualA=resid.';
members.residualV = covr.'; % variance of residuals.
members.deviance =dev;
members.deviance_null = [];
members.df_null = [];
members.dispersnfit=psi_value;
members.dispersn = 1;
members.R2 = R2;
members.R2adj = R2adj;
members.numiter = iter;
members.convergence = kvg;
members.options = options;

members.note = '';
members.date = datestr(now);
  end % compute_nonlm


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
%Compare  small NONLM versus large one
%
%   CALL     [pvalue] = compare(object2)
%
%	  The standard hypothesis test of a larger linear regression 
%	  model against a smaller one. The standard F-test is used.
%	  The output is the p-value, the residuals from the smaller 
%	  model, and the residuals from the larger model.
%
%	  See also regglm
    %error(nargchk(1,1,nargin))
    narginchk(1,1)
    try
    fn = object2.fieldnames();
    if any(~ismember({'family','link','dispersnfit','deviance','X','df','numvar'},fn));
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
      warning('WAFO:REGNONLM','Small model not included in large model, result is rubbish!')
    end
    
    pmq = abs(nL-ns);
    disp(' ')
    disp('                       Analysis of Deviance')
    if true %options.estdispersn   
      localstat = abs(devL-devs)/disprsn/pmq;
      localpvalue = 1-cdff(localstat,pmq,dfL);
      disp('Model    DF      Residual deviance      F-stat        Pr(>F)')
    else
      localstat = abs(deviL-devs)/disprsn;
      localpvalue = 1-cdfchi2(localstat,pmq);
       disp('Model    DF      Residual deviance      Chi2-stat        Pr(>Chi2)')
    end
    
   
    fprintf('Small    %d       %12.4f       %12.4f    %12.4f \n',dfs,devs,localstat,localpvalue)
    fprintf('Full     %d       %12.4f \n',dfL,devL)
    disp(' ')
  end

 function anode()
    disp(' ')
    disp('                       Analysis of Deviance')
    if true %options.estdispersn
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
     if options.estdispersn && ( strcmpi(family,'binomial') || strcmpi(family,'poisson'))
        warning('WAFO:REGNONLM','using F-test with a %s distribution is inappropriate',family)
%       elseif ~options.estdispersn
%         warning('WAFO:REGGLM','using F-test with a fixed dispersion is inappropriate')
      end
  end
  function summary()
    
    disp('Call:')
    fprintf('nonlm(formula = y ~ f(x,p), family = %s)\n',members.family)
    disp(' ')
    disp('Deviance Residuals:')
    disp('    Min       1Q         Median       3Q        Max  ')
    fprintf('%2.4f     ',percentile(members.residual,[0 0.25 0.5 0.75 1]))
    disp(' ')
    disp(' Coefficients:')
    if true %options.estdispersn
      disp('            Estimate      Std. Error     t value       Pr(>|t|)')
    else
      disp('            Estimate      Std. Error     z value       Pr(>|z|)')
    end
    almat = [members.params(:) members.params_std(:),members.params_tstat(:),members.params_pvalue(:)];
    
    for ix1 = 1:members.numvar
      fprintf('p_%d          %2.4f        %2.4f        %2.4f        %2.4f\n',ix1,almat(ix1,:))
    end
    disp(' ')
    fprintf('(Dispersion parameter for %s family taken to be %2.2f)\n',members.family,members.dispersnfit)
    disp(' ')

    fprintf('Residual deviance: %2.4f  on %d  degrees of freedom\n',members.deviance,members.df)



end % summary

function [y,ylo,yup] = predict(Xnew,varargin)
%NONLM/PREDICT Predict from a fitted NONLM object
%
%  CALL [y,ylo,yup] = predict(Xnew,options)
%
%  y        = predicted value
%  ylo,yup  = 100(1-alpha)% confidence interval for y
%  
%  Xnew     =  new covariate
%  options  = options struct defining the calculation
%         .alpha : confidence coefficient (default 0.05)
%      
%
% Example
%  y=[2 0 3 1 5 5 6 9 5 9].';  % number of successes in n(i) trials
%  n=[10 10 10 10 10 10 10 10 10 10].';
%  x=(1:10)';  % Covariate
%  b = glm([y n],x,'binomial','link','loglog');
%  plot(x, y./n,'o',x,b.mu'./n,'-','LineWidth',2), shg
%  [y1,ylo,yup] = b.predict(x); % alternative call
%  hold on, plot(x,y1,x,[ylo,yup],'r')
% 
% See also regnonlm

opts = struct('alpha',0.05,'nsim',100,'cimean',true);
if nargin==2 && strcmpi(Xnew,'defaults')
  y = opts;
  return
end
%error(nargchk(0,inf,nargin))
narginchk(0,inf)
opts = parseoptions(opts,varargin{:});


if nargin<1 || isempty(Xnew)
  Xnew = members.X;
else
  %n = size(Xnew,1);
  notnans = ~(any(~isfinite(Xnew), 2) );
  Xnew = Xnew(notnans,:);
end
[n,p] = size(Xnew);

try
  y = feval(F,Xnew,members.params(:));
catch
  error('Number of covariates must match the number of regression coefficients')
end
if nargout>1
  rparams = rndnormnd(members.params,members.params_cov,opts.nsim);
  ry = zeros(n,opts.nsim);
  for ix1 = 1:opts.nsim
    ry(:,ix1) = feval(F,Xnew,rparams(:,ix1));
  end
  yloup = percentile(ry.',[opts.alpha/2 1-opts.alpha/2]).';

  ylo = min(yloup,[],2);
  yup = max(yloup,[],2);
end
end % predict

end %


function prt=dfdp(x,f,p,dp,func)
% numerical partial derivatives (Jacobian) df/dp for use with leasqr
% --------INPUT VARIABLES---------
% x=vec or matrix of indep var(used as arg to func) x=[x0 x1 ....]
% f=func(x,p) vector initialsed by user before each call to dfdp
% p= vec of current parameter values
% dp= fractional increment of p for numerical derivatives
%      dp(j)>0 central differences calculated
%      dp(j)<0 one sided differences calculated
%      dp(j)=0 sets corresponding partials to zero; i.e. holds p(j) fixed
% func=string naming the function (.m) file
%      e.g. to calc Jacobian for function expsum prt=dfdp(x,f,p,dp,'expsum')
%----------OUTPUT VARIABLES-------
% prt= Jacobian Matrix prt(i,j)=df(i)/dp(j)
%================================

% Copyright (C) 1992-1994 Richard Shrager
% Copyright (C) 1992-1994 Arthur Jutan
% Copyright (C) 1992-1994 Ray Muzic
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA



m=size(x,1); if (m==1), m=size(x,2); end  %# PAK: in case #cols > #rows
n=length(p);      %dimensions
ps=p; prt=zeros(m,n);del=zeros(n,1);       % initialise Jacobian to Zero
for j=1:n
      del(j)=dp(j) .*p(j);    %cal delx=fract(dp)*param value(p)
      if p(j)==0
           del(j)=dp(j);     %if param=0 delx=fraction
      end
      p(j)=ps(j) + del(j);
      if del(j)~=0, f1=feval(func,x,p);
           if dp(j) < 0, prt(:,j)=(f1-f)./del(j);
           else
                p(j)=ps(j)- del(j);
                prt(:,j)=(f1-feval(func,x,p))./(2 .*del(j));
           end
      end
      p(j)=ps(j);     %restore p(j)
end
end

function members1 = mkmemberstruct()


fn = {'family','link','options','numvar','numobs','df','df_null','params',...
  'params_ci','params_cov','params_std','params_corr','params_tstat',...
  'params_pvalue','mu','eta','X','Y','W','sample_size','residual',...
  'residualP','residualD','residualA','residualT','residualV','deviance','deviance_null',...
  'dispersnfit','dispersn','R2','R2adj','numiter','convergence','note','date'};
fn(2,:) = {[]};
members1 = struct(fn{:});
end

%      .params    : estimated model parameters.
%      .paramsci  : 100(1-alpha)% confidence interval for model parameters
%      .df        : degrees of freedom for error.
%      .tstat     : t statistics for model's estimated parameters.
%      .pvalue    : p value for model's estimated parameters.
%      .stderrors : standard errors for estimated parameters
%      .corr      : correlation matrix for estimated parameters.
%      .mu        : fitted values for the model, i.e, F(x,phat.params).
%      .residual  : residual for the model (Y-phat.mu).
%      .residualT : standardized residual
%      .deviance  : deviance for the model.
%      .dispersnfit: The estimated error variance
%      .R2        : coefficient of multiple determination.