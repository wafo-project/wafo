function [methods_] = regglm(Y,X,varargin) 
%REGGLM Generalized Linear Model regression
%
% CALL  model = regglm(Y,X,options) 
%
%   model = fitted GLM model object with methods
%      .compare() : Compare small GLM object versus large one
%      .predict() : Predict from a fitted GLM object
%      .summary() : Display summary of fitted GLM object.
%
%       Y = response variable or dependent variable (n x 1 matrix). In the 
%           Binomial models Y is a n x 2 matrix where the first column holds 
%           the sample values and the second column holds the samples size.
%       X = covariate matrix or independent variables (n x p matrix). n is  
%           the number of observations and p is the number of covariates or
%           independent variables.
% options = struct defining performance of GLM with fieldnames: 
%  .family : distribution for the model (default 'normal'). Options are
%           'binomial','gamma','invnorm','normal' and 'poisson'. 
%  .link   : link function for the model given as numeric value, text string
%            or cell array of three function handles {L LD LI}, that define the 
%            link (L), the derivative of the link (LD), and the inverse link (LI).
%            A Numeric value define the exponent of the power-function.
%            Valid option for text string is: 
%            'identity', 'inverse', 'log', 'loglog', 'sqrt','cloglog','logit','probit'. 
%            The default link function depends on the distribution. Default for
%            'binomial','gamma','invnorm','normal','poisson' is
%            'logit','inverse',-2,'identity','log', respectively.
%            For Binomial models acceptable links are: cloglog, logit, loglog
%            and probit.
% .constant: if TRUE include constant term in the model otherwise omit it.
% .offset  : vector or offset values (default 0)
% .weights : vector/scalar of weigths, i.e. inverse variance at each Y(i) (default 1)
% .estdispersn : If FALSE use theoretical dispersion. (default for poisson and binomial)
%            Otherwise estimate dispersion parameter. (default for other distributions) 
% .maxiter : maximum number of iterations for the IRLS algorithm.
% .accuracy: accuracy in convergence for the IRLS algorithm. (default 0.001)
% .mustart : Start value for E(Y)              (default Y)
% .alpha   : Confidence coefficent             (default 0.05)
% .deletecolinear : If true delete colinear covarites (default)
%  
% REGGLM provide a framework for fitting  a wide class of generalized
% linear models (GLM) to data. In a GLM, each outcome of the dependent
% variables, Y, is assumed to be generated from one of the distributions
% from the exponential family ,such as Binomial-, Poisson-, Normal-, Gamma-
% and Inverse Gaussian- models, where the mean, MU, of the distribution
% depends on the independent variables, X,  through:   
%
%    E(Y) = MU = LI(X*B)   and    L(MU) = ETA = X*B
%
% where LI, L and B is the inverse link-function, link-function and the
% unknown parameters, respectively. ETA is used to denote a linear predictor.
% The parameters are fitted by iteratively  reweighted least squares algorithm
% (IRLS).
%
% NaN and infinite values in Y and X are removed before calculation begins.
%
% Examples
%
% y=[2 0 3 1 5 5 6 9 5 9].';  % number of successes in n(i) trials
% n=[10 10 10 10 10 10 10 10 10 10].';
% x=(1:10)';  % Covariate
% b = regglm([y n],x,'family','binomial','link','logit');
% [y1,ylo,yup] = b.predict(x); 
% b.display() % members and methods
% b.get()     % return members
% b.summary()
% plot(x, y./n,'o',x,y1,'-', x,[ylo,yup],'r:', 'LineWidth',2)
%
% x=(1:10)';  % Covariate
% y = [0 1 0 0 1 0 0 0 1 1]';
% n = ones(size(y));
% b = regglm([y n],[],'family','binomial');
% b1 = regglm(y ,[],'family','binomial'); % same thing
% plot(x, y./n,'o',x,b.get('mu')'./n,'-','LineWidth',2)
%
% y=[1 1 2 1 3 2 3 2 3 3]';
% n = repmat(3,size(y));
% b = regglm([y n],x,'family','binomial'); figure(gcf+1)
% plot(x, y./n,'o',x,b.get('mu')'./n,'-','LineWidth',2)
%
%  x=(1:10)';  % Covariate
%  y= log(x+randn(10,1));
%  b1 = regglm(y,x,'link','log');
%  b2 = regglm(y,[x x.^2 exp(x)],'link','log');
%  b1.compare(b2) % Test if small modell is better than large model
%
% See also regglm>predict, regglm>glmlink, reglm, reglogit, regnonlm

% Copyright (C) 2003, 2007  Stefano Favaro, Per A. Brodtkorb
% 
%     REGGLM is free software; you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation; either version 3 of the License, or
%     (at your option) any later version.
% 
%     REGGLM is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
    


%% References    
% [1] Dobson, A.J. (1990), "An Introduction to Generalized Linear Models", 
%     CRC Press.
% 
% [2] Gill, J. (2001), "Generalized Linear Models: a unified approach", 
%     (quantitative applications in the social sciences, 134). Sage 
%     Pubblications.
% 
% [3] McCullagh, P., and J.A. Nelder (1990), "Generalized Linear Models", 
%     CRC Press.
%      
% [4] Myers, R.M., Montgomery, D.C. and Vining, G.G. (2001), "Generalized 
%     Linear Models with applications in engineering and the science", 
%     John Wiley and Sons.
%
% http://en.wikipedia.org/wiki/Generalized_linear_model#Count_data


% History
% By pab dec 2007
% Based on the functions in GLMBOX written in SciLab by Stefano Favaro (sfavaro@katamail.com)



options = struct('family','normal','link',[],'estdispersn',[],'constant',true,'offset',0,...
  'weights',1,'mustart',[],'maxiter',500,'accuracy',0.001,'alpha',0.05,'deletecolinear',true);
if nargin==1 && strcmpi(Y,'defaults')
  methods_ = options;
  return
end
error(nargchk(2,inf,nargin))
options = parseoptions(options,varargin{:});
members = mkmemberstruct();
members.options = options;
members.dispersn = 1;
members.date = datestr(now);

check_options()
check_XY()
compute_glm()

methods_.display = @display_;
methods_.fieldnames = @fieldnames_;
methods_.get = @get_;
%methods_.set = @set_;
methods_.glmlink = @glmlink;
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
    if (~isnumeric(options.weights) || any(options.weights<0))
      error('Negative weights not allowed.')
    end
    
    validDistributions = {'binomial','gamma','invnorm','normal','poisson'};
    % and corresponding default link
    defaultLinks = {'logit','inverse',-2,'identity','log'};
    idx =  find(strncmpi(options.family,validDistributions,1));
    switch numel(idx)
      case 0
        error('Unsupported family: %s',options.family)
      case 1
        options.family = validDistributions{idx};
      otherwise
        error('Unknown error?')
    end
    if isempty(options.link)
      options.link = defaultLinks{idx};
    end
    if isempty(options.estdispersn)
      switch options.family
        case {'poisson','binomial'}
          options.estdispersn = false;
        otherwise
          options.estdispersn = true;
      end
    end
    validLinksBinomial = {'logit','cloglog','probit','loglog'};
    if strncmpi(options.family,'binomial',1)
      if ischar(options.link) && ~(any(strcmpi(options.link,validLinksBinomial)))
        error('Error in input parameters: link function not correct in this model.')
      end
    elseif ischar(options.link) && any(strcmpi(options.link,validLinksBinomial))
      error('Error in input parameters: link function not usual in this model.')
      
    end
  end % check_options

  function check_XY()
    if isempty(X)
      X = zeros(size(Y,1),0);
    else
      [n1,p]=size(X);
      if p>n1
        X=X';
      end
    end
    if ((size(Y,1)~=size(X,1)))
      error('Error in input parameters: the number of rows in Y must equal the number of rows in X.')
    end

    % INPUT PARAMETERS SET UP
    n = size(Y,1);

    if options.constant==1
      X=[ones(n,1),X];
    end
    offset = options.offset;
    weights = options.weights;

    notnans = ~(any(~isfinite([Y, X]), 2) | ~isfinite(offset) | ~isfinite(weights) | weights ==0);
    Y = Y(notnans,:);
    X = X(notnans,:);
    if numel(offset)==n
      offset = offset(notnans);
      options.offset = offset;
    end
    if numel(weights)==n
      weights = weights(notnans);
      options.weights = weights;
    end
    if numel(options.mustart)>1
      options.mustart = options.mustart(notnans);
    end

    %[n,p]=size(X);


    % Make sure X is full rank
    s = svd(X);
    tol = max(size(X)) * eps(max(s));
    ix = find(s<=tol);
    if any(ix) && options.deletecolinear
      X(:,ix) = [];
      %p = size(X,2);
      txt = sprintf(' %d,',ix);
      txt(end) = '';
      warning('WAFO:REGGLM','Covariate matrix is singular. Removing column(s):%s',txt)
    end



    members.sample_size = 1;
    if strncmpi(options.family,'binomial',1)
      if((size(Y,2)~=2))
        if ~all(Y==0 |Y==1)
          error('Error in input parameters: Y must be a two columns matrix: the first column for the sample value and the second column for the sample size.')
        end
      else
        members.sample_size=Y(:,2);
        if any(members.sample_size<=0)
          error('Sample size must be positive')
        end
      end
      Y=Y(:,1)./members.sample_size;
      if any(Y < 0 | 1<Y)
        error('Y must be positive but less or equal to sample size.')
      end
    else
     
      if((size(Y,2)~=1))
        error('Error in input parameters: Y must be a one column matrix.')
      end
      switch options.family
        case {'gamma' , 'invnorm'}
          if any(Y<=0)
            warning('WAFO:GLM','Y must be positive for the %s distribution',family)
          end
        case 'poisson'
          if any(Y<0)
            warning('WAFO:GLM','Y must be positive for the Poisson distribution')
          end
      end
    end
  end % check_XY

  function mu = mustart()
    if isempty(options.mustart)
      linkLogOrReci = {'log','inverse'};
      if strncmpi(options.family,'binomial',1)
        validLinksBinomial = {'logit','cloglog','probit','loglog'};
        if ischar(options.link)&& any(strcmpi(options.link,validLinksBinomial))
          mu=(members.sample_size.*Y+0.5)./(members.sample_size+1);
        elseif ischar(options.link)&& any(strcmpi(options.link,linkLogOrReci))
          mu=Y+0.5*(Y==0);
        else
          mu=Y;
        end
      elseif  ischar(options.link)&& any(strcmpi(options.link,linkLogOrReci))
        mu=Y+(Y==0);
      else
        mu=Y;
      end
    else
      mu = options.mustart;
    end
  end % mustart

  function compute_glm()

   
    [n,p] = size(X);
    
    [mufun,stdfun] = moment_funs(options.family,members.sample_size);
    
   
    clink = glmlink(options.link,1); %members.sample_size);

    flink = clink{1};%function handle to link function
    dlink = clink{2};%function handle to derivative of link function
    ilink = clink{3};%function handle to inverse link function

    mu = mustart();
    eta = flink(mu);
    BetaV = zeros(p,options.maxiter);
    cont  = 0;
    weights = options.weights;
    offset = options.offset;

    % ITERATIVELY REWEIGHTED LEAST SQUARE ALGORITHM (IRLS)

    for i=1:options.maxiter
    
      deta = dlink(mu);
      %deta = clink{2}(mu);
      %W = spdiags(weights./(deta.^2.*psi_value.*varianceff),0,n,n);
      W = spdiags(sqrt(weights)./(abs(deta).*stdfun(mu)),0,n,n);
      XTW = X.'*W;
      %XTWX = (X'*(W*X));
      %R = (X'*W*X)^(-1);
      R = pinv(XTW*XTW.');
      z = eta-offset+((Y-mu).*deta);
      BetaV(:,i) = R*XTW*(W*z);
      eta  = X*BetaV(:,i)+offset;
      mu = mufun(ilink(eta));
      
      %mu = clink{3}(eta);
  
      if i>1 && ~any((abs(BetaV(:,i)-BetaV(:,i-1)))>options.accuracy*max(sqrt(eps),BetaV(:,i-1)))
        cont=cont+1;
      end
      if cont>20
        break
      end
    end
    Y = Y.*members.sample_size;
    mu = mu.*members.sample_size;
    [dev,res,resP,resD,resA] = glmdev(options.family,Y,mu,members.sample_size,weights);
    members.deviance = dev;
    members.residual = res.';
    members.residualP = resP.';
    members.residualD = resD.';
    members.residualA = resA.';

    if options.constant
      wtdmu  = sum(weights .* Y)/sum(weights)/numel(Y);
    else
      wtdmu = ilink(offset);
    end
    members.deviance_null  = glmdev(options.family,Y,wtdmu,members.sample_size,weights);
    members.df_null =  n - options.constant;
    members.df = max(n-p,0);

    if options.estdispersn %dispersion_parameter=='mean_deviance'
      psi_value = members.deviance/members.df; % better estimate for psi_value?
    else
      psi_value = 1;
    end
    %   else
    %   if numel(weights)==1
    %     weights = weights(ones(n,1));
    %   end
    %   t=(mu>0);
    %   switch family
    %     case 'gamma'
    %       psi_value = sum(weights(t).*((Y(t)-mu(t)).^2)./((mu(t)).^2))/df;
    %     case 'normal'
    %       psi_value = sum(weights(t).*(Y(t)-mu(t)).^2)/df;
    %     case 'invnorm'
    %       psi_value=sum(weights(t).*(((Y(t)-mu(t)).^2)./(mu(t).^3)))/df;
    %     case 'poisson'
    %       psi_value= sum(weights(t).*((Y(t)-mu(t)).^2)./((mu(t))))/df;
    %     case 'binomial'
    %       psi_value= sum(weights(t).*((Y(t)-mu(t)).^2)./(mu(t)-(mu(t).^2)))/df;
    %   end
    %   end
    %end


    members.dispersnfit=psi_value;
    %BetaV(:,i:options.maxiter)=[];
    pcov = R*psi_value;
    pstd = sqrt(diag(pcov));
    members.params    = BetaV(:,i-1).';
    members.params_cov = pcov;
    members.params_std = pstd.';
    members.params_corr = pcov./(pstd*pstd.');
    members.params_tstat = members.params./members.params_std;

    if options.estdispersn %dispersion_parameter=='mean_deviance'
      members.params_pvalue=2.*cdft(-abs(members.params_tstat),members.df);
      pcrit = -pstd.'*invt(options.alpha/2,members.df);
    else
      members.params_pvalue=2.*cdfnorm(-abs(members.params_tstat));
      pcrit = -pstd.'*invnorm(options.alpha/2);
    end
    members.params_ci = [members.params+pcrit;members.params-pcrit];


    [n,p]=size(X);
    members.numvar = p;
    members.numobs = n;
    %Calculate R^2 (Ref Draper & Smith p.46)
    r=corrcoef([Y,mu(:)]);
    members.R2=r(1,2).^2; % % coefficient of determination
    members.R2adj =  max(1 - (1-members.R2)* (n-1)/(n-p-1),0); % adjusted coefficient of determination
    members.mu    = mu.';
    members.eta   = eta.';
    members.X     = X;
    members.W     = W;
    members.numiter = i-20;
    members.convergence = i<options.maxiter;
    if rcond(XTW*XTW.')<1e-14
      txt1 = 'Covariate matrix is illconditioned!';
    else
      txt1 = '';
    end
    members.options = options;
    members.family = options.family;
    members.link = options.link;
    members.note = txt1;
    members.date = datestr(now);


  end % compute_glm

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
  %Compare  small GLM versus large one
%
%   CALL     [pvalue] = compare(object2)
%
%	  The standard hypothesis test of a larger linear regression 
%	  model against a smaller one. The standard F-test is used.
%	  The output is the p-value, the residuals from the smaller 
%	  model, and the residuals from the larger model.
%
%	  See also fitls  
    error(nargchk(1,1,nargin))
    try
    fn = object2.fieldnames();
    if any(~ismember({'options','dispersnfit','deviance','X','df','numvar'},fn));
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
    if any(any((As-Al*(Al\As))>500*eps)) || ~strcmpi(members2.options.family,members.options.family) || ~strcmpi(members2.options.link,members.options.link)
      warning('WAFO:REGGLM','Small model not included in large model, result is rubbish!')
    end
    
    pmq = abs(nL-ns);
    disp(' ')
    disp('                       Analysis of Deviance')
    if options.estdispersn   
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
     if options.estdispersn && ( strcmpi(options.family,'binomial') || strcmpi(options.family,'poisson'))
        warning('WAFO:REGGLM','using F-test with a %s distribution is inappropriate',options.family)
%       elseif ~options.estdispersn
%         warning('WAFO:REGGLM','using F-test with a fixed dispersion is inappropriate')
      end
    
    
  end

  function anode()
    disp(' ')
    disp('                       Analysis of Deviance')
    if options.estdispersn
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
     if options.estdispersn && ( strcmpi(options.family,'binomial') || strcmpi(options.family,'poisson'))
        warning('WAFO:REGGLM','using F-test with a %s distribution is inappropriate',options.family)
%       elseif ~options.estdispersn
%         warning('WAFO:REGGLM','using F-test with a fixed dispersion is inappropriate')
      end
  end
  function summary()
    txtlink = members.options.link;
    if ~ischar(txtlink)
      txtlink = func2str(options.link{1});
    end
disp('Call:')
fprintf('glm(formula = %s(y) ~ x, family = %s)\n',txtlink,options.family)
disp(' ')
disp('Deviance Residuals:')
disp('    Min       1Q         Median       3Q        Max  ')
fprintf('%2.4f     ',percentile(members.residualD,[0 0.25 0.5 0.75 1]))
disp(' ')
disp(' Coefficients:')
if options.estdispersn
  disp('            Estimate      Std. Error     t value       Pr(>|t|)')
else
  disp('            Estimate      Std. Error     z value       Pr(>|z|)')
end
almat = [members.params(:) members.params_std(:),members.params_tstat(:),members.params_pvalue(:)];
if options.constant
  fprintf('(Intercept)  %2.4f        %2.4f        %2.4f        %2.4f\n',almat(1,:))
end
for ix1 = 1:members.numvar-options.constant
  fprintf('x_%d          %2.4f        %2.4f        %2.4f        %2.4f\n',ix1,almat(ix1+options.constant,:))
end
disp(' ')
fprintf('(Dispersion parameter for %s family taken to be %2.2f)\n',options.family,members.dispersnfit)
disp(' ')
if options.constant
  fprintf('    Null deviance: %2.4f  on %d  degrees of freedom\n',members.deviance_null,members.df_null)
end
fprintf('Residual deviance: %2.4f  on %d  degrees of freedom\n',members.deviance,members.df)

anode()

end % summary

function [y,ylo,yup] = predict(Xnew,varargin)
%REGGLM/PREDICT Predict from a fitted GLM object
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
%
% Example
%  y=[2 0 3 1 5 5 6 9 5 9].';  % number of successes in n(i) trials
%  n=[10 10 10 10 10 10 10 10 10 10].';
%  x=(1:10)';  % Covariate
%  b = glm([y n],x,'binomial','link','loglog');
%  plot(x, y./n,'o',x,b.mu'./n,'-','LineWidth',2), shg
%  [y1,ylo,yup] = b.predict(x); 
%  hold on, plot(x,y1,x,[ylo,yup],'r')
% 
% See also regglm

opts = struct('size',1,'alpha',0.05,'type','eta','cimean',true);
if nargin==2 && strcmpi(Xnew,'defaults')
  y = opts;
  return
end
error(nargchk(0,inf,nargin))

opts = parseoptions(opts,varargin{:});

if strcmpi(options.family,'binomial')
  ssize = opts.size(:);
else
  ssize = 1;
end


if nargin<1 || isempty(Xnew)
  Xnew = members.X;
else
  n = size(Xnew,1);
  if members.options.constant==1
    Xnew=[ones(n,1),Xnew];
  end
  notnans = ~(any(~isfinite(Xnew), 2) );
  Xnew = Xnew(notnans,:);
end
[n,p] = size(Xnew);
clink = glmlink(options.link,1);
ilink = clink{3};

  

if p ~= members.numvar
  error('Number of covariates must match the number of regression coefficients')
end
eta = (Xnew*members.params(:)+members.options.offset);

y = ssize.*ilink(eta);
if nargout>1

  pcov = members.params_cov;

  [U S V]=svd(pcov,0);
  R=(U*sqrt(S)*V'); %squareroot of pcov
  %[R,P] = genchol(pcov);
  %R = chol(pcov)
  if opts.cimean
    varxb = sum((Xnew*R).^2,2);
  else
    varxb = sum((Xnew*R).^2,2)+members.dispersnfit;
  end
  if members.options.estdispersn
    crit = -invt(opts.alpha/2, members.df);
  else
    crit = -invnorm(opts.alpha/2);
  end
  
  ecrit = crit * sqrt(varxb(:));
  yloup = [ssize.*ilink(eta-ecrit) ssize.*ilink(eta+ecrit)];
  ylo = min(yloup,[],2);
  yup = max(yloup,[],2);
end
end % predict
end % glm main



function [mufun,stdfun] = moment_funs(family,sample_size)
%MOMENT_FUNS Return function handles to mean and standard deviation

sml = realmin^(0.2);
mufun = @(mu) max(mu, sml); % mufun make sure mu is in valid range
switch family
  case 'gamma'
    stdfun = @(mu) mufun(mu);
    %stdfun = @(mu) mufun(mu)+(mu<=0);
    
  case 'invnorm'
    stdfun = @(mu) mufun(mu).^(3/2);
  case 'normal'
    mufun = @(mu) mu;
    stdfun = @(mu) ones(length(mu),1);
  case 'poisson'
    stdfun = @(mu) sqrt(mufun(mu)) + eps * (mu<=0);
    %varianceff=mu+0.00001*(mu==0);
   
  case 'binomial'
    mufun = @(mu) min(max(mu,eps),(1-eps));
    stdfun1 = @(mu) sqrt(mu).*sqrt(1-mu)./sqrt(sample_size);
    stdfun = @(mu) stdfun1(mufun(mu));
    %mu=mu+0.0001*(mu<0.0001);
    %mu=mu-0.0001*(mu>members.sample_size-0.0001);
    %varianceff=mu-((mu.^2)./members.sample_size);
end

end
function  predict()
%REGGLM>PREDICT Predict from a fitted GLM object
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
%
% Example
%  y=[2 0 3 1 5 5 6 9 5 9].';  % number of successes in n(i) trials
%  n=[10 10 10 10 10 10 10 10 10 10].';
%  x=(1:10)';  % Covariate
%  b = regglm([y n],x,'family','binomial','link','loglog');
%  plot(x, y./n,'o',x,b.get('mu')./n,'-','LineWidth',2), shg
%  [y1,ylo,yup] = b.predict(x); 
%  hold on, plot(x,y1,x,[ylo,yup],'r')
% 
% See also regglm

% Do not remove 
% This is a trick in order to get help regglm>predict to work
end

function model()
%REGGLM>MODEL 
%
% Member variables
%      .options      : input options as explained below.
%      .df           : degrees of freedom for error.
%      .params       : estimated model parameters by IRLS algorithm.
%      .params_ci    : 100(1-alpha)% confidence interval for model parameters
%      .params_tstat : t statistics for model's estimated parameters.
%      .params_pvalue: p value for model's estimated parameters.
%      .params_std   : standard errors for estimated parameters
%      .params_corr  : correlation matrix for estimated parameters.
%      .mu           : fitted values for the model.
%      .eta          : linear predictor for the model.
%      .residual     : residual for the model (Y-members.mu).
%      .dispersnfit  : The estimated error variance
%      .deviance     : deviance for the model.
%
% Methods
%      .predict   : Predict from a fitted GLM object
%      .summary   : Display summary of fitted GLM object.
%      .compare   : Compare small GLM versus large one
%      .get       : Return member variables
%
% See also regglm

% Do not remove 
% This is a trick in order to get help regglm>model to work
end

function [clnk]=glmlink(fun,m)
%GLM>GLMLINK Function handle to some common link functions
%
% CALL clnk=glmlink(fun,m)
%      clnk=glmlink(clnk,m)
% 
%  clnk = {lnk,dlnk,ilnk}  cell array of function handles defining the link
%         (lnk),  the derivative of the link (dlnk) and the inverse link (ilnk).
%  fun  = scalar or string defining the link. Options are
%         'identity','logit','log','loglog','cloglog','inverse',
%         'probit','sqrt'. The scalar value defines the exponent of the
%         power-function.
%  m    = sample size (when fun is 'logit','loglog','cloglog' or 'probit')
% 
% Example
%  clnk = glmlink(2);
%  x = linspace(0,3)
%  plot(x,clnk{1}(x))
% 
% See also glm

% Copyright (C) 2003, 2007  Stefano Favaro, Per A. Brodtkorb
% 
%     GLM is free software; you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation; either version 3 of the License, or
%     (at your option) any later version.
% 
%     GLM is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
    


% References    
% [1] Dobson, A.J. (1990), "An Introduction to Generalized Linear Models", 
%     CRC Press.
% 
% [2] Gill, J. (2001), "Generalized Linear Models: a unified approach", 
%     (quantitative applications in the social sciences, 134). Sage 
%     Pubblications.
% 
% [3] McCullagh, P., and J.A. Nelder (1990), "Generalized Linear Models", 
%     CRC Press.
%      
% [4] Myers, R.M., Montgomery, D.C. and Vining, G.G. (2001), "Generalized 
%     Linear Models with applications in engineering and the science", 
%     John Wiley and Sons.
%
% http://en.wikipedia.org/wiki/Generalized_linear_model#Count_data


% History
% By pab dec 2007
% Based on the functions in GLMBOX written in SciLab by Stefano Favaro (sfavaro@katamail.com)


if iscell(fun)
  clink = fun;
  if numel(clink)~=3 || any(~cellfun(@(x) isa(x,'function_handle'),clink))
    error('Link cellarray must contain 3 function handles)')
  end
  return
end
pw = [];
if isnumeric(fun)
  pw = fun;
  fun = 'power';
end
if strcmpi(fun,'power')
  if isempty(pw)
    error('Unknown exponent for power function. Supply exponent, please.')
  end
  switch pw
    case  0, fun = 'log';
    case  1, fun = 'identity';
    otherwise 
  end
elseif strcmpi(fun,'inverse')
  pw = -1;
  fun = 'power';
elseif strcmpi(fun,'sqrt')
  pw = 0.5;
  fun = 'power';  
end

switch fun
  case 'identity'
    lnk  = @(mu) mu;
    ilnk = @(eta) eta;
    dlnk = @(mu) ones(size(mu));
  case 'cloglog'
    
    lnk  = @(x)cloglog(x,m);
    ilnk = @(x)cloglogi(x,m);
    dlnk = @(x)cloglogd(x,m);
  case 'loglog'
%%
    lnk  = @(x)loglog(x,m);
    ilnk = @(x)loglogi(x,m);
    dlnk = @(x)loglogd(x,m);
  case 'log'
   sml = (realmin).^(1/4);
   lnk = @(mu) log(max(mu,sml));
   lb  = log(sml);  ub = -lb;
   ilnk = @(eta) exp(min(max(eta,lb),ub));
   dlnk = @(mu) 1./max(mu,sml);
%     lnk  = @log1;
%     ilnk = @logi;
%     dlnk = @logd;
  case 'logit'
    lnk  = @(x)logit(x,m);
    ilnk = @(x)logiti(x,m);
    dlnk = @(x)logitd(x,m);
  case 'power'
    lnk  = @(x)power1(x,pw);
    ilnk = @(x)poweri(x,pw);
    dlnk = @(x)powerd(x,pw);
   case 'probit'
    lnk  = @(x)probit(x,m);
    ilnk = @(x)probiti(x,m);
    dlnk = @(x)probitd(x,m);
    case 'inverse'
    lnk  = @inverse1;
    ilnk = @inversei;
    dlnk = @inversed;
  otherwise
    error('Unknown link function.')
end
clnk = {lnk,dlnk,ilnk};
end
function p = chkprb(p,m,tol)
if nargin<3
  tol = eps;
end

p = min(max(tol,p),m*(1-tol));
%p=p-(p>m-tol)*tol;
%p=p+(p<tol)*tol;
end


%% probit
function eta=probit(mu,m)
   mu = chkprb(mu,m);
   eta=-sqrt(2)*erfcinv(2*(mu./m));
end
function mu = probiti(eta,m)
   lb = -sqrt(2)*erfcinv(2*eps);
   ub = -lb;
   mu=m.*erfc(-min(max(eta,lb),ub)/sqrt(2))/2; 
end
function deta = probitd(mu,m)
   mu = chkprb(mu,m);
   deta=(sqrt(2*pi)./m).*exp((erfcinv((2*mu./m))).^2);
end

%% loglog
function eta = loglog(mu,m)
   mu = chkprb(mu,m);
   eta=log(-log(mu./m));
end
function mu = loglogi(eta,m)
  persistent lb ub
  if isempty(lb) || isempty(ub)
    lb = log(-log1p(-eps));
    ub = log(-log(eps));
  end
   mu=m.*(exp( -exp(min(max(eta,lb),ub))));
end
function deta = loglogd(mu,m)
   mu = chkprb(mu,m);
   deta=1./(mu.*log((mu./m)));
end


%% cloglog
function eta = cloglog(mu,m)
   mu = chkprb(mu,m);
   eta=log(-log1p(-mu./m));
end
function mu = cloglogi(eta,m)
  persistent lb ub
  if isempty(lb) || isempty(ub)
    lb = log(-log1p(-eps));
    ub = log(-log(eps));
  end
  mu=m.*(-expm1( -exp(min(max(eta,lb),ub))));
end
function deta = cloglogd(mu,m)
   mu = chkprb(mu,m);
   deta=1./((mu-m).*log1p(-(mu./m)));
end


% %% log function
% function eta = log1(mu)
% %LOG    Log function.
%    sml = (realmin).^(1/4);
%    eta=log(max(mu,sml));
% end
% function mu = logi(eta)
% %LOGI    inverse Log function.
%    lb = log(realmin)/4;
%    ub = -lb;
%    mu=exp(min(max(eta,lb),ub));
% end
% function deta = logd(mu)
% %LOG    derivative of Log function.
%    sml = (realmin).^(1/4);
%    deta = 1./max(mu,sml);
% end

%% Logit function
 
function z = logit(p,m)
%LOGIT    Logit function.
p = chkprb(p,m);
z = log(p./(m-p));
end
function p = logiti(z,m)
%LOGITI Compute the inverse of logit function.
  lb = log(realmin)/4;
  ub = -lb;
  p = m ./ (1+exp(-min(max(z,lb),ub)));
end

function dz = logitd(p,m)
%LOGITD Compute the derivative of logit function.    
  p = chkprb(p,m);
  dz=m./(p .*(m-p));
end

%% Power
function eta = power1(mu,pw)
%POWER Computes power function
  if pw~=1
   sml = realmin^(1/max(4*abs(pw),4));
   lb = sml;
   ub = 1/sml;
   mu = min(max(mu,lb),ub);
  end
  eta=mu.^pw;
end
function mu = poweri(eta,pw)
%POWERI Computes inverse power function
  if pw~=1
    sml = realmin^(1/max(4*abs(pw),4));
    lb = sml;
    ub = 1/sml;
    eta = min(max(eta,lb),ub);
  end
  
   mu = eta.^(1/pw);
end
function deta = powerd(mu,pw)
%POWERD Computes derivative of power function
  if pw~=1
   sml = realmin^(1/max(4*abs(pw),4));
   lb = sml;
   ub = 1/sml;
   mu = min(max(mu,lb),ub);
  end
 deta=pw.*(mu.^(pw-1));
end

%% reciprocal
function eta = inverse1(mu)
   mu=mu+(mu==0);
   eta=1./mu;
end
function mu = inversei(eta)
   eta=eta+(eta==0);
   mu=1./eta;
   mu=mu.*(mu>0)+(mu<=0);
end 
function deta = inversed(mu)
   mu=mu+(mu==0);
   deta=-1./(mu.^2);
end


function [dev,res,resP,resD,resA] = glmdev(distr,Y,mu,ssize,weights)
%REGGLM>GLMDEV Computes deviance for model as well as the residuals
%
%  CALL [dev,res,resP,resD,resA] = glmdev(distr,Y,mu,ssize,weights)
%
%   dev  = deviance for model
%   res  = residual Y-mu
%   resP = Pearson residual.
%   resD = Deviance residual. 
%   resA = Anscombe residual.
%
%  distr = distributions
%  Y       = vector of responses
%  mu      = fitted model
%  ssize   = sample size
%  weights = observation weights
%

res = Y-mu;
switch distr
  case 'gamma'
    mu=mu+0.00000001*(mu==0);
    yy=Y+(Y==0).*mu;
    dev=2*sum(weights.*(-log(yy./mu)+(Y-mu)./mu));
    if nargout>2
      resP=(Y-mu)./(mu);
      resD=(sign(Y-mu)).*sqrt(2*((-log(yy./mu)+(Y-mu)./mu)));
      resA=3*(Y.^(1/3)-mu.^(1/3))./mu.^(1/3);
    end
  case 'normal'
    if nargout>2
      resP=res;
      resD=res; %(sign(res)).*sqrt((res).^2);
      resA=res;
    end
    dev=sum(weights.*(res).^2);
  case 'invnorm'
    if nargout>2
      resP=(res)./(((mu.^3)).^(0.5));
      resD=(sign(res)).*sqrt(((res).^2)./((mu.^2).*Y));
      resA=(log(Y)-log(mu))./mu;
    end
    dev=(sum(((weights.*(res).^2)./((mu.^2).*Y))));
  case 'poisson'
    yy=Y+0.5*(Y==0);
    mumu=mu+0.5.*(mu==0);
    if nargout>2
      resP=(res)./(mu.^0.5);
      resD=(sign(res)).*sqrt(2.*(Y.*log(yy./mumu)-(Y-mumu)));
      resA=((3/2)*(((Y).^(2/3))-((mu).^(2/3))))./(mu.^(1/6));
    end
    dev=2*sum(weights.*(Y.*log(yy./mumu)-(Y-mumu)));
  case 'binomial'
    mu=(mu)+((mu)==0);
    yfix=(Y)+((Y)==0).*mu;
    mmmu=(ssize-mu)+(ssize==mu);
    yy=ssize-(Y);yy=yy+(yy==0).*(mmmu);
    if nargout>2
      resP=(Y-mu)./(((mu-(mu.*(mu./ssize)))).^0.5);
      resD=(sign(Y-mu)).*sqrt(2.*((yfix.*log(yfix./mu))+(yy.*log(yy./mmmu))));
      
      beta_value =  beta(2/3,2/3);
      betainc_valueY = betainc(Y./ssize,2/3,2/3);
      betainc_valuemu = betainc(mu./ssize,2/3,2/3);
      resA=(beta_value).*(betainc_valueY-betainc_valuemu)./(max(2.2204e-016,((mu./ssize).*(1-(mu./ssize)))).^(1/6)./sqrt(ssize));
    end
    dev=2*sum(weights.*((yfix.*log(yfix./mu))+(yy.*log(yy./mmmu))));
end
end % glmdev

function members1 = mkmemberstruct()

fn = {'family','link','options','numvar','numobs','df','df_null','params',...
  'params_ci','params_cov','params_std','params_corr','params_tstat',...
  'params_pvalue','mu','eta','X','Y','W','sample_size','residual',...
  'residualP','residualD','residualA','deviance','deviance_null',...
  'dispersnfit','dispersn','R2','R2adj','numiter','convergence','note','date'};
fn(2,:) = {[]};
members1 = struct(fn{:});
end



% Old call kept just in case        
% switch family
%   case 'gamma'
%     mu=members.mu+0.00000001*(members.mu==0);yy=Y+(Y==0).*mu;
%     members.residualP=(Y-members.mu)./(members.mu);
%     members.residualD=(sign(Y-members.mu)).*sqrt(2*((-log(yy./mu)+(Y-mu)./mu)));
%     members.residualA=3*(Y.^(1/3)-members.mu.^(1/3))./members.mu.^(1/3); 
%     members.deviance=2*sum(weights.*(-log(yy./mu)+(Y-mu)./mu));
%   case 'normal'
%     members.residualP=Y-members.mu;
%     members.residualD=(sign(Y-members.mu)).*sqrt((Y-members.mu).^2);
%     members.residualA=Y-members.mu; 
%     members.deviance=sum(weights.*(Y-members.mu).^2);
%   
%   case 'invnorm'
%     members.residualP=(Y-members.mu)./(((members.mu.^3)).^(0.5));
%     members.residualD=(sign(Y-members.mu)).*sqrt(((Y-members.mu).^2)./((members.mu.^2).*Y));
%     members.residualA=(log(Y)-log(members.mu))./members.mu;
% 
%     members.deviance=(sum(((weights.*(Y-members.mu).^2)./((members.mu.^2).*Y))));
%  
%  
%   case 'poisson'
%     yy=Y+0.5*(Y==0);
%     mumu=members.mu+0.5.*(members.mu==0);
%     members.residualP=(Y-members.mu)./(members.mu.^0.5);
%     members.residualD=(sign(Y-members.mu)).*sqrt(2.*(Y.*log(yy./mumu)-(Y-mumu)));
%     members.residualA=((3/2)*(((Y).^(2/3))-((members.mu).^(2/3))))./(members.mu.^(1/6));
%     
%     members.deviance=2*sum(weights.*(Y.*log(yy./mumu)-(Y-mumu)));
% 
%   case 'binomial'
%   mu=(members.mu)+((members.mu)==0);
%   yfix=(Y)+((Y)==0).*mu;
%   mmmu=(members.sample_size-mu)+(members.sample_size==mu);
%   yy=members.sample_size-(Y);yy=yy+(yy==0).*(mmmu);
% 
%   members.residualP=(Y-members.mu)./(((members.mu-(members.mu.*(members.mu./members.sample_size)))).^0.5);
%   members.residualD=(sign(Y-members.mu)).*sqrt(2.*((yfix.*log(yfix./mu))+(yy.*log(yy./mmmu))));
% 
%   beta_value =  beta(2/3,2/3);
%   betainc_valueY = betainc(Y./members.sample_size,2/3,2/3);
%   betainc_valuemu = betainc(members.mu./members.sample_size,2/3,2/3);
%   members.residualA=(beta_value).*(betainc_valueY-betainc_valuemu)./(max(2.2204e-016,((members.mu./members.sample_size).*(1-(members.mu./members.sample_size)))).^(1/6)./sqrt(members.sample_size));
% 
%   members.deviance=2*sum(weights.*((yfix.*log(yfix./mu))+(yy.*log(yy./mmmu))));
%  
% end        

% 
%      .family : distribution for the model.
%      .link : link function for the model.
%      .numvar : number of covariates (number of columns for X matrix). 
%      .numobs : number of observations (number of rows for Y matrix).
%      .params : estimated parameters by IRLS algorithm.
%      .df     : degrees of freedom for error.
%      .tstat  : t statistics for model's estimated parameters.
%      .pvalue : p value for model's estimated parameters.
%      .stderrors : standard errors for estimated parameters
%      .corr   : correlation matrix for estimated parameters.
%      .mu     : fitted values for the model.
%      .eta    : linear predictor for the model.
%      .residual  : residual for the model (Y-members.mu).
%      .residualP : Pearson residual.
%      .residualD : Deviance residual. 
%      .residualA : Anscombe residual.
%      .deviance  : deviance for the model.
%      .dispersnfit: 
%      .dispersn   : fixed or estimated dispersion parameter for the model.
%      .numiter   : number of iterations to reach convergence.
