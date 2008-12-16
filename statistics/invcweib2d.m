function [x, ind]= invcweib2d(x1,p,varargin)
%INVCWEIB2D Inverse of the conditional 2D weibull cdf of X2 given X1.
%
% CALL: X2 =  invcweib2d(X1 ,P,param,rtol)  
% 
%   X2    = inverse of the conditional cdfweib2d given X1 at the
%           probabilities in P.
%   X1    = conditional variable 
%   A1, A2 = scale parameters    
%   B1, B2 = shape parameters
%      C12 = interaction parameter between X1 and X2
%     phat = Distribution parameter struct
%            as returned from FITWEIB2D.  
%  options = struct with fieldnames:
%     .logp    : if TRUE, probability, p, returned as log(p).
%     .condon  : 0 it returns the the regular cdf of X1 and X2 (default)
%                1 it returns the conditional cdf given X1
%                2 it returns the conditional cdf given X2
%     .releps  : Requested relative error (Default 1e-3) 
%     .iter_max: Maximum number of iterations (default 500)
%
%   The size of X is the common size of the input arguments X1 and P. A scalar input  
%   functions as a constant matrix of the same size as the other inputs.    
%
%   INVCWEIB2D uses Newton's method to converge to the solution.
%
% Example
%  x = linspace(eps,6,200); x2 = 3;
%  phat = {2 2  3 2.5 .8};
%  F = cdfweib2d(x,x2,phat{:},'condon',2);
%  x1 = invcweib2d(x2,F,phat{[3:4 1:2 5]}); 
%  semilogy(abs(x-x1)./x+eps)
%
% See also  pdfweib2d, cdfweib2d, rndweib2d, fitweib2d, momweib2d

% Tested on: matlab 5.2
%History
% revised pab 13.11.2000, minor fixes
% By Per A. Brodtkorb 19.11.98

error(nargchk(3,15,nargin))
Np = 5;
options = struct('covariance',[],'alpha',0.05,...
  'lowertail',true,'logp',false,'abseps',1e-90,'releps',1e-5,'max_iter',500); % default options
[params,options] = parsestatsinput(Np,options,varargin{:});
if numel(options)>1
  error('Multidimensional struct of distribution parameter not allowed!')
end

[a1,b1,a2 b2,c12] = deal(params{:});

if isempty(a1)||isempty(b1)||isempty(a2)||isempty(b2)||isempty(c12)
  error('Requires either 7 input arguments or that input argument 3 is FDATA.'); 
end

[icode  x1 p a1 b1 a2 b2 c12 ] = iscomnsize(x1,p,a1,b1,a2, b2, c12);
if ~icode 
    error('Requires non-scalar arguments to match in size.');
end

%   Initialize X to zero.
x = zeros(size(p));

 ok =(0<=p & p<=1 &  a1 > 0 & b1 > 0 & a2 > 0 & b2 > 0 &  abs(c12)<1);
k = find(p<0 | p>1 | ~ok);
if any(k),
    tmp  = NaN;
    x(k) = tmp(ones(size(k))); 
end

% The inverse cdf of 0 is 0, and the inverse cdf of 1 is inf.  
k0 = find(p == 0 & ok);
if any(k0),
    x(k0) = zeros(size(k0)); 
end

k1 = find((p == 1| isinf(x1) )& ok );
if any(k1), 
  tmp = Inf;
    x(k1) = tmp(ones(size(k1))); 
end

% Newton's Method
% Permit no more than maxcount interations.
maxcount = options.max_iter;
count = 0;

k = find(~isinf(x1) & p > 0  &  p < 1 & ok);

if ~any(k),return,end

pk = p(k);

% Supply a starting guess for the iteration.
cvar=linspace(0,max(x1(k)),20); % make sure that the comp. of mn and v do 
                                % not consume valuable time
[mn v ]=momweib2d(params{:},'condon',1,'cvar',cvar); %slow
mn = interp1(cvar,mn,x1(k),'linear'); %fast
v = interp1(cvar,v,x1(k),'linear'); %fast

switch 2, %method
  case 1,%Use a lognormal distribution. 
  temp = log(v + mn .^ 2); 
  mu = 2 * log(mn) - 0.5 * temp;
  sa = -2 * log(mn) + temp;
  xk = exp(invnorm(pk,mu,sa.^2));
case 2, %   Use a normal distribution. probably the fastest choice
  xk=abs(invnorm(pk,mn,v));
  %xk((xk<=0))=1e-3;
    if any(isnan(xk))
      disp(['Warning: NaNs   ' num2str(sum(isnan(xk)))])
    end
case 3, % use weibull distribution slowest choice
	%   beta=fsolve('gamma(1+1./x).^2./gamma(1+2./x)-P1.^2./(P2+P1.^2)',parm(4).*ones(size(mn)),[],[],mn,v);
	%   alpha=mn./gamma(1+1./beta);
	%   xk=invweib(pk,alpha,beta); 
end
%x(k)=xk;
%return
h = ones(size(pk)); 

% Break out of the iteration loop for three reasons:
%  1) the last update is very small (compared to x)
%  2) the last update is very small (compared to sqrt(eps))
%  3) There are more than 100 iterations. 

opts.condon = 1;
opts.releps = options.releps;

k2=find((abs(h) > sqrt(eps)*abs(xk))  &  abs(h) > sqrt(eps));
while(any(k2) && count < maxcount), 
                                 
    count = count + 1;
    h(k2)  = (cdfweib2d(x1(k(k2)),xk(k2),params{:},opts) - pk(k2)) ./ pdfweib2d(x1(k(k2)),xk(k2),params{:},opts);
   
    xnew = xk(k2) - h(k2);
    % Make sure that the current xnew>0.
    ix = find(xnew <= 0);
    if any(ix),
        xnew(ix) = xk(k2(ix)) / 10;
        h(k2(ix)) = xk(k2(ix))-xnew(ix);
    end
    xk(k2) = xnew;
   % disp(sprintf('Iteration %d.  Number of points left:  %d', count ,numel(k2)))
    %if any(isnan(xk)),  disp(['Warning: values out of range   '...
    %  num2str(sum(isnan(xk)))]),   end
    
    % if not converged decrease the relative tolerance in the calculation of cdf
    if length(k2)<length(k)/4 && count>maxcount/4, opts.releps=opts.releps/10;end
    k2=find((abs(h) > options.releps*abs(xk))  &  abs(h) > options.abseps);   
end


% Store the converged value in the correct place
x(k) = xk;

if count == maxcount, 
    disp('Warning: INVCWEIB2D did not converge.');
    str = ['The last steps were all less than:  ' num2str(max(abs(h(k2))))];
    disp(str)
    ind=k(k2);
  else
    ind=[];
end
