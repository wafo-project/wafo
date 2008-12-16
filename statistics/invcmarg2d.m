function [x, ind]= invcmarg2d(x1,p,phat)
%INVCMARG2D Inverse of the conditional cdf of X2 given X1.
% 
%   CALL: x2 =  invcmarg2d(x1,P,phat)  
%
%  x2 = the inverse of cdfmarg2d given x1 and P
%  x1 = quantiles
%  P  = 0..1
%  phat = parameter structure (see fitmarg2d)
%
%   The size of X2 is the common size of the input arguments X1 and P. 
%   A scalar input functions as a constant matrix of the same size as the
%   other inputs.    
%
%   INVCMARG2D uses Newton's method to converge to the solution.
%
% See also   cdfmarg2d,  pdfmarg2d, fitmarg2d, rndmarg2d

%   References:
%      Plackett, R. L. (1965) "A class of bivariate distributions."
%                                J. Am. Stat. Assoc. 60. 516-22
%      [1]  Michel K. Ochi,
%       OCEAN TECHNOLOGY series 6
%      "OCEAN WAVES, The stochastic approach", Cambridge
%      1998 pp. 133-134.


%  tested on: matlab 5.2
% history
% revised pab 8.11.1999
%  - updated header info
%  - changed phat from vectro to structure
% By Per A. Brodtkorb 01.02.99

error(nargchk(3,3,nargin))
[icode  x1 p ] = iscomnsize(x1,p);

if ~icode 
    error('Requires non-scalar arguments to match in size.');
end

%VDIST=lower(phat.dist{1});
%HDIST= lower(phat.dist{2});

%   Initialize X to zero.
x = zeros(size(p));
%size(x),size(x1)
ok=(p > 0  &  p < 1);
k = find(~ok);
if any(k),
    tmp  = NaN;
    x(k) = tmp(ones(size(k))); 
end

% The inverse cdf of 0 is 0, and the inverse cdf of 1 is inf.  
%k0 = find(p == 0);
%if any(k0),
%    x(k0) = zeros(size(k0)); 
%end

k1 = find((p == 1| isinf(x1) ));
if any(k1), 
  tmp = Inf;
    x(k1) = tmp(ones(size(k1))); 
end

% Newton's Method
% Permit no more than count_limit interations.
count_limit = 150;
count = 0;

k = find(~isinf(x1) & ok );

pk = p(k);

% Supply a starting guess for the iteration.
cvar=linspace(eps,max(x1(k)),15); % make sure that the comp. of mn and v do 
                                % not consume valuable time
[mn v ]=mommarg2d(phat,'condon',1,'cvar',cvar); %slow
mn = interp1(cvar,mn,x1(k),'linear'); %fast
v = interp1(cvar,v,x1(k),'linear'); %fast

switch 2, %method
  case 1,%Use a lognormal distribution. 
  temp = log(v + mn .^ 2); 
  mu = 2 * log(mn) - 0.5 * temp;
  sigma = -2 * log(mn) + temp;
  xk = (invlognorm(pk,mu,sigma.^2));
  %xk = exp(invnorm(pk,mu,sigma));
case 2, %   Use a normal distribution. probably the fastest choice
  xk=abs(invnorm(pk,mn,v));
  %xk((xk<=0))=1e-3;
    if any(isnan(xk))
      disp(['Warning: NaNs   ' num2str(sum(isnan(xk)))])
    end

end
%x(k)=xk;
%return
h = ones(size(pk)); 

% Break out of the iteration loop for three reasons:
%  1) the last update is very small (compared to x)
%  2) the last update is very small (compared to sqrt(eps))
%  3) There are more than 100 iterations. This should NEVER happen.
 
k2=find((abs(h) > sqrt(eps)*abs(xk))  &  abs(h) > sqrt(eps));
while(any(k2) && count < count_limit), 
                                 
    count = count + 1;
    h(k2)  = (cdfmarg2d(x1(k(k2)),xk(k2),phat,'condon',1) - pk(k2)) ./ pdfmarg2d(x1(k(k2)),xk(k2),phat,'condon',1);
   % if any(isnan(h(k2))),      h(isnan(h))=-0.01;    end
    xnew = xk(k2) - h(k2);
    % Make sure that xnew>0
    ix = find(xnew <= 0);
    if any(ix),
        xnew(ix) = xk(k2(ix)) / 10;
        h(k2(ix)) = xk(k2(ix))-xnew(ix);
    end
    xk(k2) = xnew;
    %disp(sprintf('Iteration  %d,  Number of points left:  %d', count,length(k2))),
    %if any(isnan(xk)),  disp(['Warning: values out of range   ' num2str(sum(isnan(xk)))]),   end
    k2=find((abs(h) > sqrt(eps)*abs(xk))  &  abs(h) > sqrt(eps));   
end


% Store the converged value in the correct place
x(k) = xk;

if count == count_limit, 
    warning('WAFO:WSTATS:INVCMARG2D','INVCMARG2D did not converge.');
    outstr = ['The last steps were all less than:  ' num2str(max(abs(h(k2))))];
    disp(outstr)
    ind=k(k2);
  else
    ind=[];
end





