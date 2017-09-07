function [y ,abserr] = hypgf(a,b,c,x,varargin)
%HYPGF  Hypergeometric function F(a,b,c,x) 
%
% CALL:   [y ,abserr] = hypgf(a,b,c,x,abseps,releps) 
%
%       y = F(a,b,c,x)
%  abserr = absolute error estimate
% a,b,c,x = input parameters
% abseps  = requested absolute error
% releps  = requested relative error  
%
% HYPGF calculates one solution to Gauss's hypergeometric differential
% equation:
%
%    x*(1-x)Y''(x)+[c-(a+b+1)*x]*Y'(x)-a*b*Y(x) = 0
% where
%   F(a,b,c,x) = Y1(x) = 1 + a*b*x/c + a*(a+1)*b*(b+1)*x^2/(c*(c+1))+.... 
%
%   
% Many elementary functions are special cases of F(a,b,c,x):
% 1/(1-x) = F(1,1,1,x) = F(1,b,b,x) = F(a,1,a,x)
% (1+x)^n = F(-n,b,b,-x)
% atan(x) = x*F(.5,1,1.5,-x^2)
% asin(x) = x*F(.5,.5,1.5,x^2)
% log(x)  = x*F(1,1,2,-x)
% log(1+x)-log(1-x) = 2*x*F(.5,1,1.5,x^2)
%  
%  NOTE: only real x, abs(x) < 1 and c~=0,-1,-2,... are allowed.
%
% Examples:
% x = linspace(-.99,.99)';
% [Sn1,err1] = hypgf(1,1,1,x);
% plot(x,abs(Sn1-1./(1-x))+eps,'b',x,err1+eps,'r'),set(gca,'yscale','log');
% [Sn2,err2] = hypgf(.5,.5,1.5,x.^2);
% plot(x,abs(x.*Sn2-asin(x))+eps,'b',...
%      x,abs(x.*err2)+eps,'r');
% set(gca,'yscale','log');
%
% close all;



% Reference:
% Kreyszig, Erwin (1988)
% Advanced engineering mathematics
% John Wiley & Sons, sixth edition, pp 204.

% Revised pab 30.08.2004
% -added dea3 for more stable error estimation
% revised pab 18.04.2001
% - updated help header
% - added example
% - added reference
% By Per A. Brodtkorb 17.11.98
[csize,  x, a, b, c] = comnsize(x,a,b ,c);
if any(isnan(csize))
  error('Requires non-scalar arguments to match in size.');
end
  
[y,abserr] = gethgf(a,b,c,x,varargin{:});
 
return

function  [fsum, err]=gethgf(a,b,c,x,absEps,relEps,Kmax)
%  error(nargchk(4,7,nargin))
narginchk(4,7)
  if (nargin<7||isempty(Kmax))
    Kmax  = 10000;  
  end
  Kmin = 2;
  if (nargin<6 || isempty(relEps)),
    relEps   = eps*10;
  end
  if (nargin<5||isempty(absEps))
    absEps   = 0;%eps*10;
  end
  useDEA = 1;
  %abseps = max(abseps,eps*10); 
   
  fsum  = zeros(size(x));   
  delta = fsum;
  err   = fsum;
  
  ok    = ~((round(c)==c & c<=0) | abs(x)>1);
  k1=find(~ok);
  if any(k1), 
    warning('WAFO:HYPGF','Illegal input: c = 0,-1,-2,... or abs(x)>1')
    fsum(k1) = NaN;
    err(k1)  = NaN;
  end
  k0=find(ok & abs(x)==1);
  if any(k0)
    cmab = c(k0)-a(k0)-b(k0);
    fsum(k0) = exp(gammaln(c(k0))+gammaln(cmab)-...
		   gammaln(c(k0)-a(k0))-gammaln(c(k0)-b(k0)));
    err(k0) = eps;
    k00 = find(real(cmab)<=0);
    if any(k00)
      err(k0(k00)) = nan;
      fsum(k0(k00)) = nan;
    end
  end
  k=find(ok & abs(x)<1);
  if any(k),
    delta(k) = ones(size(k));
    fsum(k)  = delta(k);
    
    k1 = k;
    E = cell(1,3);
    E{3} = fsum(k);
    converge = 'n';
    for  ix=0:Kmax-1,
      delta(k1) = delta(k1).*((a(k1)+ix)./(ix+1)).*((b(k1)+ix)./(c(k1)+ ...
						  ix)).*x(k1);
      fsum(k1) = fsum(k1)+delta(k1);
      
      E(1:2) = E(2:3);
      E{3}   = fsum(k1);
      
      if ix>Kmin
	if useDEA,
	  [Sn, err(k1)] = dea3(E{:});
	  k00 = find((abs(err(k1))) <= max(absEps,abs(relEps.*fsum(k1))));
	  if any(k00)
	    fsum(k1(k00)) = Sn(k00);
	  end
	  if (ix==Kmax-1)
	    fsum(k1) = Sn;
	  end
	  k0 = (find((abs(err(k1))) > max(absEps,abs(relEps.*fsum(k1)))));
	   if any(k0),% compute more terms
	     %nk=length(k0);%# of values we have to compute again
	     E{2} = E{2}(k0);
	     E{3} = E{3}(k0);
	   else
	     converge='y';
	     break;
	   end
	else
	  err(k1) = 10*abs(delta(k1));
	   k0 = (find((abs(err(k1))) > max(absEps,abs(relEps.* ...
							fsum(k1)))));
	   if any(k0),% compute more terms
	     %nk=length(k0);%# of values we have to compute again
	   else
	     converge='y';
	     break;
	   end
	end
	k1 = k1(k0);
      end
    end
    if ~strncmpi(converge,'y',1)
      disp(sprintf('#%d values did not converge',length(k1)))
    end
  end
  %ix
  return
function [result,abserr]  =dea3(E0,E1,E2)
%DEA3 Extrapolate a slowly convergent sequence
%
%  CALL: [result,abserr]  =dea3(E0,E1,E2)
%
%  result   = extrapolated value
%  abserr   = absolute error estimate  
%  E0,E1,E2 = 3 values of a convergent sequence to extrapolate
%  
% DEA3  attempts to extrapolate nonlinearly to a better estimate of the
%      sequence's limiting value, thus improving the rate of
%       convergence. Routine is based on the epsilon algorithm
%            of P. Wynn. 
% 
% Example % integrate sin(x) from 0 to pi/2
%   for I = 5:7
%      NPARTS = 2.^I;
%      x = linspace(0,pi/2,2.^I+1);
%      Ei(I-4) = trapz(x,sin(x))  
%   end
%   [En, err] = dea3(Ei(1),Ei(2),Ei(3))  
%   truErr = Ei-1  
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% -*- Mode: Matlab -*- %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% dea3.m --- 
%% Author          : Per Andreas Brodtkorb
%% Created On      : Wed Aug 18 12:40:23 2004
%% Last Modified By: Per Andreas Brodtkorb
%% Last Modified On: Wed Aug 18 12:51:49 2004
%% Update Count    : 3
%% Status          : Unknown, Use with caution!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% REFERENCES  "Acceleration de la convergence en analyse numerique",
%                 C. Brezinski, "Lecture Notes in Math.", vol. 584,
%                 Springer-Verlag, New York, 1977.
  
  
%error(nargchk(3,3,nargin))
narginchk(3,3)
ten = 10.0d0;
one = 1.0d0;
small  = eps; %spacing(one)
delta2 = E2 - E1;
delta1 = E1 - E0;
err2   = abs(delta2);
err1   = abs(delta1);
tol2   = max(abs(E2),abs(E1)) * small;
tol1   = max(abs(E1),abs(E0)) * small;

result = zeros(size(E0));
abserr = result;

k0 = find( ( err1 <= tol1 ) | (err2 <= tol2));
if any(k0)
  %C           IF E0, E1 AND E2 ARE EQUAL TO WITHIN MACHINE
  %C           ACCURACY, CONVERGENCE IS ASSUMED.
  result(k0) = E2(k0);
  abserr(k0) = err1(k0) + err2(k0) + E2(k0)*small*ten;
end
k1 = find(~(( err1 <= tol1 ) | (err2 <= tol2)));
if any(k1)
  ss = one./delta2(k1) - one./delta1(k1);
  k2 = k1(((abs(ss.*E1(k1)) <= 1.0d-3)));
  if any(k2)
    result(k2) = E2(k2);
    abserr(k2) = err1(k2) + err2(k2) + E2(k2)*small*ten;
  end
  k3 = find(~(abs(ss.*E1(k1)) <= 1.0d-3));
  if any(k3)
    k33 = k1(k3);
    result(k33) = E1(k33) + one./ss(k3);
    abserr(k33) = err1(k33) + err2(k33) + abs(result(k33)-E2(k33));
  end
end
return
   %   end subroutine dea3
   
   
% function  [fsum, delta]=gethgf2(a,b,c,z)
%    Kmax=10000;   fac=1;   tmps=0;
%    tmp   = (gamma(a).*gamma(b))./(fac.*gamma(c));
%    delta=tmp;
%    tol=eps;
%    for  ix=1:Kmax
%      fac   = fac*ix;
%      tmp   = (gamma(a+ix).*gamma(b+ix))./(fac.*gamma(c+ix));
%      delta = tmp*(z.^ix); 
% 
%      %                              don't bother with more terms than we use
%      if (max(abs(delta)) <= tol) , break ,end
%      %                            update
%      tmps  = tmps + delta 
%    
%    end
%    if ix>=Kmax, disp('Some values did not converge'), end
%     fsum=(gamma(c)./(gamma(a).*gamma(b))).*tmps;
%      delta=(gamma(c)./(gamma(a).*gamma(b))).*delta;
%    return

   
    

  


