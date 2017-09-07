function [x,xlo,xup] = invbeta(F,varargin)
%INVBETA  Inverse of the Beta distribution function
%
% CALL:  x = invbeta(F,a,b,options)
%        [x,xlo,xup] = invbeta(F,phat,options)
%
%      x   = inverse cdf for the Beta distribution evaluated at F.
%  xlo,xup = 100*(1-ALPHA)% confidence interval for X.
%        F = lower or upper tail probability 
%     a, b = distribution parameters
%     phat = Distribution parameter struct
%            as returned from FITBETA.  
%  options = struct with fieldnames:
%         .lowertail: if TRUE (default), F = Prob[X <= x],
%                     otherwise, F = Prob[X > x].
%         .logp     : if TRUE, probability, p, input as log(p).
%         .alpha    : Confidence coefficent        (default 0.05)
%         .abseps   : Requested absolute error     (default 1e-90)
%         .releps   : Requested relative error     (default sqrt(eps))
%         .max_iter : Maximum number of iterations (default 100)
%
% INVBETA computes the inverse Beta distribution function via a
% combination of Newton method and a bisection search.
% The Beta PDF is defined by:
% 
%    f = x^(a-1)*(1-x)^(b-1)/H(a,b)    0<= x <= 1, a>0, b>0
%
% NOTE: To compute accurate uppertail quantiles use the identity
%      1 - INVBETA(F,A,B) = INVBETA(1-F,B,A).
%
% Example:
%   a=1;b=2;    
%   opt = {'lowertail',false,'logp',false};
%   F0 = [logspace(-300,-1) linspace(0.11,0.5)];
%   x  = invbeta(F0,a,b,opt{:});
%   F  = cdfbeta(x,a,b,opt{:});
%   semilogy(abs(F-F0)./F0+eps); % relative error
%
%   a=1;b=2;
%   x0 = [logspace(log10(realmin),-1) linspace(0.2,0.5)];
%   F  = cdfbeta(x0,a,b);
%   x  = invbeta(F,a,b);
%   semilogy(abs(x-x0)./x0+eps); % relative error
%
%   close all;
%
% See also pdfbeta, cdfbeta, rndbeta, fitbeta, mombeta

% Copyright (C) 1995 Anders Holtsberg, 2000 Per A. Brodtkorb
%
% This program is free software; you can redistribute it and/or modify
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


% tested on matlab 5.3
%History:
% revised pab sept 2007
% -improved the Newton-Raphson method
% revised pab nov2005
% -commented out disp statement
%revised pab 29.10.2000
% adapted from stixbox
% -added nargchk, comnsize
% -refined the Newton-Raphson method
%       Anders Holtsberg, 27-07-95

options = struct('proflog',false,'alpha',0.05,...
  'lowertail',true,'logp',false,'abseps',1e-200,'releps',sqrt(eps),'max_iter',100); % default options
if (nargin==1 && nargout <= 1 && isequal(F,'defaults'))
  x = options; 
  return
end
%error(nargchk(2,inf,nargin))
narginchk(2,inf)
Np = 2;
[params,options,tmp,phat] = parsestatsinput(Np,options,varargin{:});
if numel(options)>1
  error('Multidimensional struct of distribution parameter not allowed!')
end

[a,b] = deal(params{:});

if options.logp
  if options.lowertail
    logF = F;
  else
    %F = -expm1(R);
    logF = log1p(-exp(F));
  end
else
  if options.lowertail
    logF = log(F);
  else
    logF = log1p(-F);
  end
end

[icode logF,a,b] = iscomnsize(logF,a,b);
if ~icode,
  error('F, a and b must be of common size or scalar');
end

x = zeros(size(logF));

ok = (a>0 & b>0);

k = find(-inf<logF & logF<0 & ok);
if any(k)
 
  bk = b(k); %min(b(k),100000);
  if any(bk>10^5), 
    warning('WAFO:INVBETA','b is too large'),
  end
  ak = a(k);
  logFk = logF(k);
  
  
  % new call
  %   swap1 = (logFk>log(0.5) & ak>bk);
  %   swap2 = (logFk<log(0.5) & ak<bk);
  %   swapi = swap1 | swap2;
%   swapi = logFk>log(0.5);
%   if any(swapi(:))
%     % using the identity 1-B(x,a,b) = B(1-x,b,a)
%     logFk(swapi) = log(-expm1(logFk(swapi)));
% %     logFk(swap1) = log(-expm1(logFk(swap1)));
% %     logFk(swap2) = log1p(-exp(logFk(swap2)));
%     tmpa = ak(swapi);
%     ak(swapi) = bk(swapi);
%     bk(swapi) = tmpa;
%   end
  xstart = max(min(ak ./ (ak+bk),1-sqrt(eps)),sqrt(eps));
  
  k1 = find(ak==1);
  if any(k1)
    xstart(k1) = min(max(-expm1(log(-expm1(logFk(k1)))./bk(k1)),realmin*1000),1-eps);
  end
  k1 = find(bk==1);
  if any(k1)
    xstart(k1) = min(max(exp(logFk(k1)./ak(k1)),realmin*1000),1-eps);
  end

  
  %xnew = newtonx(xstart,logFk,ak,bk,options);
  [xnew] = newtonlogx(xstart,logFk,ak,bk,options);
  %  xnew(swapi) = 1-xnew(swapi);
  
  x(k)=xnew;
end

k2=find(logF==0 & ok);
if any(k2)
  x(k2)=1;
end

k3=find(~ok);
if any(k3)
  x(k3)=NaN;
end


if nargout >= 2  % Compute confidence bounds.
  if isempty(phat)
    error('Must have distribution struct!')
  end
  pcov = phat.covariance;
  %pcov = options.covariance;
  alpha = options.alpha;
  if options.proflog
    error('Confidence interval usin proflog not implemented!')
%     xlo = x;
%     xup = x;
%     for ix =1:numel(x)
%       [Lp,CI] = proflog(phat,'i',1,'x',x(ix),'link',@lnkbeta,'alpha',alpha);
%       xlo(ix) = CI(1);
%       xup(ix) = CI(2);
%     end
  else

   h = 1e-5;
  % Approximate the variance of X on the logit scale
  if options.lowertail
    logitp = log(x./(1-x));
    dBda = h./(betainc(x,a+h,b)-F);
    dBdb = h./(betainc(x,a,b+h)-F);
  else
    logitp = -log(x./(1-x));
    dBda = h./(betainc(x,b,a+h)-F);
    dBdb = h./(betainc(x,b+h,a)-F);
  end
 
  
  dL   = 1 ./ (x.*(1-x)+realmin); % derivative of logit(X) w.r.t. X
  dLda = dBda.* dL;      % dlogitp/da = dF/da * dlogitp/dF
  dLdb = dBdb.* dL;    % dlogitp/db = dF/db * dlogitp/dF
  varLogitp = pcov(1,1).*dLda.^2 + 2.*pcov(1,2).*dLda.*dLdb + pcov(2,2).*dLdb.^2;
    
  if any(varLogitp(:) < 0)
    error('Covariance must be a positive semi-definite matrix.');
  end
    
  % Use a normal approximation on the logit scale, then transform back to
  % the original CDF scale
  zcrit = -invnorm(alpha/2) * sqrt(varLogitp);
  zlo = logitp - zcrit;
  zup = logitp + zcrit;
  
  xlo = 1 ./ (1 + exp(-zlo));
  xup = 1 ./ (1 + exp(-zup));
  end
end



function [xnew, didnt] = newtonlogx(xstart,logPrbk,ak,bk,options)
 %NEWTONLOGX Newton-Raphson iteration using logarithmic scale on x and Prb

  
  %tail = 'lower';
  xnew = xstart;
  uk   = log(xnew);
  du     = ones(size(ak));
  dx     = du;
  prb = du;
  umin = log(realmin);
  iterNo = 0;
  max_count = options.max_iter; %500;
 
  ix        = find(uk(:));
  while (any(ix) && iterNo<max_count)
    
    xi     = xnew(ix);
    prb(ix) = cdfbeta(xi,ak(ix),bk(ix));
    logprb  = max(log(prb(ix)),-1000);
    
    du(ix) = (logprb - logPrbk(ix)).*prb(ix)./max(xi.*pdfbeta(xi,ak(ix),bk(ix)),realmin);
    
    %du(ix)   = max(-3, min(3, du(ix)));
    %uk(ix)   = uk(ix)-du(ix);
    ui   = uk(ix)-du(ix);
    %     % Make sure that the current guess is larger than UMIN and less than 0
    %     % by going 6/10 and 9/10 of the way towards UMIN and 0, respectively
     uk(ix) = ui + 0.1*(du(ix) - 9*ui).*(ui>0) + 0.1*(4*du(ix)-6*(ui - umin)).*(ui<umin) ;
     
    
    xnew(ix) = exp(uk(ix));
    dx(ix)   = xnew(ix)-xi;
    if iterNo<5
      ix = find(abs(log(prb)-logPrbk)> options.releps.*abs(logPrbk) | (abs(dx) > options.releps*abs(xnew)));
    else
      ix = find((abs(dx) > options.releps*abs(xnew))  &  abs(dx) > options.abseps );
    end 
    iterNo = iterNo+1;  
    %disp(['Iteration ',num2str(iterNo),'  Number of points left:  ' num2str(length(ix)) ]),
  end
  
  
  problems = abs(log(prb)-logPrbk)>sqrt(options.releps).*abs(logPrbk);
  %problems = (abs(dF)> sqrt(options.releps).*F );% & abs(dF)> options.abseps);
  didnt = [];
  if (iterNo == max_count || any(problems(:))),    
    F = exp(logPrbk);
    dF = prb-F;
    if any(problems(:)),    
      tozero = problems & F<=dF;
      xnew(tozero) = 0;
      toone = problems & (1-F<=-dF);
      xnew(toone) = 1;
      if max(min(F(problems),1-F(problems)))<eps
        reason = 'The tail probabilities are too small or too large for the given parameters A and B!';
      else
        reason = 'The parameters A or B are too small or too large for the given tail probability!';
      end
    else
      reason = '';
    end
    
    badcdf = find(problems(:) |  abs(dx(:))>sqrt(options.releps));    
    didnt = unique([ix(:);badcdf]);
    numDidnt = length(didnt);
    didnt = didnt(1);
    badA = ak(didnt);
    badB = bk(didnt);
    badF = F(didnt);
    baddF = dF(didnt);
    %badX = x(didnt);
    outstr = sprintf('a = %g, b=%g, F = %g dF = %g\n',[badA(:),badB(:),badF(:),baddF(:)].');
    warning('WAFO:INVBETA','invbeta did not converge for %d values:\n%s\n%s',numDidnt,reason,outstr);    
  end
  
  
%%  function X = binv(A,B,VAPP,P)
% C***********************************************************************
% C
% C PURPOSE: This function uses a series expansion method to compute
% C          percentage points for the Beta distribution.
% C
% C ARGUMENTS:
% C     A,B  -  The parameters of the Beta distribution.
% C          A and B must both be positive real numbers of the form
% C          (integer)/2.
% C     VAPP -  The initial approximation to the Beta percentile.
% C          VAPP is a real number between 0.0 and 1.0.
% C     P    -  The probability for which the inverse Beta percentile is
% C          to be evaluated.
% C          P is a real number between 0.0 and 1.0.
% C
% C EXTERNAL FUNCTION CALLED:
% C     BETA  -  BETA(A,B,X) returns the value of the density function for
% C          a Beta(A,B) random variable.
% C          A and B are elements in  {N/2 | N is a positive integer}.
% C          0 <= X <= 1.
% C
% C     DBETAI - DBETAI(X,A,B) returns the probability that a random
% C          variable from a Beta distribution having parameters (A,B)
% C          will be less than or equal to X --- P{B(A,B) <= X} = DBETAI.
% C          A and B are positive real numbers.
% C
% C***********************************************************************


% %       PARAMETER               ( ERROR = 1.0D-8, ERRAPP = 1.0D-3 )
% %       DOUBLE PRECISION          D(2:20,0:18)
%        V = VAPP;
%        VHOLD = 0.0D0
%        LOOPCT = 2
%  10    if ((DABS((V-VHOLD)/V).>=ERRAPP).AND.(LOOPCT~=0)) %THEN
%           VHOLD = V;
%           LOOPCT = LOOPCT - 1;
% %C        (          USE DBETAI TO FIND  F(V) = PROB{ BETA(A,B) <= V }. )
% %C        (          AND THEN COMPUTE    Q = (P - F(V))/f(V).           )
%           Q = (P-cdfbeta(V,A,B))/pdfbeta(V,A,B);
% %C        (                       LET D(N,K) = C(N,K)*Q**(N+K-1)/(N-1)! )
%           T = 1.0D0 - V;
%           S1 = Q*(B-1.0D0)/T;
%           S2 = Q*(1.0D0-A)/V;
%           D(2,0) = S1 + S2;
%           TAIL = D(2,0)*Q/2.0D0;
%           V = V + Q + TAIL;
%           K = 3;
%  20       IF ((DABS(TAIL/V)>ERROR) && (K<=20)) % THEN
% %C           (                                   FIRST FIND  D(2,K-2).  )
%              S1 = Q*(DBLE(K)-2.0D0)*S1/T;
%              S2 = Q*(2.0D0-DBLE(K))*S2/V;
%              D(2,K-2) = S1 + S2;
% %C           (  NOW FIND  D(3,K-3), D(4,K-4), D(5,K-5), ... , D(K-1,1). )
%              for I=3:K-1
%                  SUM1 = D(2,0)*D(I-1,K-I);
%                  BCOEFF = 1.0D0
%                  for J = 1:K-I
%                      BCOEFF = (BCOEFF*(K-I-J+1))/(J);
%                      SUM1 = SUM1 + BCOEFF*D(2,J)*D(I-1,K-I-J);
%                  end % 30   CONTINUE
%                  D(I,K-I) = SUM1 + D(I-1,K-I+1)/DBLE(I-1)
%              end % 40              CONTINUE
% %C           ( AND THEN COMPUTE D(K,0) AND USE IT TO EXPAND THE SERIES. )
%              D(K,0) = D(2,0)*D(K-1,0) + D(K-1,1)/(K-1);
%              TAIL = D(K,0)*Q/(K);
%              V = V + TAIL;
% %C           (                           CHECK FOR A DIVERGENT SERIES.  )
%              IF ((V <= 0.0D0) .OR. (V >= 1.0D0))  THEN
%                 %PRINT *,'SERIES IN BINV DIVERGES'
%                 X = -1.0D0
%                 return
%                 %GO TO 50
%              END % IF
%              K = K+1
%              GO TO 20
%           END %IF
%           GO TO 10
%        END %IF
%  50    BINV = V
%        % END BINV

% function xnew = newtonx(xstart,logPrbk,ak,bk,options)
% %NEWTONX Newton-Raphson iteration using logarithmic scale on Prb
% 
%   %tail = 'upper';
%   xnew = xstart;
%   dx   = ones(size(ak));
%   prb = dx;
%   
%   iterNo    = 0;
%   max_count = options.max_iter; %500;
%  
%   ix        = find(xnew(:));
%   while (any(ix) && iterNo<max_count)
%     
%     xi = xnew(ix);
%     prb(ix) = cdfbeta(xi,ak(ix),bk(ix));
%     logprb = max(log(prb(ix)),-1000);
% 
%     
%     dx(ix) = -(logprb - logPrbk(ix)).*prb(ix)./max(pdfbeta(xi,ak(ix),bk(ix)),realmin);
%     
%     xnew(ix)   = max(xi/10,min(xi+(1-xi)*9/10, xi - dx(ix)));
%     dx(ix) = xi-xnew(ix);
%    
%     ix = find((abs(dx) > options.releps*abs(xnew))  &  abs(dx) > options.abseps );
%     
%     iterNo = iterNo+1; 
%     disp(['Iteration ',num2str(iterNo),'  Number of points left:  ' num2str(length(ix)) ]),
%   end
%    
%   problems = abs(log(prb)-logPrbk)>sqrt(options.releps).*abs(logPrbk);
%  
%   if (iterNo == max_count)|| any(problems(:)),    
%     F = exp(logPrbk);
%     dF = prb-F;
%     if any(problems(:)),    
%       tozero = problems & F<=dF;
%       xnew(tozero) = 0;
%       toone = problems & (1-F<=-dF);
%       xnew(toone) = 1;
%       reason = 'The parameters A or B are too small or too large for the given tail probability!';
%     else
%       reason = '';
%     end
%      badcdf = find(problems(:) | (isfinite(ak(:)) & abs(dx(:))>sqrt(options.releps)));
%     didnt = unique([ix(:);badcdf]);
%     numDidnt = length(didnt);
%     didnt = didnt(1);
% 
%     badA = ak(didnt);
%     badB = bk(didnt);
%     badF = F(didnt) ; %-expm1(logPrbk(didnt));
%     baddF = dF(didnt);
%     %badX = x(didnt);
%     %outstr = sprintf('a = %g, F = %g\n',[badA(:),badF(:)].');
%     %warning('WAFO:INVBETA','invbeta did not converge for \n%s',outstr);
%     outstr = sprintf('a = %g, b=%g, F = %g dF = %g\n',[badA(:),badB(:),badF(:),baddF(:)].');
%     warning('WAFO:INVBETA','invbeta did not converge for %d values:\n%s\n%s',numDidnt,reason,outstr);    
%   end