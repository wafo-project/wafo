function [x,xlo,xup] = invgam(F,varargin)
%INVGAM Inverse of the Gamma distribution function
%
% CALL:  x = invgam(F,a,b,options)
%       [x,xlo,xup] = invgam(F,phat,options)
%
%        x = inverse cdf for the Gamma distribution evaluated at F
%  xlo,xup = 100*(1-alpha) % confidence bounds of x.
%        a = parameter, a>0
%        b = parameter, b>0 (default b=1)
%     phat = Distribution parameter struct
%            as returned from FITGAM.  
%  options = struct with fieldnames:
%         .lowertail: if TRUE (default), F = Prob[X <= x],
%                     otherwise, F = Prob[X > x].
%         .logp     : if TRUE, probability, p, input as log(p).
%         .alpha    : Confidence coefficent        (default 0.05)
%         .abseps   : Requested absolute error     (default 1e-90)
%         .releps   : Requested relative error     (default sqrt(eps))
%         .max_iter : Maximum number of iterations (default 500)
%
% INVGAM computes the inverse Gamma distribution function via a
% combination of Newton method and a bisection search.
%
% Example:
%   a=1;b=1;    
%   opt = {'lowertail',false,'logp',false}
%   F0 = [logspace(log10(realmin),-1) linspace(0.2,1-1e-3) logspace(log10(1-sqrt(eps)),log1p(-eps)/log10(10))];
%   %F0 = [logspace(-300,-1) linspace(0.11,0.5)];
%   x  = invgam(F0,a,b,opt{:});
%   F  = cdfgam(x,a,b,opt{:});
%   semilogy(abs(F-F0)./F0+eps), shg % relative error
% 
%   opt = {'lowertail',false,'logp',true}
%   x0 = [logspace(-90,log10(a/10)) linspace(a/5,3*a) logspace(log10(10*a),log10(50*a))];
%   F  = cdfgam(x0,a,b,opt{:});
%   x  = invgam(F,a,b,opt{:});
%   N0 = length(x0);
%   semilogy(1:N0,abs(x-x0)./x0+eps),shg % relative error
%
% See also  pdfgam, cdfgam, fitgam, rndgam, momgam


% Reference: Johnson, Kotz and Balakrishnan (1994)
% "Continuous Univariate Distributions, vol. 1", p. 494 ff
% Wiley

% Tested on; Matlab 5.3
% History: 
% adapted from stixbox ms 26.06.2000
% Revised by jr 01-August-2000
% - Added approximation for higher degrees of freedom (see References)
% added b parameter ms 23.08.2000
% revised pab 23.10.2000
%  - added comnsize, nargchk
% revised pab 4nov2005
%  -improved the starting guess for x.
% revised pab sept 2007
% -improved the Newton-Raphson method

%     This program is free software; you can redistribute it and/or modify
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


error(nargchk(2,inf,nargin))
options = struct('proflog',false,'alpha',0.05,...
  'lowertail',true,'logp',false,'abseps',1e-90,'releps',sqrt(eps),'max_iter',500); % default options

Np = 2;
[params,options,tmp,phat] = parsestatsinput(Np,options,varargin{:});
if numel(options)>1
  error('Multidimensional struct of distribution parameter not allowed!')
end

[a,b] = deal(params{:});
if isempty(b), b=1;end

if options.logp
  F(F>0) = NaN;
  logPrb = F;
  
  interior = log(realmin)<=logPrb & logPrb<log1p(-realmin);
else
  F(F<0 | 1<F) = NaN;
  logPrb = log(F);
  interior = (0<F& F<1); 
end

  
[iscmn,F,a,b,logPrb] = iscomnsize(F,a,b,logPrb);
if ~iscmn
  error('F, a and b must be of common size or scalar.')
end

q = zeros(size(a));

ok = (~isnan(logPrb) & a>0 & b>0);

k2=find(((options.lowertail & logPrb>=log1p(-realmin)) | (~options.lowertail & logPrb<=log(realmin)))     &ok);
if any(k2)
  q(k2) = inf;
end


k3=find(~ok);
if any(k3)
  q(k3)=nan;
end


k1 = find(interior & ok);
if any(k1)
  ak      = a(k1);
  xnew    = aprx_invgamma(F(k1),ak,options);
  if true
  logPrbk = logPrb(k1);
  smlprb = logPrbk<=log(0.5);
  logPrbk(~smlprb) = log(-expm1(logPrbk(~smlprb)));
  if options.lowertail
    if any(smlprb(:))
      xnew(smlprb) = newtonlogx(xnew(smlprb),logPrbk(smlprb),ak(smlprb),options);
    end
    if any(~smlprb(:))
      xnew(~smlprb) = newtonx(xnew(~smlprb),logPrbk(~smlprb),ak(~smlprb),options);
    end
  else
    if any(smlprb(:))
      xnew(smlprb) = newtonx(xnew(smlprb),logPrbk(smlprb),ak(smlprb),options);
    end
    if any(~smlprb(:))
      xnew(~smlprb) = newtonlogx(xnew(~smlprb),logPrbk(~smlprb),ak(~smlprb),options);
    end
  end
  end
   q(k1) = xnew;   
end

x = q.*b;

if nargout >= 2
  if isempty(phat)
    error('Must have distribution struct!')
  end
  
  alpha = options.alpha;
  if options.proflog
    error('Confidence interval usin proflog not implemented!')
    xlo = x;
    xup = x;
    for ix =1:numel(x)
      [Lp,CI] = proflog(phat,'i',1,'x',x(ix),'link',@lnkgam,'alpha',alpha);
      xlo(ix) = CI(1);
      xup(ix) = CI(2);
    end
  else
    % Compute confidence bounds using the delta method.
    % dFda = dFdx*dxda = dFdx*dqda*b =>  dqda = -dFda/dFdx
    dqda = -dgammainc(q,a,tail)./pdfgam(q,a);
    
    pcov = phat.covariance;
    
    % Approximate the variance of x=q*b on the log scale.
    %    dlogx/da = dlogx/dq * dqda = dqda/q
    %    dlogx/db = 1/b
    logx = log(x);
    varlogx = pcov(1,1).*(dqda./q).^2 + 2.*pcov(1,2).*(dqda./x) + pcov(2,2)./(b.^2);
    if any(varlogx(:) < 0)
      error('Covariance must be a positive semi-definite matrix.');
    end
    zcrit = -invnorm(alpha/2)* sqrt(varlogx);
    
    % Convert back to original scale
    xlo = exp(logx - zcrit);
    xup = exp(logx + zcrit);
  end
end


function xnew = newtonx(xstart,logPrbk,ak,options)
%NEWTONX Newton-Raphson iteration on the upper tail

  tail = 'upper';
  xnew = xstart;
  dx   = ones(size(ak));
  
  iterNo    = 0;
  max_count = options.max_iter; %500;
 
  ix        = find(xnew(:));
  while (any(ix) && iterNo<max_count)
    
    xi = xnew(ix);
    logprb = max(gammaincln(xi,ak(ix),tail),-1000);
    prb = exp(logprb);
    
    dx(ix) = -(logprb - logPrbk(ix)).*prb./max(pdfgam(xi,ak(ix)),realmin);
    
    xnew(ix)   = max(xi/10,min(10*xi, xi - dx(ix)));
    dx(ix) = xi-xnew(ix);
   
    ix = find((abs(dx) > options.releps*abs(xnew))  &  abs(dx) > options.abseps );
    
    iterNo = iterNo+1; 
    %disp(['Iteration ',num2str(iterNo),...
    %'  Number of points left:  ' num2str(length(ix)) ]),
  end
   
 
  if (iterNo == max_count), 
     badcdf = find(isfinite(ak(:)) & abs(dx(:))>sqrt(options.releps));
    didnt = unique([ix(:);badcdf]);
    %didnt = didnt(1);
    badA = ak(didnt);
    badF = -expm1(logPrbk(didnt));
    %badX = x(didnt);
    outstr = sprintf('a = %g, F = %g\n',[badA(:),badF(:)].');

    warning('WAFO:INVGAM','invgam did not converge for \n%s',outstr);
%     for iy = 1:length(didnt)
%       [x(didnt(iy)),flag,val] = fzero(sprintf('(cdfgam(x,%g)-%g)',badA(iy),badF(iy)),badX(iy));
%     end
    
  end

 function xnew = newtonlogx(xstart,logPrbk,ak,options)
 %NEWTONLOGX Newton-Raphson iteration on the lower tail

  
  tail = 'lower';
  xnew = xstart;
  uk   = log(xnew);
  du     = ones(size(ak));
  dx     = du;
  logprb = du;
  
  iterNo = 0;
  max_count = options.max_iter; %500;
 
  ix        = find(uk(:));
  while (any(ix) && iterNo<max_count)
    
    xi     = xnew(ix);
    logprb(ix) = max(gammaincln(xi,ak(ix),tail),-1000);
    prb    = exp(logprb(ix));
    
    du(ix) = (logprb(ix) - logPrbk(ix)).*prb./max(xi.*pdfgam(xi,ak(ix)),realmin);
    
    du(ix)   = max(-3, min(3, du(ix)));
    uk(ix)   = uk(ix)-du(ix);
    xnew(ix) = exp(uk(ix));
    dx(ix)   = xnew(ix)-xi;
    
    ix = find((abs(dx) > options.releps*abs(xnew))  &  abs(dx) > options.abseps );
    
    iterNo = iterNo+1;  
    %disp(['Iteration ',num2str(iterNo),'  Number of points left:  ' num2str(length(ix)) ]),
  end
  prb = exp(logprb);
  F = exp(logPrbk);
  dF = prb-F;
  problems = (abs(dF)> sqrt(options.releps).*F & abs(dF)> options.abseps);
  
  if (iterNo == max_count || any(problems(:))),    
    if any(problems(:)),
      tozero = problems & F<=dF;
      xnew(tozero) = 0;
      reason = 'The parameter A is too small for the given lower tail probability!';
    else
      reason = '';
    end
    
     badcdf = find(problems(:) | isfinite(ak(:)) & abs(dx(:))>sqrt(options.releps));
    didnt = unique([ix(:);badcdf]);
    numDidnt = length(didnt);
    didnt = didnt(1);
    badA = ak(didnt);
    badF = F(didnt);
    baddF = dF(didnt);
    %badX = x(didnt);
    outstr = sprintf('a = %g, F = %g dF = %g\n',[badA(:),badF(:),baddF(:)].');
    warning('WAFO:INVGAM','invgam did not converge for %d values:\n%s\n%s',numDidnt,reason,outstr);    
  end
  
  



  function q = aprx_invgamma(F,a,options)
  %  APRX_INVGAMMA Approximate inverse of the gamma distribution 
  %
  %  CALL 
  %
  
  
  if nargin>2
    z = invnorm(F,0,1,options);
  else
    z = invnorm(F);
  end
%     q = normq2gamq(z,a);
%   
% function q = normq2gamq(z,a)
  
  % Z = quantiles from a normal distribution
      
  % % This approximation is from Johnson et al, p. 348, Eq. (17.33).
  % is OK for a>1 and quite good for a>100
  sa = sqrt(a);
  q  = a.*(1 - 1./(9.*a)+z./(3.*sa)  ).^3 ;
  
  % For small q and a, gammainc(q,a) ~= q^a,
  sml = (q<0.025) | (q<0.9 & a<=0.3);
  if any(sml(:))
    if ~isscalar(a), a  = a(sml); end
    q(sml) = max(cdfnorm(z(sml)).^(1./a),100*realmin);
  end
    
   

 
return

% % Find lowest possible probability
% 
% ginv2 = @(F,a) a.*(1-(1./(9.*a))+invnorm(F)./(3.*sqrt(a))).^3;
% 
% xlo =  max(a.*(1-(1./(9.*a))-37./(3.*sqrt(a))).^3,realmin)
% 
% myfun1 = @(x1) gammaincln((x1),ak(1),tail);
%     
%     [x0,y0] = fplot(myfun1,[xlo ak(1)/10]);
% 
% 
% ak = logspace(0,6,256);
% xmax = logspace(0,6,257);
% [Xmax,AK] = meshgrid(xmax,ak);
% logg = gammaincln(Xmax,AK);
% logg(logg<log(realmin)) = log(realmin);
% f = wdata(logg,{xmax,ak});
% f.labels = {'X','A'}
% v = f.contourf; fcolorbar(v)
% set(gca,'xscale','log')
% set(gca,'yscale','log') 
% 
% 
% scalet = 'linear';
% a = logspace(-10,1,101);x = logspace(-10,log(0.025)); [A,X] = meshgrid(a,x);F = gammainc(X,A), 
% xn1 = F.^(1./A); 
% z = invnorm(F);
% xn2 = max(A.*(1 - 1./(9.*A)+z./(3.*sqrt(A))  ).^3,realmin);
% loga  = log(A);
% sigma2 = (log1p(A) - loga);
% mu    = loga - 0.5 .* sigma2;
% xn3 = exp(mu + sqrt(sigma2).*z);
% v = -5:1;
% figure(1),cs = contourf(log10(a),x,log10(abs(X-xn1)./X),v); fcolorbar(v),vline(log10(0.3)),set(gca,'yscale',scalet)
% figure(2),cs = contourf(log10(a),x,log10(abs(X-xn2)./X),v); fcolorbar(v),vline(log10(0.3)),set(gca,'yscale',scalet)
% figure(3),cs = contourf(log10(a),x,log10(abs(X-xn3)./X),v); fcolorbar(v),vline(log10(0.3)),set(gca,'yscale',scalet)
% 
% 
% 
% scalet = 'log';
% F0 = [logspace(-3,-1) linspace(0.11,0.5)];
% a = 100000;x = invgam2(F0,a); [A,X] = meshgrid(a,x);F = gammainc(X,A), 
% xn1 = F.^(1./A); 
% z = invnorm(F);
% xn2 = max(A.*(1 - 1./(9.*A)+z./(3.*sqrt(A))  ).^3,realmin);
% loga  = log(A);
% sigma2 = (log1p(A) - loga);
% mu    = loga - 0.5 .* sigma2;
% %mu = loga;
% %sigma2 = 1./(a);
% xn3 = exp(mu + sqrt(sigma2).*z);
% v = -5:1;
% figure(1),plot(x,(abs(X-xn1))); set(gca,'yscale',scalet)
% figure(2),plot(x,(abs(X-xn2))); set(gca,'yscale',scalet)
% figure(3),plot(x,(abs(X-xn3))); set(gca,'yscale',scalet)
% 
