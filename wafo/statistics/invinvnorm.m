function [x,xlo,xup] = invinvnorm(F,varargin)
%INVINVNORM Inverse of the Inverse Gaussian distribution function
%
% CALL:  x = invinvnorm(F,m,l,options)
%
%        x = inverse cdf for the Inverse Gaussian distribution evaluated at F
%  xlo,xup = 100*(1-alpha) % confidence bounds of x.
%     m,l  = parameters (see pdfinvnorm)
%     phat = Distribution parameter struct
%            as returned from FITINVNORM.  
%  options = struct with fieldnames:
%         .lowertail: if TRUE (default), F = Prob[X <= x],
%                     otherwise, F = Prob[X > x].
%         .logp     : if TRUE, probability, p, returned as log(p).
%         .alpha    : Confidence coefficent    (default 0.05)
%         .abseps   : Requested absolute error     (default 1e-90)
%         .releps   : Requested relative error     (default sqrt(eps))
%         .max_iter : Maximum number of iterations (default 500)
%
% Example:
%   a=1;b=1;    
%   opt = {'lowertail',false,'logp',false};
%   F0 = [logspace(log10(realmin),-1) linspace(0.2,1-1e-3),...
%         logspace(log10(1-sqrt(eps)),log1p(-eps)/log10(10))];
%   %F0 = [logspace(-300,-1) linspace(0.11,0.5)];
%   x  = invinvnorm(F0,a,b,opt{:});
%   F  = cdfinvnorm(x,a,b,opt{:});
%   semilogy(abs(F-F0)./F0+eps); % relative error
%
%   close all;
%
% See also pdfinvnorm, cdfinvnorm, rndinvnorm, fitinvnorm, mominvnorm

% Reference: 
% Cohen & Whittle, (1988) "Parameter Estimation in Reliability
% and Life Span Models", p. 259 ff, Marcel Dekker.


% Tested on; Matlab 5.3
% History: 
% adapted from stixbox ms 14.08.2000
% ad hoc solutions to improve convergence added ms 15.08.2000
% revised pab 25.10.2000
% - added nargchk + comnsize
% - changed the ad hoc solutions a bit

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


error(nargchk(3,9,nargin))
options = struct('covariance',[],'alpha',0.05,...
  'lowertail',true,'logp',false,'abseps',1e-90,'releps',sqrt(eps),'max_iter',500); % default options
Np = 2;
[params,options] = parsestatsinput(Np,options,varargin{:});
if numel(options)>1
  error('Multidimensional struct of distribution parameter not allowed!')
end

[m,l] = deal(params{:});
m(m<=0) = nan;
l(l<=0) = nan;

if options.logp
  F(F>0) = NaN;
  logPrb = F;
  
  interior = log(realmin)<=logPrb & logPrb<log1p(-realmin);
else
  F(F<0 | 1<F) = NaN;
  logPrb = log(F);
  interior = (0<F& F<1); 
end


[icode, F, m, l,logPrb] = iscomnsize (F,m, l,logPrb);
if ~icode 
  error ('F, m and l must be of common size or scalar');
end

x=zeros(size(F));
ok=(F>=0 & F<=1 &(m>0&(l>0)));

k=find(interior & ok);
if any(k)
  mk=m(k);lk=l(k);  Fk=F(k);
  % Need  a better starting guess here for 
  %  1) l<1 and F close to 1
  %  2) l>5 and F close to 1 or 0
  % Supply a starting guess
  [mk,v]=mominvnorm(mk,lk);
  temp = log(v + mk .^ 2); 
  mu = 2 * log(mk) - 0.5 * temp;
  sa = -2 * log(mk) + temp;
  xk = invlognorm(Fk,mu,sa,options);

  
  %
  
  %  x(k)= newton(xk,Fk,mk,lk,v,options);
  
  logPrbk = logPrb(k);
  smlprb = logPrbk<=log(0.5);
  logPrbk(~smlprb) = log(-expm1(logPrbk(~smlprb)));
  if options.lowertail
    if any(smlprb(:))
      xk(smlprb) = newtonlogx(xk(smlprb),logPrbk(smlprb),mk(smlprb),lk(smlprb),options);
    end
    if any(~smlprb(:))
      xk(~smlprb) = newtonx(xk(~smlprb),logPrbk(~smlprb),mk(~smlprb),lk(~smlprb),options);
    end
  else
    if any(smlprb(:))
      xk(smlprb) = newtonx(xk(smlprb),logPrbk(smlprb),mk(smlprb),lk(smlprb),options);
    end
    if any(~smlprb(:))
      xk(~smlprb) = newtonlogx(xk(~smlprb),logPrbk(~smlprb),mk(~smlprb),lk(~smlprb),options);
    end
  end
  x(k) = xk;
end

k1=find(F==1&ok);
if any(k1)
  tmp=Inf;
  x(k1)=tmp(ones(size(k1)));
end

k2=find(~ok);
if any(k1)
  tmp=NaN;
  x(k2)=tmp(ones(size(k2)));
end


if nargout>1
% TODO % Implement  xlo and xup
 warning('WAFO:INVINVNORM','xlo and xup not implemented yet')
 xlo = nan;
 xup = xlo;
end

function xk = newton(xk,Fk,mk,lk,v,options)
 dx = ones(size(xk));
  options.logp = false;
  ix =find(xk);
  count=1;  max_count=options.max_iter;
  mstep=2*sqrt(v); % maximum step in each iteration
  while (any(ix)&&(count<max_count));
    count=count+1;
    xi=xk(ix);mi=mk(ix);li=lk(ix);
    dx(ix) = (cdfinvnorm(xi,mi,li,options) - Fk(ix))./pdfinvnorm(xi,mi,li);
    dx(ix) = min(abs(dx(ix)),mstep(ix)/count).*sign(dx(ix)); 
    xi = xi - dx(ix);
    % Make sure that the current guess is larger than zero.
    xk(ix) = xi + 0.5*(dx(ix) - xi) .* (xi<=0);
    
    ix=find((abs(dx) > options.releps*abs(xk))  &  abs(dx) > options.abseps);
    disp(['Iteration ',num2str(count),'  Number of points left:  ' num2str(length(ix)) ]),
  end

  if (count==max_count)
    disp('Warning: INVINVNORM did not converge!')
    disp(['The last steps were all less than: ' num2str(max(abs(dx(ix))))])
 
    badcdf = find(isfinite(mk(:)) & abs(dx(:))>sqrt(options.releps));
    didnt = unique([ix(:);badcdf]);
    %didnt = didnt(1);
    badA = mk(didnt);
    badF = Fk(didnt);
    %badX = x(didnt);
    outstr = sprintf('a = %g, F = %g\n',[badA(:),badF(:)].');
     numDidnt = numel(didnt);
    warning('WAFO:INVINVNORM','invinvnorm did not converge for %d values:\n%s\n%s',numDidnt,outstr);    
 
%     for iy = 1:length(didnt)
%       [x(didnt(iy)),flag,val] = fzero(sprintf('(cdfgam(x,%g)-%g)',badA(iy),badF(iy)),badX(iy));
%     end
    
  end

function xnew = newtonx(xstart,logPrbk,mk,lk,options)
%NEWTONX Newton-Raphson iteration on the upper tail
 options.lowertail = false;
 options.logp = true;
  xnew = xstart;
  dx   = ones(size(mk));
  
  iterNo    = 0;
  max_count = options.max_iter; %500;
 
  ix        = find(xnew(:));
  while (any(ix) && iterNo<max_count)
    
    xi = xnew(ix);
    logprb = max(cdfinvnorm(xi,mk(ix),lk(ix),options),-1000);
    prb = exp(logprb);
    
    dx(ix) = -(logprb - logPrbk(ix)).*prb./max(pdfinvnorm(xi,mk(ix),lk(ix)),realmin);
    if iterNo>100
      dx(ix) = dx(ix)*10/iterNo;
    end
    xnew(ix)   = max(xi/10,min(10*xi, xi - dx(ix)));
    dx(ix) = xi-xnew(ix);
   
    ix = find((abs(dx) > options.releps*abs(xnew))  &  abs(dx) > options.abseps );
    
    iterNo = iterNo+1; 
    %disp(['Iteration ',num2str(iterNo),'  Number of points left:  ' num2str(length(ix)) ]),
  end
   
 
  if (iterNo == max_count), 
     badcdf = find(isfinite(mk(:)) & abs(dx(:))>sqrt(options.releps));
    didnt = unique([ix(:);badcdf]);
    %didnt = didnt(1);
    badA = mk(didnt);
    badF = -expm1(logPrbk(didnt));
    %badX = x(didnt);
    outstr = sprintf('a = %g, F = %g\n',[badA(:),badF(:)].');

    numDidnt = numel(didnt);
    warning('WAFO:INVINVNORM','invinvnorm did not converge for %d values:\n%s\n%s',numDidnt,outstr);    
%     for iy = 1:length(didnt)
%       [x(didnt(iy)),flag,val] = fzero(sprintf('(cdfgam(x,%g)-%g)',badA(iy),badF(iy)),badX(iy));
%     end
    
  end

 function xnew = newtonlogx(xstart,logPrbk,mk,lk,options)
 %NEWTONLOGX Newton-Raphson iteration on the lower tail

  
  options.lowertail = true;
  options.logp = true;
  xnew = xstart;
  uk   = log(xnew);
  du     = ones(size(mk));
  dx     = du;
  logprb = du;
  
  iterNo = 0;
  max_count = options.max_iter; %500;
 
  ix        = find(uk(:));
  while (any(ix) && iterNo<max_count)
    
    xi     = xnew(ix);
    logprb(ix) = max(cdfinvnorm(xi,mk(ix),lk(ix),options),-1000);
    prb    = exp(logprb(ix));
    
    du(ix) = (logprb(ix) - logPrbk(ix)).*prb./max(xi.*pdfinvnorm(xi,mk(ix),lk(ix)),realmin);
    
    
    du(ix)   = max(-3, min(3, du(ix)));
    uk(ix)   = uk(ix)-du(ix);
    xnew(ix) = exp(uk(ix));
    dx(ix)   = xnew(ix)-xi;
    
    ix = find((abs(dx) > options.releps*abs(xnew))  &  abs(dx) > options.abseps );
    
    iterNo = iterNo+1;  
   % disp(['Iteration ',num2str(iterNo),'  Number of points left:  ' num2str(length(ix)) ]),
  end
  prb = exp(logprb);
  F = exp(logPrbk);
  dF = prb-F;
  problems = (abs(dF)> sqrt(options.releps).*F & abs(dF)> options.abseps);
  
  if (iterNo == max_count || any(problems(:))),    
    if any(problems(:)),
      tozero = problems & F<=dF;
      xnew(tozero) = 0;
      reason = 'The parameter m is too large for the given lower tail probability!';
    else
      reason = '';
    end
    
     badcdf = find(problems(:) | isfinite(mk(:)) & abs(dx(:))>sqrt(options.releps));
    didnt = unique([ix(:);badcdf]);
    numDidnt = length(didnt);
    didnt = didnt(1);
    badA = mk(didnt);
    badF = F(didnt);
    baddF = dF(didnt);
    %badX = x(didnt);
    outstr = sprintf('a = %g, F = %g dF = %g\n',[badA(:),badF(:),baddF(:)].');
    warning('WAFO:INVINVNORM','invinvnorm did not converge for %d values:\n%s\n%s',numDidnt,reason,outstr);    
  end
