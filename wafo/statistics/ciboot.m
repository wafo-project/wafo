function  [ci, y] = ciboot(x,fun,method,alpha,B,varargin)
%CIBOOT   Bootstrap confidence interval.
%
%	  CALL ci = ciboot(x,fun,method,alpha,B)
%	  
%   ci      = 100(1-alpha)% Confidence intervals for FUN(X) where every 
%              row is of the form [LeftLimit, PointEstimate, RightLimit]. 
%   x       = data
%   fun     = function handle or name of function
%   method  = Bootstrap method
%           1.  Normal approximation (std is bootstrap).
%           2.  Simple bootstrap principle (bad, don't use).
%           3.  Studentized, std is computed via jackknife (If T is
%             'mean' this done the fast way via the routine TESTMEAN1B).
%           4.  Studentized, std is 30 samples' bootstrap.
%           5.  Efron's  method.              (default)
%           6.  Efron's  method with bias correction (BC).  
%   alpha   = Confidence coefficent           (default 0.05)
%   B       = The number of bootstrap resamples 
%             (default B = 1000 for method 3 and fun='mean'
%                      B = 500 for method 1,2,5,6
%                      B = 200 otherwise)
%
% CIBOOT Compute a 100(1-alpha)% confidence interval for fun(X) based on a 
% bootstrap technique. Often FUN(X) is a number but it may also be a 
% vector or even a matrix. Every row of the result CI is of the form 
%
%      [LeftLimit, PointEstimate, RightLimit]
%
%  and the corresponding element of T(X) is found by noting 
%  that t = T(X); t = t(:); is used in the routine. 
%
%  Example 
%    abstol = 1e-5;
%    reltol = 0.5;
%    X =  [0.636620   0.995173   0.879845   0.292531   0.721456;...
%          0.190022   0.806725   0.448761   0.370546   0.079624].';
%    C = cov(X); 
%    method = 6;
%    ci = ciboot(X,'cov', method);
%    assert(C, [ 0.072462,   0.036949;
%                0.036949   0.078306], abstol);
%    assert(C(:), ci(:,2), -abstol)
%    assert(ci, [ 0.0144   0.0724   0.1345;
%                -0.0205   0.0369   0.0760;
%                -0.0205   0.0369   0.0760;
%                 0.0167   0.0783   0.1449], -reltol);
%
%  See also stdboot and stdjack.
 
%       Anders Holtsberg, 1994, 1998
%       Copyright (c) Anders Holtsberg
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



if nargin < 5, B = []; end
if nargin < 4 || isempty(alpha), alpha = 0.05; end
if nargin < 3 || isempty(method), method = 5; end

if min(size(x)) == 1
   x = x(:);
end


%[n,nx] = size(x);

% === 1 ================================================

if method == 1 
   if isempty(B), B = 500; end
   s = feval(@stdboot,x,fun,B,varargin{:});
   s = s(:);
   t0 = feval(fun,x,varargin{:});
   t0 = t0(:);
   z = -invnorm(alpha/2);
   ci = [t0-z*s, t0, t0+z*s];
   return;
end

% === 2 5 6 ==============================================

if method == 2 || method == 5 || method == 6
   if isempty(B), B = 500; end
   [s,y] = feval(@stdboot,x,fun,B,varargin{:});
   t0 = feval(fun,x,varargin{:});
   t0 = t0(:);
   if method == 2 || method == 5
      q = percentile(y',[alpha/2, 1-alpha/2]',1);
      if method == 2
         ci = [2*t0-q(2,:)', t0, 2*t0-q(1,:)'];
      else
         ci = [q(1,:)', t0, q(2,:)'];
      end
   else
      t00 = t0(:,ones(size(y,2),1));
      J = ((y<t00)+(y<=t00)).'/2;
      z0 = invnorm(mean(J));
      beta = cdfnorm(invnorm([alpha/2,1-alpha/2]')*...
             ones(1,length(z0))+2*ones(2,1)*z0);
      q = zeros(2,length(z0));
      for i=1:length(z0)
         q(:,i) = percentile(y(i,:),beta(:,i),1);
      end
      ci = [q(1,:)', t0, q(2,:)'];   
   end
   return
end

% === 3 'mean' ===========================================

if ischar(fun)
  funname = fun;
elseif isa(fun,'function_handle')
  funname = func2str(fun);
else
  funname = '';
end
if method == 3 && strcmp(funname,'mean')
   if isempty(B), B = 1000; end
   [dummy1,ci,dummy2,y] = testmean1boot(x,alpha,B);
   return
end

% === 3 4 ================================================

if method == 3 || method == 4
   if isempty(B), B = 200; end
   
   if method == 3 
      fun2 = @(x) stdjack(x,fun,varargin{:});
   else
      fun2 = @(x) stdboot(x,fun,30,varargin{:});
   end
   xb = x;
   t0 = feval(fun,xb,varargin{:});
   s0 = feval(@stdboot,xb,fun,200,varargin{:});
   t0 = t0(:);
   s0 = s0(:);
   y = zeros(length(t0(:)),B);
   tic

   fprintf(' B-i\n')
   for i = 1:B
      xb = rndboot(x);
      tb = feval(fun,xb,varargin{:});
      sb = feval(fun2,xb);
      tb = tb(:);
      sb = sb(:);
      y(:,i) = (tb-t0) ./ sb;
      if toc > 1, tic, fprintf('%4.0f \n',B-i), end
   end
   q = percentile(y',[alpha/2, 1-alpha/2]',1);
   ci = [t0-s0.*q(2,:)', t0, t0-s0.*q(1,:)'];
   return
end