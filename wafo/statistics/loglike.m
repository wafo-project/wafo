function [LL,pcov,H]=loglike(phat,data,varargin)
%LOGLIKE Negative Log-likelihood function.
% 
% CALL:  [L,pcov,H] = loglike(phat,data,opt1,opt2,...,optN,dist);
%
%       L    = -sum(log(f(phat|data))), i.e., the negative log-likelihood 
%              function with parameters phat given the data.
%       pcov = Asymptotic covariance matrix of phat (if phat is estimated by
%              a maximum likelihood or maximum product spacing method).
%       H    = Hessian matrix ( double derivative of loglikelihood function)
%       phat = [A1,A2,...,Am], vector of distribution parameters.
%       data = x1, or cellarray of vectors, ie., {x1,x2,..xN} 
%              where each xi are column vectors of data points.
%  opt1,opt2,
%  ...,optN  = additional input to distribution
%       dist = function handle or string defining the PDF to use.
%           
%  LOGLIKE is a utility function for maximum likelihood estimation. 
%  This works on any PDF having the following calling syntax:
%
%       f = pdftest(x1,x2,..,xN,A1,A2,..,Am,opt1,opt2,...,optN)
%
% where x1,x2,...,xN contain the points to be evaluated and A1,A2... are
% the distribution parameters. 
%
% Example: MLE and asymptotic covariance of phat:
%   R = rndweib(1,3,0,100,1);                      
%   phat0 = [1.3 2.5];                             %initial guess
%   phat = fminsearch(@(x)loglike(x,R,'pdfweib') ,phat0);
%   [L, cov] = loglike(phat,R,'pdfweib');
%   [phat2] = fitweib(R);                     % compare with fitweib
%
% See also logps, proflog


% Copyright (C) 2000 Per A. Brodtkorb
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



%Tested on: matlab 5.3
% History:
% revised pab May 2007
% -added xmin
% -made sure LL is a number if x is NaN or infinite
% -dist can now also be a function handle.
% by pab 31.10.2000


error(nargchk(3,inf,nargin))

params = num2cell(phat(:).',1);
logp = true;
xmin = realmin; % trick to avoid log of zero
if logp
  myfun = @(x) max(x,100*log(realmin));
else
  myfun = @(x) log(x+xmin);
end
options  = cat(2,varargin(1:end-1),{'logp',logp}); % cell array of vectors with data points
dist   = varargin{end};
switch class(dist)
  case {'char','function_handle'} % OK
  otherwise
  error('Distribution is unspecified')
end
if isnumeric(data)
  data = {data};
end
d = size(data{1});
if any(d(2:end)>1)
  for ii=1:length(data)
    data{ii}  = data{ii}(:); %% make sure it is a vector.
  end
end

%xmin = eps;
x0 = feval(dist,data{:},params{:},options{:});
if any(~isfinite(x0))
  LL = -sum(myfun(x0(isfinite(x0)))) + 100*log(realmax)*sum(~isfinite(x0));
else
  LL = -sum(myfun(x0)); %negative log likelihood function
end


if nargout >1
  Nd     = length(x0);
  np     = length(params);
  
  % pab 07.01.2001: Always choose the stepsize h so that
  % it is an exactly representable number.
  % This is important when calculating numerical derivatives and is
  % accomplished by the following.
  delta  = donothing(eps^0.4+2)-2;
  delta2 = delta^2;
  
  diffmethod = 4;
  switch  diffmethod,
    case 1,
      
    % Approximate 1/(nE( (d L(x|theta)/dtheta)^2)) with
    %             1/(d L(theta|x)/dtheta)^2  
    % This is a bad approximation especially for the off diagonal elements
    [pcov,H] = estcov1;
    case 2,
      % Approximate 1/(nE( (d L(x|theta)/dtheta)^2)) with
      %             1/(d^2 L(theta|x)/dtheta^2)
      %  using central differences
      % This is usually a much better estimate than case 1 and slightly
      % faster than case 3.
      
      [pcov,H] = estcov2;
      if any(diag(pcov)<0)
         warning('WAFO:loglike','Estimated Covariance may not be accurate!')
        [pcov,H] = estcov1;
      end
    
    case 3,
      % Approximate 1/(nE( (d L(x|theta)/dtheta)^2)) with
      %             1/(d^2 L(theta|x)/dtheta^2)
      % using differentiation matrices
      % This is the same as case 2 when N=3
      [pcov, H] = estcov3;
    case 4,
      % Approximate the derivate using the function Hessian
      %  This is much more accurate than the methods above.
      fun = @(phat1)myloglike(phat1, data, myfun, dist, options);
      [pcov, H] = estcov4(fun, phat);
      return
      if any(det(pcov)<0)
        [pcov,H] = estcov2;
        if any(det(pcov)<0)
           warning('WAFO:loglike','Estimated Covariance may not be accurate!')
          [pcov,H] = estcov1;
        end
      end
  end
end

  function [pcov,H] = estcov1
    % Approximate 1/(nE( (d L(x|theta)/dtheta)^2)) with
    %             1/(d L(theta|x)/dtheta)^2
    % This is a bad approximation especially for the off diagonal elements
    xp     = zeros(Nd,np);
    sparam = params;
    for ix = 1:np,
      sparam{ix}= params{ix}+delta;
      xp(:,ix) = feval(dist,data{:},sparam{:},options{:});
      sparam{ix}= params{ix};
    end
    J = (myfun(xp) - repmat(myfun(x0),1,np))./delta;
    [Q,R]= qr(J,0);
    Rinv = R\eye(np);
     
    pcov = Rinv'*Rinv;
    H = -pinv(pcov);
  end % estcov1

  function [pcov,H] = estcov2
    % Approximate 1/(nE( (d L(x|theta)/dtheta)^2)) with
    %             1/(d^2 L(theta|x)/dtheta^2) 
    %  using central differences
    % This is usually a much better estimate than case 1 and slightly
    % faster than case 3.
    H = zeros(np);             % Hessian matrix
    for ix=1:np,
      sparam = params;
      sparam{ix}= params{ix}+delta;
      x  = feval(dist,data{:},sparam{:},options{:}); 
      fp = sum(myfun(x));
      sparam{ix} = params{ix}-delta;
      x  = feval(dist,data{:},sparam{:},options{:}); 
      fm = sum(myfun(x));
      H(ix,ix) = (fp+2*LL+fm)/delta2;
      for iy=ix+1:np,
        sparam{ix} = params{ix}+delta;
        sparam{iy} = params{iy}+delta;
        x   = feval(dist,data{:},sparam{:},options{:});
        fpp = sum(myfun(x));
        sparam{iy} = params{iy}-delta;
        x   = feval(dist,data{:},sparam{:},options{:});
        fpm = sum(myfun(x));
        sparam{ix} = params{ix}-delta;
        x   = feval(dist,data{:},sparam{:},options{:});
        fmm = sum(myfun(x));
        sparam{iy} = params{iy}+delta;
        x   = feval(dist,data{:},sparam{:},options{:});
        fmp = sum(myfun(x));
        H(ix,iy) = ((fpp+fmm)-(fmp+fpm))/(4*delta2);
        H(iy,ix) = H(ix,iy);
        sparam{iy} = params{iy};
      end
    end
    
    % invert the Hessian matrix (i.e. invert the observed information number)
    pcov = -pinv(H); 
  end % estcov2

  function [pcov,H] = estcov3
    % Approximate 1/(nE( (d L(x|theta)/dtheta)^2)) with
    %             1/(d^2 L(theta|x)/dtheta^2) 
    % using differentiation matrices
    % This is the same as case 2 when N=3
    if 1, % 
      xn =[     1;     0;    -1];
      D(:,:,1) =[...
	    1.5000   -2.0000    0.5000;...
	    0.5000    0.0000   -0.5000;...
	    -0.5000    2.0000   -1.5000];
      D(:,:,2) =[...
	    1.0000   -2.0000    1.0000;...
	    1.0000   -2.0000    1.0000;...
	    1.0000   -2.0000    1.0000];
    else % If you have the differentiation matrix suite toolbox you may 
      %use this 
      % By increasing N better accuracy might be expected
      %N=3; % NB!: N must be odd
      %[xn D]=chebdif(N,2);  % construct differentiation matrix
    end
    N=length(xn);
    xn=xn.';
    % Construct differentiation matrices
    D11 = kron(D(:,:,1),D(:,:,1));  
    %D20 = kron(D(:,:,2),eye(N));
    %D02 = kron(eye(N),D(:,:,2));
    H   = zeros(np);           % Hessian matrix
    LL2 = zeros(N,N);          % Log likelihood evaluated at phat and in
                               % the vicinity of phat
			       
    N2 = (N+1)/2;              % The middle indices
    LL2(N2,N2) = -LL; % = sum(log(x))
    for ix=1:np,
      sparam = params;
      for iy = [1:N2-1 N2+1:N];
        sparam{ix}= params{ix}+xn(iy)*delta;
        x = feval(dist,data{:},sparam{:},options{:});
        %sparam{ix} = params{ix}+xn(ones(Nd,1),iy)*delta;
        %x = feval(dist,repmat(data,[1 length(iy)]),sparam{:});
        LL2(iy,N2) = sum(myfun(x)).';
      end
      %sparam=params;
      H(ix,ix) = (D(N2,:,2)*LL2(:,N2))/delta2;
      for iy=ix+1:np,
        for iz=1:N,
          sparam{ix} = params{ix}+xn(iz)*delta;
          for iw=[1:N2-1 N2+1:N];
            sparam{iy}= params{iy}+xn(iw)*delta;
            x = feval(dist,data{:},sparam{:},options{:});
            %sparam{iy} = params{iy}+xn(ones(Nd,1),iw)*delta;
            %x = feval(dist,repmat(data,[1,length(iw)]),sparam{:});
            LL2(iz,iw) = sum(myfun(x));
          end
        end
        H(ix,iy) = D11((N^2+1)/2,:)*LL2(:)/delta2;
        H(iy,ix) = H(ix,iy);
      end
    end
    % invert the Hessian matrix (i.e. invert the observed information number)
    pcov = -pinv(H); 
  end

  function [pcov,H,err] = estcov4(fun, phat),
    [H,err] = hessian(fun, phat,'RombergTerms',4);
    pcov = -pinv(H);
  end
  
end % loglike
function L = myloglike(phat1, data, myfun, dist, options)
    sparam2 = num2cell(phat1,1);
    x  = feval(dist,data{:},sparam2{:},options{:}); 
    L = sum(myfun(x));
 end
  
function y=donothing(x)
  y=x;
end