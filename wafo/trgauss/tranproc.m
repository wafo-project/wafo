function y = tranproc(x,ff)
%TRANPROC Transforms process X and up to four derivatives 
%         using the transformation f.
%
%  CALL:  y = tranproc(x,f);
%
%     x = input data matrix with 1+N columns, [X X1 ... XN], where
%         Xi  is the i'th time derivative of  X. 0<=N<=4.
%     f = [x,f(x)], transform function, y = f(x).
%     y = output data matrix with 1+N columns, [Y Y1 ...YN], of
%         transformed data, where Y = f(X) and  Yi is the i'th time
%         derivative of Y = f(X).
%
% By the basic rules of derivation:
%   Y1 = f'(X)*X1
%   Y2 = f''(X)*X1^2 + f'(X)*X2
%   Y3 = f'''(X)*X1^3 + f'(X)*X3 + 3*f''(X)*X1*X2
%   Y4 = f''''(X)*X1^4 + f'(X)*X4 + 6*f'''(X)*X1^2*X2 
%      + f''(X)*(3*X2^2 + 4*X1*X3) 
%
% The differentiation of  f  is performed numerically with a central 
% difference method with linear extrapolation towards the beginning 
% and end of f, respectively. 
%
% Example: % Derivative of g and the transformed Gaussian model.
%  x = linspace(-6,6,501)';
%  g = hermitetr(x);
%  gder = tranproc([g(:,1) ones(size(g,1),1)],g);
%  gder(:,1) = g(:,1);  
%  plot(g(:,1),[g(:,2),gder(:,2)])
%  plot(g(:,1),pdfnorm(g(:,2)).*gder(:,2),g(:,1),pdfnorm(g(:,1)))
%  legend('Transformed model','Gaussian model')
%
% See also  trangood.

% Tested on: matlab 5.1
% history:
% revised pab 09.01.2001
% -added check on hn to make sure the spacing is not too dense
% -added nmax
% -updated help header
% by ???

%error(nargchk(2,2,nargin))
narginchk(2,2)
N    = size(x,2)-1; % N = number of derivatives
nmax = ceil((max(ff(:,1))-min(ff(:,1)))*10^(7/max(N,1)));
f    = trangood(ff,size(ff,1),min(x(:,1)),max(x(:,1)),nmax);

n  = size(f,1);
y  = x;
xu = 1+(n-1)*(x(:,1)-f(1,1))/(f(n,1)-f(1,1));

fi = floor(xu);

i  = find(fi==n);
fi(i) = fi(i)-1;

xu = xu-fi;
y(:,1) = f(fi,2)+(f(fi+1,2)-f(fi,2)).*xu;


if N>0
  hn = f(2,1)-f(1,1);
  if hn^N<sqrt(eps)
    disp('Numerical problems may occur for the derivatives in tranproc.')
    warning('WAFO:TRANPROC','The sampling of the transformation may be too small.')
  end
  % Transform X with the derivatives of  f.
  fxder = zeros(size(x,1),N);
  fder  = f;
  for k=1:N, % Derivation of f(x) using a difference method.
    n = size(fder,1);
    %fder = [(fder(1:n-1,1)+fder(2:n,1))/2 diff(fder(:,2))./diff(fder(:,1))];
    fder = [(fder(1:n-1,1)+fder(2:n,1))/2 diff(fder(:,2))/hn];
    fxder(:,k) = tranproc(x(:,1),fder);
  end;
   %(-fder(ix+2,2)+8*fder(ix+1,2) - ...
%	      8*fder(ix-1,2)+fder(ix-2,2))./(12*hn);
  % Calculate the transforms of the derivatives of X.
  % First time derivative of y: y1 = f'(x)*x1
  y(:,1+1)=fxder(:,1).*x(:,1+1); 
  if N>1
    % Second time derivative of y: 
    %             y2 = f''(x)*x1.^2+f'(x)*x2
    y(:,1+2)=fxder(:,2).*x(:,1+1).^2 + fxder(:,1).*x(:,1+2);
    if N>2
      % Third time derivative of y: 
      %      y3 = f'''(x)*x1.^3+f'(x)*x3 +3*f''(x)*x1*x2
      y(:,1+3)=fxder(:,3).*x(:,1+1).^3 + fxder(:,1).*x(:,1+3) + ...
	  3*fxder(:,2).*x(:,1+1).*x(:,1+2);
      if N>3
	 % Fourth time derivative of y: 
	 %    y4 = f''''(x)*x1.^4+f'(x)*x4
	 %    +6*f'''(x)*x1^2*x2+f''(x)*(3*x2^2+4x1*x3) 
	y(:,1+4)=fxder(:,4).*x(:,1+1).^4 + fxder(:,1).*x(:,1+4) + ...
	    6*fxder(:,3).*x(:,1+1).^2.*x(:,1+2) + ...
	    fxder(:,2).*(3*x(:,1+2).^2+4*x(:,1+1).*x(:,1+3));
	if N>4
	  warning('WAFO:TRANPROC',['Transformation of derivatives of order>4 not supported' ...
		' in tranproc.'])
	end
      end
    end
  end
end



