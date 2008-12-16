function [g, t0] = hermitefun(coef,x,ma,sa,gdef)
%HERMITEFUN Calculates the transformation by a Hermite polynomial.
%            Assumption: a Gaussian process, Y, is related to the
%                      non-Gaussian process, X, by Y = g(X). 
%
% CALL:  [g,test] = hermitefun(coef,x,ma,sa,gdef)
%
%    g     = [x g(x)] a two column matrix with the transformation g(x).
%    test  = int (g(x)-x)^2 dx  where int. limits is given by X. This
%           is a measure of departure of the data from the Gaussian model.
%   coef   = [c3 c4], vector with polynomial coefficients, see below.
%   x      =  vector of x values (default linspace(-5,5,513)'*sa+ma)
%   ma, sa = mean and standard deviation of the process, respectively.
%            (default ma=0,sa=1)
%   gdef   = integer defining the transformation (default 1)
%
%  If gdef<0 hardening model (kurtosis < 3)
%     g(x) =  xn - c3(xn^2-1) - c4*(xn^3-3*xn) 
%  where 
%    xn = (x-ma)/sa
% 
%  If gdef>=0 softening model (kurtosis >= 3)
%     G(y) = mean + K*sa*[ y + c3(y^2-1) + c4*(y^3-3*y) ]
%   where
%     y  = g(x) = G^-1(x)
%     K  = 1/sqrt(1+2*c3^2+6*c4^2)
%
% See also  hermitetr, ochitr, dat2tr


% Tested on: matlab 5.3
% by pab 
% - default x is now levels([-5 5 513])*sa+ma -> 
%  better to have the discretization
%  represented with exact numbers, especially when calculating
%  derivatives of the transformation numerically.

error(nargchk(1,5,nargin))
if nargin<5||isempty(gdef), gdef =1 ; end
if nargin<4||isempty(sa), sa = 1;end
if nargin<3||isempty(ma), ma = 0; end
if nargin<2||isempty(x),  x  = linspace(-5*sa+ma,5*sa+ma,513)'; end


c3 = coef(1); c4 = coef(2);
if ~isreal(c4)||~isreal(c3)
 error('Unable to calculate the polynomial')
end
xn = (x-ma)/sa;

if gdef<0
  p = [-c4 -c3 1+3*c4 c3];
  g = [x polyval(p,xn)];
else
  g = [x zeros(length(xn),1)];
  K = 1/sqrt(1+2*c3^2+6*c4^2);
  p = K*[c4 c3 1-3*c4 -c3];%+[ 0 0 0 ma]; % polynomial coefficients
end



% Strip leading zeros and throw away.
inz= find(abs(p)>eps);
p = p(inz(1):end);
% Check if it is a strictly increasing function.
k  = length(p);
dp = (k-1:-1:1).*p(1:k-1);  % Derivative
r  = roots(dp);             % Find roots of the derivative
r(abs(imag(r))>eps)=[]; % Keep only real roots

if any(r)
  % Check if it is possible to invert the polynomial for the given
  % interval
  if gdef<0
    x00 = r;
  else
    x00 = sa*polyval(p,r)+ma;
  end
  txt1 = 'The polynomial is not a strictly increasing function.';
  txt2 = ['The derivative of g(x) is infinite at x = ' num2str(x00')];
  warning('WAFO:HERMITEFUN','%s \n %s ',txt1,txt2)
  
  txt2 = ['for the given interval x = ' num2str(x([1 end]).',3)];
  
  if any(x(1)<= x00 & x00 <= x(end)), 
    cdef = 1; 
  else
    cdef = sum( xor(x00 <= x(1) , x00 <= x(end)));
  end
  if mod(cdef,2),
    txt1 = 'Unable to invert the polynomial';
    error('WAFO:HERMITEFUN','%s \n %s ',txt1,txt2)
  end
  txt1='However, successfully inverted the polynomial ';
  disp(sprintf('%s \n %s ',txt1,txt2))
end


% Inverting the polynomial
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
switch k*(gdef>=0)
  case 0, % g(x) ready
  case 2, % Linear
    g(:,2) = xn;             % (x-p(2))/p(1);
  case 3, % Quadratic: Solve c3*u^2+u-c3 = (x-ma)/(K*sa) = xn/K
    u =  .5*(-1-sqrt(1+4*c3*(c3+xn/K)))/c3;
    if 0,%y0>y(1)
      g(:,2) = u/c3;             % choose the largest solution
    else
      g(:,2) = -(c3+xn/K)./u/c3; % choose the smallest solution
    end
  case 4, % Cubic: solve u + c3*(u.^2-1) + c4*(u.^3-3*u)= (x-ma)/(K*sa)	  
    b = 1/(3*c4); 
    x0 = c3*b ;  
    % substitue u = z-x0  and divide by c4 => z^3 + 3*c*z+2*q0  = 0
    c  = b-1-x0^2;
    q0 = x0^3-1.5*b*(x0+xn/K); % 
    if any(r) && cdef ==0 % Three real roots
      d        = sqrt(-c);
      theta1   = acos(-q0./d^3)/3;
      th2      = [0 -2*pi/3 2*pi/3];
      x1       = abs(2*d*cos(theta1(ceil(length(xn)/2)) + th2)-x0);
      [tmp ix] = min(x1);   % choose the smallest solution
      g(:,2)   = 2*d*cos(theta1 + th2(ix))-x0;
    else                 %Only one real root exist
      q1 = sqrt((q0).^2+c^3);
      % Find the real root of the monic polynomial
      A0 = (q1-q0).^(1/3);
      B0 = -(q1+q0).^(1/3);
      g(:,2) = A0+B0-x0;   % real root
      %% The other complex roots are given by
      %x= -(A0+B0)/2+(A0-B0)*sqrt(3)/2-x0;
      %x=-(A0+B0)/2+(A0-B0)*sqrt(-3)/2-x0;
    end
 case 6,
   % old call
   g(:,2) = xn;
   g(:,1) = sa*polyval(p,xn)+ma;
   %g(:,1) = sa*K*(xn + c3*(xn.^2-1) + c4*(xn.^3-3*xn))+ma;
   
   % interpolate so that g(:,1) have equidistant spacing 
   g(:,2)=smooth(g(:,1),g(:,2),1,x,1);
   g(:,1)=x;
end

if nargout>1,
  t0 = trapz(xn,(xn-g(:,2)).^2);
end

