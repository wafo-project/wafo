function [dyda,y,d2yda] = dgammainc(x,a,tail)
%DGAMMAINC Incomplete gamma function with derivatives.
%
%  CALL [dyda,y,d2yda] = dgammainc(x,a,tail)
%
%   DYDA  = First derivative of the incomplete gamma function, with respect
%           to its second argument, A.
%   Y     = Incomplete gamma function itself.
%   D2YDA = Second derivative of the incomplete gamma function with respect
%           to its second argument, A.
%   X,A   = arguments
%   tail  = 'lower' (default) or 'upper'
% 
%   The incomplete gamma function is defined as:
%
%     gammainc(x,a) = 1 ./ gamma(a) .*
%        integral from 0 to x of t^(a-1).*exp(-t) dt
%
%  The upper incomplete gamma function is defined as
%   1 - gammainc(x,a).
%
% Example
%  a  = linspace(0.1,3);
%  x1 = 1;
%  plot(a,dgammainc(x1,a))
%
%  x  = linspace(0,3);
%  [X,A] = meshgrid(x,a);
%  contour(x,a,dgammainc(X,A))
% 
% See also gamma, gammaln, gammainc, gammaincln.


if nargin < 3
  tail = 'lower';
end
% pab 07.01.2001: Always choose the stepsize h so that
% it is an exactly representable number.
% This is important when calculating numerical derivatives and is
% accomplished by the following.
h = donothing(max(1e-6*a,1e-6)+2)-2;
%h = donothing(max(eps*a,1e-6)+2)-2;
logya   = gammaincln(x,a,tail);
logyaph = gammaincln(x,a+h,tail);
y = exp(logya);
dyda = y.*(logyaph-logya)./h;

k = find(x==0);
if any(k)
  dyda(k) = 0;
end

if nargout>2
  if all(a>h)
    logyamh = gammaincln(x,a-h,tail);
    d2yda = y.*(  (logyaph-2*logya+logyamh)./(h.^2) +(dyda./y).^2);
  else
    logyap2h = gammaincln(x,a+2*h,tail);
    d2yda = y.*(  (logya-2*logyaph+logyap2h)./(h.^2) +(dyda./y).^2);
  end
  if any(k)
    d2yda(k) = 0;
  end
end

function y=donothing(x)
  y=x;
return