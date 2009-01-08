function y = gammaincln(x,a,tail)
% GAMMAINCLN Logarithm of incomplete gamma function.
%
%  CALL y = gammaincln(x,a,tail)
%
%  tail = 'lower' (default) or 'upper'.
%
% Y = GAMMAINCLN(X,A) computes the natural logarithm of the incomplete gamma
%     function for each element of X and A.  GAMMAINCLN is defined as
%  
%         LOG(GAMMAINC(X,A))
%  
%     and is obtained from GAMMAINC(X,A) compensating for the roundoff when
%    GAMMAINC is close to one.  For large X, GAMMAINCLN is approximately 
%   -exp(-X), whereas log(GAMMAINC) can be zero.
%
% Example
%   x = linspace(0,500);
%   logg = log(gammainc(x,1));
%   logg2 = gammaincln(x,1);
%   semilogx(x,logg,x,logg2,'.')
%   
% See also gammainc, gammaln

if nargin<3||isempty(tail)
  tail = 'lower';
end
tails = {'lower','upper'};

switch lower(tail)
  case 'upper'
    %problems = (x<a);
    [tail2,tail1] = deal(tails{:});
  case 'lower'
    %problems = (x>a);
    [tail1,tail2] = deal(tails{:});
  otherwise
    error('unknown tail!')
end
y = log(gammainc(x,a,tail1));
problems = -0.1<y;
if any(problems(:))
  if ~isscalar(a), a = a(problems);end
  if ~isscalar(x), x = x(problems);end
  y(problems) = log1p(-gammainc(x,a,tail2));
end

