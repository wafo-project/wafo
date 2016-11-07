function r= polyshft(p,a,b)
% POLYSHFT Polynomial coefficient shift
%
% CALL:  px = polyshft(py,a,b)
%
%   px,py = polynomial coefficients for the X and Y variable, respectively
%   a,b   = lower and upper limits 
%
% POLYSHFT shift the polynomial coefficients by a variable shift:
%
%   Y = 2*(X-.5*(b+a)/(b-a)
%
% i.e., the interval -1 <= Y <= 1 is mapped to the interval a <= X <= b
% 
% Example:
%  py = [1 0];
%  px = polyshft(py,0,5);
%  assert(polyval(px,[0 2.5 5]), [-1 0 1 ], eps);  % This is the same as the line below
%  assert(polyval(py,[-1 0 1 ]), [-1 0 1 ], eps);
% 
% See also: polyishft

% Reference 
% William H. Press, Saul Teukolsky, 
% William T. Wetterling and Brian P. Flannery (1997)
% "Numerical recipes in Fortran 77", Vol. 1, pp 184-194

% History
% by pab 2000

error(nargchk(3,3,nargin));
if (a==-1) && (b ==1),
  r = p;
  return
end
p = p(:).';
n = length(p);
%f = [ 1 -0.5*(a+b) ];
f = [2/(b-a) -(a+b)/(b-a)];
r = p(1);
for i = 1:n-1
  r = conv( r, f );      % Polynomial multiplication
  r(i+1) = r(i+1) + p(i+1);
end

%const = 2/(b-a);
%Then rescale by the factor const...
%r(1:n-1) = const.^(n-1:-1:1).*r(1:n-1);
return








