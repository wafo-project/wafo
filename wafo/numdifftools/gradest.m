function [grad,err,finaldelta] = gradest(fun,x0,varargin)
% GRADEST Gradient vector of an analytical function of n variables
% 
% CALL: [grad,err,finaldelta] = gradest(fun,x0)
%
%  grad = first partial derivatives of fun evaluated at x0.    Size 1 x N
%  err  = error estimates corresponding to each value in grad. Size 1 x N
%  finaldelta = vector of final step sizes chosen for each partial derivative.
%  fun  = analytical function to differentiate. fun must
%        be a function of the vector or array x0.
%  x0   = vector location at which to differentiate fun
%        If x0 is an nxm array, then fun is assumed to be
%        a function of N = n*m variables. 
%
% GRADEST estimate first partial derivatives of fun evaluated at x0.
% GRADEST uses derivest to provide both derivative estimates
% and error estimates. fun needs not be vectorized.
%
% Example:
%  [grad,err] = gradest(@(x) sum(x.^2),[1 2 3]); 
%  assert(grad, [ 2,4, 6], 1e-12);
%  assert(err < 1e-12);
%
%  %At [x,y] = [1,1], compute the numerical gradient
%  %of the function sin(x-y) + y*exp(x)
%
%  z = @(xy) sin(-diff(xy)) + xy(2)*exp(xy(1));
%
%  [grad2,err2 ] = gradest(z,[1 1]);
%  assert(grad2, [3.71828182845911,   1.71828182845906], 1e-12);
%  assert(err2 < 1e-12);
%
%  %At the global minimizer (1,1) of the Rosenbrock function,
%  %compute the gradient. It should be essentially zero.
%
%  rosen = @(x) (1-x(1)).^2 + 105*(x(2)-x(1).^2).^2;
%  [grad3,err3] = gradest(rosen,[1 1]);
%  assert(grad3, [0, 0], 1e-12);
%  assert(err3<1e-12);
%
% See also derivest, gradient


% Author: John D'Errico
% e-mail: woodchips@rochester.rr.com
% Release: 1.0
% Release date: 2/9/2007
% 
% Revised pab 2007
% -updated help header to wafo style

% get the size of x0 so we can reshape
% later.
%sx = size(x0);


options.DerivativeOrder = 1;
options.MethodOrder = 2;
options.Style = 'central';
options.RombergTerms = 2;
options.FixedStep = [];
options.MaxStep = 100;
options.StepRatio = 2;
options.NominalStep = [];
options.Vectorized = 'no';

if nargin==1 && isequal(fun,'defaults')
  grad  = options;
  return
end
if nargin>2
  options = parseoptions(options,varargin{:});
end
% total number of derivatives we will need to take
nx = numel(x0);


grad = zeros(1,nx);
err = grad;
finaldelta = grad;
for ind = 1:nx
  [grad(ind),err(ind),finaldelta(ind)] = derivest( ...
    @(xi) fun(swapelement(x0,ind,xi)), ...
    x0(ind),options);
end

end % mainline function end

% =======================================
%      sub-functions
% =======================================
function vec = swapelement(vec,ind,val)
% swaps val as element ind, into the vector vec
vec(ind) = val;

end % sub-function end


