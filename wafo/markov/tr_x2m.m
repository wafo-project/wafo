function M = tr_x2m(X,known)
% TR_X2M Transform X-vector to Model-structure.
% 
% CALL: M = tr_x2m(X,known)
%
% See also tr_m2x

M.x0  = X(1:2)';
M.s   = exp(X(3));
M.lam = exp(X(4));
return
