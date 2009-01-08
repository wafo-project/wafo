function [X] = tr_m2x(M)
% TR_M2X Transform Model-structure to X-vector.
%
% CALL: X = tr_m2x(M)
%
% See also tr_x2m,  estsmctp

X = [M.x0'; log(M.s); log(M.lam)];
return