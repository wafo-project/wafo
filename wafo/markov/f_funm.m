function [F,Fh] = f_funm(param,M)
%F_FUNM  Calculate min-max matrix for Model-structure. 
%
% CALL:  [F,Fh] = f_funm(param,M);
%
% Auxiliary function used by ESTSMCTP.

[F,Fh] = mktestmat(param,M.x0,M.s,M.lam);

