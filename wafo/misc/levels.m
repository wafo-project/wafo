function l=levels(param)
%LEVELS Calculates discrete levels given the parameter matrix.
%
%  CALL: ui = levels(param);
%
%        ui    = the discrete levels.
%        param = the parameter matrix = [u_min u_max num_levels].
% Example:% 
%
% param=[2 3.2 4]; levels(param)

l=linspace(param(1),param(2),param(3));
