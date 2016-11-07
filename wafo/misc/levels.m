function l=levels(param)
%LEVELS Calculates discrete levels given the parameter matrix.
%
%  CALL: ui = levels(param);
%
%        ui    = the discrete levels.
%        param = the parameter matrix = [u_min u_max num_levels].
%
% Example:
%
% param = [1, 2, 3]; 
% assert(levels(param), [1, 1.5, 2])
%

l = linspace(param(1), param(2), param(3));

%!test assert(levels([1,2,3]), [1,1.5, 2])