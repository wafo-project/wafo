function [r,xmax,xmin] = range(x)
%RANGE  Computes the range between the maximum and minimum values. 
% 
% CALL: [r,xmax,xmin] = range(x) 
%
%   r =  max(x) - min(x), i.e. a scalar or vector giving the range of x
%   x = vector or matrix
%
% Example:
%   x=rndexp(5,1,10)
%   [r,xmax,xmin] = range(x);
%
% See also  max, min, iqrange

% By pab 06.02.2000
xmin = min(x);
xmax = max(x);
r = xmax - xmin;
