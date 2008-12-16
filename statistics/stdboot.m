function [s,y] = stdboot(varargin)
%STDBOOT  Bootstrap estimate of the standard deviation of a parameter.
%
%	  s = stdboot(X,fun,B)
%	  
%         The function is equal to sqrt(diag(covboot(X,fun)))
%
%	  See also std, stdjack, covboot, and ciboot.

%       Anders Holtsberg, 02-03-95
%       Copyright (c) Anders Holtsberg



[C,y] = feval(@covboot,varargin{:});
s = sqrt(diag(C));
