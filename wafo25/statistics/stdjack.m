function [s,y] = stdjack(varargin)
%STDJACK  Jackknife estimate of the standard deviation of a parameter 
%	  estimate.
%
%	  s = stdjack(X,fun)
%	  
%	  The function is equal to sqrt(diag(covjack(X,'T')))

%       Anders Holtsberg, 28-02-95
%       Copyright (c) Anders Holtsberg


[C,y] = feval(@covjack,varargin{:});
s = sqrt(diag(C));