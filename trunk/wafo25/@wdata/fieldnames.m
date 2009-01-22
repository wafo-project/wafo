function out = fieldnames(this,varargin)
%FIELDNAMES  Return a cell array of WDATA property names.
%
%  CALL Names = fieldnames(wdata)
% 
%  returns a cell array of strings containing 
%   the names of the properties of WDATA.



%out = fieldnames(struct(this),varargin{:});
out = builtin('fieldnames',this,varargin{:});