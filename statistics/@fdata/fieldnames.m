function out = fieldnames(this,varargin)
%FIELDNAMES  Return a cell array of property names.
%
%  CALL Names = fieldnames(thisClass)
% 
%  returns a cell array of strings containing 
%   the names of the properties of thisClass.



out = fieldnames(struct(this),varargin{:});
