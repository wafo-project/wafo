function a = ismatlab()
%ISMATLAB True if function run in an Matlab environment
%
% CALL a = ismatlab()
% 
%   a = 1 if Matlab
%       0 otherwise
%
% ISMATLAB is useful in situations where functions needs to take special
% for running in matlab environment.
%
% Example
%  assert(ismatlab, exist('matlabroot','builtin')~=0)
%
% See also isoctave

a = exist('matlabroot','builtin')~=0;