function a = isoctave()
%ISOCTAVE True if function run in an Octave environment
%
% CALL a = isoctave()
% 
%   a = 1 if Octave
%       0 otherwise
%
% ISOCTAVE is useful in situations where functions needs to take special
% for running in octave environment.
%
% Example
%  isoctave
%
% See also ismatlab
a = exist('octave_core_file_name','builtin')~=0;