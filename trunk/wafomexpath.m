function p=wafomexpath
%WAFOMEXPATH Returns the path to mex-executables for the WAFO Toolbox
%
% CALL:  Str = wafomexpath 
% 
%   Str  =  a string with the path to the directory containing 
%           the mex-executables for the WAFO Toolbox,
%           including a terminating slash or backslash, depending
%           on the system on wich Matlab is run.
%           The mex files must be located in a subdirectory to
%           the WAFOEXEPATH. The subdirectory name is given by
%           "version('-release')".
% Example
% wafomexpath
% 
% See also  wafoexepath, waforoot, wafopath.



% Tested on: Matlab 5.3
% History:
% by pab 2007

release = version('-release');

p = [ fullfile(wafoexepath, release), filesep]; 
