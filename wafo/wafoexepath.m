function p=wafoexepath
%WAFOEXEPATH Returns the path to executables for the WAFO Toolbox
%
% CALL:  Str = wafoexepath 
% 
%   Str  =  a string with the path to the directory containing 
%           the executables for the WAFO Toolbox,
%           including a terminating slash or backslash, depending
%           on the system on wich Matlab is run.
%           The executable files must be located in a subdirectory to
%           the WAFO/exec directory, with the name given by
%           "lower(computer)".
% Example
% wafoexepath
%
% See also  wafomexpath, waforoot, wafopath.


%          NOTE: User no longer has to edit this file 
%          It is  automatic.

% Tested on: Matlab 5.3
% History:
% by by pab 11.08.99
%   new routine 

p =  [fullfile(waforoot, 'exec', lower(computer)) filesep]; 
