% STARTUP example startup file to put in your matlab path
%
% You should install the WAFO toolbox according to the description given
% in the install.m file. Then modify this file so that the root directory of
% WAFO defined in this file complies with the installation on your machine. 
%   
% See also  install

% add the WAFO root directory to the search path
% modify the line below if this does not comply with 
% your machine

addpath(fullfile(matlabroot,'toolbox','contrib','wafo','')) 

% adds the subdirectories of WAFO to the search path
switch 3
  case 1, wafopath
  case 2, initwafo;             % Initiate WAFO toolbox paths
  case 3, initwafo('full',1,0); % Full init., put WAFO paths last in search path
end
