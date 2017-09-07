%README Readme file for module SIMTOOLS in WAFO Toolbox 
%
% 
% Version 2.5.1   18-Jun-2017
%
% This file describes new features, bug fixes, and changes in this version 
% of the WAFO Toolbox.
%
% FIXES
% ~~~~~
%  seamovie     - new dimension standards, option to save movie as avi
%
% NEW FUNCTIONS
% ~~~~~~~~~~~~~
%  spec2field   - spectral simulation of 3D Gaussian space-time field
%  spec2waves   - spectral simulation of 2D Gaussian space-time wave
%
%-----------------------------------------------------------------
% Version 2.1.1   07-Sep-2005
%
%  FIXES
%  ~~~~~
%  spec2sdat    - small fix
%  seamovie     - fixed a bug: in matlab7: "Mv =[]; Mv(j) = getframe;"
%                 does not work, now fixed.
%  seasim       - added fwaitbar and removed disp statements
%  spec2linspec - changed the default constant controlling its
%                 performance. Can be improved further.
%
%
%  NEW FUNCTIONS
%  ~~~~~~~~~~~~~
%  rfm2dtp      - Reconstructs a sequence of turning points from a rainflow matrix.
%                 Replaces the function rfc2load.
%
%----------------------------------------------------------------
% Version 2.1.0   07-Apr-2004
%
% This file describes new features, bug fixes, and changes in this version 
% of the WAFO Toolbox.
%
%
%  FIXES
%  ~~~~~
%  disufq1.c   - is renamed to disufq.c
%  spec2nlsdat - changed call from disifq1 to disufq
%  derivate.m  - is no longer needed and hence deleted 
%  smctpsim    - Corrections concerning initial state of regime process
%                and simulation of non-Markovian regime process.
%
%----------------------------------------------------------------
% Version 2.0.5   14-Jun-2003
%
% This file describes new features, bug fixes, and changes in this version 
% of the WAFO Toolbox.
%
%
%  FIXES
%  ~~~~~
%  smctpsim    - 'RFM' option didn't work cause x was not initiated. Now fixed!
%  
