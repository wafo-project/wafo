%README Readme file for module MISC in WAFO Toolbox 
% 
% Version 2.5.1 07-March-2011
% 
%  RENAME
%  ~~~~~~
%  Function  findpeaks.m  renamed to  wfindpeaks.m  
%  Function  smooth.m     renamed to  cssmooth.n  
%       to avoid name conflict with MATLAB toolboxes
%
% Version 2007a
% 
% csort removed because the built in sort is faster.
% 
% Version 2.1.1   07-Nov-2004
%
%  FIXES
%  ~~~~~
%  hypgf    - added dea3 for more stable error estimation
%  gaussq   - added the possibility of using a function handle.
%  gaussq2d - added the possibility of using a function handle.
%  ccquad   - added the possibility of using a function handle.
%
%  NEW FUNCTIONS
%  ~~~~~~~~~~~~~
%  choices.m added because matlab 7 does not provide this function anymore.
%  
%--------------------------------------------------------------
% Version 2.1.0   07-Apr-2004
% 
% No changes in this module compared to previous version. 
%
% Version 2.0.5   11-Nov-2003
%
% This file describes new features, bug fixes, and changes in this version 
% of the WAFO Toolbox.
%
%
%  FIXES
%  ~~~~~
%  colorbarf - is no longer needed and hence deleted 
%  mkcont    - is no longer needed and hence deleted  
%  h1line.m  - is no longer needed and hence deleted 
%  distchk.m - is no longer needed and hence deleted 
%  gaussq    - added the possibility of using a inline function.
%  gaussq2d  - replaced call to distchk with comnsize 
%  convlv    - fixed a bug when nr=nl=0 and when nh~=1
%  
%  NEW FUNCTIONS
%  ~~~~~~~~~~~~~
%
%  parseoptions - Create or alter a OPTIONS structure.
%  geth1line    - replaces h1line
%  mkcontents   - replaces mkcont
%  fcolorbar    - Display color bar with discrete color axis for filled contour plot
%  fwaitbar     - Fast display of wait bar.
%  genchol      - Generalized Cholesky factorization
