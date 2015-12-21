%README New features, bug fixes, and changes in KDETOOLS.
% Version 2.5.1   31-Dec-2010
%  See wafo code
%
%----------------------------
% Version 2.1.1   15-Sep-2005
%
% Bug fixes, some new and/or updated routines.
%
%  FIXES
% ~~~~~~~~
%  deriv     - made it general for derivatives of any order
%  qlevels   - updated calls from normpdf to wnormpdf in the help-header
%  gridcount - fixed a bug for d>2;  
%  hldpi, hldpi2 - added nargchk + minor fix
%  kdefft - is removed because it is obsolete (i.e. it is a special case
%           of kdebin). (removed all dependence on this function)
%  kdebin, kde, kdefun - updated to accept KDE options structure (also
%                       updated all functions dependent on these).
%
%  NEW FUNCTIONS
%   ~~~~~~~~~~~~~
%  kde1dgui  - Graphical user interface to kde function.
%  kde2dgui  - Graphical user intergace to kde function in two dimensions.
%  kdeoptset - Create or alter a KDE options structure
%  kdedemo1  - Demonstrate the smoothing parameter impact on KDE
%  kdedemo2  - Demonstrate the difference between transformation- and ordinary-KDE
%
%--------------------------------------------------------------------
% Version 2.1.0   7-Apr-2004
%
% Bug fixes, some new and/or updated routines.
%
%   FIXES
% ~~~~~~~~
%  Replaced (+)-signs with todo comments
%  Replaced the binning code in some functions with a call to a separate
%  function, gridcount.
%
%   hbcv2     - fixed some bugs 
%   hldpi2fft - fixed some bugs, added call to gridcount   
%   hste,hstt - added inc, maxit,abseps,releps as inputs  
%   ssample   - changed ind generation to avoid dependence on stats-toolbox 
%   hns       - fixed a bug when D>1   
%   hboot     - fixed a bug in default value for hvec 
%   hldpi2    - added missing psim function
%
%   NEW FUNCTIONS
%   ~~~~~~~~~~~~~
%   gridcount   - Histogram using linear binning.
%   kernelstats - Return 2'nd order moment of kernel pdf
%   ffndgrid    - Fast 'n' Furious N-D data gridding
%
%--------------------------------------------------------------------
% Version 2.0.5   3-Jun-2003
%
% Bug fixes, some new and/or updated routines.
%
%--------------------------------------------------------------------
% Version 2.0.4   11-Apr-2001
%
% This file describes new features, bug fixes, and changes in this version 
% of the WAFO/kdetools Toolbox.
% 
%  List of changes:
%  
%  NOTE:
%  Items marked with:
%        (N)  : needs a new  name 
%         *   : not yet available
%         +   : needs more work
%         **  : have changed in a way which might affect your code!
% 
%  
%   FIXES
%   ~~~~~
%    Some but too small to mention.
% 
% 
%   ENHANCEMENTS
%   ~~~~~~~~~~~~
%
% - Added the possibility that the parameter L2 is a cellarray of parametric
%   or non-parametric transformations. This means that a empirical
%   transformation estimated by CDF2TR may be used to transform the data.
%   The following functions are affected by this change:
%     KDEBIN, KDE, KDEFUN.
%
%   NEW FUNCTIONS
%   ~~~~~~~~~~~~~
%
%   - None
