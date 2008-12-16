%README New features, bug fixes, and changes in ONEDIM.
%
%  Version 2.1.1   07-Nov-2004
%
%   NEW FUNCTIONS
%   ~~~~~~~~~~~~~
%   findextrema - Finds indices to minima and maxima of data
%
%------------------------------------------------------------------
% Version 2.1.0   07-Apr-2004
%
% FIXES
% ~~~~~ 
%  dat2tr, lc2tr, lc2tr2, cdf2tr, mctrtest, troptset, trplot, lomaxcdf - is moved
%             to trgauss from onedim.
%
%-------------------------------------------------------------------
% Version 2.0.5   14-Jun-2003
%
%   FIXES
%   ~~~~~
%
% dat2steep - removed disp statement and replaced with call to waitbar. 
% findrfc   - fixed a bug in the help header  
% troptset  - moved some code into parseoptions for easier maintainence  
%
% 
%--------------------------------------------------------------------
% Version 2.0.4   11-Apr-2001
% 
% This file describes new features, bug fixes, and changes in this version 
% of the WAFO/onedim Toolbox.
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
% - (**) All the options governing the estimation of the transformations are
% now moved into a TRANSFORM OPTIONS structure for easier transfering the
% options to subfunctions. The following functions are affected by this:
%    DAT2TR, LC2TR, MCTRTEST, RECONSTRUCT
%
% - MM2LC is vectorized even further  speed up things
% 
%   NEW FUNCTIONS
%   ~~~~~~~~~~~~~
%
%  LC2TR2   - Estimate transformation from observed crossing intensity, alternativ inversion.
%  CDF2TR   - Estimate transformation from observed Cumulative Distributions Function
%  TROPTSET - Create or alter TRANSFORM OPTIONS structure.
 
