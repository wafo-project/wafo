%README Readme file for module SPEC in WAFO Toolbox 
% Version 2.5.3   24-May-2017
%
% Several updates since 2004
%   to make routines work with Matlab 9.1 and older
%
% Version 2.1.1   14-nov-2004
%
%  FIXES
%  ~~~~~
%  spec2cov2 - rate is sometimes zero -> spec2cov2 crash: made an ad-hoc
%              solution.
%  getjonswappeakedness - Made sure gamma is 1 if Tp/sqrt(Hm0) > 5.1429
%  rotspec   - moved from private directory to wafo/spec
%  spec2spec - Moved code to rotspec + small cosmetics fixes   
%              implemented freq -> enc.
%  dir2enc   - BUG's removals plus major changes.
%  mccormick - replaced fmin with fminbnd   
%
%-----------------------------------------------------------
% Version 2.1.0   14-sep-2004
%
% This file describes new features, bug fixes, and changes in this version 
% of the WAFO Toolbox.
%
% List of changes:
%
%  FIXES
%  ~~~~~
%  rotspec -  -Removed old unused code and added example and rotateGrid to input   
%  spec2spec - Moved code to rotspec + small cosmetics fixes + some fixes
%              in the private functions dir2enc, spa2time, time2spa. 
%  dir2k1d   - removed since it is no longer needed
%
% Version 2.0.6   14-Jan-2004
%
% This file describes new features, bug fixes, and changes in this version 
% of the WAFO Toolbox.
%
% List of changes:
%
%  FIXES
%  ~~~~~
%  specq2lc - eww, and kwaves added as subfunctions, moved it to trgauss
%  eww      - is no longer needed in spec and hence is deleted
%  kwaves   - is no longer needed in spec and hence is deleted
%  cov2spec - Increased the precision in calculations and deleted L from input, added fast instead
%  spec2cov  - streamlined some code, updated help header, fixed bug in example  
%  plotspec - replaced call to colorbarf with fcolorbar and 
%              replaced call to findzlevel with clevels  
% 
%  NEW FUNCTIONS
%  ~~~~~~~~~~~~~
%  spec2cov2 - based on code from spec2XXpdf programmes
% getjonswappeakedness - Peakedness factor Gamma given Hm0 and Tp for JONSWAP
