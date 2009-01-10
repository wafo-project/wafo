%README New features, bug fixes, and changes in TRGAUSS.
% Version 2.1.1   07-Sep-2005
%
% FIXES
% ~~~~~ 
%  pdfplot  - fixed some bugs
%  mctrtest - changed order of input
%  dat2tr   - changed the restriction on parameters used in the call to
%            hermitetr
%  rind     - removed unused code + compiled rindalan24 again.
% 
%
% NEW FUNCTIONS
% ~~~~~~~~~~~~~
% chitwo2lc_sorm  - SORM-approximation of the crossing intensity for the 
%                   noncentral Chi^2 process
% chitwo2lc_sp    - Saddlepoint approximation of the crossing intensity for the 
%                   noncentral Chi^2 process 
%
%--------------------------------------------------------------------
% Version 2.1.0   07-Apr-2004
%
% FIXES
% ~~~~~ 
%  All the spec2XXX programmes are modified by removing tic, toc
%  statements and replace some code with call to spec2cov2.
%  spec2mmtpdf, spec2thpdf, rind : are totally rewritten to allow rind-options
%                           structure controlling the accuracy of integration.
%  trangood   - replaced call to sort with unique  
%  initdata   - is no longer needed and hence deleted
%  trraylpdf  - is located in module WAVEMODELS and hence deleted from
%               this directory
%  spec2vhpdf - is no longer needed and hence deleted
%  ochifun    - is moved to private directory
%  hermitefun - is moved to private directory
%  specq2lc - eww, and kwaves added as subfunctions, moved it from spec to trgauss
%  eww      - is no longer needed and hence is deleted
%  kwaves   - is no longer needed and hence is deleted
%  dat2tr, lc2tr, lc2tr2, cdf2tr, mctrtest, troptset, trplot, lomaxcdf - is moved
%             to trgauss from onedim.
%
% New routines
%  rindoptset  - Create or alter RIND OPTIONS structure.
%  initoptions - Initializes RIND options according to speed.
%  th2vhpdf    - Transform joint T-H density to V-H density (replaced
%                spec2vhpdf)
%  bvnormcdf   - Bivariate Normal cumulative distribution function
%  bvnormprb   - Bivariate normal probability  
%  mvnormpcprb - Multivariate Normal probabilities with product correlation
%  mvnormprb   - Computes multivariate normal probability by Genz' algorithm.
%  mvnortpcprb - Multivariate normal or student T probability
%                with product correlation structure.
% 
%--------------------------------------------------------------------
% Version 2.0.5   14-Jun-2003
%
% Bug fixes, some new and/or updated routines.
% 
%--------------------------------------------------------------------
% Version 2.0.4   11-Apr-2001
% 
% This file describes new features, bug fixes, and changes in this version 
% of the WAFO/trgauss Toolbox.
%
% List of changes:
% 
% NOTE:
% Items marked with:
%       (N)  : needs a new  name 
%        *   : not yet available
%        +   : needs more work
%        **  : have changed in a way which might affect your code!
%
% 
