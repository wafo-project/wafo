%README New features, bug fixes, and changes in CYCLES.
%
% Version 2.1.1 06-Sep-2005
%
% Bug fixes, some new and updated routines     
% 
%  UPDATED FUNCTIONS
%  ~~~~~~~~~~~~~
%  CMATPLOT        Plots a cycle matrix, e.g. a rainflow matrix.
%  RFMEXTRAPOLATE  Extrapolates a rainflow matrix.
%  TP2ARFC         Calculates asymmetric rainflow cycles from turning points.
%  TP2RFC          Finds the rainflow cycles from the sequence of turning points.
%
%  NEW FUNCTIONS
%  ~~~~~~~~~~~~~
%  TPEXTRAPOLATE   Extrapolates a sequence of turning points.
%
%--------------------------------------------------------------------
%
% Version 2.1.0   07-Apr-2004
%
% Bug fixes and updated routines.
% 
%  UPDATED FUNCTIONS
%  ~~~~~~~~~~~~~
%  EXTRALC   Extrapolate level crossing spectrum
%  LSPLOT    Plot load spectra
%

%--------------------------------------------------------------------
%
% Version 2.0.5   14-Jun-2003
%
% Bug fixes, some new and/or updated routines.
% 
%  NEW FUNCTIONS
%  ~~~~~~~~~~~~~
%  DAT2RFM   Calculates the rainflow matrix from a time signal.
%  LSPLOT    Plot load spectra
%
%--------------------------------------------------------------------
% Version 2.0.4   11-Apr-2001
% 
% This file describes new features, bug fixes, and changes in this version 
% of the WAFO/cycles Toolbox.
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
%  FIXES
%  ~~~~~
%   Some minor fixes.
%
%
%  ENHANCEMENTS
%  ~~~~~~~~~~~~
%  - Revised routines for rainflow counting (tp2rfc, tp2arfc, ...). 
%    They now supports the Cloormann/Seeger method (same as the AFNOR 
%    recommendation and the ASTM standard) for rainflow counting, as 
%    well as the Rychlik definition. The only difference between the 
%    methods are the treatment of the so-called residual. The routines 
%    are more versatile and handles either one column or two column signals. 
%    It is also possible to extract rainflow cycles with time. 
%
%  - Rainflow counting with side information (tp2rfm_sid, tp2arfm_sid). 
%    Performs rainflow matrix counting of a load signal which is attached 
%    with a side information signal (e.g. steering angle). Hence, it is 
%    possible to split the rainflow matrix according to the side information. 
%    The result is a number of rainflow matrices, each representing different 
%    values of the side information. See [1]
%
%  - Extrapolation of rainflow matrices and level crossings.
%    The extreme RFM (see 'lc2rfmextreme') is an approximation of the RFM for 
%    large amplitudes, based on extreme values theory for extreme level 
%    crossings. The routines 'cmat2extralc' and 'extralc' extrapolates the 
%    level crossing spectrum for high and for low levels, by using the 
%    Generalized Pareto Distribution (GPD), which is motivated by extreme 
%    value theory. The routine 'rfmextrapolate' takes a measured rainflow 
%    matrix as input and constructs an extrapolation by using the extreme 
%    RFM for large amplitudes and soothing elsewhere. See [2]
%
%  References:
%
%  1.  P. Johannesson (1999):
%      Rainflow Analysis of Switching Markov Loads.
%      PhD thesis, Mathematical Statistics, Centre for Mathematical Sciences,
%      Lund Institute of Technology.
%
%  2.  Johannesson, P., and Thomas, J-.J. (2000): 
%      Extrapolation of Rainflow Matrices. 
%      Preprint 2000:82, Mathematical statistics, Chalmers, pp. 18. 
%
%  Download: <http://www.math.chalmers.se/~parj/Publications/>
%
%  NEW FUNCTIONS
%  ~~~~~~~~~~~~~
% Cycle counting (Rainflow, min-max cycles) & Crossings. [cycles]
%   tp2arfc4p   - Asymmetric RFC and residual from TP (used by tp2arfc).
%   res2arfc    - Asymmetric rainflow cycles for a residual.
%
% Discrete loads & Cycle matrices (Rainflow matrix). [cycles]
%   dtp2arfm4p  - Asymmetric RFM and residual from discrete TP (used by dtp2arfm).
%   dtp2rfm_sid - RFM from discrete turning points with side information.
%   dtp2arfm_sid- Asymmetric RFM from discrete TP with side information.
%   cmatresamp  - Resamples a cycle matrix, random resampling.
%   cmatcombine - Combines two cycle matrices.
%
% Extrapolation & Smoothing of RFM/CMAT/LC. [cycles]
%   rfmextrapolate - Extrapolates a rainflow matrix.
%   cmat2extralc   - Extrapolate level crossing spectrum from cycle matrix.
%   wgpdfit_mld    - Routine for ML estimation of GPD with discrete data.
%   lc2rfmextreme  - Compute extreme RFM from level crossings.
%   extralc        - Extrapolate a level crossing spectrum.

