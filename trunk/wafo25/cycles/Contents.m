% Module CYCLES in WAFO Toolbox. 
% Version 2.1.1   06-Sep-2005 
%
%
% This module contains routines for cycle counting, discretization, 
% and crossings.
% 
% The routines are described in the contents of 'Fatigue in WAFO'.
% Try 'help fatigue'.
%
%   Readme           - Readme file for module CYCLES in WAFO Toolbox.
%
%   cc2amp           - Calculates the amplitudes from a cycle count.
%   cc2cmat          - Calculates the cycle count matrix from a cycle count.
%   cc2dcc           - Discretize a cycle count.
%   cc2lc            - Calculates the number of upcrossings from a cycle count  
%   ccplot           - Plots a cycle count as a point process in the plane.
%   cmat2amp         - Calculates a histogram of amplitudes from a cycle matrix.
%   cmat2extralc     - Extrapolate level crossing spectrum
%   cmat2lc          - Calculates the level crossings from a cycle matrix.
%   cmat2nt          - Calculates a counting distribution from a cycle matrix.
%   cmat2rmcmat      - Converts a cycle matrix from Range-Mean format to min-max format.
%   cmatcombine      - Combines two cycle matrices.
%   cmatplot         - Plots a cycle matrix, e.g. a rainflow matrix.
%   cmatresamp       - Resamples a cycle matrix.
%   cocc             - Plots cycles as points together with isolines of a cycle matrix.
%   dat2dtp          - The sequence of discretized turning points from a signal.
%   dat2rfm          - Calculates the rainflow matrix from a time signal.
%   dcc2cmat         - Calculates the cycle matrix for a discrete cycle count.
%   dtp2arfm         - Calculates asymmetric RFM from discrete turning points.
%   dtp2arfm4p       - Calculates asymmetric RFM from discrete turning points (4-point).
%   dtp2arfm_sid     - Asymmetric RFM from discrete TP with side information.
%   dtp2rfm          - Calculates rainflow matrix from discrete turning points.
%   dtp2rfm_sid      - Rainflow matrix from discrete turning points with side information.
%   extralc          - Extrapolate level crossing spectrum
%   lc2rfmextreme    - Compute extreme RFM from level crossings.
%   lsplot           - Plot load spectra.
%   nt2cmat          - Calculates a cycle matrix from a counting distribution.
%   nt2lc            - Calculates the level crossings from a cycle matrix.
%   res2arfc         - Calculates asymmetric rainflow cycles for a residual.
%   rfcfilter        - Rainflow filter a signal.
%   rfmextrapolate   - Extrapolates a rainflow matrix.
%   rmcmat2cmat      - Converts a cycle matrix from Range-Mean format to min-max format.
%   smoothcmat       - Smooth a cycle matrix using (adaptive) kernel smoothing
%   smoothcmat_hnorm - Bandwidth selection for kernel smoothing of a cycle matrix. 
%   tp2arfc          - Calculates asymmetric rainflow cycles from turning points.
%   tp2arfc4p        - Calculates asymmetric rainflow cycles from turning points (4-point).
%   tp2lc            - Calculates the number of upcrossings from the turning points.
%   tp2mm            - Calculates min2Max and Max2min cycles from a sequence of turning points
%   tp2rfc           - Finds the rainflow cycles from the sequence of turning points.
%   tpextrapolate    - Extrapolates a sequence of turning points.
%   fitgenpar_mld.m  - Parameter estimates for GPD data
