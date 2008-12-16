% Fatigue in WAFO (Wave Analysis for Fatigue and Oceanography).
% Version 2.1.1   28-Sep-2005 
%
% Routines for Markov models, cycle counting, damage, and fatigue.
%
%
% Cycle counting (Rainflow, min-max cycles) & Crossings. [WAFO/cycles]
%   dat2tp      - The turning points from data, optionally rainflowfiltered. [onedim]
%   dat2rfm     - Calculates the rainflow matrix from a time signal.
%   rfcfilter   - Rainflow filter turning points.
%   tp2lc       - Level crossings from turning points.
%   tp2mm       - min2Max (and Max2min) cycles from turning points.
%   tp2rfc      - Rainflow cycles from turning points.
%   tp2arfc     - Asymmetric rainflow cycles from TP.
%   tp2arfc4p   - Asymmetric RFC and residual from TP (used by tp2arfc).
%   res2arfc    - Asymmetric rainflow cycles for a residual.
%   cc2lc       - Level crossings from a cycle count.
%   cc2amp      - Amplitude histogram from a cycle count.
%   findrfc     - Finds indices to rainflow cycles of a sequence of TP. [onedim]
%   findcross   - Finds the indices of level v up- and down- crossings of a vector. [onedim]
%
% Discrete loads & Cycle matrices (Rainflow matrix). [WAFO/cycles]
%   dat2dtp     - Discrete TP from load.
%   cc2dcc      - Discretization of a cycle count.
%   dtp2rfm     - Rainflow matrix from discrete TP.
%   dtp2arfm    - Asymmetric rainflow matrix from discrete TP.
%   dtp2arfm4p  - Asymmetric RFM and residual from discrete TP (used by dtp2arfm).
%   dtp2rfm_sid - RFM from discrete turning points with side information.
%   dtp2arfm_sid- Asymmetric RFM from discrete TP with side information.
%   cc2cmat     - Estimates cycle matrix from a cycle count. 
%   dcc2cmat    - Histogram matrix from discrete class indices.
%   cmat2nt     - Counting distribution from cycle matrix.
%   nt2cmat     - Cycle matrix from counting distribution.
%   cmat2lc     - Level crossings from cycle matrix.
%   nt2lc       - Level crossings from counting distribution.
%   cmat2amp    - Amplitude histogram from cycle matrix.
%   cmat2rmcmat - Convert cycle matrix from min-max to range-mean.
%   rmcmat2cmat - Convert cycle matrix from range-mean to min-max.
%   cmatresamp  - Resamples a cycle matrix, random resampling.
%   cmatcombine - Combines two cycle matrices.
%
% Extrapolation & Smoothing of RFM/CMAT/LC. [WAFO/cycles]
%   rfmextrapolate - Extrapolates a rainflow matrix.
%   tpextrapolate  - Extrapolates a sequence of turning points.
%   cmat2extralc   - Extrapolate level crossing spectrum from cycle matrix.
%   fitgenpar_mld  - Routine for ML estimation of GPD with discrete data.
%   lc2rfmextreme  - Compute extreme RFM from level crossings.
%   extralc        - Extrapolate a level crossing spectrum.
%   smoothcmat     - Smooth a cycle matrix.
%
% Plotting. [WAFO/cycles]
%   plotcc      - Plot cycle count.
%   cmatplot    - Plot cycle matrix.
%   cocc        - Plot a cycle count with cycle matrix isolines.
%   plotlc      - Plot level crossing intensity/spectrum. [onedim]
%   lsplot      - Plot load spectra.
%
% Markov model. [WAFO/markov]
%   mat2tmat    - Convert matrix to transition matrix.
%   mc2stat     - Stationary distribution for a Markov chain.
%   mc2rfm      - RFM for a Markov chain.
%   mctp2stat   - Stationary distribution for a MCTP.
%   mctp2reverse- Calculate time-reversed MCTP.
%   mctp2rfm    - RFM for a Markov chain of turning points.
%   mctp2arfm   - Asymmetric rainflow matrix for a MCTP.
%   arfm2mctp   - Markov matrix given an asymmetric rainflow matrix. 
%   mktestmat   - Makes test matrices for min-Max (and Max-min) matrices.
%
% Switching Markov model. [WAFO/markov]
%   smc2rfm     - RFM for a Switching Markov chain.
%   smctp2stat  - Stationary distribution for a SMCTP.
%   smctp2joint - Calculate joint Markov chain for a SMCTP.
%   smctp2rfm   - RFM for a Switching Markov chain of turning points.
%   smctp2arfm  - Asymmetric rainflow matrix for a SMCTP.
%   plothmm     - Plot Hidden Markov Model (HMM).
%
% Estimate a Switching Markov model. [WAFO/markov]
%   estsmctp    - Estimate SMCTP model from an observed RFM.
%   f_smctp     - Help routine to estsmtp.
%   estmc       - Estimate Transition matrix of MC from signal.
% 
%   loglcmat    - log-Likelihood of cycle matrix.
%   chi2cmat    - Chi-square distance of cycle matrix.
%   hdcmat      - Hellinger distance of cycle matrix.
%   klcmat      - Kullback-Leibler distance of cycle matrix.
%
%   tr_p2x      - Transform P-matrix to X-vector (Used by 'estsmctp')
%   tr_x2p      - Transform X-vector to P-matrix. (Used by 'estsmctp')
%   tr_m2x      - Transform Model-structure to X-vector. (Used by 'estsmctp')
%   tr_x2m      - Transform X-vector to Model-structure. (Used by 'estsmctp')
%   scalemat    - Scale and translate a cycle matrix. (Used by 'estsmctp')
%   f_funm      - Calculate min-max matrix for Model-structure. (Used by 'estsmctp')
%
% Fatigue & damage. [WAFO/damage]
%   cc2dam      - Calculate total damage from a cycle count.
%   cmat2dam    - Calculate total damage from a cycle matrix.
%   cmat2dmat   - Calculate damage matrix from cycle matrix.
%   lc2dplus    - Upper bound for total damage from level crossings.
%   plotsn      - Plots S-N data and estimates parameters.
%   sphdam      - Calculates spherical damage for a 3-D load.
%   ftf         - Calculates fatigue failure time distribution.
%   damint      - Calculates damage intensity from counting distribution.
%   down2cc     - Calculates the most damaging cycle count given crossings.
%   roadspec    - Road spectrum.
%
% Simulation of random loads [WAFO/wsim]
%   lc2sdat     - Simulates a process with given irregularity factor and crossing spectrum.
%   rfm2dtp     - Reconstructs a sequence of turning points from a rainflow matrix. 
%   mcsim       - Simulates a discrete Markov chain
%   mctpsim     - Simulates a discrete Markov chain of turning points
%   smcsim      - Simulates a switching Markov chain
%   smctpsim    - Simulates a switching Markov chain of turning points
%   sarmasim    - Simulates a switching AR- or ARMA-process
%
% Demos. [WAFO/wdemos]
%   rfcdemo1       - Demo for switching AR(1)-processes.
%   rfcdemo2       - Demo for calculation of the RFM for a switching MCTP model.
%   itmkurs        - Initiate paths for Demo Load and Fatigue Analysis. 
%   democc         - Demo of rainflow and min-max cycle definitions.
%   democc_markmax - plots load and marks a maximum. (Used by democc) 
%   democc_rfcdef  - illustrates the definition of rainflow cycles. (Used by democc) 
%   democc_mmdef   - illustrates the definition of min-max cycles. (Used by democc) 
%   democc_tpdef   - illustrates the definition of turning points. (Used by democc) 
%   democc_plotmat - plots RFC and min-max counts. (Used by democc) 

