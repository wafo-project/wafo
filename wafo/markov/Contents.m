% Module MARKOV in WAFO Toolbox. 
% Version 2.5.3   30-05-2017 
% 
%
% This module contains routines for Markov loads, 
% switching Markov loads, and their connection to rainflow cycles.
% 
% The routines are described in the contents of 'Fatigue in WAFO'.
% Try 'help fatigue'.
%
%   Readme       - Readme file for module MARKOV in WAFO Toolbox.
%
%   arfm2mctp    - Calculates the markov matrix given an asymmetric rainflow matrix. 
%   chi2cmat     - Chi-square distance of cycle matrix.
%   estmc        - estimates transition matrix P from a time series of a Markov chain.
%   estsmctp     - Estimate SMCTP model from an observed rainflow matrix.
%   f_funm       - Calculate min-max matrix for Model-structure. 
%   f_smctp      - Auxiliary function used by ESTSMCTP
%   hdcmat       - Hellinger distance of cycle matrix.
%   hmmplot      - plots a Hidden Markov Model.
%   klcmat       - Kullback-Leibler distance of cycle matrix.
%   loglcmat     - log-Likelihood of cycle matrix.
%   mat2tmat     - Converts a matrix to a transition matrix.
%   mc2rfm       - Calculates the rainflow matrix/intensity for a Markov chain.
%   mc2stat      - Calculates the stationary distribution for a Markov chain.
%   mctp2arfm    - Calculates asymmetric rainflow matrix for a MCTP.
%   mctp2reverse - Calculates the time-reversed MCTP for a SMCTP.
%   mctp2rfm     - Calculates the rainflow matrix for a MCTP.
%   mctp2stat    - Calculates the stationary distribution for a MCTP.
%   mktestmat    - Makes test matrices for min-max (and max-min) matrices.
%   plothmm      - plots a Hidden Markov Model.
%   scalemat     - Scale and translate a cycle matrix.
%   smc2rfm      - Calculates rainflow matrix/intensity for a switching Markov chain.
%   smctp2arfm   - Calculates the asymmetric rainflow matrix for a SMCTP.
%   smctp2joint  - Calculates the joint MCTP for a SMCTP.
%   smctp2rfm    - Calculates the rainflow matrix for a SMCTP.
%   smctp2stat   - Stationary distributions for a switching MCTP.
%   tr_m2x       - Transform Model-structure to X-vector.
%   tr_p2x       - Transform P-matrix to X-vector
%   tr_x2m       - Transform X-vector to Model-structure.
%   tr_x2p       - Transforms a vector X to a transition matrix P.
