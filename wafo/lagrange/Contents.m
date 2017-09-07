% Module LAGRANGE in WAFO toolbox
% Version 1.2.2  30-05-2017
%
%
% LAGRANGE contains routines for generation and analysis of random 
%   1st and 2nd order Lagrange waves
%
%   WafoL_Tutorial     - Tutorial in pdf-format for Lagrange module
%
%   dat2crossind       - Finds indices to level v down and/or upcrossings from data
%   disper2            - Dispersion relation with possible mean flow
%   genoptset          - Creates or alters 3D generation options structure
%   ldat2lslope        - Extracts slopes at level crossings in Lagrange model
%   ldat2lwav          - Finds time/space Lagrange process from simulated components
%   ldat2lwav3D        - Generates Lagrange 3D wave process from simulated components
%   looptest           - Simulates 2D Lagrange waves to estimate folding rate
%   lwav2frontback     - Gives front/back crest periods/wavelength of wave data
%   pdfnorm2d          - Bivariate Gaussian distribution  
%   seamovieL          - Makes a movie of a 2D or 3D simulated sea structure 
%   simoptset          - Creates or alters simulation options structure
%   spec2lasym         - Simulates asymmetry measures for Lagrange waves from spectrum
%   spec2lcov          - Calculates auto- and cross-covariance functions 
%   spec2ldat          - Simulates w and x components of 2D Lagrange wave
%   spec2ldat3D        - Spectral simulation of components in 3D Lagrangian sea 
%   spec2ldat3DM       - Particle trajectory simulation according to Marc Prevosto
%   spec2ldat3DP       - Parallel version of spec2ldat3DM for trajectory simulation
%   spec2lseries       - Spectral simulation of time series in 3D Lagrangian sea 
%   spec2slcomp        - Compares 2nd order Stokes and 1st order Lagrange time waves 
%   spec2slopedat      - Simulates Lagrange waves and extracts slopes at crossings 
%   spec2slopedat3D    - Simulates values and slopes in 3D Lagrange field 
%   spec2spaceslopecdf - Computes cdf for slope at crossings of space waves 
%   spec2spaceslopepdf - Computes pdf for slope at crossings of space waves 
%   spec2steepdat      - Simulates Lagrange waves and extracts steepness and slopes
%   spec2timeslopecdf  - Computes cdf for slopes at crossings of time waves 
%   spec2timeslopepdf  - Computes pdf for slopes at crossings of time waves 
%   WafoLCh1           - Script with commands for WafoL tutorial, Chapter 1
%   WafoLCh2           - Script with commands for WafoL tutorial, Chapter 2
%   WafoLCh3           - Script with commands for WafoL tutorial, Chapter 3
%   wav2slope          - Extracts slopes at up- and downcrossings after smoothing
