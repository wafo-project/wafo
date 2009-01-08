% WAFO Toolbox
% Version 2.2.1   07-Aug-2008 
% 
%   Wave Analysis for Fatigue and Oceanography (WAFO) Toolbox
% 
% The WAFO toolbox is mainly developed at 
%   Mathematical statistics, Lund University, Sweden
%   <http://www.maths.lth.se/matstat/>
%
% WAFO consists of several modules (directories) with short descriptions below.
% The modules SPEC, TRGAUSS, WAVEMODELS, and MULTIDIM are mainly for 
% oceanographic applications.
% The modules CYCLES, MARKOV, and DAMAGE are mainly for fatigue problems. 
% The contents file for each module is shown by typing 'help module-name' 
% Type 'help fatigue' for a presentation of all routines related to fatigue.
%
% The paths to the modules are initiated by the function 'initwafo'.
%
% ONEDIM     - Data analysis of time series. Example: extraction of 
%              turning points, estimation of spectrum, covariance function.
% SPEC       - Computation of spectral moments and covariance functions. 
%              Ex: common spectra implemented, directional spectra, 
%              bandwidth measures
% TRGAUSS    - Modelling with linear, Gaussian waves. Ex: exact 
%              distributions for wave characteristics, transformed
%              Gaussian processes
% WAVEMODELS - Models for distributions of wave characteristics found in 
%              the literature. Ex: parametric models for breaking 
%              limited wave heights.
% MULTIDIM   - Multi-dimensional time series analysis.  (Under construction)
% CYCLES     - Cycle counting, discretization, and crossings. Ex: Rainflow 
%              cycles and matrix, discretization of loads.
% MARKOV     - Routines for Markov loads, switching Markov loads, and 
%              their connection to rainflow cycles.
% DAMAGE     - Calculation of damage. Ex: Damage of a rainflow count or 
%              matrix, damage matrix, S-N plot.
% SIMTOOLS   - Simulation of random processes. Ex: spectral simulation, 
%              simulation of discrete Markov chains, switching Markov
%              chains, harmonic oscillator
% STATISTICS - Statistical tools and extreme-value distributions.
%              Ex: generation of random numbers, estimation of parameters,
%              evaluation of pdf and cdf
% KDETOOLS   - Kernel-density estimation.
% MISC       - Miscellaneous routines. Ex: numerical integration, smoothing 
%              spline, binomial coefficient, water density.
% WDEMOS     - WAFO demos. 
% DOCS       - Documentation of toolbox, definitions. An overview is given 
%              in the routine wafomenu. 
% DATA       - Measurements from marine applications.
% PAPERS     - Commands that generate figures in selected scientific 
%              publications.
% SOURCE     - Fortran and C files. Information on compilation.
% EXEC       - Executable files (cf. SOURCE), pre-compiled for Solaris, 
%              Alpha-Dec or Windows. 
% 
% WAFO homepage: <http://www.maths.lth.se/matstat/wafo/>
% On the WAFO home page you will find:
%  - The WAFO Tutorial
%  - New versions of WAFO to download.
%  - Reported bugs.
%  - List of publications related to WAFO.

% By jr 00.05.26
% Revised 00.07.10
% Revised 00.07.13
% Revised by PJ 10-Apr-2001
% Revised by PJ 14-Jun-2003 
