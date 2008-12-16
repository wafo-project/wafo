% Module ONEDIM in WAFO Toolbox. 
% Version 2.2.1  05-May-2008 
% 
%   Readme       - Readme file for module ONEDIM in WAFO Toolbox.
%
% Misc
%   detrendma    - Removes a trend from data using a moving average
%   findoutliers - Finds the indices to spurious points in a timeseries
%   reconstruct  - reconstruct the spurious/missing points of timeseries
%
% Plotting
%   covplot      - Plots the auto covariance function (ACF) 1D or 2D. 
%   plotlc       - Plots level-crossing spectrum (lc) 
%   spwaveplot   - Plots specified waves from a timeseries
%   waveplot     - Plots the surface elevation of timeseries.
%
% Autocovariance and spectrum estimation
%   dat2cor      - Estimate auto correlation function from data.
%   dat2cov      - Estimate auto covariance function from data.
%   dat2spec     - Estimate one-sided spectral density from data.
%   dat2spec2    - Estimate one-sided spectral density, version 2. 
%   welch_cpsd   - Estimate cross power spectral density using Welch's method
%   welch_psd    - Estimate power spectral density using Welch's method
%   welch_psd2   - Estimate power spectrum using Welch's method alternative 
%
% Data extraction
%   dat2crossind - Finds indices to level v down and/or upcrossings from data
%   dat2lc       - Extracts level-crossing spectrum from data, 
%   dat2steep    - Extracts waveheights and steepnesses from data.
%   dat2tc       - Extracts troughs and crests from data.
%   dat2tp       - Extracts turning points from data,
%   dat2wa       - Extracts sequence of wavelengths from data.
%   ecross       - Extracts exact level v crossings 
%   findcross    - Finds indices to level v up and downcrossings of a vector
%   findextrema  - Finds indices to minima and maxima of data
%   findrfc      - Finds indices to rainflow cycles of a sequence of TP.
%   mm2lc        - Extracts level-crossing spectrum from min2Max cycles.   
%
%
%   dat2midind   - Finds indices to midpoints between a min and Max and Max and min.
%   hs2sign      - Calculates a ratio-significant value of a histogram.
%   vc2sign      - Calculates ratio-significant value of the input vector.
% 
%   parzen       - returns the N-point Parzen window in a column vector.
%   hanning      - returns the N-point Hanning window in a column vector.
