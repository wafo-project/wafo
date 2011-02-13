% Module KDETOOLS in WAFO Toolbox. 
% Version 2.5.2   07-Feb-2011 
% 
%   Readme      - Readme file for module KDETOOLS in WAFO Toolbox.
% 
% Misc
%   chi2gof     - CHI Squared Goodness Of Fit test.
%   chi2gof2    - CHI Squared Goodness Of Fit test.
%   evalpdf     - evaluates a PDF struct by interpolation
%   fftce       - Circulant Embedding of a vector or matrix
%   kdedemo1    - Demonstrate the smoothing parameter impact on KDE
%   kdedemo2    - Demonstrate the difference between transformation- and ordinary-KDE
%   kdeoptset   - Create or alter KDE OPTIONS structure.
%   qlevels     - Calculates quantile levels which encloses P% of PDF
%   qlevels2    - Calculates quantile levels which encloses P% of data
%   vsph        - calculates volume of  d-dimensional sphere with radius r0
%
% Bandwidth Selectors
%   hbcv        - Biased Cross-Validation estimate of smoothing parameter.
%   hbcv2       - Biased Cross-Validation smoothing parameter for 2D data.
%   hboot       - Bootstrap cross-validation estimate of smoothing parameter.
%   hldpi       - L-stage Direct Plug-In estimate of smoothing parameter.
%   hldpi2      - L-stage DPI estimate of smoothing parameter for 2D data
%   hldpi2fft   - L-stage DPI estimate of smoothing parameter for 2D data
%   hlscv       - Least Squares Cross-Validation estimate of smoothing parameter
%   hmns        - Multivariate Normal Scale Estimate of Smoothing Parameter.
%   hns         - Normal Scale Estimate of Smoothing Parameter.
%   hos         - Oversmoothing Parameter.
%   hscv        - Smoothed cross-validation estimate of smoothing parameter.
%   hste        - 2-Stage Solve the Equation estimate of smoothing parameter.
%   hstt        - Scott-Tapia-Thompson estimate of smoothing parameter.
%
% Kernel density estimators and data-gridders
%   bincount    - 1-dimensional Bin Count
%   ffndgrid    - Fast 'n' Furious N-D data gridding.
%   gridcount   - D-dimensional histogram using linear binning.
%   kde         - Kernel Density Estimator.
%   kdebin      - Binned Kernel Density Estimator.
%   kdefun      - Kernel Density Estimator.
%   kde1dgui    - GUI to Kernel Density Estimator.
%   kde2dgui    - GUI to Kernel Density Estimator in two dimensions.
%
% Kernel functions
%   deriv       - 4th, 6th, 8th and 10th derivatives of the kernel function.
%   deriv2      - High order partial derivatives of the Gaussian kernel.
%   kernelstats - Return 2'nd order moment of kernel pdf
%   mkernel     - Multivariate Kernel Function.
%   mkernel2    - Multivariate Kernel Function, alternative version.
%
% Random number generators
%   mkernelrnd  - Random numbers from the Multivariate Kernel Function.
%   sample      - Random sampling of points from a data-set
%   ssample     - Random sampling from a smoothed empirical distribution
%
%   pdfstat     - Mean and variance for the PDF 2D distribution
