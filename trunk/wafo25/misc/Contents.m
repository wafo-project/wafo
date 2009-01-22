% Module MISC in WAFO Toolbox. 
% Version 2.1.1   15-Sep-2005 
% 
%   readme       - Readme file for module MISC in WAFO Toolbox 
% 
% Functions
%   betaloge     - Natural Logarithm of beta function.
%   binom        - Calculates the binomial coefficient n!/((n-k)!*k!)
%   findpeaks    - find peaks of vector or matrix possibly rainflow filtered
%   genchol      - Generalized Cholesky factorization
%   hypgf        - Hypergeometric function F(a,b,c,x) 
%   sinc         - Sin(pi*x)/(pi*x) function.
%   smooth       - Calculates a smoothing spline.
%   savgol       - Savitzky-Golay filter coefficients.
%   stirlerr     - Computes  log(n!) - log( sqrt(2*pi*n)*(n/exp(1))^n )
%   convlv       - Convolves real data set with a response function. 
%
% Oceanographic constants
%   getshipchar  - Estimate ship characteristics from value of one ship-property
%   gravity      - returns the constant acceleration of gravity 
%   wdensity     - Returns the water density 
%
% Integration
%   ccquad       - Numerical integration using a Clenshaw-Curtis quadrature.
%   gaussq       - Numerically evaluate integral, Gauss quadrature.
%   gaussq2d     - Numerically evaluate2D integral, Gauss quadrature.
%   simpson      - Numerical integration with the Simpson method
%
%   grule        - Computes nodes and weights for Gauss-Legendre quadrature.
%   qrule        - Computes nodes and weights for Gaussian quadratures.   
%   qrule2d      - Computes nodes and weights for Gaussian quadratures 
%   hrule        - Computes nodes and weights for generalized Hermite quadrature.
%   jrule        - Computes nodes and weights for Gauss-Jacobi quadrature.
%   lrule        - Computes nodes and weights for generalized Laguerre  quadrature. 
%
% Documentation utilities
%   geth1line    - Extracts the first comment line (the H1 line) of a m-file
%   mkcontents   - Makes Contents file in current working directory.
%   wafoversion  - Wave Analysis for Fatigue and Oceanography version.
%
% Misc
%   comnsize     - Calculates common size of all non-scalar arguments.
%   iscomnsize   - True if all non-scalar arguments are of common size.
%   isoctave     - True if function run in an Octave environment
%   ismatlab     - True if function run in an Matlab environment
%   keep         - Keeps workspace variables of your choice and clears the rest.
%   levels       - Calculates discrete levels given the parameter matrix.
%   loaddata     - Loads a matrix from a text file.
%   parent       - Lists parents of a Custom object class
%   parseoptions - Create or alter a OPTIONS structure.
%   where        - List all locations of one or more Matlab functions.

