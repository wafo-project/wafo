%DATASTRUCTURES of spectrum, covariance function and density (pdf) in WAFO
%
% To represent spectra, covariance functions and probability density functions
% in WAFO, the MATLAB datatype 'structured array' is used. Here follows a list
% of the fields in the struct representing S, cvf and pdf, respectively.
%
% Spectrum structure
% ~~~~~~~~~~~~~~~~~~
%  Requisite fields:
%   .type  String: 'freq', 'dir', 'k2d', k1d', 'encdir' or 'enc'.
%   .S     Spectrum values (size=[nf 1] or [np nf]).
%   .w OR .f OR .k Frequency/wave number lag, length nf.
%   .tr    Transformation function (default [] (none)). 
%   .h     Water depth (default inf).
%   .norm  Normalization flag, Logical 1 if S is normalized, 0 if not
%   .note  Memorandum string.
%   .date  Date and time of creation or change.
%  Type-specific fields:
%   .k2    Second dim. wave number lag, if .type='k2d' or 'rotk2d', length np.
%   .theta Angular lags, if .type='dir', 'rotdir' or 'encdir', length np.
%   .v     Ship speed, if .type = 'enc' or 'encdir'.
%   .phi   angle of rotation of the coordinate system
%          (counter-clocwise) e.g. azymuth of a ship.
%
% See also  createspec, plotspec
%
% Covariance function (cvf) structure
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   .R Covariance function values. Size [ny nx nt], all singleton dim. removed.
%   .x     Lag of first space dimension, length nx.
%   .y     Lag of second space dimension, length ny.
%   .t     Time lag, length nt.
%   .h     Water depth.
%   .tr    Transformation function.
%   .type  'enc', 'rot' or 'none'.
%   .v     Ship speed, if .type='enc'
%   .phi   Rotation of coordinate system, e.g.  direction of ship 
%   .norm  Normalization flag, Logical 1 if autocorrelation, 0 if covariance.
%   .Rx ... .Rtttt   Obvious derivatives of .R.
%   .note  Memorandum string.
%   .date  Date and time of creation or change.
%
% See also  createcov, spec2cov, cov2spec, covplot
%
% Probability density function (pdf) structure
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Describing a density of n variables:
%   .f     Probability density function values (n-dimensional matrix)
%   .x     Cell array of vectors defining grid for variables (n cells)
%   .labx  Cell array of label strings for the variables (n cells)
%   .title Title string                           
%   .note  Memorandum string.
%
% See also  createpdf, pdfplot

% History: 
% revised by IR 03.04.2001
% revised pab 21.01.2000 
%  - spellchecked the file 
%  - added norm to spec
%  - changed reference to specplot to plotspec
% revised by es 20.10.1999 (cell arrays in pdf-struct)
%          by es 13.10.1999
more on
help datastructures
more off

