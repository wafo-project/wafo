%README New features, bug fixes, and changes in MARKOV.
%
% Version 2.1.1 07-Sep-2005
%
% Some updated routines.
%
% UPDATED FUNCTIONS
% ~~~~~~~~~~~~~~~
% ESTSMCTP Adapted to Matlab 7 by replacing fmins by fminsearch
%
%
%--------------------------------------------------------------------
% Version 2.1.0   07-Apr-2004
% 
% No known bug fixes or updated routines.
%--------------------------------------------------------------------
% Version 2.0.5   14-Jun-2003
%
% Bug fixes, some new and/or updated routines.
%
%--------------------------------------------------------------------
% Version 2.0.4   11-Apr-2001
% 
% This file describes new features, bug fixes, and changes in this version 
% of the WAFO/markov Toolbox.
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
%  - Calculation of the expected Asymmetric RainFlow Matrix (ARFM).
%    Routines for MCTP and SMCTP. See [1]
%
%  - Inversion of an expected Asymmetric RainFlow Matrix to a Markov matrix.
%    See [1]
%
%  References:
%
%  1.  P. Johannesson (1999):
%      Rainflow Analysis of Switching Markov Loads.
%      PhD thesis, Mathematical Statistics, Centre for Mathematical Sciences,
%      Lund Institute of Technology.
%
%  Download: <http://www.math.chalmers.se/~parj/Publications/>
%
%  NEW FUNCTIONS
%  ~~~~~~~~~~~~~
%
% Markov model. [markov]
%   mctp2arfm   - Asymmetric rainflow matrix for a MCTP.
%   arfm2mctp   - Markov matrix given an asymmetric rainflow matrix. 
%
% Switching Markov model. [markov]
%   smctp2arfm  - Asymmetric rainflow matrix for a SMCTP.

