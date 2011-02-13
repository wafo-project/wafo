function options = troptset(varargin)
%TROPTSET Create or alter TRANSFORM OPTIONS structure.
%
%  CALL:  options = troptset(funcname,opts1,opts2,...,par1,val1,par2,val2,...);
%
%   options    = transformation options structure in which the named 
%                parameters have the specified values.  
%   funcname   = string giving the name of the function for which default
%                values for the options structure should be extracted.
%                Options are 'dat2tr', 'lc2tr', 'reconstruct'.
%   opts1,
%   opts2..    = options structures
%   par1,par2..= strings identifying the parameter to alter
%   val1,val2..= corresponding values the parameters are altered to.
%   
%   TROPTSET combines the default options for a function given by FUNCNAME
%   with new options structures (OPTS1,OPTS2,...) and/or with the named
%   parameters (PAR1,PAR2,...) with the corresponding values (VAL1, VAL2,...).
%   The parameters are set in the same order as the input arguments.
%   Any parameters with non-empty values of the options struct overwrite
%   the corresponding old parameters. 
%
%   The input arguments can be given in any order with one exception:
%   PARx and VALx must be given in pairs in that order.
%   Any unspecified parameters for PARx are set to []. 
%   Parameters with value [] indicate the use of default values for that
%   parameter when OPTIONS is passed to the function. It is sufficient to
%   type only the 2 first characters to uniquely identify the parameter
%   or function name.  Upper case letters for parameter names and values
%   that are strings are ignored. If an invalid string is provided, the
%   default is used.
%   
%   TROPTSET with no input arguments and no output arguments displays all 
%   parameter names and their possible values.
%
%   TROPTSET with no input arguments creates an options structure
%   OPTIONS where all the fields are set to [].
%
%   
% TROPTSET PARAMETERS
% ChkDer    - 'off' or 0: No check on the derivative of the transform.
%             'on'  or 1: Check if transform have positive derivative 
% Csm, Gsm  - Defines the smoothing of the crossing intensity 
%             and the transformation g, respectively. Valid values must 
%             be 0<=Csm,Gsm<=1. (default Csm=0.9, Gsm=0.05)
%             Smaller values gives smoother functions.
% Crossdef  - Crossing definition used in the crossing spectrum:
%             'u'   or 1: only upcrossings
%             'uM'  or 2: upcrossings and Maxima (default)
%             'umM' or 3: upcrossings, minima, and Maxima.
%             'um'  or 4: upcrossings and minima.
% Plotflag  - 'off'   or 0: No plotting (Default)
%             'final' or 1: Plot final result
%             'iter'  or 2: Monitor the development.
% Delay     - Delay time for each plot when PLOTFLAG=='iter'.
% Param     - Vector which defines the region of variation of the data x.
%             (default [-5 5 501]). 
% LinExtrap - 'off' or 0: uses a regular smoothing spline. 
%             'on'  or 1: use a smoothing spline with a constraint on the
%                         ends to ensure linear extrapolation outside the
%                         range of the data. (default)
% Cvar      - Variances for the the crossing intensity. (default  1) 
% Gvar      - Variances for the empirical transformation, g. (default  1) 
% Ne        - Number of extremes (maxima & minima) to remove from the
%             estimation of the transformation. This makes the
%             estimation more robust against outliers. (default 7)
% Ntr       - Maximum length of empirical crossing intensity or CDF.
%             The empirical crossing intensity or CDF is interpolated
%             linearly  before smoothing if their lengths exceeds Ntr.
%             A reasonable NTR will significantly speed up the
%             estimation for long time series without loosing any
%             accuracy. NTR should be chosen greater than
%             PARAM(3). (default 1000)
% multip    - 0 the data in columns belong to the same seastate (default).
%             1 the data in columns are from separate seastates.
%
% Examples:
%  troptset('lc2tr')
%  troptset('lc2tr','csm',.99)
%  troptset('csm',.99,'lc2tr')  % is the same as  troptset('lc2tr')
%
% See also  dat2tr, lc2tr, cdftr, reconstruct

% History
% revised pab 21Nov2003
%  -moved some code into parseoptions for easier maintainence    
% by pab 20.12.2000
% based on MATLAB's optimset


% Print out possible values of properties.
if (nargin == 0) && (nargout == 0)
  disp('          ChkDer: [ off | {on} ]')
  disp('             Csm: [ scalar between 0 and 1 {0.95} ]')
  disp('             Gsm: [ scalar between 0 and 1 {0.05} ]')
  disp('        Crossdef: [ u |{uM}| umM | um ]') 
  disp('           Param: [ vector {-5 5 513} ]')
  disp('       LinExtrap: [ off | {on} ]')
  disp('            Cvar: [ positive scalar or vector {1}]')
  disp('            Gvar: [ positive scalar or vector {1}]')
  disp('              Ne: [ positive integer {7}]')
  disp('             Ntr: [ positive integer {2000}]')
  disp('        Plotflag: [ off | {final} | iter ]')
  disp('           Delay: [ positive scalar {0} ]')
  return;
end

% Initialization
% Legal functions names
fnames = strvcat('dat2tr','lc2tr','cdf2tr','reconstruct'); 

% Legal parameter names
names = {'chkder','csm','gsm','crossdef','param','linextrap','cvar', ...
	 'gvar','ne','ntr','plotflag','delay','multip'};
%default values
defaultVals = {'on',0.95,0.05,'uM',[-5 5 513],'on',1,1,7,2000,'final',0,0};

options = cell2struct(defaultVals,names,2);
options = parseoptions(fnames,options,varargin{:});

return

