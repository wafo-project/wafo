function options = rindoptset(varargin)
%RINDOPTSET Create or alter RIND OPTIONS structure.
%
%  CALL:  options = rindoptset(funcname,opts1,opts2,..,par1,val1,par2,val2,..);
%
%   options    = transformation options structure in which the named 
%                parameters have the specified values.  
%   funcname   = string giving the name of the function for which default
%                values for the options structure should be extracted.
%                Options are 'rind', 'spec2mmtpdf', 'spec2thpdf'.
%   opts1,
%   opts2..    = options structures
%   par1,par2..= strings identifying the parameter to alter
%   val1,val2..= corresponding values the parameters are altered to.
%   
%   RINDOPTSET combines the default options for a function given by FUNCNAME
%   with new options structures (OPTS1,OPTS2,...) and/or with the named
%   parameters (PAR1,PAR2,...) with the corresponding values (VAL1, VAL2,...).
%   The parameters are set in the same order as the input arguments.
%   Any parameters with non-empty values of the options struct overwrite
%   the corresponding old parameters. 
%   The input arguments can be given in any order with one exception:
%   PARx and VALx must be given in pairs in that order.
%   Any unspecified parameters for PARx are set to []. 
%   Parameters with value [] indicate to use the default value for that
%   parameter when OPTIONS is passed to the function. It is sufficient to
%   type only the 2 first characters to uniquely identify the parameter
%   or function name.  Upper case letters for parameter names and values
%   that are strings are ignored. If an invalid string is provided, the
%   default is used.
%   
%   RINDOPTSET with no input arguments and no output arguments displays all 
%   parameter names and their possible values.
%
%   RINDOPTSET with no input arguments creates an options structure
%   OPTIONS where all the fields are set to their default values.
%
%   
% RINDOPTSET PARAMETERS
%  METHOD  = INTEGER defining the integration method
%            0 Integrate by Gauss-Legendre quadrature  (Podgorski et al. 1999)
%            1 Integrate by SADAPT for Ndim<9 and by KRBVRC otherwise 
%            2 Integrate by SADAPT for Ndim<20 and by KRBVRC otherwise 
%            3 Integrate by KRBVRC by Genz (1993) (Fast Ndim<101) (default)
%            4 Integrate by KROBOV by Genz (1992) (Fast Ndim<101)
%            5 Integrate by RCRUDE by Genz (1992) (Slow Ndim<1001)
%            6 Integrate by SOBNIED               (Fast Ndim<1041)
%            7 Integrate by DKBVRC by Genz (2003) (Fast Ndim<1001)
%            
%  XCSCALE = REAL to scale the conditinal probability density, i.e.,
%            f_{Xc} = exp(-0.5*Xc*inv(Sxc)*Xc + XcScale) (default XcScale =0)
%  ABSEPS  = REAL absolute error tolerance.       (default 0)
%  RELEPS  = REAL relative error tolerance.       (default 1e-3)
%  COVEPS  = REAL error tolerance in Cholesky factorization (default 1e-13)
%  MAXPTS  = INTEGER, maximum number of function values allowed. This 
%            parameter can be used to limit the time. A sensible 
%            strategy is to start with MAXPTS = 1000*N, and then
%            increase MAXPTS if ERROR is too large.    
%            (Only for METHOD~=0) (default 40000) 
%  MINPTS  = INTEGER, minimum number of function values allowed.
%            (Only for METHOD~=0) (default 0)
%  SEED    = INTEGER, seed to the random generator used in the integrations
%            (Only for METHOD~=0)(default floor(rand*1e9))
%  NIT     = INTEGER, maximum number of Xt variables to integrate
%            This parameter can be used to limit the time. 
%            If NIT is less than the rank of the covariance matrix,
%            the returned result is a upper bound for the true value
%            of the integral.  (default 1000)
%  XCUTOFF = REAL cut off value where the marginal normal
%            distribution is truncated. (Depends on requested
%            accuracy. A value between 4 and 5 is reasonable.)
%  XSPLIT  = parameters controlling performance of quadrature
%             integration:
%             if Hup>=xCutOff AND Hlo<-XSPLIT OR
%                Hup>=XSPLIT AND Hlo<=-xCutOff then
%             do a different integration to increase speed
%             in rind2 and rindnit. This give slightly different 
%            results
%            if XSPILT>=xCutOff => do the same integration allways
%            (Only for METHOD==0)(default XSPLIT = 1.5)   
%  QUADNO  = Quadrature formulae number used in integration of Xd
%            variables. This number implicitly determines number of nodes
%            used.  (Only for METHOD==0)
%  SPEED   = defines accuracy of calculations by choosing different 
%            parameters, possible values: 1,2...,9 (9 fastest,  default []).
%            If ~isempty(SPEED) the parameters, ABSEPS, RELEPS, COVEPS,
%            XCUTOFF, MAXPTS and QUADNO will be set according to
%            INITOPTIONS.
% NC1C2    = number of times to use the regression equation to restrict
%            integration area. Nc1c2 = 1,2 is recommended. (default 2) 
%            (note: works only for method >0)
% Examples:
%  rindoptset('rind')
%  rindoptset('speed',5)
%
% See also  rind, rindoptset>initoptions

% History
% revised pab may 2007
% -moved function initoptions to here as a subfunction
% Revised pab July 2004
% - added NC1C2 
% by pab 20.05.2003%  NEW FUNCTIONS
%  ~~~~~~~~~~~~~
% based on MATLAB's optimset


% Print out possible values of properties.
% if (nargin == 0) && (nargout == 0)
%   help rindoptset
%   return;
% end

% Initialization
% Legal functions names
fnames = strvcat('rind','spec2mmtpdf','spec2thpdf','spec2tpdf2'); 
% Legal parameter names
names  = {'method','xcscale',...
	  'abseps','releps','coveps',...
	  'maxpts','minpts',...
	  'seed','nit','xcutoff',...
	  'xsplit','quadno', ...
	  'speed','nc1c2'};     
vals = {3,0,0,1e-3,1e-10,...
	40000,...
	0,...
	floor(rand*1e9),...
	1000,...
	[],...
	1.5,...
	[] ,...
	[],...
	2};

% Initialize options with default values
options = cell2struct(vals,names,2);
options = parseoptions(fnames,options,varargin{:});

if ~isempty(options.speed) && options.speed>0
  options = initoptions(options.speed,options);
end
return


function options = initoptions(speed,options)
%INITOPTIONS Initializes RIND options according to speed.
%                                   
% CALL options = initoptions(speed,options)    
%
%     speed = integer defining accuracy of calculations. 
%             Valid numbers:  1,2,...,10 
%            (1=slowest and most accurate,10=fastest, but less accuracy)
%   options = rind-options structure, see RINDOPTSET fort details.
%
% RIND OPTIONS parameters initialized according to speed:
% SPEED    = Integer defining accuracy of calculations. 
% ABSEPS   = Absolute error tolerance.
% RELEPS   = Relative error tolerance.
% COVEPS   = Error tolerance in Cholesky factorization.
% XCUTOFF  = Truncation limit of the normal CDF 
% MAXPTS   = Maximum number of function values allowed.
% QUADNO   = Quadrature formulae used in integration of Xd(i)
%            implicitly determining # nodes 
%
%    
% See also  rindoptset, rind
  
%  error(nargchk(1,2,nargin))
  narginchk(1,2)
  if nargin<2||isempty(options)
    options = rindoptset;
  end
  options.speed =  min(max(speed,1),13);
  if isempty(speed)
    return
  end
 
  MAXPTS = 10000;
  options.quadno = (1:3)+ (10-min(options.speed,9)) + (options.speed==1);
  switch options.speed
    case {11,12,13},
      abseps = 1d-1;
    case (10)
      abseps = 1d-2;
    case {7,8,9}
      abseps = 1d-2;
    case {4,5,6}
      MAXPTS = 20000;
      abseps = 1d-3;
    case {1,2,3}
      MAXPTS = 30000;
      abseps = 1d-4;
  end
  if options.speed<12
    TMP = max(abs(11-abs(options.speed)),1);
    TMP = mod(TMP+1,3)+1;
    options.coveps = abseps*((1.d-1)^TMP);
  elseif options.speed<13
    options.coveps = 0.1;
  else
    options.coveps = 0.5;
  end
  
  options.releps = min(abseps ,1.d-2);
  options.maxpts = MAXPTS;    
  
  if (options.method==0) 
     % This gives approximately the same accuracy as when using 
     % RINDDND and RINDNIT    
     
     %    xCutOff= MIN(MAX(xCutOff+0.5d0,4.d0),5.d0)
     abseps = abseps*1.d-1;
  end
  options.abseps  = abseps;
  truncError      = 0.05 * max(0,options.abseps);
  options.xcutoff = max(min(abs(invnorm(truncError)),7),1.2);
  options.abseps  = max(options.abseps - truncError,0);
  
 return 
 