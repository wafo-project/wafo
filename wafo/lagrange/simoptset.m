function options  = simoptset(varargin)
%SIMOPTSET Creates or alters simulation options structure
%
%CALL:  options = simoptset(funcname,opts1,opts2,...,par1,val1,par2,val2,...);
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
%   SIMOPTSET combines the default options for a function given by FUNCNAME
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
%   SIMOPTSET with no input arguments and no output arguments displays all 
%   parameter names and their possible values.
%
%   SIMOPTSET with no input arguments creates an options structure
%   OPTIONS where all the fields are set to []. ???
%
% See also troptset, genoptset

% Tested on Matlab 8.1, 8.6
% History
% Expanded by GL May 22, 2015, to allow field 'u'
%   If  u  is non-empty, = [u1 u2 Nu], then the u-vector will be  
%       u = linspace(u1,u2,Nu)
%   If  u  is empty then the u-vector will be  
%       u = linspace(0,(Nu-1)*du,Nu)
% Expanded by GL Sept 22, 2014, with 3D simulation options
% Expanded by GL June 28, 2013, to choose timefft on/off
%   removed lp and ltheta
%   set .iseed dafault to 'shuffle'
% Started by GL Jan 18, 2007
% based on troptset by PAB

% Print out possible values of properties
if (nargin==0) && (nargout==0)
  disp('    Nt: [ positive integer or empty {[]} ]')
  disp('    Nu: [ positive integer or empty {[]} ]')
  disp('    Nv: [ positive integer or empty {[]} ]')
  disp('    dt: [ real or empty {[]} ]')
  disp('    du: [ positive real or empty {[]} ]')
  disp('    dv: [ positive real or empty {[]} ]')
  disp('    u: [ [real real positive integer] or empty ]')
  disp('    lalpha: [ real array {0} ]')   
  disp('    lbeta: [ real {0} ]')   
  disp('    ffttype: [ ffttime | {fftspace, ffttwodim} ]')
  disp('    iseed: [ shuffle (sets random seed) | {positive integer (int32)} ]')
  disp('    plotflag: [ 0 | {1,2,3} ]')    
  return
end

% Legal function names
%fnames = strvcat('spec2ldat','spec2ldat3D') 
fnames = char('spec2ldat','spec2ldat3D');

% Legal parameter names
names = {'Nt','Nu','Nv','dt','du','dv','u',...
    'lalpha','lbeta','ffttype','iseed','plotflag'};

% Default values
defaultVals = {2048,2048,[],0.5,0.5,[],[],0,0,'ffttime','shuffle',0};

options = cell2struct(defaultVals,names,2);
options = parseoptions(fnames,options,varargin{:});

if isnumeric(options.iseed),
    options.iseed=int32(options.iseed);
end
return
