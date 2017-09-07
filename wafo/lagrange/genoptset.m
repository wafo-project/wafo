function options3D  = genoptset(varargin)
%GENOPTSET Creates or alters 3D generation options structure
%
%CALL: options3D = genoptset(funcname,opts1,opts2,...,par1,val1,par2,val2,...);
%
%   options3D    = options structure in which the named 
%                parameters have the specified values.  
%   funcname   = string giving the name of the function for which default
%                values for the options structure should be extracted.
%
%   par1,par2..= strings identifying the parameter to alter
%   val1,val2..= corresponding values the parameters are altered to.
%   
%   GENOPTSET sets options for transformation of ldat3D (W,X,Y) to 
%   3D Lagrange fields or 2D timeseries
%
%   GENOPTSET with no input arguments and no output arguments displays all 
%   parameter names and their possible values.
%
% See also troptset and simoptset

% Tested on Matlab 8.1, 8.6, 9.1
% History
% Created by GL March 29 2015 for generation of 3D Lagrange fields
% Expansion of simoptset
% simoptset Expanded by GL Sept 22, 2014 with 3D simulation options
% simoptset Expanded by GL June 28, 2013 to choose timefft on/off
% Started by GL Jan 18, 2007
% based on troptset by PAB

% Print out possible values of properties
if (nargin==0) && (nargout==0)
  disp(' type:      [ moviedata | { field, timeseries } ]')
  disp(' t0:        [ real time point for single field or empty ]')
  disp(' PP:        [ 2 x n  array of real coordinates for  n  time series or empty ]')
  disp(' start:     [ [x,y]  real coordinates for lover left corner of fields or empty ]')
  disp(' end:       [ [x,y]  real coordinates for upper right corner of fields or empty ]')  
  disp(' rate:      [ positive integer for interpolation of time series or empty ]')  
  disp(' plotflag:  [ on | { off } ]') 
  return
end

% Legal parameter names
names = {'type','t0','PP','start','end','rate','plotflag'};

% Default values
defaultVals = {'moviedata',[],[],[10 10],[],1,'on'};

options3D = cell2struct(defaultVals,names,2);
options3D = parseoptions(options3D,varargin{:});

return
