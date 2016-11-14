function options = specoptset(varargin)
%SPECOPTSET Create or alter SPECTRUM OPTIONS structure.
%
%  CALL:  options = specoptset(funcname,opts1,opts2,...,par1,val1,par2,val2,...);
%
%   options    = transformation options structure in which the named 
%                parameters have the specified values.  
%   funcname   = string giving the name of the function for which default
%                values for the options structure should be extracted.
%                Options are 'dat2dspec', 'dat2spec'.
%   opts1,
%   opts2..    = options structures
%   par1,par2..= strings identifying the parameter to alter
%   val1,val2..= corresponding values the parameters are altered to.
%   
%   SPECOPTSET combines the default options for a function given by FUNCNAME
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
%   SPECOPTSET with no input arguments and no output arguments displays all 
%   parameter names and their possible values.
%
%   SPECOPTSET with no input arguments creates an options structure
%   OPTIONS where all the fields are set to [].
%
%   
% SPECOPTSET PARAMETERS
%     nharm     = number of harmonics in Fourier coefficients (Fcof)     (default 2)
%   gravity     = acceleration of gravity                                (default see gravity)
%  wdensity     = water density,                                         (default see wdensity)
%      bet      = 1, theta given as directions toward which waves travel (default)
%                -1, theta given as directions from which waves come 
%      igam     = 1, if z is measured + upward from mean water level     (default)
%                 2, if z is measured + downward from mean water level
%                 3, if z is measured + upward from sea floor
% x-axisdir     = angle clockwise from true north to + x-axis in degrees (default 90)
% y-axisdir     = angle clockwise from true north to + y-axis in degrees (default 0)
%  plotflag     = 0 no plotting (default)
%                 1 plots the spectrum, S,     
%                 2 plot  plots spectral density and the directional spreading
%  dflag        = specifies a detrending performed on the signal before estimation.
%                 'mean','linear' or 'ma' (= moving average)             (default 'mean')   
%  ftype        = frequency type, 'w' or 'f'  (default 'w')
% Plotflag      = 'off'   or 0: No plotting (Default)
%                 'final' or 1: Plot final result
%                 'iter'  or 2: Monitor the development.
% Delay         = Delay time for each plot when PLOTFLAG=='iter'.
% message       = Level of screeen display.
% maxIter       = maximum number of iterations used in estimate. various
%                 effects for differnt methods
% maxCoef       = maximum spectral coefficient
% coefAbsTol    = acceptance tolerance for spectral coefficents.
% errorTol      = acceptance tolerance for the directional spreading.
% minModelOrder = Minimum model order (BDM, EMEM)
% maxModelOrder = Maximum model order (BDM, EMEM)
% relax         = start relaxation parameter
%
% Examples:
%  true_opt = struct('dflag','mean','ftype', 'w',...
%         'thtype','r','noverlap', 0,'window', [], 'message', 2000,...
%         'plotflag', 'off','delay', 0,'nharm', 2,'gravity', 9.80629386676700,...
%         'wdensity', 1027.78633943154,'bet', 1,'igam', 1,'x_axisdir', 90,...
%         'y_axisdir', 0, 'maxiter', 25, 'coefabstol',  0.01,...
%         'errortol',  0.01,'minmodelorder',  1,'maxmodelorder', [],...
%         'relax',  1, 'maxcoef', []);
%  assert(specoptset('dat2dspec'), true_opt, 1e-10);
%  true_opt.ftype='f';
%  true_opt.thtype='degrees';
%  assert(specoptset('dat2dspec','ftype','f','thtype','degrees'), ...
%         true_opt, 1e-10);
%
% See also  dat2spec, dat2dspec

% History
% revised pab, Dec 2003
%  replaced code with call to parseoptions  
% by pab 22.10.2002
% based on the idea of MATLAB's optimset


% Print out possible values of properties.
if (nargin == 0) && (nargout == 0)
%  disp('          method: [ {''cov''} | ''psd'' | ''psdo'' | ''pmem'' | ''pburg'' ] (dat2spec)')
%  disp('          method: [ {''mlm''} | ''imlm'' | ''mem'' | ''emem'' | ''fem'' ] (dat2dspec)')
  disp('           dflag: [ {''mean''} | ''linear'' | ''ma'' (= moving average)')
  disp('           ftype: [ {''w''} | ''f'' ]')
  disp('          thtype: [ {''radians''} | ''degrees'' ]') 
%  disp('            nfft: [ {128} ]')
  disp('        noverlap: [ {0} ]')
  disp('          window: [ vector {hanning(nfft)}]')
  disp('         gravity: [ positive scalar or vector {gravity}]')
  disp('        wdensity: [ positive integer {wdensity}]')
  disp('             bet: [ {1} | 1 ]' )
  disp('            igam: [ {1} | 2 |3 ]')
  disp('       x_axisdir: [ {90} ]')
  disp('       y_axisdir: [ {0} ]')
  disp('         message: [ positive integer {2000}]')
  disp('        Plotflag: [ off | {final} | iter ]')
  disp('           Delay: [ positive scalar {0} ]')
  disp('         maxIter: [ {30} ]')
  disp('      coefAbsTol: [ {0.01} ]')
  disp('         maxcoef: [ [] ]')
  disp('        errorTol: [ {0.01} ]')
  disp('   minModelOrder: [ {1} ]')
  disp('   maxModelOrder: [ {[]} ]')
  disp('           relax: [ {1} ]')
  
  return;
end

% Initialization
% Legal functions names
fnames = strvcat('dat2spec','dat2dspec'); 
% Legal parameter names
names  = lower({'dflag','ftype','thtype','noverlap','window',...
		'message','plotflag','delay', ...
		'nharm','gravity','wdensity',...
		'bet', 'igam','x_axisdir','y_axisdir',...
		'maxIter','coefAbsTol','errorTol',...
		'minModelOrder','maxModelOrder','relax','maxCoef'});     

%default values
defaultVals = {'mean','w','radians',0,[],....
	       2000,'final',0,...
	       2,gravity,wdensity,...
	       1,1,90, 0,...
	       30, 0.01,0.01,...
	       1 , [], 1, []};

options = cell2struct(defaultVals,names,2);
options = parseoptions(fnames,options,varargin{:});
return