function options = kdeoptset(varargin)
%KDEOPTSET Create or alter KDE OPTIONS structure.
%
%  CALL:  options = kdeoptset(funcname,opts1,opts2,..,par1,val1,par2,val2,..);
%
%   options    = transformation options structure in which the named 
%                parameters have the specified values.  
%   funcname   = string giving the name of the function for which default
%                values for the options structure should be extracted.
%                Options are 'kdebin', 'kde', 'fftkde'.
%   opts1,
%   opts2..    = options structures
%   par1,par2..= strings identifying the parameter to alter
%   val1,val2..= corresponding values the parameters are altered to.
%   
%   KDEOPTSET combines the default options for a function given by FUNCNAME
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
%   KDEOPTSET with no input arguments and no output arguments displays all 
%   parameter names and their possible values.
%
%   KDEOPTSET with no input arguments creates an options structure
%   OPTIONS where all the fields are set to their default values.
%
%   
% KDEOPTSET PARAMETERS
%   kernel = String defining the kernel function. Options are:
%            'epanechnikov'  - Epanechnikov kernel. (default)
%            'biweight'      - Bi-weight kernel.
%            'triweight'     - Tri-weight kernel.  
%            'triangular'    - Triangular kernel.
%            'gaussian'      - Gaussian kernel
%            'rectangular'   - Rectangular kernel. 
%            'laplace'       - Laplace kernel.
%            'logistic'      - Logistic kernel.
%   hs     = smooting parameter vector/matrix.
%            (default compute from data using hsMethod)
% hsMethod = string defining the method to compute the  smooting
%            parameter vector/matrix, hs, if not explicitly given.
%            Valid options are one of the methods:
%            'hns', 'hos', 'hldpi', 'hbcv', 'hste', 'hboot', 'hlscv',
%            'hscv', 'hstt'
%            (default 'hns')  
%   alpha  = sensitivity parameter               (default 0 regular KDE)
%            A good choice might be alpha = 0.5 ( or 1/D)
%            alpha = 0      Regular  KDE (hs is constant)
%            0 < alpha <= 1 Adaptive KDE (Make hs change)  
%   L2     = vector of transformation parameters (default 1 no transformation)
%            or a cellarray of transformations parameters
%            or non-parametric transformations.
%            t(xi;L2) = xi^L2*sign(L2)   for L2(i) ~= 0
%            t(xi;L2) = log(xi)          for L2(i) == 0 
%            If single value of L2 is given then the transformation
%            is the same in all directions.
%   inc    = scalar defining the dimension of kde (default 128)
%            (A value below 50 is very fast to compute but may give some
%            inaccuracies.  Values between 100 and 500 give very accurate 
%            results)  (Only used in function kdebin.)
%  
%  Note that only the first 4 letters of the kernel name is needed. 
%  If L2~=1 KDE transforms the data before estimation. The final estimate
%  is obtained by transforming back by a simple change of variables.
%  Beaware of spurious spikes close to the edges when L2~=1.
%  These spikes are due to numerical problems close to the edges.
%
% (1) If a single non-zero value of hs is given then smoothing is the 
%     same in all directions.
% (2) If hs = [H1, H2] then smoothing is by H1 in the 1'st direction and H2
%     in the 2'nd direction.  If hs = 0 or hs = [] is specified, automatic 
%     values are used computed with hsMethod. 
% (3) If hs = [H1, H3;H3 ,H2] then the smoothing can be in orientations
%     different from those of the co-ordinate directions.  
%     Example: HS=(.3*eye(D))*sqrtm(cov(data))
%     Note: det(HS)>0
% 
%  If hs is too large we over smooth introducing a large bias
%  in the estimated density. If hs is too small we under smooth
%  introducing spurious bumps.
%
%  Also note that HNS only gives a optimal value when the distribution of 
%  the (transformed) DATA is Gaussian. If the distribution is asymmetric, 
%  multimodal or have long tails the HNS may return a too large smoothing
%  parameter  i.e. KDE may oversmooth and thus give a large bias.  To
%  remedy this try:
%
% (1) different values for HS and check by eye, i.e. 
%     for 2D data count the points lying within the contours to roughly
%     check that the calculated contours approximately encloses the
%     corresponding percentage of the data.  
% (2) with alpha > 0, i.e., adaptive KDE (Good choice for skew distributions)
% (3) with L2~=1, i.e., transformation KDE. A reasonable value of L2 is
%     found when a normal-plot of the transformed DATA is approximately
%     linear. Beaware that due to numerical problems the back
%     transformation may result in spurious spikes close to where Xi is
%     zero. These spikes should be removed!  
%  
% Examples:
%  assert(kdeoptset('kdebin'), struct('kernel', 'epan', 'hs', [], ...
%         'hsMethod', 'hns', 'alpha', 0, 'L2', 1, 'inc', 128));
%
% See also  kde, kdebin

% History
% by pab 04.02.2005%  NEW FUNCTIONS
%  ~~~~~~~~~~~~~
% based on MATLAB's optimset


% Print out possible values of properties.
if (nargin == 0) && (nargout == 0)
  help kdeoptset
  return;
end

% Initialization
% Legal functions names
fnames = strvcat('kde','kdebin'); 
% Legal parameter names
names  = {'kernel','hs','hsMethod','alpha','L2','inc'};     
vals = {'epan',[],'hns',0,1,128};

% Initialize options with default values
options = cell2struct(vals,names,2);
options = parseoptions(fnames,options,varargin{:});
if ~isempty(options.alpha)
  options.alpha=min(abs(options.alpha),1);
end
%if isnumeric(options.hs)
%  options.hsMethod = [];
%end
return
 
