function [F,h] = cc2cmat(param,cc,ddef,method,h,NOsubzero,alpha)
% CC2CMAT Calculates the cycle count matrix from a cycle count.
% using (0) Histogram, (1) Kernel smoothing, (2) Kernel smoothing.
%
% CALL:  [F,h] = cc2cmat(param,cc,ddef,method,h,NOsubzero,alpha);
%
% Input:
%   param     = Parameter vector, [a b n], defines the grid.
%   cc        = Cycle count with minima in column 1 and maxima in column 2. [nx2]
%   ddef      =  1: causes peaks to be projected upwards and troughs 
%                   downwards to the closest discrete level (default).
%             =  0: causes peaks and troughs to be projected
%                   the closest discrete level.
%             = -1: causes peaks to be projected downwards and the 
%                   troughs upwards to the closest discrete level.
%   method    =  0: Histogram. (Default)
%                1: Kernel estimator (constant bandwidth). 
%                2: Adaptiv kernel estimator (local bandwidth). 
%   h         = Bandwidth (Optional, Default='automatic choice')
%   NOsubzero = Number of subdiagonals that are set to zero
%               (Optional, Default = 0, only the diagonal is zero)
%   alpha     = Parameter for method (2) (Optional, Default=0.5).
%               A number between 0 and 1.
%               alpha=0 implies constant bandwidth (method 1).
%               alpha=1 implies most varying bandwidth.
%
% Output:
%   F       = Estimated cycle matrix.
%   h       = Selected bandwidth.
%
% See also  dcc2cmat, cc2dcc, smoothcmat

% Tested  on Matlab 5.3
%
% History:
% Correction by PJ 14-Feb-2000
%   Changed 'smthcmat' to 'smoothcmat'
% Revised by PJ  01-Nov-1999
%   updated for WAFO
% Created by PJ (Pär Johannesson) 1997
%   from 'Toolbox: Rainflow Cycles for Switching Processes V.1.0'

% Check input arguments

ni = nargin;
no = nargout;
error(nargchk(2,7,ni));

if ni<3, ddef=[]; end
if ni<4, method=0; end
if ni<5, h=[]; end
if ni<6, NOsubzero=[]; end
if ni<7, alpha=[]; end

if method < 0 || method >2
  error('Input argument "method" should be 0, 1 or 2');
end

u = levels(param); % Discretization levels

n = param(3);     % Size of matrix
N = length(cc);   % Total number of cycles

% Compute Histogram

dcc = cc2dcc(param,cc,ddef);
F = dcc2cmat(dcc,n);

% Smooth by using Kernel estimator ?

if method >= 1   
  [F,h] = smoothcmat(F,method,h,NOsubzero,alpha);
end
