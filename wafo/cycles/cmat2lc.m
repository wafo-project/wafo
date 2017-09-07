function lc = cmat2lc(param,F)
%CMAT2LC Calculates the level crossings from a cycle matrix.
%
% CALL:  lc = cmat2lc(param,F);
% 
% Input: 
%   param = Parameter vector, [a b n], defines the discretization.
%   F     = Cycle matrix (e.g. rainflow matrix) [nxn]
% Output:
%   lc    = a two column matrix with levels and number of upcrossings.
%
% Example:
%  x = load('sea.dat'); 
%  TP = dat2tp(x); 
%  RFC = tp2rfc(TP); 
%  param = [-2 2 151]; 
%  F = cc2cmat(param, RFC);
%  lc = cmat2lc(param, F);
%  plot(lc(:,1), lc(:,2));
%
%  close all;
%
% See also  cc2cmat

% Tested on Matlab 6.0
%
% History:
% Revised by jr 01-Apr-2001
% - Example added
% - Updated help 
% Created by PJ (P�r Johannesson) 14-Jan-2000

% Check input arguments
ni = nargin;
no = nargout;
%error(nargchk(2,2,ni));
narginchk(2,2)
lc = [levels(param)' diag(cmat2nt(F))];
