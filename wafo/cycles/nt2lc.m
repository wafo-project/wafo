function lc = nt2lc(param,NT)
%NT2LC Calculates the level crossings from a cycle matrix.
%
% CALL:  lc = nt2lc(param,NT);
%
% Input:
%   param = Parameter vector, [a b n], defines the discretization.
%   NT    = Coutning distribution. [nxn]
% Output;
%   lc    = a two column matrix with levels and number of upcrossings.
%
% Example: 
%   F0 = round(triu(rand(50),1)*10);
%   NT = cmat2nt(F0);
%   param = [-1 1 50];
%   lc = nt2lc(param,NT);
%   plotlc(lc);
%
%   close all;
%
% See also  cmat2nt, cc2cmat

% Tested on Matlab 6.0
%
% History:
% Created by PJ (P�r Johannesson) 18-May-2000

% Check input arguments
ni = nargin;
no = nargout;
error(nargchk(2,2,ni));

lc = [levels(param)' diag(NT)];
