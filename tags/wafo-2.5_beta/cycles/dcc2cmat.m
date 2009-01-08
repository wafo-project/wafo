function F = dcc2cmat(dcc,n)
% DCC2CMAT  Calculates the cycle matrix for a discrete cycle count.
%
% CALL:  F = dcc2cmat(dcc,n);
%
%   F      = Cycle matrix
%   dcc    = a two column matrix with a discrete cycle count.
%   n      = Number of discrete levels.
%
% The discrete cycle count takes values from 1 to n.
% 
% A cycle count is transformed into a discrete cycle count by
% using the function CC2DCC.
%
% See also  cc2cmat, cc2dcc, cmatplot

% Tested on Matlab 5.3
%
% History:
% Revised by PJ 12-Jan-2000
%   Corrected help text
% Revised by PJ (Pär Johannesson) 01-Nov-1999
%   updated for WAFO
% Copied from WAT Ver. 1.2

F=zeros(n);
k=length(dcc);
for i=1:k,
  F(dcc(i,1),dcc(i,2)) = F(dcc(i,1),dcc(i,2))+1;
end

