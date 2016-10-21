function ccplot(varargin)
%CCPLOT Plots a cycle count as a point process in the plane.
%
% CALL: ccplot(cc,psize)
%
%   cc    = a two column matrix with cycles.
%   psize = point size, (optional, default value is 12).
%
% See also  tp2mm, tp2rfc, tp2tc

% Tested  on Matlab  5.3
%
% History:
% Revised by PJ 26-Jul-2000
%   Now works when cc is 4 column matrix.
% Revised by PJ (Pär Johannesson) 18-May-2000
%   When input cc is a cell-array, 
%   then each cell i plotted in a subplot. 
% Revised by PJ (Pär Johannesson) 01-Nov-1999
%   updated for WAFO
% Copied from WAT Ver. 1.2

% Check input arguments

plotcc(varargin{:})
