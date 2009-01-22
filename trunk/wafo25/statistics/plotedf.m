function [varargout] = plotedf(z,varargin)
%PLOTEDF Plot Empirical Distribution Function  
%          and optionally compare it with distribution G.
%
%  CALL:  F = plotedf(X,G,plotflag,sym);
%
%        F  = empirical distribution of X, two column matrix.
%        X  = data vector.
%        G  = cdf, two column matrix or FDATA object as 
%            returned from the FITXXX functions(optional).
%  plotflag = 0  no plotting
%             1 plot cdf F(x) (default)  
%             See pdfplot for more details on plotoptions 
%       sym = {s1,s2} cell array or comma separated list of plot 
%             symbols for F and G, respectively.
%             (default {'b','r--'})
% 
% NOTE:  SYM can be given anywhere after X
%
% Example:
%   R = rndgev(.8,1,11,1,100);
%   x = linspace(5,15,200);
%   plotedf(R,[x;cdfgev(x,.8,1,11)]','g')
%
% See also  edf, cumtrapz, plotedfcnd


% Copyright (C) 2000 WAFO-group
%
% This program is free software; you can redistribute it and/or modify
%     it under the terms of the GNU Lesser General Public License as published by
%     the Free Software Foundation; either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU Lesser General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.


% Tested on: Matlab 5.3, 5.2, 5.1
% History:
% revised pab 25.10.2000
%  - made call to plotedfcnd instead -> making maintainence easier.
% revised pab 07.03.2000
% - enabled so that f may be empty while plotflag is given 
% modified by Per A. Brodtkorb 10.11.98
% to accept both pdf and cdf's. Also enabled new plotting features,
% plotting of  probability of exceedances on a semilogy paper ....
% revised ms 13.06.2000
% - pdf usage removed (can't distinguish between a cdf and an increasing pdf)
% - changed axis labelling in figures
% - changed to staircaseplot when plotflag=2,3
% - revised header info
% - moved the conditional version to edfcnd


error(nargchk(1,5,nargin))
[varargout{1:nargout}] = plotedfcnd(z,-inf,varargin{:});









