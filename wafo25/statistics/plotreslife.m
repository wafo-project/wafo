function mrl=plotreslife(data,varargin)
%PLOTRESLIFE Plot Mean Residual Life (mean excess vs thresholds)
%
% CALL  [mrl,u] = plotreslife(data,options)
%
%  mrl    = Mean residual life values, i.e., mean excesses over thresholds, u.
%  u      = vector of thresholds
% data    = vector of data.
% options = options structure defining estimation of mrl. 
%        .Nmin : Minimum number of extremes to include. (Default Nmin = 3).
%        .umin : Minimum threshold (default min(data))
%        .umax : Maximum threshold (default max(data))
%        .Nu   : number of threshold values (default min(length(data),100))
%        .alpha: Confidence coefficient (default 0.05)
%
% PLOTRESLIFE displays mean excesses over thresholds. If the data comes from a
% generalized Pareto distribution then the plot is a linear function of u.
%
% Example
%  opt = plotreslife('defaults'); % Get default options
%  opt.Nu = 20;
%  R = rndgenpar(0.2,2,2,100,1);
%  plotreslife(R,opt);
%  
% See also fitgenpar, reslife


%
%     This program is free software; you can redistribute it and/or modify
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


mrl = reslife(data,varargin{:});

if ~isstruct(mrl)
  plot(mrl,'.')
end
