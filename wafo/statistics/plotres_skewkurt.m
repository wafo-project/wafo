function lmr=plotres_skewkurt(varargin)
%PLOTRES_SKEWKURT Plot residual L-kurtosis to L-skewness ratio  vs thresholds
%
% CALL  lmr = plotres_skewkurt(data,options)
%
% lmr    = L-moment ratio for excesses over thresholds, u.
%
% data    = vector of data.
% options = options structure defining estimation of mrl. 
%        .Nmin : Minimum number of extremes to include. (Default Nmin = 3).
%        .umin : Minimum threshold (default min(data))
%        .umax : Maximum threshold (default max(data))
%        .Nu   : number of threshold values (default min(length(data),100))
%        .alpha: Confidence coefficient (default 0.05)
%        .plotflag:
%
% RES_SKEWKURT estimate L-moment ratio for excesses over thresholds. The 
% purpose of LMR is to determine the threshold where the upper tail of the 
% data can be approximated with the generalized Pareto distribution (GPD).
% The GPD is appropriate for the tail, if the LMR is one as function of
% threshold, u. Theoretically in the GPD model the ratio
%
%   KU*(5+SK)/(SK*(1+5*SK)) = 1
%
% where KU and SK are L-kurtosis and L-skewness, respectively.
%
% Example
%  opt = res_skewkurt('defaults'); % Get default options
%  opt.Nu = 20;
%  xn = load('sea.dat');
%  Ie = findpot(xn,0,5);
%  mrl = plotres_skewkurt(xn(Ie,2),opt);
%  
% See also reslife, fitgenpar, fitgenparrange, disprsnidx

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


lmr = res_skewkurt(varargin{:});
if ~isstruct(lmr)
  plot(lmr,'.'), hold on,
  plot(lmr.args,ones(size(lmr.args)),'k:'),hold off
  m = mean(lmr.data);
  s = std(lmr.data);
  ax = axis;
  ax(3) = m-1.5*s;
  ax(4) = m+1.5*s;
  axis(ax)
end
