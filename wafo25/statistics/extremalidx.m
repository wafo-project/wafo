function  [ei, Tc] = extremalidx(data)
%EXTREMALIDX Extremal Index measuring the dependence of data
%
%  CALL [EI,Tc] = extremalidx(data)
%
%  EI   = Extremal index
%  Tc   = minimum distance between clusters.
%  data = twocolumn matrix with sampled times in the first column
%         and values the second columns. Sampling frequency must be 1Hz.
%
%   EXTREMALIDX estimate the Extremal Index (EI) which is a measure of
%   independence. The EI is one if the data are independent and less than
%   one if there are some dependence. The extremal index can also be
%   intepreted as the reciprocal of the mean cluster size.
%
%
% Example
%  xn = load('sea.dat');
%  Ie = findpot(xn,0,5);
%  di = disprsnidx(xn(Ie,:),'Tb', 100);
%  plot(di) % a threshold around 1 seems appropriate.
%  Ie10 = findpot(xn,1.2,0);
%  [ei10,Tmin10] = extremalidx(Ie10);
%  Tmin = xn(Tmin10,1);
%  IeTmin = findpot(xn,1,Tmin);
%  ei = extremalidx(IeTmin); % should be closer to one
%
% See also reslife, fitgenparrange, disprsnidx, findpot, decluster


% Reference
% Christopher A. T. Ferro, Johan Segers (2003)
% Inference for clusters of extreme values
% Journal of the Royal Statistical Society: Series B (Statistical Methodology) 65 (2), 545–556.
% doi:10.1111/1467-9868.00401

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


[n m]  = size(data);


if min(m,n)==1,  % only times are given
 data  = data(:);
end
ei = [];
Tc = [];
switch numel(data)
  case 0, return
  case 1,ei = 1;  return
end

Ti = diff(sort(data(:,1))); % interexceedance times
ei = lcl_extrml_idx(Ti);


if nargout>1
 if ei==1
   Tc = min(Ti);
 else
   N = length(Ti)+1;
   C = floor(N*ei)+1;
   sTi = -sort(-Ti);
   Tc = sTi(min(C,N-1)); % declustering time
 end
end

function ei = lcl_extrml_idx(Ti)
Tmax = max(Ti); 
if Tmax<=1,
  ei = 0;
elseif Tmax<=2
  ei = min(1,2*mean(Ti).^2/mean(Ti.^2));
else
  ei = min(1,2*mean(Ti-1).^2/mean((Ti-1).*(Ti-2)));
end

